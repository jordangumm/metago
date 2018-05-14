import sys
import os
import random
import string

import click
from pyflow import WorkflowRunner, CommandTaskRunner, TaskManager, WorkflowRunnerThreadSharedData, QCaller
from pyflow import *
from subprocess import call


class FluxQCaller(QCaller):
    def run(self):
        GlobalSync.subprocessControl.acquire()
        try :
            tmp_proc = subprocess.Popen(' '.join(self.cmd), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
            self.lock.acquire()
            try:
                self.proc = tmp_proc
                # handle the case where Popen was taking its good sweet time and a killProc() was sent in the meantime:
                if self.is_kill_attempt: self.killProc()
            finally:
                self.lock.release()

            if self.is_kill_attempt: return

            for line in self.proc.stdout :
                self.results.outList.append(line)
            self.results.retval = self.proc.wait()
        finally:
            GlobalSync.subprocessControl.release()
        self.results.isComplete = True


class FluxRunMode(object):
    data = { "local" : ModeInfo(defaultCores=1,
                                defaultMemMbPerCore=siteConfig.defaultHostMemMbPerCore,
                                defaultIsRetry=False),
             "sge"   : ModeInfo(defaultCores=getSGEJobsDefault(),
                                defaultMemMbPerCore="unlimited",
                                defaultIsRetry=True),
             "flux"   : ModeInfo(defaultCores=getSGEJobsDefault(),
                                defaultMemMbPerCore="unlimited",
                                defaultIsRetry=True) }


class FluxRetryParam(RetryParam):
    def _finalize(self):
        """
        decide whether to turn retry off based on retry and run modes:
        """
        if (self._retry_mode == "nonlocal") and \
                (not FluxRunMode.data[self._run_mode].defaultIsRetry) :
            self.max = 0
        else :
            self.max = int(self._retry_max)


class FluxWorkflowRunner(WorkflowRunner):
    def _cdata(self) :
        # We're doing this convoluted setup only to avoid having a
        # ctor for ease of use by the client. See what pyFlow goes
        # through for you client code??
        #
        try:
            return self._constantData
        except AttributeError:
            self._constantData = FluxWorkflowRunnerThreadSharedData()
            return self._constantData

    def _startTaskManager(self) :
        # Start a new task manager if one isn't already running. If it is running
        # provide a hint that a new task has just been added to the workflow.
        #
        if (self._tman is not None) and (self._tman.isAlive()) :
            self._tdag.isFinishedEvent.set()
            return
        if not self._cdata().isTaskManagerException :
            self._tman = FluxTaskManager(self._cdata(), self._tdag)
            self._tman.start()


class FluxWorkflowRunnerThreadSharedData(WorkflowRunnerThreadSharedData):
    @staticmethod
    def _validateFixParam(param):
        """
        validate and refine raw run() parameters for use by workflow
        """

        param.mailTo = setzer(param.mailTo)
        param.schedulerArgList = lister(param.schedulerArgList)
        if param.successMsg is not None :
            if not isString(param.successMsg) :
                raise Exception("successMsg argument to WorkflowRunner.run() is not a string")

        # create combined task retry settings manager:
        param.retry=FluxRetryParam(param.mode,
                               param.retryMax,
                               param.retryWait,
                               param.retryWindow,
                               param.retryMode)

        # setup resource parameters
        if param.nCores is None :
            param.nCores = FluxRunMode.data[param.mode].defaultCores

        # ignore total available memory settings in non-local modes:
        if param.mode != "local" :
            param.memMb = "unlimited"

        if param.mode in ("sge", "flux") :
            if siteConfig.maxSGEJobs != "unlimited" :
                if ((param.nCores == "unlimited") or
                    (int(param.nCores) > int(siteConfig.maxSGEJobs))) :
                    param.nCores = int(siteConfig.maxSGEJobs)

        if param.nCores != "unlimited" :
            param.nCores = int(param.nCores)
            if param.nCores < 1 :
                raise Exception("Invalid run mode nCores argument: %s. Value must be 'unlimited' or an integer no less than 1" % (param.nCores))

        if param.memMb is None :
            if param.nCores == "unlimited" :
                param.memMb = "unlimited"
            mpc = FluxRunMode.data[param.mode].defaultMemMbPerCore
            if mpc == "unlimited" :
                param.memMb = "unlimited"
            else :
                param.memMb = mpc * param.nCores
        elif param.memMb != "unlimited" :
            param.memMb = int(param.memMb)
            if param.memMb < 1 :
                raise Exception("Invalid run mode memMb argument: %s. Value must be 'unlimited' or an integer no less than 1" % (param.memMb))

        # verify/normalize input settings:
        if param.mode not in FluxRunMode.data.keys() :
            raise Exception("Invalid mode argument '%s'. Accepted modes are {%s}." \
                            % (param.mode, ",".join(FluxRunMode.data.keys())))

        if param.mode == "sge" :
            # TODO not-portable to windows (but is this a moot point -- all of sge mode is non-portable, no?):
            def checkSgeProg(prog) :
                proc = subprocess.Popen(("which", prog), stdout=open(os.devnull, "w"), shell=False)
                retval = proc.wait()
                if retval != 0 : raise Exception("Run mode is sge, but no %s in path" % (prog))
            checkSgeProg("qsub")
            checkSgeProg("qstat")

        if param.mode == "flux":
            def checkSgeProg(prog) :
                proc = subprocess.Popen(("which", prog), stdout=open(os.devnull, "w"), shell=False)
                retval = proc.wait()
                if retval != 0 : raise Exception("Run mode is sge, but no %s in path" % (prog))
            checkSgeProg("qsub")
            checkSgeProg("qstat")

        stateDir = os.path.join(param.dataDir, "state")
        if param.isContinue == "Auto" :
            param.isContinue = os.path.exists(stateDir)

        if param.isContinue :
            if not os.path.exists(stateDir) :
                raise Exception("Cannot continue run without providing a pyflow dataDir containing previous state.: '%s'" % (stateDir))

        for email in param.mailTo :
            if not verifyEmailAddy(email):
                raise Exception("Invalid email address: '%s'" % (email))


class FluxTaskManager(TaskManager):
    """
    This class runs on a separate thread from workflowRunner,
    launching jobs based on the current state of the TaskDAG
    """
    def _getCommandTaskRunner(self, task) :
        """
        assist launch of a command-task
        """

        # shortcuts:
        payload = task.payload
        param = self._cdata.param

        if payload.cmd.cmd is None :
            # Note these should have been marked off by the TaskManager already:
            raise Exception("Attempting to launch checkpoint task: %s" % (task.fullLabel()))

        isForcedLocal = ((param.mode != "local") and (payload.isForceLocal))

        # mark task resources as occupied:
        if not isForcedLocal :
            if self.freeCores != "unlimited" :
                if (self.freeCores < payload.nCores) :
                    raise Exception("Not enough free cores to launch task")
                self.freeCores -= payload.nCores

            if self.freeMemMb != "unlimited" :
                if (self.freeMemMb < payload.memMb) :
                    raise Exception("Not enough free memory to launch task")
                self.freeMemMb -= payload.memMb

        if payload.mutex is not None :
            self.taskMutexState[payload.mutex] = True

        TaskRunner = None
        if param.mode == "local" or payload.isForceLocal or payload.isCmdMakePath :
            TaskRunner = LocalTaskRunner
        elif param.mode == "sge" :
            TaskRunner = SGETaskRunner
        elif param.mode == "flux" :
            TaskRunner = FluxTaskRunner
        else :
            raise Exception("Can't support mode: '%s'" % (param.mode))

        #
        # TODO: find less hacky way to handle make tasks:
        #
        taskRetry = payload.retry

        if payload.isCmdMakePath :
            taskRetry = copy.deepcopy(payload.retry)
            taskRetry.window = 0

            if param.mode == "local" or payload.isForceLocal :
                launchCmdList = ["make", "-j", str(payload.nCores)]
            elif param.mode == "sge" :
                launchCmdList = siteConfig.getSgeMakePrefix(payload.nCores, payload.memMb, param.schedulerArgList)
            elif param.mode == "flux":
                launchCmdList = siteConfig.getSgeMakePrefix(payload.nCores, payload.memMb, param.schedulerArgList) 
            else :
                raise Exception("Can't support mode: '%s'" % (param.mode))

            launchCmdList.extend(["-C", payload.cmd.cmd])
            payload.launchCmd = Command(launchCmdList, payload.cmd.cwd, payload.cmd.env)

        #
        # each commandTaskRunner requires a unique tmp dir to write
        # wrapper signals to. TaskRunner will create this directory -- it does not bother to destroy it right now:
        #

        # split the task id into two parts to keep from adding too many files to one directory:
        tmpDirId1 = "%03i" % ((int(task.id) / 1000))
        tmpDirId2 = "%03i" % ((int(task.id) % 1000))
        taskRunnerTmpDir = os.path.join(self._cdata.wrapperLogDir, tmpDirId1, tmpDirId2)

        return TaskRunner(task.runStatus, self._cdata.getRunid(),
                          task.fullLabel(), payload.launchCmd,
                          payload.nCores, payload.memMb,
                          taskRetry, param.isDryRun,
                          self._cdata.taskStdoutFile,
                          self._cdata.taskStderrFile,
                          taskRunnerTmpDir,
                          param.schedulerArgList,
                          self._cdata.flowLog,
                          task.setRunstate)

    @lockMethod
    def harvestTasks(self) :
        """
        Check the set of running tasks to see if they've completed and update
        Node status accordingly:
        """
        notrunning = set()
        for task in self.runningTasks.keys() :
            if self.stopped() : break
            trun = self.runningTasks[task]
            if not task.runStatus.isComplete.isSet() :
                if trun.isAlive() : continue
                # if not complete and thread is dead then we don't know what happened, very bad!:
                task.errorstate = 1
                task.errorMessage = "Thread: '%s', has stopped without a traceable cause" % (trun.getName())
            else :
                task.errorstate = task.runStatus.errorCode
                task.errorMessage = task.runStatus.errorMessage

            if task.errorstate == 0 :
                task.setRunstate("complete")
            else:
                task.setRunstate("error")

            notrunning.add(task)

            if not task.isError() :
                self._infoLog("Completed %s: '%s' launched from %s" % (task.payload.desc(), task.fullLabel(), namespaceLabel(task.namespace)))
            else:
                msg = task.getTaskErrorMsg()

                if self._cdata.isTaskSubmissionActive() :
                    # if this is the first error in the workflow, then
                    # we elaborate a bit on the workflow's response to
                    # the error. We also send any email-notifications
                    # for the first error only:
                    msg.extend(["Shutting down task submission. Waiting for remaining tasks to complete."])

                self._errorLog(msg)
                if self._cdata.isTaskSubmissionActive() :
                    self._cdata.emailNotification(msg, self._flowLog)

                # Be sure to send notifications *before* setting error
                # bits, because the WorkflowRunner may decide to
                # immediately shutdown all tasks and pyflow threads on
                # the first error:
                self._cdata.setTaskError(task)

        # recover task resources:
        for task in notrunning :
            self._removeTaskFromRunningSet(task)


class FluxTaskRunner(SGETaskRunner):


    def getFullCmd(self):
        # qsub options
        #
        qsubCmd = ['echo', '"{} {}"'.format(sys.executable, ' '.join(self.wrapperCmd)), '|',
                 "qsub",
                 "-V",  # import environment variables from shell
                 "-q", 'fluxm',
                 "-o", self.wrapFile,
                 "-e", self.wrapFile]

        qsubCmd.extend(self.schedulerArgList)

        return tuple(qsubCmd)

    def runOnce(self, retInfo) :

        def qcallWithTimeouts(cmd, maxQcallAttempt=1) :
            maxQcallWait = 180
            qcall = None
            for i in range(maxQcallAttempt) :
                qcall = FluxQCaller(cmd,self.infoLog)
                qcall.start()
                qcall.join(maxQcallWait)
                if not qcall.isAlive() : break
                self.infoLog("Trial %i of sge command has timed out. Killing process for cmd '%s'" % ((i + 1), cmd))
                qcall.killProc()
                self.infoLog("Finished attempting to kill sge command")

            return qcall.results

        # 1) call qsub, check for errors and retrieve taskId:
        #
        if os.path.isfile(self.wrapFile): os.remove(self.wrapFile)

        # write extra info, just in case we need it for post-mortem debug:
        qsubFile = os.path.join(os.path.dirname(self.wrapFile), "qsub.args.txt")
        if os.path.isfile(qsubFile): os.remove(qsubFile)
        qsubfp = open(qsubFile, "w")
        for arg in self.getFullCmd() :
            qsubfp.write(arg + "\n")
        qsubfp.close()

        results = qcallWithTimeouts(self.getFullCmd())

        isQsubError = False
        self.jobId = None
        if len(results.outList) != 1 :
            isQsubError = True
        else :
            w = results.outList[0].split('.')
            if (len(w) > 3) and (w[2] == "arc-ts") and (w[3] == "umich") :
                self.setNewJobId(int(w[0]))
            else :
                isQsubError = True

        if not results.isComplete :
            self._killJob()  # just in case...
            retInfo.taskExitMsg = ["Job submission failure -- qsub command timed-out"]
            return retInfo

        if isQsubError or (self.jobId is None):
            retInfo.taskExitMsg = ["Unexpected qsub output. Logging %i line(s) of qsub output below:" % (len(results.outList)) ]
            retInfo.taskExitMsg.extend([ "[qsub-out] " + line for line in results.outList ])
            return retInfo

        if results.retval != 0 :
            retInfo.retval = results.retval
            retInfo.taskExitMsg = ["Job submission failure -- qsub returned exit code: %i" % (retInfo.retval)]
            return retInfo

        # No qsub errors detected and a flux job_number is acquired -- success!
        self.infoLog("Task submitted to flux queue with job_number: %i" % (self.jobId))


        # 2) poll jobId until flux indicates it's not running or queued:
        #
        queueStatus = Bunch(isQueued=True, runStartTimeStamp=None)

        def checkWrapFileRunStart(result) :
            """
            check wrapper file for a line indicating that it has transitioned from queued to
            running state. Allow for NFS delay or incomplete file
            """
            if not os.path.isfile(self.wrapFile) : return
            for line in open(self.wrapFile) :
                w = line.strip().split()
                if (len(w) < 6) or (w[4] != "[wrapperSignal]") :
                    # this could be incomplete flush to the signal file, so
                    # don't treat it as error:
                    return
                if w[5] == "taskStart" :
                    result.runStartTimeStamp = timeStrToTimeStamp(w[0].strip('[]'))
                    result.isQueued = False
                    return


        # exponential polling times -- make small jobs responsive but give sge a break on long runs...
        ewaiter = ExpWaiter(5, 1.7, 60)

        pollCmd = ("/bin/bash", "--noprofile", "-o", "pipefail", "-c", "qstat | grep {}".format(self.jobId, self.jobId))
        while not self.stopped():
            results = qcallWithTimeouts(pollCmd, 6)
            isQstatError = False
            job_completed = False
            try:
                tokens = results.outList[0].split()
            except:
                continue

            if len(tokens) != 6: continue
            if tokens[4] == 'C': break # task is complete

            if isQstatError :
                if not results.isComplete :
                    retInfo.taskExitMsg = ["The qstat command for sge job_number %i has timed out for all attempted retries" % (self.jobId)]
                    self._killJob()
                else :
                    retInfo.taskExitMsg = ["Unexpected qstat output or task has entered sge error state. Sge job_number: %i" % (self.jobId)]
                    retInfo.taskExitMsg.extend(["Logging %i line(s) of qstat output below:" % (len(results.outList)) ])
                    retInfo.taskExitMsg.extend([ "[qstat-out] " + line for line in results.outList ])
                    # self._killJob() # leave the job there so the user can better diagnose whetever unexpected pattern has occurred
                return retInfo

            # also check to see if job has transitioned from queued to running state:
            if queueStatus.isQueued :
                checkWrapFileRunStart(queueStatus)
                if not queueStatus.isQueued :
                    self.setRunstate("running", queueStatus.runStartTimeStamp)

            ewaiter.wait()

        if self.stopped() :
            # self._killJob() # no need, job should already have been killed at the stop() call...
            return retInfo

        lastJobId = self.jobId

        # if we've correctly communicated with SGE, then its roll is done here
        # if a job kill is required for any of the error states above, it needs to be
        # added before this point:
        self.jobId = None

        wrapResult = self.getWrapFileResult()

        if wrapResult.taskExitCode is None :
            retInfo.taskExitMsg = ["Sge job_number: '%s'" % (lastJobId)]
            retInfo.taskExitMsg.extend(self.getWrapperErrorMsg())
            retInfo.retval = 1
            return retInfo
        elif wrapResult.taskExitCode != 0 :
            retInfo.taskExitMsg = self.getExitMsg()

        retInfo.retval = wrapResult.taskExitCode
        retInfo.isAllowRetry = True

        # success! (for sge & taskWrapper, but maybe not for the task...)
        return retInfo
