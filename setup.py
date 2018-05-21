from setuptools import setup, find_packages

setup(
    name='metaGO',
    version='0.1',
    packages=find_packages(),
    scripts=['bin/metago'],

    install_requires=[],

    author='Jordan Gumm',
    author_email='jordan@variantanalytics.com',
    description='metagenomic cli for automated analysis'
)
