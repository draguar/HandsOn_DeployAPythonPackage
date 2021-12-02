#!/usr/bin/env python

import setuptools
from distutils.core import setup

setup(name='HandsOnDeployAPackage',
      version='1.7',
      description='Hands-on: deploy a Python package',
      author='Guilom',
      author_email='guillaume.vanel@protonmail.com',
      url='https://github.com/draguar/HandsOn_DeployAPythonPackage',
      py_modules=['ourCode'],
      license_files = '../LICENSE',
      install_requires = ['pandas>=1.3.4', 'matplotlib>=3.5.0'])
