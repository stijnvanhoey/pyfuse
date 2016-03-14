#!/usr/bin/env python

from setuptools import setup

setup(name='pyfuse',
      version='1.0',
      description='Python package for creating and running lumped hydrological models',
      url='https://github.com/stijnvanhoey/pyfuse',
      author='S. Van Hoey',
      author_email='stijnvanhoey@gmail.com',
      packages=['pyfuse'],
      license='BSD 3-clause New or Revised License',
      install_requires=['matplotlib', 'numpy', 'scipy'],
      keywords='hydrological modelling, fuse, Clark, model structure',)
