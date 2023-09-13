from os import path
from setuptools import setup, find_packages

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'requirements.txt')) as requirements_file:
    # Parse requirements.txt, ignoring any commented-out lines.
    requirements = [line for line in requirements_file.read().splitlines()
                    if not line.startswith('#')]

setup(name='a5py',
      version='0.0',
      description='ASCOT5 python library',
      url='http://github.com/ascot/python/a5py',
      license='LGPL',
      packages=find_packages(),
      zip_safe=False,
      install_requires=requirements,
      scripts=[
        'bin/a5removegroup',
        'bin/a5copygroup',
        'bin/a5editoptions',
        'bin/a5combine',
        'bin/a5ascot4input',
        'bin/a5setactive',
        'bin/a5gui',
        'bin/a5ls',
        'bin/a5makecompatible',
        'bin/a5doxygen',
        'bin/test_ascot.py'
      ],
      package_data={'a5py/ascotpy': ['lib.so']},
      include_package_data=True)
