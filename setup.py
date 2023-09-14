from os import path
from setuptools import setup, find_packages

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'requirements.txt')) as requirements_file:
    # Parse requirements.txt, ignoring any commented-out lines.
    requirements = [line for line in requirements_file.read().splitlines()
                    if not line.startswith('#')]

setup(name='a5py',
      version='5.5',
      description='ASCOT5',
      url='https://github.com/ascot4fusion/ascot5',
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
        'bin/a5makecompatible',
      ],)
