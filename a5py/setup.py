from setuptools import setup, find_packages

setup(name='a5py',
      version='0.0',
      description='ASCOT5 python library',
      url='http://github.com/ascot/python/a5py',
      license='LGPL',
      packages=find_packages(),
      zip_safe=False,
      install_requires=[
          'numpy',
          'h5py',
          'prompt_toolkit'
      ],
      scripts=[
        'bin/a5removeruns',
        'bin/a5editoptions',
        'bin/a5setoptions',
        'bin/a5separatemarkers',
        'bin/a5poincare',
        'bin/a5combine',
        'bin/a5ascot4input',
        'bin/a5copyinput',
        'bin/a5continuerun',
        'bin/a5setactive',
        'bin/a5GUI',
        'bin/a5ls',
        'bin/a5vol',
        'bin/a5doxygen'
      ])
