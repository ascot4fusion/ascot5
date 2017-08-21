from setuptools import setup, find_packages

setup(name='a5py',
      version='0.0',
      description='ASCOT5 python library',
      url='http://github.com/ascot/python/a5py',
      license='LGPL',
      packages=find_packages(),
      zip_safe=False,
      scripts=[
        'bin/a5removeruns'
      ])


