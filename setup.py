import os
import subprocess
import sys

from setuptools import find_packages, setup
from setuptools.command.install import install


class InstallCommand(install):
    """Install dependencies for PETL""" 
    
    description = 'Install dependencies for PETL (pyCCF, JAVELIN, PLIKE)'
    user_options = install.user_options + [
        ('ujav=', None, 'if True, will install locally'),
        ('fcompile=', None, 'Fortran compiler used for JAVELIN'),
        ('plike=', None, 'if True, will install PLIKE'),
    ]
    
    def initialize_options(self):
        install.initialize_options(self)
        self.ujav = True
        self.fcompile = None
        self.plike = False
        
    def finalize_options(self):
        install.finalize_options(self)
        
    def run(self):
        """Run command"""
        
        command = ['sh', 'build_dep.sh']
        
        if self.ujav:
            command.append('-u true')
        else:
            command.append('-u false')
            
        if self.fcompile is not None:
            command.append('-f %s' % fcompile )
            
        if self.plike:
            command.append('-p true')
        else:
            command.append('-p false')
            
        command.append('>> build_dep.log')
            
        subprocess.check_call(command)
        install.run(self)





NAME = 'petl'
VERSION = '0.1.0'
URL = 'https://github.com/Zstone19/petl'

AUTHOR = 'Zachary Stone'
EMAIL = 'stone28@illinois.edu'
PACKAGES = find_packages(where='src')

#Get requirements.txt
f = open('requirements.txt', 'r')
REQUIREMENTS = []

for line in f.readlines():
    REQUIREMENTS.append( line.strip('\n\r')  )
f.close()



setup(
    name=NAME,
    version=VERSION,
    url='https://github.com/Zstone19/petl',
    author=AUTHOR,
    author_email=EMAIL,
    maintainer=AUTHOR,
    maintainer_email=EMAIL,
    license='MIT',
    install_requires=REQUIREMENTS,
    python_requires='>=3.8,<3.11',
    packages=PACKAGES,
    package_dir={'':'src'},
    include_package_data=True,
    cmdclass={'install': InstallCommand,}
)