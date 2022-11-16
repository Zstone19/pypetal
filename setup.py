import os
import subprocess

from setuptools import find_packages, setup

import distutils.cmd
import distutils.log

class InstallDependencies(distutils.cmd.Command):
    """Install dependencies for PETL""" 
    
    description = 'Install dependencies for PETL'
    user_options = [
        ('user=', 'u', 'if True, will install locally'),
        ('fcompile=', 'f', 'Fortran compiler used for JAVELIN'),
        ('plike=', 'p', 'if True, will install PLIKE'),
    ]   

    def initialize_options(self):
        self.user = False
        self.fcompile = None
        self.plike = False
        
    def finalize_options(self):
        self.cwd = os.getcwd()
        
    def run(self):
        """Run command"""
        
        command = ['sh', 'build_dep.sh']
        
        if self.user:
            command.append('-u true')
        else:
            command.append('-u false')
            
        if self.fcompile is not None:
            command.append('-f %s' % self.fcompile )
            
        if self.plike:
            command.append('-p true')
        else:
            command.append('-p false')
            
        subprocess.check_call(command)





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
    cmdclass={'install_dep': InstallDependencies,}
    
)