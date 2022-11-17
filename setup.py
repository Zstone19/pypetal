from setuptools import setup


NAME = 'petl'
VERSION = '0.1.0'
URL = 'https://github.com/Zstone19/petl'

AUTHOR = 'Zachary Stone'
EMAIL = 'stone28@illinois.edu'

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
    packages=['petl'],
    package_dir={'petl':'./src/petl'},
    include_package_data=True,
    cmdclass={'install': InstallCommand,}
)