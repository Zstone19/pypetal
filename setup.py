from setuptools import setup

NAME = 'pypetal'
VERSION = '0.1.0'
URL = 'https://github.com/Zstone19/pypetal'

AUTHOR = 'Zachary Stone'
EMAIL = 'stone28@illinois.edu'

#Get requirements.txt
REQUIREMENTS = []
with open('requirements.txt', 'r', encoding='UTF-8') as f:
    for line in f.readlines():
        REQUIREMENTS.append( line.strip('\n\r')  )


subdirs = ['', '.drw_rej', '.fromfile', '.pyccf', '.pyroa', '.pyzdcf', '.utils', '.weighting']

setup(
    name=NAME,
    version=VERSION,
    url='https://github.com/Zstone19/pypetal',
    author=AUTHOR,
    author_email=EMAIL,
    maintainer=AUTHOR,
    maintainer_email=EMAIL,
    license='MIT',
    install_requires=REQUIREMENTS,
    python_requires='>=3.8,<3.11',
    packages=['pypetal' + x for x in subdirs],
    package_dir={'pypetal':'./src/pypetal'},
    include_package_data=True
)
