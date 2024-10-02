from setuptools import setup, find_packages

setup(
    name='ArchaeoPyDating',  
    version='0.2',  
    packages=find_packages(),  
    install_requires=[                    
        'numpy',                          
        'matplotlib',                     
        'scipy',                          
        'future',                         
        'pyshp',                          
    ],
    description='Python tool for archaeomagnetic dating',
    author='Mario Serrano',  
    author_email='marioser@ucm.es',  
    url='https://github.com/Mariossb/ArchaeoPyDating',
    package_data={
        'ArchaeoPyDating': ['curves/**/*', 'external/*'],
    }
)