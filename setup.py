from setuptools import setup, find_packages

# Import the version number
from cansnpscorer import __version__

setup(
    name='CanSNPCreator',
    version=__version__,
    description='A tool for generating and evaluating CanSNPs',
    author='Andreas Sjödin',
    author_email='your.email@example.com',
    packages=find_packages(),  # Automatically discover and include all packages in the package directory
    install_requires=[
        # Add your project dependencies here
        # For example: 'numpy>=1.18.5',
    ],
    entry_points={
        'console_scripts': [
            'CanSNPCreator=cansnpscorer.cansnp_scorer:main',
        ],
    },
)
