from setuptools import setup, find_packages

setup(
    name='CanSNPScorer',
    version='0.1.0',
    description='A tool for evaluating CanSNPs',
    author='Your Name',
    author_email='your.email@example.com',
    packages=find_packages(),  # Automatically discover and include all packages in the package directory
    install_requires=[
        # Add your project dependencies here
        # For example: 'numpy>=1.18.5',
    ],
    entry_points={
        'console_scripts': [
            'cansnp-evaluator=cansnpscorer.cansnp_evaluator:main',
            'cansnp-extractor=cansnpscorer.cansnp_extractor:main',
        ],
    },
)

