#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()
print(readme)

requirements = ['pycolt>=0.3', 'cclib>=1.7', 'numpy>=1.17', 'openbabel>=3.1', 'orbkit>=1.1']

setup_requirements = ['pytest-runner', ]

test_requirements = ['pytest>=3', ]

setup(
    author="Felix Plasser",
    python_requires='>=3.6, !=3.0.*, !=3.1.*, !=3.2.*, !=3.3.*, !=3.4.*, !=3.5*',
    classifiers=[
        'Intended Audience :: Developers',
        'License :: OSI Approved :: GNU General Public License 3.0',
        'Natural Language :: English',
        'Programming Language :: Python :: >=3.6',
    ],
    description="Theoretical Density, Orbital Relaxation and Exciton analysis",
    install_requires=requirements,
    license="GNU General Public License 3.0",
    long_description=readme,
    include_package_data=True,
    keywords='theodore',
    name='theodore',
    packages=find_packages(include=['theodore', 'theodore.*']),
    entry_points = {
        'console_scripts': ['theodore=theodore:run',]
        },
    setup_requires=setup_requirements,
    tests_require=test_requirements,
    url='https://github.com/felixplasser/theodore-qc',
    version='2.4.5',
    zip_safe=False,
)
