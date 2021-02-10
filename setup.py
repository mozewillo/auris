#!/usr/bin/env python
# -*- coding: utf-8 -*-

import setuptools
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()
setuptools.setup(
    name="auris",
    version="1.0",
    author="Aleksandra Możwiłło",
    author_email="alexmozwillo@gmail.com",
    description="Package for identifying common mutations and analyzing "
                "genomes for the existence of specific substrings.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github',
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
    ],
    python_requires='>=3.8.0',
    install_requires=['numpy', 'pandas', 'regex', 'matplotlib>=3', 'argparse'],
)
