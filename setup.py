from setuptools import setup, find_packages, Extension
import numpy, sys
import re

setup(
    name="DEBRaVE",
    version="0.0.1",
    description="DEBRaVE: Detached Eclipse Binary Radial Velocity Extractor",
    url="https://github.com/ganesh1076/DEBRaVE",
    author="",
    author_email="",
    packages=find_packages(),
    keywords="Astronomy Binaries Radial Velocity Spectroscopy TODCOR"
)