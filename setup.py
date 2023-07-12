from setuptools import setup, find_packages, Extension
import numpy, sys
import re

# auto-updating version code stolen from RadVel
def get_property(prop, project):
    result = re.search(
        r'{}\s*=\s*[\'"]([^\'"]*)[\'"]'.format(prop),
        open(project + "/__init__.py").read(),
    )
    return result.group(1)

def get_requires():
    reqs = []
    for line in open("requirements.txt", "r").readlines():
        reqs.append(line)
    return reqs

setup(
    name="DEBRaVE",
    version=get_property("__version__", "orbitize"),
    description="DEBRaVE: Detached Eclipse Binary Radial Velocity Extractor",
    url="https://github.com/ganesh1076/DEBRaVE",
    author="",
    author_email="",
    packages=find_packages(),
    keywords="Astronomy Binaries Radial Velocity Spectroscopy TODCOR",
    install_requires=get_requires(),
)