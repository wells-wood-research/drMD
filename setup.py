## Package metadata and dependencies live in pyproject.toml (managed with uv).
## This shim only keeps setuptools package discovery for legacy builds.
from setuptools import setup, find_packages

setup(packages=find_packages())
