from setuptools import setup, find_packages

setup(
    name="sccmec_typer",
    version="0.1.0",
    packages=find_packages(),
    scripts=['bin/sccmec_typer.py'],
    install_requires=[
        "pandas",
    ],
    author="Your Name",
    description="A minimap2-based tool for SCCmec typing",
)
