import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="wgregseq",
    version="0.0.1",
    author="Tom Roeschinger, Scott Saunders, Manuel Razo-Meja",
    author_email="troeschi@caltech.edu",
    description="This repository contains all active research materials for the Reg-Seq project.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/tomroesch/Reg-Seq2",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)