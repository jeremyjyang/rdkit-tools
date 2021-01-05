import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="rdkwrapper",
    version="0.0.9",
    author="Jeremy Yang",
    author_email="jeremyjyang@gmail.com",
    description="Wrappers and helpers for RDKit.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/unmtransinfo/rdkit-tools",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
