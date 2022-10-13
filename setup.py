import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="rdktools",
    version="0.9.5",
    author="Jeremy Yang",
    author_email="jeremyjyang@gmail.com",
    description="Tools and helpers for RDKit.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/jeremyjyang/rdkit-tools",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.8',
)
