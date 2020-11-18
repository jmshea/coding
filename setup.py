import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="binext",
    version="0.94",
    author="John M. Shea",
    author_email="jshea@ece.ufl.edu",
    description="Opinionated module for working with binary Galois extension fields",
    long_description=long_description,
    long_description_content_type="text/markdown",
    keywords="finite fields binary extension galois minimal polynomials conjugates arithmetic",
    url="https://github.com/jmshea/coding/tree/main/binext",
    packages=['binext'],
    install_requires=[
        'numpy',
        'pandas'
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent"
    ],
    test_suite="tests"
)
