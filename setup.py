import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="astroutils", # Replace with your own username
    version="1.0.0",
    author="Héctor Salas O.",
    author_email="hector.salas.o@gmail.com",
    description="A small pakage with scripts for astronomy",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/hsalas/astroutils",
    packages=setuptools.find_packages(),
    install_requires=['aplpy', 'astropy', 'readline', 'numpy', 'matplotlib'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
