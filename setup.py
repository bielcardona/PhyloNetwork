from setuptools import setup, find_packages
setup(name="PhyloNetwork",
    version="1.0b1",
    author="Gabriel Cardona, David Sanchez",
    author_email="bielcardona@gmail.com, dscharles@gmail.com",
    license="BSD",
    keywords="phylogenetic trees networks",
    packages=['pyphylonetwork'],
    package_dir = {'pyphylonetwork': 'src/pyphylonetwork'},
    install_requires = [
        'networkx',
        'pyparsing',
        'numpy'
    ]
)
