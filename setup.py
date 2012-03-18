from setuptools import setup, find_packages
setup(name="phylonetwork",
    version="1.0b2",
    author="Gabriel Cardona, David Sanchez",
    author_email="bielcardona@gmail.com, dscharles@gmail.com",
    license="BSD",
    keywords="phylogenetic trees networks",
    packages=['phylonetwork'],
    package_dir = {'phylonetwork': 'src/phylonetwork'},
    install_requires = [
        'networkx',
        'pyparsing',
        'numpy'
    ]
)
