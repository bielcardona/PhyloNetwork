from setuptools import setup
import versioneer
from os import path

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md')) as f:
    long_description = f.read()

setup(name="phylonetwork",
      version=versioneer.get_version(),
      cmdclass=versioneer.get_cmdclass(),
      author="Gabriel Cardona",
      author_email="bielcardona@gmail.com",
      license="BSD",
      keywords="phylogenetic trees networks",
      packages=['phylonetwork'],
      # package_dir = {'phylonetwork': 'src/phylonetwork'},
      python_requires='>=3.7',
      install_requires=[
          'networkx>=2',
          'pyparsing',
          'numpy',
          'cached_property'
      ],
      url='https://github.com/bielcardona/PhyloNetworks',
      long_description=long_description,
      long_description_content_type='text/markdown'
      )
