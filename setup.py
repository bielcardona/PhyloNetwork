from setuptools import setup
import versioneer

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
          'networkx>=2,<3',
          'pyparsing',
          'numpy'
      ]
      )
