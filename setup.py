# -*- coding: utf-8 -*-
import pathlib
from setuptools import setup
import haplotype_plot

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()

setup(
    name='haplotype_plot',
    version=haplotype_plot.version,
    python_requires='>=3',
    long_description=README,
    long_description_content_type="text/markdown",
    packages=['haplotype_plot', 'haplotype_plot.tests'],
    url='https://github.com/neobernad/haplotype_plot',
    license='MIT',
    author='José Antonio Bernabé-Díaz',
    author_email='joseantonio.bernabe1@um.com',
    entry_points={
        'console_scripts': [
            'haplotype_plot = haplotype_plot.main:main',
        ]
    },
    description='Generates haplotype plots from VCF files.',
    install_requires=['scikit-allel', 'h5py', 'PyVCF', 'numpy', 'matplotlib', 'seaborn', 'pathlib', 'Cython'],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Topic :: Software Development :: Build Tools',
        'Topic :: Scientific/Engineering',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: Unix',
        'Operating System :: MacOS',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.6'
    ]
)
