# -*- coding: utf-8 -*-
from setuptools import setup

import haplotype_plot

setup(
    name='haplotype_plot',
    version=haplotype_plot.version,
    packages=['haplotype_plot', 'haplotype_plot.tests'],
    url='https://github.com/neobernad/haplotype_plot',
    license='MIT',
    author='José Antonio Bernabé Díaz',
    author_email='joseantonio.bernabe1@gmail.com',
    description='Generates haplotype plots from VCF or H5 files.', install_requires=['scikit-allel', 'h5py', 'PyVCF',
                                                                                     'numpy', 'matplotlib', 'seaborn']
)
