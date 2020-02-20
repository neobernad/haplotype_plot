# -*- coding: utf-8 -*-
from setuptools import setup

import haplotyper

setup(
    name='haplotype_plot',
    version=haplotyper.version,
    packages=['haplotyper', 'haplotyper.tests'],
    url='https://github.com/neobernad/haplotype_plot',
    license='MIT',
    author='José Antonio Bernabé Díaz',
    author_email='joseantonio.bernabe1@gmail.com',
    description='Generates haplotype plots from VCF or H5 files.', install_requires=['scikit-allel', 'h5py']
)
