# -*- coding: utf-8 -*-
import logging
import allel

logger = logging.getLogger()


def vcf_to_hdf5(vcf_path: str, hdf5_path: str):
    """ Conversion from VCF or VCF.gz to HDF5 file format.

    Parameters:
        vcf_path (str): Input path to the VCF file.
        hdf5_path (str): Output path to the HDF5 file.

    """
    logger.debug("Converting VCF file '{vcf_path}'".format(vcf_path=vcf_path))
    allel.vcf_to_hdf5(vcf_path, hdf5_path, fields='*', overwrite=True)
    logger.debug("HDF5 file stored in '{hdf5_path}'".format(hdf5_path=hdf5_path))
