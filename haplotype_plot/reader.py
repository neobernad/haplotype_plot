# -*- coding: utf-8 -*-
import logging
import vcf  # from PyVCF
import allel
import os

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


def get_samples(vcf_path: str) -> list:
    """ Reads the sample list from a VCF file.

        Parameters:
            vcf_path (str): Input path to the VCF file.
        Returns:
            A list of strings (samples names)
    """
    sample_list = list()
    vcf_file = open(vcf_path, 'r')
    vcf_reader = vcf.Reader(vcf_file)
    sample_list = vcf_reader.samples
    vcf_file.close()
    return sample_list


def vcf_to_hdf5(vcf_path: str, hdf5_path: str):
    """ Conversion from VCF or VCF.gz to HDF5 file format.

    Parameters:
        vcf_path (str): Input path to the VCF file.
        hdf5_path (str): Output path to the HDF5 file.

    """
    # sample_list = get_samples(vcf_path)
    logger.debug("Converting VCF file '{vcf_path}'".format(vcf_path=vcf_path))
    if os.path.exists(hdf5_path):
        logger.debug("File '{hdf5_path}' already exists. I will remove it.".format(hdf5_path=hdf5_path))
        os.remove(hdf5_path)
    allel.vcf_to_hdf5(vcf_path, hdf5_path, fields='*', overwrite=False)
    logger.debug("HDF5 file stored in '{hdf5_path}'".format(hdf5_path=hdf5_path))
