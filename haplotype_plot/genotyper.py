# -*- coding: utf-8 -*-
import logging
import allel
import h5py
import os
import os.path
import haplotype_plot.constants as constants
import haplotype_plot.reader as reader
import haplotype_plot.filter as strainer
import haplotype_plot.haplotyper as haplotyper
import numpy as np

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


def debug_hdf5(hdf5_path: str):
    callset = h5py.File(hdf5_path, mode='r')
    logger.debug("Callset keys: {keys}".format(keys=callset.keys()))
    logger.debug("Callset keys of {key}: '{values}'".format(key=constants.HDF5_CALLDATA_KEY,
                                                            values=callset[constants.HDF5_CALLDATA_KEY].keys()))
    logger.debug("Callset keys of {key}: '{values}'".format(key=constants.HDF5_VARIANTS_KEY,
                                                            values=callset[constants.HDF5_VARIANTS_KEY].keys()))


def get_genotypes_n_variants(callset: h5py.File) -> (allel.GenotypeChunkedArray, allel.VariantChunkedTable):
    """ Returns a the genotypes and variants tables from an hdf5.

    Parameters:
        hdf5_path (str): Input path to the HDF5 file.
    Returns:
        Tuple(allel.GenotypeChunkedArray, allel.VariantChunkedTable)
            - allel.GenotypeChunkedArray as the genotypes table
            - allel.VariantChunkedTable as the variants table with columns 'CHROM' and 'POS'
    """
    logger.debug("Retrieving genotypes via '{key}' key".format(key=constants.HDF5_CALLDATA_GENOTYPE_KEY))
    genotypes = allel.GenotypeChunkedArray(callset['calldata/GT'])
    logger.debug("Retrieving variants via '{key}' key'".format(key=constants.HDF5_VARIANTS_KEY))
    variants = allel.VariantChunkedTable(callset["variants"], names=['CHROM', 'POS'])
    return genotypes, variants


def _has_chromosome(variants: allel.VariantChunkedTable, chrom: str):
    """ Returns whether a chromosome is listed in a VariantChunkedTable.

        Parameters:
            variants (VariantChunkedTable): Variants presents in the HDF5 file.
            chrom (str): Chromosome to look for in the variants table.
        Returns:
            bool: True if chromosome 'chrom' is present in the table 'variants'
    """
    np_array = strainer.variants_filter_by_chrom(variants, chrom)
    return np.count_nonzero(np_array)


def get_sample_index(sample_list: list, sample: str) -> int:
    """ Returns the index of a sample in a sample list.

        Parameters:
            sample_list (list: str): List of samples (strings) present in the VCF file.
            sample (str): Sample name in the 'sample_list' to be selected as parental genotype.
        Returns:
            integer: The index of input sample in the sample list.
    """
    try:
        parental_sample_index = sample_list.index(sample)
    except ValueError:
        msg = "Sample '{sample}' is not found in VCF sample list '{sample_list}'".format(
            sample=sample,
            sample_list=sample_list
        )
        logger.error(msg)
        raise ValueError(msg)
    return parental_sample_index


def _sort_genotypes(genotypes: allel.GenotypeChunkedArray,
                    parental_sample_index: int) -> allel.GenotypeChunkedArray:
    """ Sorts an allel.GenotypeChunkedArray object, placing the sample selected as parental
        'parental_sample_index' as the first genotype listed in 'genotypes'.

        Parameters:
            genotypes (allel.GenotypeChunkedArray): GenotypesChunkedArray object.
            parental_sample_index (int): Index in the 'sample_list' of the parental sample
        Returns:
            The allel.GenotypeChunkedArray object but sorted.
    """
    selected_progeny = np.repeat(True, genotypes.n_samples)
    selected_progeny[parental_sample_index] = False
    selected_parental = np.logical_not(selected_progeny)

    genotype_progeny = genotypes.subset(None, selected_progeny)[:]
    genotype_parental = genotypes.subset(None, selected_parental)[:]

    genotype = genotype_parental.concatenate(genotype_progeny, axis=1)
    return genotype


def process(vcf_file_path: str, chrom: str,
            parental_sample: str,
            zygosis: haplotyper.Zygosity) -> haplotyper.HaplotypeWrapper:
    """ Returns a 'haplotyper.HaplotypeWrapper' object.

        Parameters:
            vcf_file_path (str): Input path to the VCF file.
            chrom (str): What chromosome should be considered for the haplotype process.
            parental_sample (str): Sample name that is considered as the parental one.
            zygosis (haplotyper.Zygosity): Whether the VCF has only homozygous variants or has heterozygous
        Returns:
           haplotyper.HaplotypeWrapper: Wrapper containing genotypes and variants.
    """
    vcf_file_abspath = os.path.abspath(vcf_file_path)
    vcf_path = os.path.dirname(vcf_file_abspath)
    hdf5_filename = os.path.splitext(vcf_file_abspath)[0] + constants.HDF5_EXT
    hdf5_file_path = os.path.join(vcf_path, hdf5_filename)

    reader.vcf_to_hdf5(vcf_file_path, hdf5_file_path)
    logger.debug("Loading HDF5 file '{hdf5_path}'".format(hdf5_path=hdf5_file_path))
    callset = h5py.File(hdf5_file_path, mode='r')
    sample_list = list(callset["samples"])
    parental_sample_index = get_sample_index(sample_list, parental_sample)
    genotypes, variants = get_genotypes_n_variants(callset)
    if not _has_chromosome(variants, chrom):
        msg = "Chromosome '{chrom}' not found in the VCF '{vcf_file_path}'".format(
            chrom=chrom,
            vcf_file_path=vcf_file_path
        )
        logger.error(msg)
        raise ValueError(msg)

    genotypes = _sort_genotypes(genotypes, parental_sample_index)
    genotypes_uc, variants_uc = strainer.filters_for_haplotyping(genotypes, variants, chrom)
    genotypes_uc, variants_uc = strainer.filter_phasing(genotypes_uc, variants_uc)

    haplotype_wrapper = haplotyper.HaplotypeWrapper(
        vcf_path, genotypes_uc, variants_uc, chrom, sample_list, parental_sample
    )

    if zygosis == haplotyper.Zygosity.HOM:
        haplotype_wrapper.calc_homozygous_haplotypes()
    elif zygosis == haplotyper.Zygosity.HET:
        haplotype_wrapper.calc_heterozygous_haplotypes()
    else:
        msg = "Unknown zygosis '{zygosis}' parameter".format(zygosis=zygosis)
        logger.error(msg)
        raise ValueError(msg)
    return haplotype_wrapper
