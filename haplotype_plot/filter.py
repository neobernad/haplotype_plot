# -*- coding: utf-8 -*-
import logging
import numpy as np
import allel

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


def variants_filter(variants: allel.VariantChunkedTable, filter_expression: str) -> np.ndarray:
    variants_np_array = variants[:]  # Conversion to allel.model.ndarray.VariantTable
    variant_selection = variants_np_array.eval(filter_expression)[:]
    return variant_selection


def variants_filter(variants: allel.VariantTable, filter_expression: str) -> np.ndarray:
    variant_selection = variants.eval(filter_expression)
    return variant_selection


def variants_filter_by_chrom(variants: allel.VariantChunkedTable, chrom: str) -> np.ndarray:
    filter_expression = '(CHROM == \'' + chrom + '\')'
    return variants_filter(variants, filter_expression)


def variants_filter_by_chrom(variants: allel.VariantTable, chrom: str) -> np.ndarray:
    filter_expression = '(CHROM == \'' + chrom + '\')'
    return variants_filter(variants, filter_expression)


def variants_filter_by_chrom_n_pos(variants: allel.VariantChunkedTable, chrom: str, pos: int) -> np.ndarray:
    filter_expression = '(CHROM == \'' + chrom + '\') & (POS == ' + str(pos) + ')'
    return variants_filter(variants, filter_expression)


def variants_filter_by_chrom_n_pos(variants: allel.VariantTable, chrom: str, pos: int) -> np.ndarray:
    filter_expression = '(CHROM == \'' + chrom + '\') & (POS == ' + str(pos) + ')'
    return variants_filter(variants, filter_expression)


def filters_for_haplotyping(genotypes: allel.GenotypeChunkedArray,
                            variants: allel.VariantChunkedTable,
                            chrom: str) -> (allel.GenotypeArray,
                                            allel.VariantTable):
    """ Performs a series of filters to prepare the 'genotypes' and 'variants' object
        for haplotyping.

        Parameters:
            genotypes (allel.GenotypeChunkedArray): GenotypesChunkedArray object.
            variants (allel.VariantChunkedTable): VariantChunkedTable object.
            chrom (str): What chromosome should be considered for the haplotype process.
        Returns:
            Tuple (allel.GenotypeArray, allel.VariantTable):
                - allel.GenotypeArray: GenotypeArray object
                - allel.VariantTable: VariantTable object
    """
    # Filter by chrom
    np_array_variants_in_chr = variants_filter_by_chrom(variants, chrom)
    logger.debug("There are {count_variants_in_chr} variants in chromosome {chrom}".format(
        count_variants_in_chr=np.count_nonzero(np_array_variants_in_chr),
        chrom=chrom
    ))

    # Filter by segregating SNPs
    allele_count = genotypes.count_alleles()
    np_array_log_sec = allele_count.is_segregating()
    logger.debug("There are {count_log_sec} segregating SNPs".format(
        count_log_sec=np.count_nonzero(np_array_log_sec)
    ))
    np_array_variants_to_keep = np_array_variants_in_chr & np_array_log_sec
    logger.debug("Number of variants to keep {count_variants_to_keep}".format(
        count_variants_to_keep=np.count_nonzero(np_array_variants_to_keep)
    ))
    # Subsets: perform the subset and load the results into memory uncompressed
    genotypes_uc = genotypes.subset(np_array_variants_to_keep, range(0, genotypes.n_samples))[:]
    variants_np_array = variants[:]
    variants_uc = variants_np_array.compress(np_array_variants_to_keep)
    return genotypes_uc, variants_uc


def filter_phasing(genotypes_un: allel.GenotypeArray,
                   variants_un: allel.VariantTable,
                   window_size: int = 100) -> (allel.GenotypeArray, allel.VariantTable):
    """ Filters 'genotypes_un' and 'variants_un' whether there are phasing genotypes.

        Parameters:
            genotypes_un (allel.GenotypeArray): GenotypeArray object.
            variants_un (allel.VariantTable): VariantTable object.
            window_size (int): Number of previous heterozygous sites to include when phasing each parent.
        Returns:
            Tuple (allel.GenotypeArray, allel.VariantTable):
                - allel.GenotypeArray: GenotypeArray object
                - allel.VariantTable: VariantTable object
    """
    # Phasing by transmission
    genotypes_phased = allel.phase_by_transmission(genotypes_un, window_size=window_size)
    # locate variants where all genotype calls were phased
    np_array_loc_phased_all = np.all(genotypes_phased.is_phased, axis=1)
    count_loc_phased_all = np.count_nonzero(np_array_loc_phased_all)
    if count_loc_phased_all == 0:
        # At least 1 variant could not be phased
        logger.warning("Could not find enough phased calls. Using all calls")
        return genotypes_phased[:], variants_un

    logger.debug("Found {num_phased} phased variant calls".format(num_phased=count_loc_phased_all))
    genotypes_phased_all = genotypes_phased[np_array_loc_phased_all]
    variants_un_phased_all = variants_un.compress(np_array_loc_phased_all)
    return genotypes_phased_all, variants_un_phased_all
