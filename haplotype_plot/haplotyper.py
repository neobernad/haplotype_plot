# -*- coding: utf-8 -*-
import logging
import allel
import numpy as np
import os
from enum import Enum

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class Zygosity(Enum):
    UNDEFINED = 1
    HOM = 2  # HOMOZYGOUS
    HET = 3  # HETEROZYGOUS


class HaplotypeWrapper(object):
    DEFAULT_PLOT_NAME = "haplotypes.png"

    """
    Class to wrap the following information:
        * __input_path (str): Absolute path where the input data is.
        * __genotypes (allel.GenotypeArray): Genotypes in the current analysis.
        * __variants (allel.VariantTable): Variants in the current analysis.
        * __chrom (str): Chromosome where variants/genotypes are
        * __sample_list (List[str]): List of samples (from the VCF).
        * __parent_sample (str): Name of the sample (as shown in the VCF) used as parent.
        * __zygosity (Zygosity): Whether the haplotypes are homozygous or heterozygous.
        * __parent_haplotypes (allel.HaplotypeArray): Parent haplotypes.
        * __parent_n_progeny_haplotypes (allel.HaplotypeArray): Parent haplotypes concatenated with progeny haplotypes.
    """

    def __init__(self, input_path: str, genotypes: allel.GenotypeArray,
                 variants: allel.VariantTable, chrom: str, sample_list: list,
                 parent_sample: str):
        self.__input_path: str = input_path
        self.__genotypes: allel.GenotypeArray = genotypes
        self.__variants: allel.VariantTable = variants
        self.__chrom: str = chrom
        self.__sample_list: list = sample_list
        self.__parent_sample: str = parent_sample
        self.__zygosity: Zygosity = Zygosity.UNDEFINED
        self.__parent_haplotypes: allel.HaplotypeArray = allel.HaplotypeArray(
            # Empty allel.HaplotypeArray
            np.empty((self.__genotypes.n_variants, self.__genotypes.n_samples), dtype='i1')
        )
        self.__parent_n_progeny_haplotypes: allel.HaplotypeArray = allel.HaplotypeArray(
            np.empty((self.__genotypes.n_variants, self.__genotypes.n_samples), dtype='i1')
        )

    def __str__(self):
        attrs = vars(self)
        return ', '.join("%s: %s" % item for item in attrs.items())

    @property
    def input_path(self) -> str:
        return self.__input_path

    @property
    def genotypes(self) -> allel.GenotypeArray:
        return self.__genotypes

    @property
    def variants(self) -> allel.VariantTable:
        return self.__variants

    @property
    def chrom(self) -> str:
        return self.__chrom

    @property
    def sample_list(self) -> list:
        return self.__sample_list

    @property
    def parent_sample(self) -> str:
        return self.__parent_sample

    @property
    def parent_haplotypes(self) -> allel.HaplotypeArray:
        return self.__parent_haplotypes

    @property
    def parent_n_progeny_haplotypes(self) -> allel.HaplotypeArray:
        return self.__parent_n_progeny_haplotypes

    def generate_plot_path(self) -> str:
        return os.path.join(self.input_path, self.DEFAULT_PLOT_NAME)

    def is_homozygous(self):
        return self.__zygosity == Zygosity.HOM

    def is_heterozygous(self):
        return self.__zygosity == Zygosity.HET

    def __set_homozygous(self):
        self.__zygosity = Zygosity.HOM

    def __set_heterozygous(self):
        self.__zygosity = Zygosity.HET

    def calc_homozygous_haplotypes(self):
        """ Calculates a the parent and progeny haplotypes from a given 'allel.GenotypeArray'.
            It considers ONLY ONE of the alleles due to it supposes the genotypes are homozygous.

            Parent haplotypes are found in self.__parent_haplotypes.
            Parent plus progeny haplotypes are found in self.__parent_n_progeny_haplotypes.
        """
        # Parent genotype
        genotypes_parent = self.genotypes[:, 0]
        # Convert to haplotype array
        haplotypes_parent = genotypes_parent.to_haplotypes()
        # Pull out the "left" allele (haplotypes) from the other samples in the VCF, treated as progeny
        # Here we assume the genotypes are homozygous
        haplotypes_rest_varieties = allel.HaplotypeArray(self.genotypes[:, 1:, 0])
        # Stack parent's haplotypes alongside haplotypes it transmitted to its progeny
        haplotypes_parent_n_progeny = haplotypes_parent.concatenate(haplotypes_rest_varieties, axis=1)
        self.__parent_haplotypes = haplotypes_parent
        self.__parent_n_progeny_haplotypes = haplotypes_parent_n_progeny
        self.__set_homozygous()

    def calc_heterozygous_haplotypes(self) -> (allel.HaplotypeArray, allel.HaplotypeArray):
        """ Returns a the parent and progeny haplotypes from a given 'allel.GenotypeArray'.
            It considers BOTH alleles due to it supposes the genotypes are heterozygous.

            Parent haplotypes are found in self.__parent_haplotypes.
            Parent plus progeny haplotypes are found in self.__parent_n_progeny_haplotypes.
        """
        # Parent genotype is indexed in position 0
        genotypes_parent = self.genotypes[:, 0]
        # Convert to haplotype array
        haplotypes_parent = genotypes_parent.to_haplotypes()
        # Pull out the both allele (haplotypes) from the other samples in the VCF, treated as progeny
        # Skip genotype 0 (parent)
        left_alleles = allel.HaplotypeArray(self.genotypes[:, :, 0])
        right_alleles = allel.HaplotypeArray(self.genotypes[:, :, 1])

        # Initially, copy parent alleles
        haplotypes_parent_n_progeny = allel.HaplotypeArray(haplotypes_parent, copy=True)
        wanted_variants = np.repeat(True, self.genotypes.n_variants)  # Get all variants
        wanted_sample = np.repeat(False, self.genotypes.n_samples)  # Set to False initially
        # Start at 1, we already copied the parents haplotypes in 'haplotypes_parent_n_progeny'
        for sample_index in range(1, self.genotypes.n_samples):  # Skip haplotype 0 (parent), it is already inserted
            wanted_sample[sample_index] = True
            subset_left_alleles = left_alleles.subset(wanted_variants, wanted_sample)[:]
            subset_right_alleles = right_alleles.subset(wanted_variants, wanted_sample)[:]
            wanted_sample[sample_index] = False
            haplotypes_parent_n_progeny = haplotypes_parent_n_progeny.concatenate(subset_left_alleles, axis=1)
            haplotypes_parent_n_progeny = haplotypes_parent_n_progeny.concatenate(subset_right_alleles, axis=1)

        self.__parent_haplotypes = haplotypes_parent
        self.__parent_n_progeny_haplotypes = haplotypes_parent_n_progeny
        self.__set_heterozygous()
