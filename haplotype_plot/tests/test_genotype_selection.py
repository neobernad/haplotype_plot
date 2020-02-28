import unittest
import logging
import os
import haplotype_plot.genotyper as genotyper
import haplotype_plot.conversion as converter
from haplotype_plot.haplotyper import Zygosity

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

dir_path = os.path.dirname(os.path.realpath(__file__))


class TestGenotyper(unittest.TestCase):
    vcf_one_chr_path = os.path.join(dir_path, "data/chr01.vcf")
    hdf5_one_chr_path = os.path.join(dir_path, "data/chr01.h5")
    chrom = "chr01"
    parental_sample = "SAMPLE1"
    sample_list = None

    def setUp(self) -> None:
        self.sample_list = converter.get_samples(self.vcf_one_chr_path)

    def test_genotype_n_variants(self):
        haplotype_wrapper = genotyper.process(self.vcf_one_chr_path, self.chrom,
                                              self.sample_list, self.parental_sample,
                                              Zygosity.HOM)
        assert haplotype_wrapper

    def test_genotype_filtering(self):
        haplotype_wrapper = genotyper.process(self.vcf_one_chr_path, self.chrom,
                                              self.sample_list, self.parental_sample,
                                              Zygosity.HOM)
        genotypes_uc_len = len(haplotype_wrapper.genotypes)
        variants_uc_len = len(haplotype_wrapper.variants)
        logger.debug("Length genotypes: {length}".format(length=genotypes_uc_len))
        logger.debug("Length variants: {length}".format(length=variants_uc_len))
        assert genotypes_uc_len == variants_uc_len


if __name__ == '__main__':
    unittest.main()
