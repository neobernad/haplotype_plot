import unittest
import pprint
import haplotype_plot.genotyper as genotyper
import haplotype_plot.conversion as converter


class TestGenotyper(unittest.TestCase):

    def test_genotype_n_variants(self):
        hdf5_path = "data/chr01.hdf5"
        genotypes, variants = genotyper.get_genotypes_n_variants(hdf5_path)
        print(genotypes)
        print(variants)
        assert genotypes
        assert variants

    def test_callset_keys(self):
        hdf5_path = "data/chr01.hdf5"
        genotyper.debug_hdf5(hdf5_path)

    def test_genotype_filtering(self):
        vcf_path = "data/chr01.vcf"
        chrom = "chr01"
        parental_sample = "SAMPLE1"
        sample_list = converter.get_samples(vcf_path)
        haplotype_wrapper = genotyper.process_homozygous(vcf_path, chrom, sample_list, parental_sample)
        genotypes_uc_len = len(haplotype_wrapper.genotypes)
        variants_uc_len = len(haplotype_wrapper.variants)
        print("Length genotypes: " + str(genotypes_uc_len))
        print("Length variants: " + str(variants_uc_len))
        assert genotypes_uc_len == variants_uc_len


if __name__ == '__main__':
    unittest.main()
