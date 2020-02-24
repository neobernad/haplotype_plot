import unittest
import pprint
import haplotyper.genotyper as genotyper


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


if __name__ == '__main__':
    unittest.main()
