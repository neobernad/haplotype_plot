# -*- coding: utf-8 -*-
import unittest
import os.path
import os
import haplotype_plot.reader as reader
import haplotype_plot.writer as writer

dir_path = os.path.dirname(os.path.realpath(__file__))


class TestConversion(unittest.TestCase):
    vcf_one_chr_path = os.path.join(dir_path, "data/chr01.vcf")
    hdf5_one_chr_path = os.path.join(dir_path, "data/chr01.h5")
    vcf_many_chr_path = os.path.join(dir_path, "data/chr01_02_03.vcf")
    hdf5_many_chr_path = os.path.join(dir_path, "data/chr01_02_03.h5")
    vcf_test_output = os.path.join(dir_path, "data/test.vcf")

    def test_vcf_to_hdf5(self):
        reader.vcf_to_hdf5(self.vcf_one_chr_path, self.hdf5_one_chr_path)
        assert os.path.exists(self.hdf5_one_chr_path)

        reader.vcf_to_hdf5(self.vcf_many_chr_path, self.hdf5_many_chr_path)
        assert os.path.exists(self.hdf5_many_chr_path)

    def test_get_samples(self):
        sample_list = reader.get_samples(self.vcf_one_chr_path)
        print("Samples in '{vcf_path}' are: {sample_list}".format(vcf_path=self.vcf_one_chr_path,
                                                                  sample_list=sample_list))
        assert sample_list


if __name__ == '__main__':
    unittest.main()
