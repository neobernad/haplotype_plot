# -*- coding: utf-8 -*-
import unittest
import os.path
import haplotype_plot.conversion as converter


class TestConversion(unittest.TestCase):

    def test_vcf_to_hdf5(self):
        vcf_path = "data/chr01.vcf"
        hdf5_path = "data/chr01.hdf5"
        converter.vcf_to_hdf5(vcf_path, hdf5_path)
        assert os.path.exists(hdf5_path)

        vcf_path = "data/chr01_02_03.vcf"
        hdf5_path = "data/chr01_02_03.hdf5"
        converter.vcf_to_hdf5(vcf_path, hdf5_path)
        assert os.path.exists(hdf5_path)

    def test_get_samples(self):
        vcf_path = "data/chr01.vcf"
        sample_list = converter.get_samples(vcf_path)
        print("Samples in '{vcf_path}' are: {sample_list}".format(vcf_path=vcf_path,
                                                                  sample_list=sample_list))
        assert sample_list


if __name__ == '__main__':
    unittest.main()
