import collections
import unittest
import logging
import os
import vcf  # from PyVCF
import numpy as np
import copy
import haplotype_plot.genotyper as genotyper
import haplotype_plot.reader as reader
import haplotype_plot.filter as strainer
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
        self.sample_list = reader.get_samples(self.vcf_one_chr_path)

    def test_genotype_n_variants(self):
        haplotype_wrapper = genotyper.process(self.vcf_one_chr_path, self.chrom, self.parental_sample,
                                              Zygosity.HOM)
        assert haplotype_wrapper

    def test_genotype_filtering(self):
        haplotype_wrapper = genotyper.process(self.vcf_one_chr_path, self.chrom, self.parental_sample,
                                              Zygosity.HOM)
        genotypes_uc_len = len(haplotype_wrapper.genotypes)
        variants_uc_len = len(haplotype_wrapper.variants)
        logger.debug("Length genotypes: {length}".format(length=genotypes_uc_len))
        logger.debug("Length variants: {length}".format(length=variants_uc_len))
        assert genotypes_uc_len == variants_uc_len

    # Find current chrom and pos in the VCF
    # vcf_reader.fetch is faster, but requires pysam and it will not work on windows
    # @unittest.skip("Not completed")
    def test_save(self):
        # https://github.com/jamescasbon/PyVCF/issues/82
        haplotype_wrapper = genotyper.process(self.vcf_one_chr_path, self.chrom, self.parental_sample,
                                              Zygosity.HOM)
        assert haplotype_wrapper
        output_vcf_file = os.path.splitext(os.path.abspath(self.vcf_one_chr_path))[0] + ".phased.vcf"
        print(output_vcf_file)
        fp_vcf_file = open(self.vcf_one_chr_path, 'r')
        fp_output_vcf_file = open(output_vcf_file, 'w')
        vcf_reader = vcf.Reader(fp_vcf_file)
        f_keys = list(vcf_reader.formats.keys())
        f_keys.remove("GT")
        f_keys.insert(0, "GT")
        vcf_writer = vcf.Writer(fp_output_vcf_file, vcf_reader)
        # print(haplotype_wrapper.variants[0:])
        data_variants = haplotype_wrapper.variants.copy()
        data_genotypes = haplotype_wrapper.genotypes.copy()
        len_variants = len(haplotype_wrapper.variants)
        if len_variants != len(data_genotypes.is_phased):
            raise RuntimeError("Number of phased calls mismatch with the number of variants")

        for record in vcf_reader:
            selected_element = strainer.variants_filter_by_chrom_n_pos(haplotype_wrapper.variants,
                                                                       record.CHROM, record.POS)
            data_variants_compressed = data_variants.compress(selected_element)
            if data_variants_compressed:
                new_record = copy.deepcopy(record)
                data_genotypes_compressed = data_genotypes.compress(selected_element)
                if len(data_genotypes_compressed) > 1:
                    raise RuntimeError("Multiple records returned by filter. We should retrieve only one here.")
                data_genotypes_vector = data_genotypes_compressed[0]

                for index, sample in enumerate(haplotype_wrapper.sample_list):
                    call = record.genotype(sample)
                    new_record.samples[index].data = collections.namedtuple('CallData', f_keys)
                    f_vals = [record.samples[index].data[vx] for vx in range(len(f_keys))]
                    handy_dict = dict(zip(f_keys, f_vals))
                    new_gt = '|'.join(map(str, data_genotypes_vector[index]))
                    for f in f_keys:
                        if f == 'GT':
                            handy_dict[f] = new_gt
                            continue
                        handy_dict[f] = record.samples[index][f]
                    new_vals = [handy_dict[x] for x in f_keys]
                    new_record.samples[index].data = new_record.samples[index].data._make(new_vals)
                    new_record.samples[index].gt_nums = new_gt
                    for allele_index in range(len(new_record.samples[index].gt_alleles)):
                        new_record.samples[index].gt_alleles[allele_index] = data_genotypes_vector[index][allele_index]
                vcf_writer.write_record(new_record)

        fp_vcf_file.close()
        fp_output_vcf_file.close()

    def test_callset_keys(self):
        genotyper.debug_hdf5(self.hdf5_one_chr_path)


if __name__ == '__main__':
    unittest.main()
