# -*- coding: utf-8 -*-
import collections
import copy
import logging
import os
import vcf
import allel
import haplotype_plot.filter as strainer
from haplotype_plot.haplotyper import HaplotypeWrapper

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


def save_phased_variants(haplotype_wrapper: HaplotypeWrapper,
                         input_vcf_file: str, output_vcf_file: str = None):
    """ Saves a file with the phased variants from the object HaplotypeWrapper.

        Parameters:
            haplotype_wrapper (HaplotypeWrapper): HaplotypeWrapper object.
            input_vcf_file (str): Path of the input VCF file
            output_vcf_file (str): Output VCF phased file. If 'None' specified,
            the output file is calculated based on the input_vcf_file parameter.
            For instance:
                input_vcf_file = "/path/to/vcf/file.vcf"
                Then, if 'output_vcf_file' is None, the following path is generated:
                output_vcf_file = "/path/to/vcf/file.phased.vcf"

    """
    # https://github.com/jamescasbon/PyVCF/issues/82
    if output_vcf_file is None:
        output_vcf_file = os.path.splitext(os.path.abspath(input_vcf_file))[0] + ".phased.vcf"
    logger.debug("Saving phased variants...")
    fp_vcf_file = open(input_vcf_file, 'r')
    fp_output_vcf_file = open(output_vcf_file, 'w')
    vcf_reader = vcf.Reader(fp_vcf_file)
    # Place "GT" as first element on the format keys
    f_keys = list(vcf_reader.formats.keys())  # Ordered dict -> list
    f_keys.remove("GT")
    f_keys.insert(0, "GT")
    vcf_writer = vcf.Writer(fp_output_vcf_file, vcf_reader)
    data_variants = haplotype_wrapper.variants.copy()
    data_genotypes = haplotype_wrapper.genotypes.copy()
    len_variants = len(data_variants)
    if len_variants != len(data_genotypes.is_phased):
        raise RuntimeError("Number of phased calls mismatch with the number of variants")

    try:
        for record in vcf_reader:
            selected_element = strainer.variants_filter_by_chrom_n_pos(haplotype_wrapper.variants,
                                                                       record.CHROM, record.POS)
            data_variants_compressed = data_variants.compress(selected_element)
            if data_variants_compressed:
                new_record = copy.deepcopy(record)
                data_genotypes_compressed = data_genotypes.compress(selected_element)
                if len(data_genotypes_compressed) > 1:  # We expect to retrieve only 1 row from data_genotypes.compress
                    raise RuntimeError("Multiple records returned by filter. We should retrieve only one here.")
                data_genotypes_vector = data_genotypes_compressed[0]
                for index, sample in enumerate(haplotype_wrapper.sample_list):
                    new_record.samples[index].data = collections.namedtuple('CallData', f_keys)
                    # f_vals = [record.samples[index].data[vx] for vx in range(len(f_keys))]
                    # handy_dict = dict(zip(f_keys, f_vals))
                    handy_dict = dict()
                    phased_gt = '|'.join(map(str, data_genotypes_vector[index]))
                    # Sometimes CallData fields (GT, DP) might differ with FORMAT (GT, AD)
                    sample_calldata_fields = record.samples[index].data._fields
                    for f in f_keys:
                        if f in sample_calldata_fields:
                            if f == 'GT':
                                handy_dict[f] = phased_gt
                                continue
                            handy_dict[f] = record.samples[index][f]
                        else:
                            handy_dict[f] = ""
                    new_vals = [handy_dict[x] for x in f_keys]  # values of the keys
                    new_record.samples[index].data = new_record.samples[index].data._make(new_vals)
                    new_record.samples[index].gt_nums = phased_gt
                    for allele_index in range(len(new_record.samples[index].gt_alleles)):
                        new_record.samples[index].gt_alleles[allele_index] = data_genotypes_vector[index][allele_index]
                vcf_writer.write_record(new_record)
    finally:
        fp_vcf_file.close()
        fp_output_vcf_file.close()

    logger.debug("Variants saved in '{output_vcf_file}'".format(output_vcf_file=output_vcf_file))
