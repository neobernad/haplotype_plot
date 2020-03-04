# -*- coding: utf-8 -*-
import logging
import os
from haplotype_plot.haplotyper import HaplotypeWrapper

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


def hdf5_to_vcf(input_vcf_file_path: str, haplotype_wrapper: HaplotypeWrapper):
    output_vcf_path = os.path.dirname(os.path.abspath(input_vcf_file_path))
    print(output_vcf_path)
    logger.debug("Ok")
