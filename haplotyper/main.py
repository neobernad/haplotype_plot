# -*- coding: utf-8 -*-
import haplotyper.genotyper as genotyper


def main():
    vcf_file_path = "tests/data/chr01.vcf"
    chrom = "chr01"
    parental_sample = "SAMPLE1"
    genotyper.process_haplotypes(vcf_file_path, chrom, parental_sample)


if __name__ == '__main__':
    main()
