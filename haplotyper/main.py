# -*- coding: utf-8 -*-
import haplotyper.genotyper as genotyper
import haplotyper.conversion as converter


def main():
    vcf_file_path = "tests/data/chr01.vcf"
    chrom = "chr01"
    parental_sample = "SAMPLE1"
    sample_list = converter.get_samples(vcf_file_path)
    genotypes_uc, variants_uc = genotyper.process(vcf_file_path, chrom, sample_list, parental_sample)
    parent_haplotypes, all_haplotypes = genotyper.get_haplotypes(genotypes_uc, variants_uc)


if __name__ == '__main__':
    main()
