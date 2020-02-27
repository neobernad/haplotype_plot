# -*- coding: utf-8 -*-
import haplotype_plot.genotyper as genotyper
import haplotype_plot.conversion as converter
import haplotype_plot
import argparse


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--vcf', required=True, action='store',
                        help="`path to the input VCF")
    parser.add_argument('-c', '--chr', required=True, action='store',
                        help="chromosome to plot from the VCF")
    parser.add_argument('-p', '--parental', required=True, action='store',
                        help="sample name from the VCF used as parental haplotype")
    parser.add_argument('--version', action='version',
                        version='version: {version}'.format(version=haplotype_plot.version))
    return parser.parse_args()


def main():
    args = parse_args()
    print(args)


if __name__ == '__main__':
    main()
