# -*- coding: utf-8 -*-
import haplotype_plot.genotyper as genotyper
import haplotype_plot.conversion as converter
import haplotype_plot.haplotyper as haplotyper
import haplotype_plot.plot as hplot
import haplotype_plot
import argparse


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Creates a haplotype plot from a VCF file")
    parser.add_argument('-v', '--vcf', required=True, action='store',
                        help="path to the input VCF file")
    parser.add_argument('-c', '--chr', required=True, action='store',
                        help="chromosome to plot from the VCF")
    parser.add_argument('-p', '--parental', required=True, action='store',
                        help="sample name from the VCF used as parental haplotype")
    parser.add_argument('--version', action='version',
                        version='version: {version}'.format(version=haplotype_plot.version))
    parser.add_argument('-z', '--zygosis', help="zygosis of the input VCF file",
                        choices=[haplotyper.Zygosity.HOM.name, haplotyper.Zygosity.HET.name])
    parser.add_argument("--conf",
                        metavar="KEY=VALUE",
                        nargs='+',
                        help="set a number of key-value pairs to modify the default plot configuration "
                             "(do not put spaces before or after the = sign). "
                             "Possible values are: 'start=integer', 'end=integer'. "
                             "If a value contains spaces, you should define "
                             "it with double quotes: "
                             'foo="this is a sentence".')
    return parser.parse_args()


def main():
    args = parse_args()
    sample_list = converter.get_samples(args.vcf)
    haplotype_wrapper = genotyper.process(args.vcf, args.chr,
                                          sample_list, args.parental,
                                          haplotyper.Zygosity[args.zygosis])
    plotter = hplot.Plotter(haplotype_wrapper)
    plotter.plot_haplotypes(override_conf=args.conf)


if __name__ == '__main__':
    main()