# -*- coding: utf-8 -*-
import haplotype_plot.genotyper as genotyper
import haplotype_plot.writer as writer
import haplotype_plot.haplotyper as haplotyper
import haplotype_plot.plot as hplot
import haplotype_plot
import argparse


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Creates a haplotype plot from a VCF file")
    parser.add_argument('-v', '--vcf', required=True, action='store',
                        help="path to the input sorted VCF file")
    parser.add_argument('-o', '--output', required=False, action='store',
                        help="path to the output PNG haplotype plot file")
    parser.add_argument('-c', '--chr', required=True, action='store',
                        help="chromosome to plot from the VCF")
    parser.add_argument('-p', '--parental', required=True, action='store',
                        help="sample name from the VCF used as parental haplotype")
    parser.add_argument('--phase', required=False, action='store_true',
                        help="if specified, it will phase the genotypes. "
                             "An output VCF with the phased genotypes is created")
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
    haplotype_wrapper = genotyper.process(args.vcf, args.chr, args.parental,
                                          haplotyper.Zygosity[args.zygosis])

    plotter = hplot.Plotter(haplotype_wrapper)
    plotter.plot_haplotypes(output_file=args.output, override_conf=args.conf)
    if args.phase:
        writer.save_phased_variants(haplotype_wrapper, args.vcf)


if __name__ == '__main__':
    main()
