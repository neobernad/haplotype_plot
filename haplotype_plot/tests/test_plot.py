import unittest
from typing import List
import logging
import haplotype_plot.genotyper as genotyper
import haplotype_plot.conversion as converter
import haplotype_plot.plot as hplot
import numpy as np

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class TestPlotting(unittest.TestCase):
    vcf_path = "data/chr01.vcf"
    chrom = "chr01"
    parental_sample = "SAMPLE1"
    sample_list = None

    def setUp(self) -> None:
        self.sample_list = converter.get_samples(self.vcf_path)

    def test_plot_config(self):
        plot_config = hplot.PlotConfig()
        print(plot_config)

    def test_generate_homozygous_yticks(self):
        haplotype_wrapper = genotyper.process_homozygous(self.vcf_path, self.chrom,
                                                         self.sample_list, self.parental_sample)
        plotter = hplot.Plotter(haplotype_wrapper)
        labels = plotter.get_ytickslabels()
        logger.info(labels)

    def test_generate_heterozygous_yticks(self):
        haplotype_wrapper = genotyper.process_homozygous(self.vcf_path, self.chrom,
                                                         self.sample_list, self.parental_sample)
        plotter = hplot.Plotter(haplotype_wrapper)
        labels = plotter.get_ytickslabels()
        logger.info(labels)

    def test_plot_homozygous_haplotypes(self):
        haplotype_wrapper = genotyper.process_homozygous(self.vcf_path, self.chrom,
                                                         self.sample_list, self.parental_sample)
        plotter = hplot.Plotter(haplotype_wrapper)
        ytickslabels = plotter.get_ytickslabels()

        custom_config = hplot.PlotConfig(
            title="Parental '{parent}' in '{chrom}'".format(parent=self.parental_sample, chrom=self.chrom),
            xtickslabels=plotter.get_xtickslabels(),
            ytickslabels=ytickslabels,
            start=0,
            end=1000,
            size_x=10,
            size_y=len(ytickslabels) * .2
        )

        plotter.plot_haplotypes(custom_config)

    def test_plot_heterozygous_haplotypes(self):
        haplotype_wrapper = genotyper.process_heterozygous(self.vcf_path, self.chrom,
                                                           self.sample_list, self.parental_sample)
        plotter = hplot.Plotter(haplotype_wrapper)
        plotter.plot_haplotypes()


if __name__ == '__main__':
    unittest.main()
