import unittest
from typing import List
import logging
import haplotyper.genotyper as genotyper
import haplotyper.conversion as converter
import haplotyper.plot as hplot
import numpy as np

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class TestPlotting(unittest.TestCase):
    vcf_path = "data/chr01.vcf"
    chrom = "chr01"
    parental_sample = "SAMPLE1"

    def test_plot_config(self):
        plot_config = hplot.PlotConfig()
        print(plot_config)

    def test_generate_yticks(self):
        sample_list = converter.get_samples(self.vcf_path)
        genotypes_uc, variants_uc = genotyper.process(self.vcf_path, self.chrom, sample_list, self.parental_sample)
        labels = hplot.get_ytickslabels(sample_list, self.parental_sample)
        logger.info(labels)

    def test_plot_haplotypes(self):
        plot_title = "{parent} haplotypes in chr {chrom}".format(
            parent=self.parental_sample,
            chrom = self.chrom
        )
        sample_list = converter.get_samples(self.vcf_path)
        genotypes_uc, variants_uc = genotyper.process(self.vcf_path, self.chrom, sample_list, self.parental_sample)
        parent_haplotypes, all_haplotypes = genotyper.get_haplotypes(genotypes_uc)

        ytickslabels = hplot.get_ytickslabels(sample_list, self.parental_sample)

        plot_config = hplot.PlotConfig(
            title="Parent {parent} in chr {chrom}".format(parent=self.parental_sample, chrom=self.chrom),
            xtickslabels=hplot.get_xtickslabels(variants_uc),
            ytickslabels=ytickslabels,
            start=50,
            end=60,
            size_x=10,
            size_y=len(ytickslabels) * .2
        )

        hplot.plot_haplotypes(variants_uc, parent_haplotypes,
                              all_haplotypes, sample_list,
                              self.parental_sample, plot_config)


if __name__ == '__main__':
    unittest.main()
