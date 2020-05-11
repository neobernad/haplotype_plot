import unittest
import logging
import os
import haplotype_plot.genotyper as genotyper
import haplotype_plot.reader as reader
import haplotype_plot.haplotyper as haplotyper
import haplotype_plot.plot as hplot

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

dir_path = os.path.dirname(os.path.realpath(__file__))


class TestPlotting(unittest.TestCase):
    vcf_path = os.path.join(dir_path, "data/chr01.vcf")
    chrom = "chr01"
    parental_sample = "SAMPLE4"
    sample_list = None

    def setUp(self) -> None:
        self.sample_list = reader.get_samples(self.vcf_path)

    def test_plot_config(self):
        plot_config = hplot.PlotConfig()
        logger.debug(plot_config)

    def test_generate_heterozygous_yticks(self):
        heterozygous = haplotyper.Zygosity.HET
        haplotype_wrapper = genotyper.process(self.vcf_path, self.chrom,
                                              self.parental_sample, heterozygous)
        plotter = hplot.Plotter(haplotype_wrapper)
        labels = plotter.get_ytickslabels()
        logger.debug("Parent: {parent}".format(parent=self.parental_sample))
        logger.debug(labels)

    def test_generate_homozygous_yticks(self):
        homozygous = haplotyper.Zygosity.HOM
        haplotype_wrapper = genotyper.process(self.vcf_path, self.chrom,
                                              self.parental_sample, homozygous)
        plotter = hplot.Plotter(haplotype_wrapper)
        labels = plotter.get_ytickslabels()
        logger.debug("Parent: {parent}".format(parent=self.parental_sample))
        logger.debug(labels)

    def test_plot_homozygous_haplotypes(self):
        homozygous = haplotyper.Zygosity.HOM
        haplotype_wrapper = genotyper.process(self.vcf_path, self.chrom,
                                              self.parental_sample, homozygous)
        plotter = hplot.Plotter(haplotype_wrapper)
        ytickslabels = plotter.get_ytickslabels()

        custom_config = hplot.PlotConfig(
            title="Parental '{parent}' in '{chrom}'".format(parent=self.parental_sample, chrom=self.chrom),
            xtickslabels=plotter.get_xtickslabels(),
            ytickslabels=ytickslabels,
            start=0,
            end=1000,
            size_x=10,
            size_y=len(ytickslabels) * .2,
            show=True
        )

        plotter.plot_haplotypes(custom_config)

    def test_plot_heterozygous_haplotypes(self):
        heterozygous = haplotyper.Zygosity.HET
        haplotype_wrapper = genotyper.process(self.vcf_path, self.chrom,
                                              self.parental_sample, heterozygous)
        plotter = hplot.Plotter(haplotype_wrapper)
        user_conf = list(["show=False", "xtickslabels=False", "size_y=5"])
        plotter.plot_haplotypes(override_conf=user_conf)


if __name__ == '__main__':
    unittest.main()
