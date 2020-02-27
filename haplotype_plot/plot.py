# -*- coding: utf-8 -*-
import logging
import allel
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import haplotype_plot.genotyper as genotyper
import haplotype_plot.haplotyper as haplotyper
import haplotype_plot.conversion as converter
import seaborn as sns

sns.set_style('white')
sns.set_style('ticks')

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class PlotConfig:

    def __init__(self, title: str = None, xtickslabels=None, ytickslabels=None,
                 start: int = 0, end: int = 0, size_x: float = 10.0, size_y: float = 3.0):
        self.title = title
        self.xtickslabels = xtickslabels
        self.ytickslabels = ytickslabels
        self.start = start
        self.end = end
        self.size_x = size_x
        self.size_y = size_y
        self.__check_properties_integrity()

    def __check_properties_integrity(self):
        msg = None
        if self.start > self.end:
            msg = "PlotConfig start position is greater than end position {start} > {end}".format(
                start=self.start,
                end=self.end
            )
            logger.error(msg)
            raise ValueError(msg)

        if self.size_x <= 0:
            msg = "PlotConfig X axis size {size_x} is 0 or negative".format(
                size_x=self.size_x
            )
            logger.error(msg)
            raise ValueError(msg)

        if self.size_y <= 0:
            msg = "PlotConfig Y axis size {size_y} is 0 or negative".format(
                size_y=self.size_y
            )
            logger.error(msg)
            raise ValueError(msg)

    def __str__(self):
        attrs = vars(self)
        return ', '.join("%s: %s" % item for item in attrs.items())


class Plotter(object):

    def __init__(self, haplotype_wrapper: haplotyper.HaplotypeWrapper):
        self.__haplotype_wrapper: haplotyper.HaplotypeWrapper = haplotype_wrapper

    def __str__(self):
        attrs = vars(self)
        return ', '.join("%s: %s" % item for item in attrs.items())

    def __default_config(self) -> PlotConfig:
        ytickslabels = self.get_ytickslabels()
        default_plot_config = PlotConfig(
            title="Parent {parent} in chr {chrom}".format(
                parent=self.__haplotype_wrapper.parent_sample,
                chrom=self.__haplotype_wrapper.chrom
            ),
            xtickslabels=self.get_xtickslabels(),
            ytickslabels=ytickslabels,
            start=0,
            end=0,
            size_x=10,
            size_y=len(ytickslabels) * .2
        )
        return default_plot_config

    def get_painting(self) -> np.ndarray:
        return allel.paint_transmission(self.__haplotype_wrapper.parent_haplotypes,
                                        self.__haplotype_wrapper.parent_n_progeny_haplotypes)

    def __get_heterozygous_ytickslabels(self) -> list:
        ytickslabels = list()
        for sample in self.__haplotype_wrapper.sample_list:
            ytickslabels.append(sample + "_1")
            ytickslabels.append(sample + "_2")
        return ytickslabels

    def __get_homozygous_ytickslabels(self) -> list:
        sample_np_array = np.array(self.__haplotype_wrapper.sample_list)
        num_genotypes = len(self.__haplotype_wrapper.sample_list)
        parental_sample_index = genotyper.get_sample_index(self.__haplotype_wrapper.sample_list,
                                                           self.__haplotype_wrapper.parent_sample)
        selected_progeny = np.repeat(True, num_genotypes)
        selected_progeny[parental_sample_index] = False
        ytickslabels: list = [self.__haplotype_wrapper.parent_sample + "_1",
                             self.__haplotype_wrapper.parent_sample + "_2"]
        ytickslabels.extend(sample_np_array[selected_progeny])
        return ytickslabels

    def get_ytickslabels(self) -> list:
        if self.__haplotype_wrapper.is_homozygous():
            return self.__get_homozygous_ytickslabels()
        elif self.__haplotype_wrapper.is_heterozygous():
            return self.__get_heterozygous_ytickslabels()
        else:
            msg = "HaplotypeWrapper has UNDEFINED zygosity"
            logger.error(msg)
            raise ValueError(msg)

    def get_xtickslabels(self) -> list:
        return list(self.__haplotype_wrapper.variants["POS"])

    def plot_haplotypes(self, plot_config: PlotConfig = None):
        painting = self.get_painting()
        if not plot_config:
            plot_config = self.__default_config()
        self.plot_transmission(painting, plot_config)

    def _crop_painting_dimension(self, painting: np.ndarray, plot_config: PlotConfig) -> np.ndarray:
        start: int = plot_config.start
        end: int = plot_config.end

        if end == 0:  # Show full plot
            end = len(painting)

        plot_config.xtickslabels = plot_config.xtickslabels[start:end]
        return painting[start:end]

    def plot_transmission(self, painting: np.ndarray, plot_config: PlotConfig):
        fig, ax = plt.subplots(figsize=(plot_config.size_x, plot_config.size_y))
        palette = sns.color_palette("Spectral", 10)

        # map painting codes onto colours
        cmap = mpl.colors.ListedColormap([
            'grey',  # 0 = undetermined
            palette[4],  # 1 = allele inherited from first parental haplotype
            palette[3],  # 2 = allele inherited from second parental haplotype
            palette[8],  # 3 = reference allele, also carried by both parental haplotypes
            palette[7],  # 4 = non-reference allele, also carried by both parental haplotypes
            palette[0],  # 5 = non-parental allele (i.e., Mendelian error)
            palette[9],  # 6 = either or both parental alleles missing
            'white',  # 7 = missing allele
        ])
        # Delimit painting with user start/end range
        painting = self._crop_painting_dimension(painting, plot_config)
        # Plot painting
        ax.pcolormesh(painting.T, cmap=cmap, vmin=0, vmax=7)

        # Tidy up axes
        # ax.set_xticks(np.arange(painting.shape[0]))
        # ax.set_yticks(np.arange(painting.shape[1]) + .5)
        ax.set_yticks(np.arange(len(plot_config.ytickslabels)) + .5)
        ax.set_yticklabels(plot_config.ytickslabels)
        ax.set_xticks(np.arange(len(plot_config.xtickslabels)) + .5)
        ax.set_xticklabels(plot_config.xtickslabels, rotation=45)
        ax.set_ylabel('Progeny haplotypes')
        ax.set_xlabel('Variants')
        if plot_config.title:
            ax.set_title(plot_config.title)
        plt.tight_layout()
        plt.show()
