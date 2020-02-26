# -*- coding: utf-8 -*-
import logging
import allel
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import haplotyper.genotyper as genotyper
import haplotyper.conversion as converter
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


def get_painting(parent_haplotypes: allel.HaplotypeArray,
                 all_haplotypes: allel.HaplotypeArray) -> np.ndarray:
    return allel.paint_transmission(parent_haplotypes, all_haplotypes)


def get_ytickslabels(sample_list: list, parental_sample: str) -> list:
    sample_np_array = np.array(sample_list)
    num_genotypes = len(sample_list)
    parental_sample_index = genotyper.get_sample_index(sample_list, parental_sample)
    selected_progeny = np.repeat(True, num_genotypes)
    selected_progeny[parental_sample_index] = False
    yticklabels: list = [parental_sample + "_1", parental_sample + "_2"]
    yticklabels.extend(sample_np_array[selected_progeny])
    return yticklabels


def get_xtickslabels(variants: allel.VariantChunkedTable) -> list:
    return list(variants["POS"])


def plot_haplotypes(variants_uc: allel.VariantChunkedTable, parent_haplotypes: allel.HaplotypeArray,
                    all_haplotypes: allel.HaplotypeArray, sample_list: list, parental_sample: str,
                    plot_config: PlotConfig):
    painting = get_painting(parent_haplotypes, all_haplotypes)
    plot_transmission(painting, plot_config)


def _crop_painting_dimension(painting: np.ndarray, plot_config: PlotConfig):

    # TODO: This plots works for homozygous haplotypes, create another one for heterozygous
    # SAMPLE1_1, SAMPLE1_2, SAMPLE2_1, SAMPLE2_2, etc..

    start: int = plot_config.start
    end: int = plot_config.end

    if end == 0:  # Show full plot
        end = len(painting)

    plot_config.xtickslabels = plot_config.xtickslabels[start:end]
    return painting[start:end]


def plot_transmission(painting: np.ndarray, plot_config: PlotConfig):

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
    painting = _crop_painting_dimension(painting, plot_config)
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
