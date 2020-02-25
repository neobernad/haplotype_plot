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
                 start=0, end=0):
        self.title = title
        self.xtickslabels = xtickslabels
        self.ytickslabels = ytickslabels
        self.start = start
        self.end = end

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
    painting: np.ndarray = get_painting(parent_haplotypes, all_haplotypes)
    plot_transmission(painting, plot_config)


def plot_transmission(painting: np.ndarray, plot_config: PlotConfig):
    # set figure height depending on number of haplotypes
    fig, ax = plt.subplots(figsize=(12, .2 * painting.shape[1]))
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

    # plot painting
    ax.pcolormesh(painting.T, cmap=cmap, vmin=0, vmax=7)

    # tidy up axes
    # yticklabels = ["Gal1", "Gal2", "PSU", "TN", "Ch", "Am", "Roch", "TB", "PSA"]
    # ytickslabels = range(painting.shape[1]
    # TODO: Use plot_config.ytickslabels to establish ax.set_yticks instead of using the shape
    ax.set_yticks(np.arange(painting.shape[1]) + .5)
    ax.set_xticks(np.arange(painting.shape[0]))
    ax.set_yticklabels(plot_config.ytickslabels)
    # TODO: Print the xticks as the variant number depending on the start and end range in the plot_config
    # ax.set_xticklabels(plot_config.xtickslabels)
    ax.set_ylabel('Progeny haplotypes')
    ax.set_xlabel('Variants')
    if plot_config.title:
        ax.set_title(plot_config.title)
    plt.show()
