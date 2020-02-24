# -*- coding: utf-8 -*-
import logging
import numpy as np
import allel


def variants_filter(variants: allel.VariantChunkedTable, filter_expression: str):
    variants_np_array = variants[:]  # Conversion to allel.model.ndarray.VariantTable
    variant_selection = variants_np_array.eval(filter_expression)[:]
    return np.count_nonzero(variant_selection)
