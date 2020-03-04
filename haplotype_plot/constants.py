# -*- coding: utf-8 -*-
import numpy

HDF5_VARIANTS_KEY: str = "variants"
HDF5_CALLDATA_KEY: str = "calldata"
HDF5_GENOTYPE_KEY: str = "GT"
HDF5_CALLDATA_GENOTYPE_KEY: str = HDF5_CALLDATA_KEY + "/" + HDF5_GENOTYPE_KEY
HDF5_EXT: str = ".h5"


MISSING_INT = -1
MISSING_FLOAT = float('nan')
MISSING_STR = ''
MISSING_BYTE = b''
MISSING_BOOL = False
TRUE_INT = 1
FALSE_INT = 0

GT_FIELD = '/calls/GT'
GQ_FIELD = '/calls/GQ'
ID_FIELD = '/variations/id'
ALT_FIELD = '/variations/alt'
REF_FIELD = '/variations/ref'
QUAL_FIELD = '/variations/qual'
DP_FIELD = '/calls/DP'
AO_FIELD = '/calls/AO'
RO_FIELD = '/calls/RO'
AD_FIELD = '/calls/AD'
CHROM_FIELD = '/variations/chrom'
POS_FIELD = '/variations/pos'
INFO_FIELD = '/variations/info'


class _MissingValues():
    def __init__(self):
        self._missing_values = {int: MISSING_INT,
                                'Integer': MISSING_INT,
                                float: MISSING_FLOAT,
                                'Float': MISSING_FLOAT,
                                str: MISSING_STR,
                                'String': MISSING_STR,
                                numpy.int8: MISSING_INT,
                                numpy.int16: MISSING_INT,
                                numpy.int32: MISSING_INT,
                                numpy.float16: MISSING_FLOAT,
                                numpy.float32: MISSING_FLOAT,
                                numpy.bool_: MISSING_BOOL,
                                numpy.bytes_: MISSING_BYTE,
                                bool: MISSING_BOOL}

    def __getitem__(self, dtype):
        str_dtype = str(dtype)
        if dtype in self._missing_values:
            return self._missing_values[dtype]
        elif isinstance(dtype, str):
            if 'str' in dtype:
                return MISSING_STR
            elif 'int' in dtype:
                return MISSING_INT
            elif 'float' in dtype:
                return MISSING_FLOAT
            elif dtype[0] == 'S':
                return MISSING_BYTE
            elif dtype[:2] == '|S':
                return MISSING_BYTE
        elif 'int' in str_dtype:
            return MISSING_INT
        elif 'float' in str_dtype:
            return MISSING_FLOAT
        elif 'bool' in str_dtype:
            return MISSING_BOOL
        elif str_dtype[:2] == '|S':
            return MISSING_BYTE
        elif str_dtype[:2] == '<U':
            return MISSING_STR
        else:
            raise ValueError('No missing type defined for type: ' + str(dtype))


MISSING_VALUES = _MissingValues()

# Speed is related to chunksize, so if you change snps-per-chunk check the
# performance
SNPS_PER_CHUNK = 600
VCF_FORMAT = 'VCFv4.2'
DEF_DSET_PARAMS = {
                   # Not much slower than lzf but compresses much more
                   'compression': 'gzip',
                   'shuffle': True,
                   # checksum, slower but safer
                   'fletcher32': True}
