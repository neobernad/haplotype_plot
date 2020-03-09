# Haplotype plot

Generates haplotype plots for VCF files.
**Requires Python 3 to run**.

[![Build Status](https://travis-ci.org/neobernad/haplotype_plot.svg?branch=master)](https://travis-ci.org/neobernad/haplotype_plot)
[![codecov](https://codecov.io/gh/neobernad/haplotype_plot/branch/master/graph/badge.svg)](https://codecov.io/gh/neobernad/haplotype_plot)

# Pypi

You may install haplotype plot from pip:

```bash
pip install haplotype-plot
```

# Usage

```bash
usage: main.py [-h] -v VCF [-o OUTPUT] -c CHR -p PARENTAL [--phase]
               [--version] [-z {HOM,HET}] [--conf KEY=VALUE [KEY=VALUE ...]]

Creates a haplotype plot from a VCF file

optional arguments:
  -h, --help            show this help message and exit
  -v VCF, --vcf VCF     path to the input sorted VCF file
  -o OUTPUT, --output OUTPUT
                        path to the output PNG haplotype plot file
  -c CHR, --chr CHR     chromosome to plot from the VCF
  -p PARENTAL, --parental PARENTAL
                        sample name from the VCF used as parental haplotype
  --phase               if specified, it will phase the genotypes. An output
                        VCF with the phased genotypes is created
  --version             show program's version number and exit
  -z {HOM,HET}, --zygosis {HOM,HET}
                        zygosis of the input VCF file
  --conf KEY=VALUE [KEY=VALUE ...]
                        set a number of key-value pairs to modify the default
                        plot configuration (do not put spaces before or after
                        the = sign). Possible values are: 'start=integer',
                        'end=integer'. If a value contains spaces, you should
                        define it with double quotes: foo="this is a
                        sentence".
```

For examples of use please go to the [documentation](https://neobernad.github.io/haplotype_plot/#/).

---