# Use cases

---

## Homozygous VCF

Plot the haplotypes of a VCF considering the following assumptions: their variants are homozygous (only one allele is considered in the plot), the chromosome of interest is `chr01` and the parental line being `SAMPLE1`.

```bash
python3 haplotype_plot/main.py -v "haplotype_plot/tests/data/chr01.vcf" -c "chr01" -p "SAMPLE1" -z HOM
```

![Example 1](images/example_1.png)

Both alleles from parental `SAMPLE1` are distinguished via `_1` and `_2` postfixes.

## Homozygous VCF, X axis redimension

In the use case shown above we had **too many** variants in the X axis to display, right? We can redimension this axis via `start` and `end` configuration (`--conf`) arguments:

```bash
python3 haplotype_plot/main.py  -v "tests/data/chr01.vcf" -c "chr01" -p "SAMPLE1" -z HOM --conf start=100 end=150
```

![Example 2](images/example_2.png)

The plot shows the 100th variant (1893) until the 150th one (2901).