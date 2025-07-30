# Heterozygosity

`spack load vcftools@0.1.14`

`vcftools --vcf NW_022145594.1_hetero.vcf --het --out NW_022145594.1_hetero_heterozygocity`

For each individual, gives observed and expected homozygosity.

But first, specify inversion region:

`vcftools --vcf NW_022145594.1_hetero.vcf --chr NW_022145594.1 --from-bp 12702886 --to-bp 16793794 --recode --out NW_022145594.1_hetero_onlyinv.vcf`

