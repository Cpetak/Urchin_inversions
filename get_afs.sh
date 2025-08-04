
cname=610
bcftools view -i 'ANN[*] ~ "HIGH"' ${cname}_snpeff.vcf > high_${cname}_snpeff.vcf
vcftools --vcf high_${cname}_snpeff.vcf --freq --out high_${cname}_snpeff_freq.txt
bcftools view -i 'ANN[*] ~ "LOW"' ${cname}_snpeff.vcf > low_${cname}_snpeff.vcf
vcftools --vcf low_${cname}_snpeff.vcf --freq --out low_${cname}_snpeff_freq.txt
bcftools view -i 'ANN[*] ~ "MODERATE"' ${cname}_snpeff.vcf > mod_${cname}_snpeff.vcf
vcftools --vcf mod_${cname}_snpeff.vcf --freq --out mod_${cname}_snpeff_freq.txt

cname=610_inv
bcftools view -i 'ANN[*] ~ "HIGH"' ${cname}_snpeff.vcf > high_${cname}_snpeff.vcf
vcftools --vcf high_${cname}_snpeff.vcf --freq --out high_${cname}_snpeff_freq.txt
bcftools view -i 'ANN[*] ~ "LOW"' ${cname}_snpeff.vcf > low_${cname}_snpeff.vcf
vcftools --vcf low_${cname}_snpeff.vcf --freq --out low_${cname}_snpeff_freq.txt
bcftools view -i 'ANN[*] ~ "MODERATE"' ${cname}_snpeff.vcf > mod_${cname}_snpeff.vcf
vcftools --vcf mod_${cname}_snpeff.vcf --freq --out mod_${cname}_snpeff_freq.txt
