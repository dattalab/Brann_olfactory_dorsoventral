#!/usr/bin/env bash

cd ../data/tables

VCF=mgp_REL2021_snps.vcf.gz
if [[ ! -e $VCF ]]; then
    wget https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-2112-v8-SNPs_Indels/mgp_REL2021_snps.vcf.gz
    wget https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-2112-v8-SNPs_Indels/mgp_REL2021_snps.vcf.gz.csi
fi

bcftools view -O b -R mm39_olfr.bed $VCF >mgp_REL2021_snps.olfr.vcf.gz
bcftools index mgp_REL2021_snps.olfr.vcf.gz

cat >strains.txt <<'EOF'
C57BL_6NJ
CAST_EiJ
EOF

# filter to only keep homozygous variants that differ between strains
bcftools view -S strains.txt -c 1 -Oz -o mgp_REL2021_snps.olfr.CAST.vcf.gz mgp_REL2021_snps.olfr.vcf.gz
bcftools index mgp_REL2021_snps.olfr.CAST.vcf.gz

## make snpEff database to annotate
## https://github.com/pcingola/SnpEff/issues/536#issuecomment-2189201996
# vi snpEff.config
# $(mm39.genome : Mouse)
# java -Xmx4g -jar snpEff.jar build -gtf22 -v mm39

# module load java
export MGP= # <tables_folder>
cd ~/snpEff

bcftools annotate -x "INFO,^FORMAT/GT,FORMAT/PL" \
    $MGP/mgp_REL2021_snps.olfr.CAST.vcf.gz -Oz \
    -o $MGP/snpEff_annotated/mgp_REL2021_snps.olfr.CAST.clean.vcf.gz

java -Xmx8g -jar snpEff.jar -v mm39 \
    $MGP/snpEff_annotated/mgp_REL2021_snps.olfr.CAST.clean.vcf.gz \
    >$MGP/snpEff_annotated/mgp_REL2021_snps.olfr.CAST.clean.ann.vcf.gz

# extract reads from BAM file that overlap OR genes
samtools merge -r -b f1.txt -L $META/MGP/mm39_olfr.bed -o bam_gene_subset/F1_mm39_olfr.bam -@ $SLURM_CPUS_PER_TASK