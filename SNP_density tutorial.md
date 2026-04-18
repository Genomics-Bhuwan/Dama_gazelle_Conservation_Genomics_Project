##### SNP density tutorial
- SNP density is a tool used for visualizing the genetic diversity across the genome.
- It represents the number of SNPs in a heatmap.
- It is used for comparing genetic diversity between different populations or closely related species.
- To calculate the SNP density, count the number of SNPs over a certain length of the genome, typically ranging from 100kb to 1 megabase.
- Use VCFtools with the snpden function(Danecek et al., 2011).
- Vcf as input file and a table with the snp count for reach region interval.

```bash
- I have five samples. What I need to do is split of the combined vcf to per sample vcf.
- Sometime, if a command is ran witout reporting errors, it doesnot guarantee that the output is correct. In that case, check the output manually and calcualte some metrics for input and output. You might get a clue to what went wrong. 
- To be sure, I calculated the numbers of different genotupes before split, after it and after extraction of heterozygotes sites. The number match ideally.

```

#### From Beginning from the last time
- Make vcf bialleic after all the hard and soft filtering
- Only keeping the 1-17 autosomes.
- 


#############################################################################################################
#############################################################
##### I am doing from beginning for the Dama one last time as suggested by Sergei
```bash
#!/usr/bin/env bash
# This script assumes that draw_variant_window_densities.py (MACE script) is installed globally. 

# DP1800 filter is set incorrectly and should be inverted. FILTER column in vcf file should contain a lable of the applied filter only if variant have failed to path the filter. Other filters (usual GATK) seems to be OK. 
# So good variants should contain "DP1800" only. 

#bcftools filter -i 'FILTER="DP1800"' -O z Dama_gazelle_biallelic_snps_autosomes.vcf.gz > Dama_gazelle_biallelic_snps_autosomes.filtered.vcf.gz;

mkdir -p  split_samples; 
for SAMPLE in `bcftools query -l Dama_gazelle_biallelic_snps_autosomes.filtered.vcf.gz`;
    do
    echo "Handling ${SAMPLE}..."
    mkdir -p split_samples/${SAMPLE}
    # count all genotypes from the filtered file for all samples
    bcftools view -s ${SAMPLE} -O v Dama_gazelle_biallelic_snps_autosomes.filtered.vcf.gz | grep -v "^#" | cut -f 10 | cut -f 1 -d : | sort  | uniq -c > Dama_gazelle_biallelic_snps_autosomes.filtered.${SAMPLE}.all.genotype.test;
    # get per-sample files without reference-like (0/0 and 0|0) and uncalled (./.) sites
    bcftools view -s ${SAMPLE} -c 1 -O z Dama_gazelle_biallelic_snps_autosomes.filtered.vcf.gz > split_samples/${SAMPLE}/Dama_gazelle_biallelic_snps_autosomes.filtered.${SAMPLE}.vcf.gz
    # extract heterozygous positions only
    bcftools filter -i "GT='het'" -O z split_samples/${SAMPLE}/Dama_gazelle_biallelic_snps_autosomes.filtered.${SAMPLE}.vcf.gz > split_samples/${SAMPLE}/Dama_gazelle_biallelic_snps_autosomes.filtered.${SAMPLE}.hetero.vcf.gz

	for FILE in split_samples/${SAMPLE}/*.vcf.gz; 
		do 
		zcat ${FILE} | grep -v "^#" | cut -f 10 | cut -f 1 -d : | sort  | uniq -c > ${FILE}.genotype.test; 
	    done

	draw_variant_window_densities.py -i split_samples/${SAMPLE}/Dama_gazelle_biallelic_snps_autosomes.filtered.${SAMPLE}.hetero.vcf.gz  \
	                                 -o split_samples/${SAMPLE}/Dama_gazelle_biallelic_snps_autosomes.filtered.${SAMPLE}.hetero.w1000k_s100k \
	                                 -w 1000000 \
	                                 -s 100000 \
	                                 --title "Heterozygous SNP densities (${SAMPLE})" \
	                                 -a Dama_gazelle.chr.whitelist  \
	                                 -z Dama_gazelle.chr.orderlist \
	                                 --scaffold_syn_file Dama_gazelle.chr.syn  \
	                                 --syn_file_key_column 0 \
	                                 --syn_file_value_column 1 \
	                                 --density_thresholds 0,0.1,0.5,0.75,1.0,1.25,1.5,2.0,2.5 \
	                                 --hide_track_label \
	                                 --rounded;

	draw_variant_window_densities.py -i split_samples/${SAMPLE}/Dama_gazelle_biallelic_snps_autosomes.filtered.${SAMPLE}.hetero.vcf.gz  \
	                                 -o split_samples/${SAMPLE}/Dama_gazelle_biallelic_snps_autosomes.filtered.${SAMPLE}.hetero.w100k_s10k \
	                                 -w 100000 \
	                                 -s 10000 \
	                                 --title "Heterozygous SNP densities (${SAMPLE})" \
	                                 -a Dama_gazelle.chr.whitelist  \
	                                 -z Dama_gazelle.chr.orderlist \
	                                 --scaffold_syn_file Dama_gazelle.chr.syn  \
	                                 --syn_file_key_column 0 \
	                                 --syn_file_value_column 1 \
	                                 --density_thresholds 0,0.1,0.5,0.75,1.0,1.25,1.5,2.0,2.5 \
	                                 --hide_track_label \
	                                 --rounded;
	done

```

#### Step 2. 
- Make a file named "Dama_gazelle.chr.orderlist" and keep below name in this file
  ```bash
chr1
chr2
chr3
chr4
chr5
chr6
chr7
chr8
chr9
chr10
chr11
chr12
chr13
chr14
chr15
chr16
chr17
```
- Dama_gazelle.chr.syn
```bash
1	chr1
2	chr2
3	chr3
4	chr4
5	chr5
6	chr6
7	chr7
8	chr8
9	chr9
10	chr10
11	chr11
12	chr12
13	chr13
14	chr14
15	chr15
16	chr16
17	chr17
```
- Dama_gazelle.chr.whitelist
```bash
1
2
3
4
5
6
7
8
9
10
11
12
13
14
15
16
17
```
