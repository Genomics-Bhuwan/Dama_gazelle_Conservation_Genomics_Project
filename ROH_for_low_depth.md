#### I am trying to calculate the ROH for low depth sequences from 2X-10X. 
- Steps I have to do each time before running ROHan for each samples
-  Find out teh TSTV ratio for each sample as the package ROHan needs this.
### Step 1. For Grant's gazelle.
```bash
module load biocontainers/default
module load bcftools/1.17

bcftools stats Grants_Gazelle_filtered_maf05_biallelic.vcf.gz | grep "^TSTV"
```
