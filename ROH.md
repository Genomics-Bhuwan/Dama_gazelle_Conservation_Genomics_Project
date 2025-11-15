#### ROHs.
- Runs of Homozygosity are the stretches of DNA where an individual inherits identical haplotypes from both parents.
- Occurs because of inbreeding which can occur when closely related with invidual mates.
- Runs of Homozygosity are stretches of DNA where an individual inherits identical haplotypes from both parents. They often are a result of inbreeding, which can occur when 
- **BCFtools* roh will be used to estimate the RUNS OF HOMOZYGOSITY.

#### Setting up for analyses.
- Check the depth of the BAM files.
- It is required to have at least ~15X depth in genome to confidently assess for runs of homozygosity.
- It will be done using BCFtools
- For Dama gazelle, the depth of each sample is >90X and we should be good to go with this.

                                                                                                                                                 
#                                                                                                                                                                       
# ----------------Modules------------------------- #                                                                                                                    
```bash
module load bcftools-1.15                                                                                                                                               
```

#### Run the ROH using BCFtools.
- bcftools roh: This command runs the roh plugin from bcftools to detect runs of homozygosity.
- AF-dflt 0.4: This option sets the default allele frequency to 0.4 when the frequency is missing in the input file.
- G30: This option sets the phred-scaled genotype quality threshold to 30. Genotypes below this quality threshold will be treated as missing.
```bash
bcftools roh -G30 --AF-dflt 0.4 $VCF -o ${OUTDIR}/Dama_gazelle_ROH_all.txt

```

#### Get the file with only Runs of homozygosity for all the samples. Nothing else.
- command will produce a new file that only contains the ROH regions RG lines for all give samples.
- If your bcftools roh output file contains results for all five samples, then running:
- Doesn't include other metadata or summary lines like (like AF-dflt, number of sites, etc.)
- It will include all ROH segments for all five samples in the file, not filtered per sample.
  
```bash
grep "RG" Dama_gazelle_ROH_all.txt > Dama_gazelle_ROH_all_RG.txt  
```


#### Plotting in R.

library (ggplot2)
library(tidyverse)
library(dplyr)

############################################
### Plotting bar graph of ROH categories ###
############################################

#read in data
dat_0.1 <- read.csv ("RunsofHomozygosity.csv", header=TRUE)
nrow(dat_0.1)
#7338
##Remove ROH under 100kb (100kb=100,000bp)
dat <- subset(dat, Length.bp. > 100000)
nrow(dat)
#472

##################################################
##################################################
#Let's look at the number of ROH in each sample ## 
##################################################
##################################################
dat_Atlantisia_rogersi <- subset(dat, Sample == "Atlantisia_rogersi")
nrow(dat_Atlantisia_rogersi) 
#183
dat_GuamRail_HOW_N23_0063 <- subset(dat, Sample == "GuamRail_HOW_N23_0063")
nrow(dat_GuamRail_HOW_N23_0063) 
#61
dat_GuamRail_HOW_N23_0068 <- subset(dat, Sample == "GuamRail_HOW_N23_0068")
nrow(dat_GuamRail_HOW_N23_0068) 
#127
dat_Rallus_limicola <- subset(dat, Sample == "Rallus_limicola")
nrow(dat_Rallus_limicola) 
#101

#Let's make a bar graph for length vs. number of ROH

ROH_lengthvsnumber <- read.csv("ROH_lengthvsNumber.csv", header=TRUE)

 ggplot(ROH_lengthvsnumber,aes(x=Total_length_ROH,y=NumberofROH,label=Sample, color=Sample)) + geom_point(size = 3) + theme_classic()

#What different patterns do we see? 

######################################################
######################################################
### Let's make a graph of different lengths of ROH ###
######################################################
######################################################

##################################
# Get ROH values between 0.1-1Mb
##################################
dat_0.1_1Mb <- subset(dat, Length.bp. < 1000000)
sum <- aggregate(dat_0.1_1Mb$Length.bp., by=list(Category=dat_0.1_1Mb$Sample), FUN=sum)
sum

#               Category        x
#1    Atlantisia_rogersi 42456648
#2 GuamRail_HOW_N23_0063 13222867
#3 GuamRail_HOW_N23_0068 28095326
#4       Rallus_limicola 18974775

write.csv (dat_0.1_1Mb , "samples_between_0.1_1mb.csv")
##################################
# Get ROH values between 1-5Mb
##################################
dat_1_10Mbtest <- subset(dat, Length.bp. > 1000000)
dat_1_10Mb <- subset(dat_1_10Mbtest, Length.bp. < 5000000)
sum_2_5Mb <- aggregate(dat_1_10Mb$Length.bp., by=list(Category=dat_1_10Mb$Sample), FUN=sum)
sum_2_5Mb

#               Category       x
#1 GuamRail_HOW_N23_0063 5020811
#2       Rallus_limicola 1055016


##################################
# Get ROH values between 5-10Mb
##################################
dat_10_100Mb <- subset(dat, Length.bp. > 5000000)
dat_10_100Mbfinal <- subset(dat_10_100Mb, Length.bp. < 10000000)
sum_5_10Mb <- aggregate(dat_10_100Mbfinal$Length.bp., by=list(Category=dat_10_100Mbfinal$Sample), FUN=sum)
sum_5_10Mb

#               Category        x
#1 GuamRail_HOW_N23_0063 20835957

##################################
# Get ROH values between 10-100Mb
##################################
dat_10_100Mb <- subset(dat, Length.bp. > 10000000)
dat_10_100Mbfinal <- subset(dat_10_100Mb, Length.bp. < 100000000)
sum_10_100Mb <- aggregate(dat_10_100Mbfinal$Length.bp., by=list(Category=dat_10_100Mbfinal$Sample), FUN=sum)

#no rows 

######################
#Plotting by category #
######################

datrohlength <-read.csv("ROH_byLength_Category.csv", header=TRUE)

p <- ggplot(datrohlength, aes(fill=ROHcat, y=Length, x=reorder(Sample,-Totalsize ))) + 
    geom_bar(position="stack", stat="identity", color="black") + theme_classic() + coord_flip() + scale_fill_manual(values = c( "#e6ab02", "#d95f02" ,"#7570b3")) + scale_x_discrete(labels=c("5.0e+08" = "0.5", "1.0e+09" = "1","1.5e+09" = "1.5"))+  theme(axis.text.x = element_text(color="black", size = 10), axis.text.y = element_text(color="black", size = 10))
p




######################################################
######################################################
############# Plot cummulative ROH ###################
######################################################
######################################################

dat <- read.csv("RunsofHomozygosity_forcummulative_final.csv", header=TRUE)

#Filter out ROH
dat <- subset(dat, Length.bp. > 100000)
dat <- subset(dat, Quality > 80)

#Plotting
p <- ggplot(dat, aes(x=Length.bp., y=Cummulative, color=Sample)) + theme_classic() + scale_x_log10()   +
  geom_line(size=1.2, alpha=0.75) + scale_color_manual(values=c("orangered2", "mediumseagreen", "royalblue", "orange", "mediumseagreen" ))
p
q <- p + geom_line(aes(size = Bold))  +
  scale_size_manual(values = c(0.1, 1.5))  
q 
