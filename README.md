# LDcnv
Correlation-based copy number variation detection by SNP array data


## Author
Xizhi Luo, Fei Qin, Guoshuai Cai, Feifei Xiao

## Description
The detection of copy number variants (CNVs) is identifying mean shift in genetic intensities to locate chromosomal breakpoints. Many segmentation algorithms have been developed with a strong assumption of independent observations in the genetic loci, and they assume each locus has an equal chance to be a breakpoint (i.e., boundary of CNVs). However, this assumption is violated in the genetics perspective due to the existence of correlation among genomic positions such as linkage disequilibrium (LD). To generate more accurate CNVs, we therefore proposed a novel algorithm, LDcnv, that models the CNV data with its biological characteristics relating to genetic correlation (i.e., LD). 

## Installation
Install dependent packages "modSaRa" and "modSaRa2"
```r
install.packages("devtools")
library(devtools)
install_github("FeifeiXiaoUSC/modSaRa",subdir="package")
install_github("FeifeiXiaoUSC/modSaRa2",subdir="Package")
```
Install "LDcnv"
```r
install_github("adamluo12/LDcnv")
```
## Running LDcnv
LDcnv can take two types of inputs: (i) Log R ratio (LRR) and B Allele Frequency; (ii) LRR only
### Examples
(1) Load example data
```r
library(LDcnv)
data(example.data.lrr)
data(example.data.baf)
data(example.data.map)
```
```r
> head(example.data.lrr)
     X4787234217_R01C01.Log.R.Ratio X4787234217_R02C01.Log.R.Ratio X4787234217_R04C01.Log.R.Ratio
5001                   -0.011823630                   0.0479399800                    0.009759435
5002                    0.021461580                  -0.0103044200                    0.031431820
5003                    0.022457390                   0.0897155500                    0.110280500
5004                    0.000943704                  -0.0577174600                   -0.095047650
5005                    0.150350700                   0.0009902801                   -0.051548820
5006                   -0.082549070                   0.0328653200                   -0.041928530
> head(example.data.baf)
     X4787234217_R01C01.B.Allele.Freq X4787234217_R02C01.B.Allele.Freq X4787234217_R04C01.B.Allele.Freq
5001                        1.0000000                      0.997480600                        0.9986784
5002                        0.5085036                      0.001814378                        0.0000000
5003                        0.0000000                      0.000000000                        0.0000000
5004                        1.0000000                      1.000000000                        1.0000000
5005                        0.5099864                      0.508929600                        1.0000000
5006                        0.0000000                      0.525136400                        0.0000000
> head(example.data.map)
            Name Chr Position
107988 CN_473963   1    51598
108023 CN_473964   1    51671
108028 CN_473965   1    51686
108172 CN_477984   1    52015
108604 CN_473981   1    52783
108610 CN_473982   1    52800
```
(2) GC-wave adjustment by PennCNV
```
genomic_wave.pl -adjust -gcmodel hh550.gcmodel inputfile
```

(3) Run LDcnv with LRR and BAF (The packages have to be load sequentially)
```r
library(modSaRa2)
library(LDcnv)
LDcnv_eCN(lrr = example.data.lrr,baf = example.data.baf,map = example.data.map,outname="out")
```
(4) Run LDcnv with LRR (The packages have to be load sequentially)
```r
library(modSaRa) 
library(LDcnv)
LDcnv_lrr(lrr = example.data.lrr,map = example.data.map,alpha=0.01,outname="out1")
```
(5) The output file is tab delimited and has 9 columns with rows corresponding to CNV events. The columns include sample names, chromosome, CNV start marker, CNV end marker, CNV start position, CNV end position, CNV length in b, CNV length in markers, CNV status (deletion or duplication).
```
NA12045	1	CN_482241	CN_484341	25583341	25663344	80003	35	del
NA12045	1	CN_027059	CN_513307	65313984	65324244	10260	13	dup
NA12045	1	CN_515428	CN_515441	65373150	65382483	9333	16	dup
NA12045	1	SNP_A-2181454	CN_517824	72750353	72763370	13017	10	del
NA12045	1	CN_517824	SNP_A-8516738	72763370	72811904	48534	50	dup
```

