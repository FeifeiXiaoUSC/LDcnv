"LDcnv" 
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
Load example data
```r
library(LDcnv)
data(example.data.lrr)
data(example.data.baf)
data(example.data.map)
```
Testing with eCN (The packages have to be load sequentially.)
```r
library(modSaRa2)
library(LDcnv)
LDcnv_eCN(lrr = example.data.lrr,baf = example.data.baf,map = example.data.map,outname="out")
```
Testing with lrr (The packages have to be load sequentially.)
```r
library(modSaRa) 
library(LDcnv)
LDcnv_lrr(lrr = example.data.lrr,map = example.data.map,alpha=0.01,outname="out1")
```
