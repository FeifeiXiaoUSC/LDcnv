#' LDcnv_eCN
#' This function uses both lrr and baf intensities
#' This function annotates the identified CNV using the reference map file and output the annotation of all identified CNVs. Each line of the output describes one CNV in nine columns: individual ID; chromosome ID; CNV start marker identifier; CNV start location (in base pair units); CNV end marker identifier; CNV end location (in base pair units); length of CNV (in base pair units); length of CNV(number of markers); copy number states (duplication or deletion).
#' @param lrr the matrix of the lrr intensities. Each column describes a single sample or sequence and each row describes a single marker
#' @param baf the matrix of the baf intensities. Each column describes a single sample or sequence and each row describes a single marker
#' @param map Each line of the map file describes a single marker and must contain exactly 3 columns: chromosome ID; rs# or marker identifier; position (in bp units)
#' @param alpha the significance levels for the test to accept change-points
#' @param smooth specify whether use smooth function to remove outliers of lrr intensities
#' @param thre the threshold for CNV length,default is 10
#' @param dis.thre the threshold for distance between CNVs for merging adjacent closely located CNVs,default is 5
#' @param outname name for the output file
#' @return This function generates a text file describing all detected CNVs. In addition, it also returns a list of detected change-points for all samples.
#' @return cp a list of position index for the final change-points identified by modSaRa
#' @seealso \link{modifiedSaRa} for processing the modified SaRa method
#' @examples
#' # Input the example data of SNP genotyping data from Affymatrix Human SNP Array 6.0 platform.
#' # The map file displays annotation of the markers including the chromosome and location
#' # information of each SNP or CNV marker.
#' data(example.data.lrr)
#' data(example.data.baf)
#' data(example.data.map)
#' LDcnv_eCN(lrr = example.data.lrr,baf = example.data.baf,map = example.data.map,outname="out")
#' # The following file will be generated: "out.csv"
#' # This file contains CNV output for each individual.
#' # Each line represents one CNV detected from one sample or sequence.
#' # For each line, the individual ID, start position, end position, length and state
#' # (duplication or deletion) of the CNV will be shown.
#' @export
LDcnv_eCN=function(lrr,baf,map,alpha=0.01,smooth=TRUE,thre=10,dis.thre=5,outname){
  lrr=data.matrix(lrr)
  if(smooth==TRUE){
    lrr <- smooth(lrr, R = 10, t = 2)
  }else{
    lrr=lrr
  }
  eCN <-function(lrr, baf){
    laf <- baf
    laf[laf > 0.5 & !is.na(laf)] <- 1-laf[laf > 0.5 & !is.na(laf)]
    e_CN <-  matrix(NA, dim(lrr)[1], dim(lrr)[2])
    cn.est.L <- list()
    for (i in 1:dim(lrr)[2]) {
      cn.est.L[[i]] <- list()  
      for (k in 1:length(lrr[,i])) {
        cn.est <-  Likeli.single(lrr=lrr[k,i],laf=laf[k,i])
        e_CN[k,i] <- cn.est$X
        cn.est.L[[i]][[k]] <- NaN*5
        cn.est.L[[i]][[k]]<- cn.est$L
      }
    }
    e_CN.smo <- smooth(e_CN, R=10, t=2)   
    return (list(e_CN= e_CN,laf=laf))
  }
  eCN=eCN(lrr=lrr,baf =baf )
  eCN.cal=eCN$e_CN
  laf=eCN$laf
  h=5
  h_vec=c(rep(1/h,h),rep(-1/h,h))
  empvar=vector()
  for(d in (h+1):(dim(eCN.cal)[1]-h)){
    local_cov=cov(t(eCN.cal[(d-h):(d+h-1),]))
    ecnvar=t(h_vec)%*%local_cov%*%h_vec
    empvar[d]=ecnvar
  }
  empsd1=sqrt(empvar[!is.na(empvar)])
  
  
  h=10
  h_vec=c(rep(1/h,h),rep(-1/h,h))
  empvar=vector()
  for(d in (h+1):(dim(eCN.cal)[1]-h)){
    local_cov=cov(t(eCN.cal[(d-h):(d+h-1),]))
    ecnvar=t(h_vec)%*%local_cov%*%h_vec
    empvar[d]=ecnvar
  }
  empsd2=sqrt(empvar[!is.na(empvar)])
  
  
  h=15
  h_vec=c(rep(1/h,h),rep(-1/h,h))
  empvar=vector()
  for(d in (h+1):(dim(eCN.cal)[1]-h)){
    local_cov=cov(t(eCN.cal[(d-h):(d+h-1),]))
    ecnvar=t(h_vec)%*%local_cov%*%h_vec
    empvar[d]=ecnvar
  }
  empsd3=sqrt(empvar[!is.na(empvar)])
  
  
  fInverse <-
    function(n = 10000, h= 10, hh = 2 * h, precise = 10000, emp,simT = 2000){      
      empirical = NULL
      for (i in 1 : simT){
        Y=vector()
        for(p in 1:n){
          Y[p]=rnorm(1,mean = 0,sd=emp[p]) 
        }
        LDF  =  Y
        LDF.pos      =  LDF
        LDF.neg      =  -LDF
        index.pos = localMax(y = LDF.pos, span = hh)
        pV.pos   = 1 - 2 * pnorm(LDF.pos[index.pos] / (emp[index.pos]))
        index.neg = localMax(y = LDF.neg, span = hh)
        pV.neg   = 1 - 2 * pnorm(LDF.neg[index.neg] / (emp[index.neg]))
        index <- c(index.pos, index.neg)
        pv <- c(pV.pos, pV.neg)
        pv <- pv[order(index)]
        index <- sort(index)
        len <- length(index)
        rm.list <- NULL
        for (j in 1 : (len - 1)) {
          if(index[j] >= index[j + 1] - h){
            rm.list <- c(rm.list, j)
          }
        }
        if (length(rm.list) > 0) {
          pv <- pv[-rm.list]
        }
        empirical <- c(empirical, pv)
        if (length(empirical) > 10 * precise) break
      }
      return(quantile(empirical, probs = c(0 : precise) / precise))
    }
  
  FINV=list()
  FINV[[1]]=fInverse(n=length(empsd1),h=5,precise =10000,emp = empsd1)
  FINV[[2]]=fInverse(n=length(empsd2),h=10,precise =10000,emp = empsd2)
  FINV[[3]]=fInverse(n=length(empsd3),h=15,precise =10000,emp = empsd3)
  cp <-  CNVout(e_CN = eCN.cal,lrr=lrr,laf=laf,map=map, FINV=FINV, alpha =alpha, thre=thre,dis.thre=dis.thre, outname = outname)$cp
}



