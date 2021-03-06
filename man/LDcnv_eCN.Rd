% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LDcnv_eCN.R
\name{LDcnv_eCN}
\alias{LDcnv_eCN}
\title{LDcnv_eCN
This function uses both lrr and baf intensities
This function annotates the identified CNV using the reference map file and output the annotation of all identified CNVs. Each line of the output describes one CNV in nine columns: individual ID; chromosome ID; CNV start marker identifier; CNV start location (in base pair units); CNV end marker identifier; CNV end location (in base pair units); length of CNV (in base pair units); length of CNV(number of markers); copy number states (duplication or deletion).}
\usage{
LDcnv_eCN(
  lrr,
  baf,
  map,
  alpha = 0.01,
  smooth = TRUE,
  thre = 10,
  dis.thre = 5,
  outname
)
}
\arguments{
\item{lrr}{the matrix of the lrr intensities. Each column describes a single sample or sequence and each row describes a single marker}

\item{baf}{the matrix of the baf intensities. Each column describes a single sample or sequence and each row describes a single marker}

\item{map}{Each line of the map file describes a single marker and must contain exactly 3 columns: chromosome ID; rs# or marker identifier; position (in bp units)}

\item{alpha}{the significance levels for the test to accept change-points}

\item{smooth}{specify whether use smooth function to remove outliers of lrr intensities}

\item{thre}{the threshold for CNV length,default is 10}

\item{dis.thre}{the threshold for distance between CNVs for merging adjacent closely located CNVs,default is 5}

\item{outname}{name for the output file}
}
\value{
This function generates a text file describing all detected CNVs. In addition, it also returns a list of detected change-points for all samples.

cp a list of position index for the final change-points identified by modSaRa
}
\description{
LDcnv_eCN
This function uses both lrr and baf intensities
This function annotates the identified CNV using the reference map file and output the annotation of all identified CNVs. Each line of the output describes one CNV in nine columns: individual ID; chromosome ID; CNV start marker identifier; CNV start location (in base pair units); CNV end marker identifier; CNV end location (in base pair units); length of CNV (in base pair units); length of CNV(number of markers); copy number states (duplication or deletion).
}
\examples{
# Input the example data of SNP genotyping data from Affymatrix Human SNP Array 6.0 platform.
# The map file displays annotation of the markers including the chromosome and location
# information of each SNP or CNV marker.
data(example.data.lrr)
data(example.data.baf)
data(example.data.map)
LDcnv_eCN(lrr = example.data.lrr,baf = example.data.baf,map = example.data.map,outname="out")
# The following file will be generated: "out.csv"
# This file contains CNV output for each individual.
# Each line represents one CNV detected from one sample or sequence.
# For each line, the individual ID, start position, end position, length and state
# (duplication or deletion) of the CNV will be shown.
}
\seealso{
\link{modifiedSaRa} for processing the modified SaRa method
}
