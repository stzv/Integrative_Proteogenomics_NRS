### IMPORTANT ####
## Script to trasnfor chipanse ID to human
## Author: Peter Horvatovich
## Date: September 10, 2018
##

# clear all variables
rm(list = ls(all = TRUE))

suppressMessages(library(MSnbase))
suppressMessages(library(mzR))
suppressMessages(library(stringr))
options(stringsAsFactors = FALSE)

#Specify the folder where synthetic and experiment spectra (mzML files) are stored 
synpep_mzml_path <- "SwisProt/SwissProt_ORF/150225_18_01_c.mzML" ### File from predicted
endopep_mzml_path <- "SwisProt/SwissProt_ORF/150225_18_01_c.mzML" ### File from synthetised

#set working directory which contains the PSM table

#Set the file name of the PSM table
infile1 <- "SwisProt/SwissProt_ORF/COPD8_SwissProt_DPD_normalSearchMSFragger_DPCOPD8_psm.tsv" ### File from predicted
infile2 <- "SwisProt/SwissProt_ORF/COPD8_SwissProt_DPD_normalSearchMSFragger_DPCOPD8_psm.tsv" ### File from synthetised

#double check column names "SpecFile" , "Peptide", "Charge", "scanNum", "SpecEValue", "Precursor" 
DF1 <- read.csv(infile1, header = T, sep = "\t") 
DF1$mzML = paste0("SwisProt/SwissProt_ORF/",sapply(DF1$Spectrum, function(x){unlist(strsplit(x, "\\."))[1]}), ".mzML") 
DF1$scan = as.numeric(sapply(DF1$Spectrum, function(x){unlist(strsplit(x, "\\."))[2]}))


DF2 <- read.csv(infile2, header = T, sep = "\t")
DF2$mzML <- paste0("SwisProt/SwissProt_ORF/",sapply(DF2$Spectrum, function(x){unlist(strsplit(x, "\\."))[1]}), ".mzML")
DF2$scan <- as.numeric(sapply(DF2$Spectrum, function(x){unlist(strsplit(x, "\\."))[2]}))

# a list of peptides to draw mnirror image.
# modifications on the peptides must be noted in same way as in the above PSM table.
peptide1 <- "YFYNQEEYVR"
peptide2 <- "YFYNQEEYVR"

df.pep1 <- DF1[DF1$Peptide == peptide1,]
if (nrow(df.pep1)>1){
  df.pep1 = df.pep1[which(df.pep1$Hyperscore == max(df.pep1$Hyperscore)),]
}


# choose the best scored PSM to draw mirror image
scanNum.pep1 <- as.integer(df.pep1$scan)
precMass.pep1 <- as.numeric(df.pep1$Observed.M.Z)
precCharge.pep1 <- as.integer(df.pep1$Charge)
mzml_file.pep1 <- df.pep1$mzML

df.pep2 <- DF2[DF2$Peptide == peptide2, ]

if (nrow(df.pep2) == 0) {
  sprintf("%s no corresponding PSMs of synthetic peptide ", peptide2)#; #next
  }
if (nrow(df.pep2) > 1) {
  df.pep2 <- df.pep2[which(df.pep2$Hyperscore == max(df.pep2$Hyperscore)),]
  }

scanNum.pep2 <- as.integer(df.pep2$scan)
precMass.pep2 <- as.numeric(df.pep2$Observed.M.Z)
mzml_file.pep2 <- df.pep2$mzML


spectra_file1 <- mzml_file.pep1 #endogenous peptide
spectra_file2 <- mzml_file.pep2 # synthetic peptide

rawdata1 <- openMSfile(df.pep1$mzML)
rawdata2 <- openMSfile(df.pep2$mzML)

exp_peaks1 <- peaks(rawdata1, scan = scanNum.pep1)
exp_peaks2 <- peaks(rawdata2, scan = scanNum.pep2)


pep1_spectrum = new("Spectrum2", intensity = exp_peaks1[,2], mz = exp_peaks1[,1], centroided=TRUE, precursorMz = precMass.pep1,
                     precScanNum = scanNum.pep1)

pep2_spectrum = new("Spectrum2",intensity = exp_peaks2[,2], mz = exp_peaks2[,1], centroided=TRUE, precursorMz = precMass.pep2,
                    precScanNum = scanNum.pep2)

#Peptide mirror image plot, top is endogenous, bottom is synthetic

jpeg(file = paste0("06_", peptide1, "_", peptide2, ".jpg"),
     width = 800, height = 600, units = "px")

plot(pep1_spectrum, pep2_spectrum, tolerance = 0.05, relative = FALSE,
     sequences = c(peptide1, peptide2), modifications = c(C = 57.021),
     z = seq(1, 2, by = 1), neutralLoss = NULL, peaks.cex = 0, peaks.lwd = 1)
legend("topright", "Predicted YFYNQEEYVR")
legend("bottomright", "Synthetic YFYNQEEYVR")

dev.off()