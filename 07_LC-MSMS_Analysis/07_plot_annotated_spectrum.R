#### MSMS Spectra similarity ####
# script to calculate similarity
# of MSMS spectra

# remove all variables
rm(list=ls(all=TRUE)) 

library(mzR)
library(MSnbase)
library(protViz)
library(stringr)
library(ggplot2)

options(stringsAsFactors = FALSE)

#To call functions defined in Spectra_functions.R
source("./07a_Spectra_functions.R")

#set working directory which contains the PSM table
# setwd("")

#set absolute file path where raw files are located
mzml_path = "e:\\data\\MSData\\Lund\\YZhang_validationPeptides\\mzML\\"

infile = "m:\\MSData\\MYakovleva\\DEanalysis\\SpectrumAI\\COPDLundPEAKS\\PSMTable\\psmtable.csv"

##input PSM table is tsv output from MS-GF plus search engine
##double check if the dataframe has column names "SpecFile" , "Peptide", "Charge", "ScanNum", "SpecEValue", "Precursor"
DF = read.csv(infile, header=T)
DF$SpectraFile = gsub(".raw", ".mzML", DF$SpectraFile)

DF$Sequence=gsub("[^A-Z]","",DF$Peptide)

Frag.ions.tolerance= 0.02 # unit Dalton
relative = FALSE

Spectra_list= vector(mode="list", length=length(unique(as.character(DF[,]$SpectraFile))))
names(Spectra_list)=unique(as.character(DF[,]$SpectraFile))

require(lattice)

pdf("PSM.annotated.spectra.pdf",width=12,height=7,useDingbats = FALSE)

for (i in 1:nrow(DF)){
    spectra_file=as.character(DF[i,]$SpectraFile)
    
    mzml_file=paste0(mzml_path,spectra_file)
    ScanNum=as.integer(DF[i,]$ScanNum)
    peptide=as.character(DF[i,]$Peptide)
    seq=DF[i,]$Sequence
    precMass=DF[i,]$Precursor
    precCharge=DF[i,]$Charge
    
    if (is.null(Spectra_list[[spectra_file]])){
        Spectra_list[[spectra_file]]=openMSfile(mzml_file,verbose=T)
    }
    
    exp_peaks<-as.data.frame(peaks(Spectra_list[[spectra_file]],scan=ScanNum))
    colnames(exp_peaks) = c("mz","intensity")
    
    predicted_peaks = predict_MS2_spectrum(Peptide = peptide, product_ion_charge = 1)
    match_ions = match_exp2predicted(exp_peaks, predicted_peaks, tolerance =Frag.ions.tolerance, relative = relative )
    
    spectrum_info = paste("Scan Number:",as.character(ScanNum),
    "precMass:",as.character(precMass),
    "precCharge:",as.character(precCharge),
    "Sequence:",seq)
    dummyAdd = (max(match_ions$intensity) - min(match_ions$intensity))
    
    dummyGgplot = ggplot(exp_peaks,aes(x=mz, ymax=intensity, ymin=0)) +geom_linerange() +
      geom_point(data = match_ions, aes(x=mz, y=intensity, color=type))+
      geom_text(data = match_ions, aes(x=mz, y=intensity + dummyAdd*0.075, label= ion),colour="black")+
      annotate("text", -Inf, Inf, label = spectrum_info, hjust = 0, vjust = 1)+
      ylab('Intensity')+xlab('M/Z')+ggtitle(peptide)+theme(plot.title = element_text(hjust = 0.5))
    
    print(dummyGgplot)
}

dev.off()

