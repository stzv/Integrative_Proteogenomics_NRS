###
## Script to compare endogenous and synthetic peptides MS/MS
## Author: Peter Horvatovich
## Modified by: Stepanka Zverinova
##

# clear all variables
rm(list = ls(all = TRUE))

suppressMessages(library(MSnbase))
suppressMessages(library(mzR))
suppressMessages(library(stringr))
options(stringsAsFactors = FALSE)
suppressMessages(library(gridExtra))
suppressMessages(library(grid))
suppressMessages(library(ggplot2))


# Peptide info
args <- commandArgs(trailingOnly = TRUE)
peptide <- args[1]


peptide_info <- read.csv("04_NovelPeptides_Merged.tsv",
                        header = TRUE, sep = "\t")


#### Highest Hyperscore
infile <- peptide_info[peptide_info$X.Peptide == peptide, "HighestSpectrumFile"]
dir <- str_split_fixed(infile, "/", n = Inf)
dir <- paste(dir[1 : length(dir) - 1], collapse = "/")

pep1_samples <- peptide_info[peptide_info$X.Peptide == peptide, "MSSamples"]

df1 <- read.csv(infile, header = TRUE, sep = "\t")
df1$mzML <- paste0(dir, "/",
                  sapply(df1$Spectrum, function(x) {unlist(strsplit(x, "\\."))[1]}), ".mzML", #nolint
                  sep = "")

df1$scan <- as.numeric(sapply(df1$Spectrum,
                              function(x) {unlist(strsplit(x, "\\."))[2]}))

df_pep1 <- df1[df1$Peptide == peptide, ]
pep1_psms <- length(df_pep1$Spectrum)

df_pep1 <- df_pep1[order(df_pep1$Hyperscore, decreasing = TRUE), ]
df_pep1 <- df_pep1[which(df_pep1$Hyperscore == max(df_pep1$Hyperscore)), ]

scan_num_pep1 <- as.integer(df_pep1$scan)
prec_mass_pep1 <- as.numeric(df_pep1$Observed.M.Z)
prec_charge_pep1 <- as.integer(df_pep1$Charge)

peptide_mod <- ifelse(df_pep1$Modified.Peptide == "",
                      df_pep1$Peptide,
                      df_pep1$Modified.Peptide)

rawdata1 <- openMSfile(df_pep1$mzML)

aa1 <- header(rawdata1, scan = scan_num_pep1)
rawdata1 <- openMSfile(df_pep1$mzML)
exp_peaks1 <- peaks(rawdata1, scan = scan_num_pep1)

# Print to file
log_file <- "05_Spectra_Info.tsv"
if (file.exists(log_file) == FALSE) {
    file.create(file.path(log_file), showWarnings = FALSE)
    header_str <- paste("Peptide", "Origin", "PSMs", "Sample#", 'Precursor_mz', "Precursor_Intensity",'Retention_Time', 'Charge', "Hyperscore", "NextScore", "Scan_#", "File", sep = "\t") #nolint
    cat(header_str, file = log_file, append = TRUE, sep = "\n")
}

prtstr <- paste(peptide,
                "HyperScore",
                pep1_psms,
                pep1_samples,
                aa1$precursorMZ,
                aa1$precursorIntensity,
                df_pep1$Retention,
                df_pep1$Charge,
                df_pep1$Hyperscore,
                df_pep1$Nextscore,
                scan_num_pep1,
                df_pep1$mzML,
                sep = "\t")
cat(prtstr, file = log_file, append = TRUE, sep = "\n")

# Fragment information
frag_file <- "05_Fragmentation_Info.tsv"
if (file.exists(frag_file) == FALSE) {
    file.create(file.path(frag_file), showWarnings = FALSE)
    header_str <- paste("mz", "intensity", "ion", "type", "pos", "z", "seq", "error", "peptide", sep = "\t") #nolint
    cat(header_str, file = frag_file, append = TRUE, sep = "\n")
}

pep1_spectrum <- new("Spectrum2", intensity = exp_peaks1[, 2],
                    mz = exp_peaks1[, 1], centroided = TRUE,
                    precursorMz = prec_mass_pep1,
                    precScanNum = scan_num_pep1)


fragments1 <- calculateFragments(peptide, pep1_spectrum,
                type = c("b", "y"),
                tolerance = 0.05,
                method = "closest",
                modifications = c(C = 57.02146)
                )

fragments1 <- fragments1[!(grepl("_", fragments1$ion) | grepl("\\*", fragments1$ion)), ] #nolint
fragments1$peptide <- paste("Highest Hyperscore", peptide, sep = " ")

frag1_spectrum <- new("Spectrum2", intensity = fragments1$intensity,
                    mz = fragments1$mz, centroided = TRUE,
                    precursorMz = prec_mass_pep1,
                    precScanNum = scan_num_pep1)

write.table(fragments1, file = frag_file,
            append = TRUE, quote = FALSE,
            row.names = FALSE, sep = "\t",
            col.names = FALSE)

#### Create plot
png(file = paste0("05_", peptide, ".png"), # nolint
    width = 800, height = 600, units = "px")


# Plot peaks
plot(mz(pep1_spectrum), intensity(pep1_spectrum),
    xlab = "m/z", ylab = "intensity",
    type = "h", col = "gray",
    main = peptide
    )

#Plot fragments
points(mz(frag1_spectrum), intensity(frag1_spectrum),
    type = "h", col = "black"
    )

text(fragments1$mz, fragments1$intensity,
    labels = fragments1$ion
    )

# Add text
mtext(paste("Precursor mass", prec_mass_pep1,
            sep = " "), side = 3, line = 1, adj = 0)

garbage <- dev.off()

#### Second hyperscore
infile2 <- peptide_info[peptide_info$X.Peptide == peptide, "SecondSpectrumFile"]

if (is.na(infile2)){stop("No second hyperscore")}

dir2 <- str_split_fixed(infile2, "/", n = Inf)
dir2 <- paste(dir2[1 : length(dir2) - 1], collapse = "/")

pep2_samples <- peptide_info[peptide_info$X.Peptide == peptide, "MSSamples"]

df2 <- read.csv(infile2, header = TRUE, sep = "\t")
df2$mzML <- paste0(dir, "/",
                  sapply(df2$Spectrum, function(x) {unlist(strsplit(x, "\\."))[1]}), ".mzML", #nolint
                  sep = "")

df2$scan <- as.numeric(sapply(df2$Spectrum,
                              function(x) {unlist(strsplit(x, "\\."))[2]}))


df_pep2 <- df2[df2$Peptide == peptide, ]
pep2_psms <- length(df_pep2$Spectrum)

df_pep2 <- df_pep2[order(df_pep2$Hyperscore, decreasing = TRUE), ]
df_pep2 <- df_pep2[which(df_pep2$Hyperscore == max(df_pep2$Hyperscore)), ]

scan_num_pep2 <- as.integer(df_pep2$scan)
prec_mass_pep2 <- as.numeric(df_pep2$Observed.M.Z)
prec_charge_pep2 <- as.integer(df_pep2$Charge)
ppm <- formatC(as.numeric(df_pep2$ppm_error), digits = 3, format = "f")

rawdata2 <- openMSfile(df_pep2$mzML)

aa2 <- header(rawdata2, scan = scan_num_pep1)
rawdata2 <- openMSfile(df_pep2$mzML)
exp_peaks2 <- peaks(rawdata2, scan = scan_num_pep2)

############# Print out table

prtstr <- paste(peptide,
                "Second Hyperscore",
                pep2_psms,
                pep2_samples,
                aa2$precursorMZ,
                aa2$precursorIntensity,
                df_pep2$Retention,
                df_pep2$Hyperscore,
                df_pep1$Nextscore,
                df_pep2$Charge,
                df_pep2$Calculated.Peptide.Mass,
                scan_num_pep2,
                df_pep2$mzML,
                sep = "\t") #nolint
cat(prtstr, file = log_file, append = TRUE, sep = "\n")

#### 
#
pep2_spectrum <- new("Spectrum2", intensity = exp_peaks2[, 2],
                    mz = exp_peaks2[, 1], centroided = TRUE,
                    precursorMz = prec_mass_pep2,
                    precScanNum = scan_num_pep2)

#
fragments2 <- calculateFragments(peptide, pep2_spectrum,
                type = c("b", "y"))
fragments2 <- fragments2[!(grepl("_", fragments2$ion) | grepl("\\*", fragments2$ion)), ] #nolint
fragments2$peptide <- paste("Second Hyperscore", peptide, sep = " ")

frag2_spectrum <- new("Spectrum2", intensity = fragments2$intensity,
                    mz = fragments2$mz, centroided = TRUE,
                    precursorMz = prec_mass_pep2,
                    precScanNum = scan_num_pep2)

write.table(fragments2, file = frag_file,
            append = TRUE, quote = FALSE,
            row.names = FALSE, sep = "\t",
            col.names = FALSE)

#### Create plot
png(file = paste0("05_", peptide, "HighestVsSecondHyperscore.png"), # nolint
    width = 800, height = 600, units = "px")

## Plot Peaks


# Plot peaks

plot(pep1_spectrum, pep2_spectrum,
    tolerance = 0.05, relative = FALSE,
    modifications = c(C = 57.021),
    sequences = c(peptide, peptide),
    z = seq(1, 2, by = 1),
    neutralLoss = NULL,
    peaks.cex = 0,
    peaks.lwd = 1
    )

#
garbage <- dev.off()
