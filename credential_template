# Step 0 R-package Installation 
# one-time installation
source("https://bioconductor.org/biocLite.R") # XCMS
biocLite("xcms")
install.packages("data.table") # data.table
install.packages("devtools") # devtools
devtools::install_github("pattilab/credential", force = T, ref="release/3.0") # credentialing 3.0 version
# one-time installation end


# everytime start a new R session. Load R packages
library(xcms)
library(data.table)
library(credential)

# Step 1 Set the working directory where you put your data.

  setwd("/Users/Lingjue/Desktop/CredentialingDemo")   # for Mac OS, change to your own directory where you put the mzXML files.
  # setwd("C:\Users\Mike\Desktop\CredentialingDemo")  # for Windows
  
  # There should have two folders under this directory, for example '1T1' and '1T2' representing two diffrent 12C/13C ratio combinations.

# Step 2: XCMS processing for the mzXML files in '1T1' and '1T2' separately
  # 1T1
  # 1. peak picking- by default "centWave" is used. For advanced XCMS users parameters can be adjusted as individual understanding of the
  xs_1 = xcmsSet("./1T1", method="centWave", ppm=15, peakwidth=c(10,30), snthresh=6, prefilter=c(3,100), mzdiff=0.01, integrate=1)
  # 2. feature 1st round grouping 
  xs_1= group(xs_1, bw=5, mzwid=.015, minfrac=0.5)
  # 3. retention time alignment
  xs_1 = retcor(xs_1, method="obiwarp",profStep=1)
  # 4. feature 2nd round grouping
  xs_1 = group(xs_1, bw=5, mzwid=.015, minfrac=0.5)
  # 5. fillPeaks
  xs_1 = fillPeaks(xs_1)

  # 1T2
  xs_2 = xcmsSet("./1T1", method="centWave", ppm=15, peakwidth=c(10,30), snthresh=6, prefilter=c(3,100), mzdiff=0.01, integrate=1)
  xs_2 = group(xs_2, bw=5, mzwid=.015, minfrac=0.5)
  xs_2 = retcor(xs_2, method="obiwarp",profStep=1)
  xs_2 = group(xs_2, bw=5, mzwid=.015, minfrac=0.5)
  xs_2 = fillPeaks(xs_2)
  
  
# Step 3: Construct input format for credentialing
  
  #1T1
  # format conversion
  dt_1=data.table(xs_1@groups)
  # find the peaks belonging to the same feature. Average intensities of these peaks as the feature intensity
  for(i in seq.int(dt_1[, .N])){
    idx_tmp = xs_1@groupidx[[i]]
    peaks_tmp = xs_1@peaks[idx_tmp,]
    dt_1[i,"intm":= mean(peaks_tmp[,"into"])]
  }
  
  #1T2
  dt_2=data.table(xs_2@groups)
  for(i in seq.int(dt_2[, .N])){
    idx_tmp = xs_2@groupidx[[i]]
    peaks_tmp = xs_2@peaks[idx_tmp,]
    dt_2[i,"intm":= mean(peaks_tmp[,"into"])]
  }
  
  # feature table construct
  ft_1 = dt_1[,c("mzmed","rtmed","intm")]
  ft_1[,"cc":=seq.int(ft_1[ ,.N])]
  ft_1 = ft_1[,c("cc","mzmed","rtmed","intm")]
  colnames(ft_1) <- c("cc","mz", "rt", "i")
  
  ft_2 = dt_2[,c("mzmed","rtmed","intm")]
  ft_2[,"cc":=seq.int(ft_1[ ,.N])]
  ft_2 = ft_2[,c("cc","mzmed","rtmed","intm")]
  colnames(ft_2) <- c("cc","mz", "rt", "i")
  
# Step 4: Feature credentialing in each group
  
  knots_f1 = findknots(ft_1, .zs = 1:2, ppmwid = 15, rtwid = 1, cd = 13.00335-12)
  credentials_f1 = credentialknots(knots_f1$knot, ppmwid = 15, rtwid = 1, mpc = c(12, 120), ratio = 1/1, ratio.lim = 0.1, maxnmer = 4, cd = 13.00335-12)
  ft_1 = ft_1[knots_f1$cc_knot[credentials_f1$knot_quipu[!is.na(quipu)],on="knot"],,on="cc"] 
    
  
  knots_f2 = findknots(ft_2, .zs = 1:2, ppmwid = 15, rtwid = 1, cd = 13.00335-12)
  credentials_f2 = credentialknots(knots_f2$knot, ppmwid = 15, rtwid = 1, mpc = c(12, 120), ratio = 1/2, ratio.lim = 0.1, maxnmer = 4, cd = 13.00335-12)
  ft_2 = ft_2[knots_f2$cc_knot[credentials_f2$knot_quipu[!is.na(quipu)],on="knot"],,on="cc"] 
  
  
# Step 5: Find common credentialed features between two conditions
  
  # rearrange pre-match
  dt1T1 <- ft_1[credentials_f1$quipu[,c("quipu","ratio")][!is.na(quipu)],,on="quipu"]
  dt1T2 <- ft_2[credentials_f2$quipu[,c("quipu","ratio")][!is.na(quipu)],,on="quipu"]
  
  dt1T1 <- dt1T1[order(dt1T1[,"quipu"], dt1T1[,"mz"]),]
  dt1T2 <- dt1T2[order(dt1T2[,"quipu"], dt1T2[,"mz"]),]
  
  
  match_cf = matchcredfeature(dt1T1,dt1T2,ppm=15,drt=1)
  
  
