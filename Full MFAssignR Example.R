# Here are some parameters for each of the functions for testing the code
# The same general information is shared in the vignette, but I thought it may be useful to have it
#consolidated here. Additionally, there is a raw mass list that can be used to test all the functions. Its
#name is Raw_Neg_ML and it is negative mode ESI data
##########################
#Package Install Procedure
install.packages("devtools")

setwd("Whatever directory the MFAssignR folder is in.")
devtools::install(MFAssignR)
####################
library(MFAssignR)
setwd("Your data directory")
Data <- read.csv("Your data.csv")

#########################
#Signal-to-noise estimation and check
SNcutCheck(Data, bin = 0.01) #Default settings for SNcutCheck.

SNplot(Data, cut = 1800, mass = 301.1, window.x = 50, window.y = 10) #Reasonable settings for SNplot

#####################
#Isotope prescreening
Isotope <- IsoFiltR(Data)  #Input for IsoFiltR

Mono <- Isotope[["Mono"]]
Iso <- Isotope[["Iso"]]
##############################
#CHO formula pre-assign

Assign <- MFAssignCHO(Mono, Iso, ionMode = "neg", lowMW =50, highMW = 1000, ppm_err = 3, H_Cmin = 0.3,
                      HetCut = "off", NMScut = "on", SN = 1800)   #Standard parameters for negative mode

Unambig1 <- Assign[["Unambig"]]
Ambig1 <- Assign[["Ambig"]]
Unassigned1 <- Assign[["None"]]

MSAssign <- Assign[["MSAssign"]]
Error <- Assign[["Error"]]
MSgroups <- Assign[["MSgroups"]]
VK <- Assign[["VK"]]
MSAssign
Error
MSgroups
VK
##################################
#Highlighting possible recalibrant series
check <- MFRecalList(Unambig1)

##################################
#Qualitative check of recalibrant series and mass recalibration.

Test <- MFRecalCheck(df = Unambig1,peaks = Mono, isopeaks = Iso, mode = "neg", SN = 1800, series1 = "O8_H_9",
                     series2 = "O6_H_3", series3 = "O4_H_2", series4 = "O13_H_13", series5 = "O15_H_16")

Plot <- Test[["Plot"]]
Plot      #This plot is slow to generate
Mono2 <- Test[["Mono"]]
Iso2 <- Test[["Iso"]]
List <- Test[["RecalList"]]


#############################################
#Final formula assignment
#Parameters for both positive and negative mode formula assignment. This uses the formula extension version

Assign <- MFAssign(Mono2, Iso2, ionMode = "neg", lowMW =50, highMW = 1000,  Nx = 3, Sx = 1,  ppm_err = 3, H_Cmin = 0.3,
                   HetCut = "off", DeNovo = 400, NMScut = "on", SN = 1800)


Unambig2 <- Assign[["Unambig"]]
Ambig2 <- Assign[["Ambig"]]
Unassigned2 <- Assign[["None"]]

MSAssign <- Assign[["MSAssign"]]
Error <- Assign[["Error"]]
MSgroups <- Assign[["MSgroups"]]
VK <- Assign[["VK"]]
MSAssign
Error
MSgroups
VK


#############################################
#Final formula assignment
#Parameters for both positive and negative mode formula assignment. This uses the older version for isolated
#assignment

Assign <- MFAssignAll(Mono2, Iso2, ionMode = "neg", lowMW =50, highMW = 1000,  Nx = 3, Sx = 1,  ppm_err = 3, H_Cmin = 0.3,
                   HetCut = "off", NMScut = "on", SN = 1800)

Unambig3 <- Assign[["Unambig"]]
Ambig3 <- Assign[["Ambig"]]
Unassigned3 <- Assign[["None"]]
MSAssign <- Assign[["MSAssign"]]
Error <- Assign[["Error"]]
MSgroups <- Assign[["MSgroups"]]
VK <- Assign[["VK"]]
MSAssign
Error
MSgroups
VK
