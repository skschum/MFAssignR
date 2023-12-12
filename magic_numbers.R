# Separation of isotope mass for the resolution enhanced KMD
## for 34S = 12 (obtained by trial and error)
## for 13C = 21 (obtained from Zhang et al (2018))

## KMDr difference S = -0.291 (32S) or 0.709 (34S)
## KMDr difference C = -0.496 (12C) or 0.503 (13C)
## difference is calculated by subtracting the KMDr value for the suspected polyisotopic mass from the monoisotopic mass

# allowance for the measurement error
## for the 13C: -0.4975 < KMDrDiff < - 0.494501, and 0.501501 < KMDrDiff < 0.5045
## for the 34S: -0.29349 < KMDrDiff < -0.29051, and 0.7075 < KMDrDiff < 0.70949 

# Mass matching step: data is broken into 10 overlapping clusters due to computational reasons
## difference in mass between matched observations is calculated
## observation pairs that match within +-5ppm of the theoretical mass difference of the isotope of interest 
## (1.003355 Da for 13 C and 1.995797 Da for 34 S) are retained and all other pairs are removed from further 
## consideration

# KM = mass * 14.01565/14   ## Kendrick mass
# KMD = KM-NM               ## Kendrick mass defect = Kendrick mass - nominal mass

# How does IsoFiltR() function work:

# 1) Mass matching
## match every mass in the spectrum with every other mass
## break the data into 10 overlapping clusters
## calculate difference in mass between each pair
## if it's within -/+5ppm of the theoretical MD of isotope of interest -> retain
## recombine data into one dataframe and filter duplicates

# 2) Isotopic KMD series
## calculate KM+KMD for 13C and 34S for preliminary monoisotopic and polyisotopic masses
## subtract the values from each other
## if the absoute value of the subtracted number is < 0.00149 -> pair of peaks is considered to be in series

# 3) KMDr calculation
## divides a repeating mass unit by some integer -> separates the isotope mass pairs by a consistent value