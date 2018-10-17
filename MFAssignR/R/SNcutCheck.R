#'Estimates the Signal to Noise cut
#'
#'For raw mass spectral data using the method
#'developed by Zhurov et al. Anal. Chem. (2014).
#'It uses the natural log of the intensity to create a histogram which then
#'is used to determine the S/N cut for the data being used.
#'The data being put into this function should be the raw mass list output from
#'a mass spectrum, if the noise peaks have been removed already, it will not work.
#'
#'The value reported in the console will be the estimated S/N cut and the plot
#'will demonstrate where the cut is being applied to the data. The value in the console
#'can then be multiplied by whatever value is desired in order to reach the value
#'to be used to cut the data.
#'
#'@param df - dataframe of intensity and ion mass, column 1 should be intensity, column 2 should be mass
#'@param SN - numeric value for situations where a predefined S/N value is desired, default is 0
#'@param bin - numeric value determining the binwidth of the histogram, default is 0.01
#'
#' @return S/N and log intensity histogram
#'   S/N printed in Console
#'   Histogram shows where the cut is being applied
#'
#' @examples
#' SNcutCheck(df, SN = 0, bin = 0.01)
#' SNcutCheck(df)
#'
#' @export



SNcutCheck <- function(df, SN = 0, bin = 0.01){

  df <- df[c(2,1)]

  names(df)[2] <- "mass"
  names(df)[1] <- "intensity"

Histo <- ggplot2::ggplot(df, ggplot2::aes(x = log(intensity))) + ggplot2::geom_histogram(binwidth = bin)
Histo
Data <- ggplot2::ggplot_build(Histo)$data  #Extracts the data from the histogram
Count <- Data[[1]]$count                   #Extracts the count data from the histogram
LogInt <- Data[[1]]$x                      #Extracts the x axis data from the histogram
Freqdf <- data.frame(Count, LogInt)
mode <- Freqdf[Count == max(Count),]
mode <- (mode[,2])
Valley <- Freqdf[(LogInt > mode & LogInt < (mode+2)),]  #Finds the valley between the noise max and the analyte
ValMin <- Valley[Valley$Count == min(Valley$Count),]    #Finds the minimum point of the valley
ValMin$Ave <- mean(ValMin$LogInt)                       #Average intensity of the minimum.
Int <- ValMin[,3]
Int <- mean(Int)
OutInt <- exp(Int)

if(SN == 0){
  df$Index <- "Bad"
  df$Index <- replace(df$Index, df$intensity > OutInt, "Good")
  print(OutInt)
}

if(SN != 0){
  df$Index <- "Bad"
  df$Index <- replace(df$Index, df$intensity > exp(SN), "Good")
  print(exp(SN))
}

Freq2 <- ggplot2::ggplot(df, ggplot2::aes(x = log(intensity))) +
  ggplot2::geom_histogram(ggplot2::aes(y = ..count.., color=Index, fill = Index), binwidth = bin) +
  ggplot2::labs(x = "log(Intensity)")
Freq2
}
##########################################
#'Plot of Mass Spectrum with highlighted S/N cut
#'
#'SNplot which plots the mass spectrum with the S/N cut denoted by different colors
#'for the mass spectrum peaks. This is useful for a qualitative look at the effectiveness
#'of the S/N cut being used.
#'
#'@param df - dataframe of intensity and ion mass, column 1 should be intensity, column 2 should be mass
#'@param cut - numeric value of the intensity cut value being investigated
#'@param mass - numeric value setting a centerpoint to look at the mass spectrum
#'@param window.x - numeric value setting the +/- range around the mass centerpoint, default is 0.5
#'@param window.y - numeric value setting the y axis value for the plot, determined by multiplying the cut by this value
#'
#' @return S/N cut colored mass spectrum
#'
#' Zoomed Mass spectrum which shows where the cut is being applied
#'
#' @examples
#' SNplot(df, cut = 1000, mass = 300, window.x = 1, window.y = 10)
#'
#' @export



SNplot <- function(df, cut, mass, window.x = 0.5, window.y = 10){ #plots a data set displaying the SN cut around a specific mass

  df <- df[c(2,1)]
  names(df)[2] <- "mass"
  names(df)[1] <- "Abundance"
  df$Index <- "Bad"
  df$Index <- replace(df$Index, df$Abundance > cut, "Good")
 SNplot <- ggplot2::ggplot(df, ggplot2::aes_string(x = "mass", xend = "mass", y = 0, yend = "Abundance")) +
    ggplot2::geom_segment(ggplot2::aes(color = Index), size = 0.65, alpha = 1) +
    ggplot2::geom_hline(yintercept = cut, linetype = "solid", size = 0.1) +
    ggplot2::coord_cartesian(xlim = c(mass-window.x, mass+window.x), ylim = c(0, cut*window.y))
  print(SNplot)
}




