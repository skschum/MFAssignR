#'Estimates the Signal to Noise cut
#'
#'For raw mass spectral data using a method
#'based on the separation of Kendrick mass
#'defect of intense analyte peaks and low
#'intensity noise peaks.
#'
#'When the KMD is calculated for all peaks
#'in a raw mass spectrum, the more intense
#'analyte peaks form "islands" surrounded
#'by a sea of low intensity noise peaks.
#'using the mathematical limits of the KMD
#'calculation the slope of a KMD plot can be
#'calculated and used in conjunction with the
#'KMD limits of chemically feasible molecular
#'formulas to isolate regions of the plot that
#'contain only noise peaks. These noise peaks can
#'then be averaged to determine the signal to noise
#'of the spectrum.
#'
#'The upper chemical limit of the KMD plot can be
#'estimated by calculating the KMD for non-hydrogenated
#'condensed carbon formulas (C10, C20, C30, etc.).
#'This provides an x,y intercept of (0,0) when the slope
#'(0.001232) of the KMD plot is used. This provides a
#'reasonable lower bound for isolating noise peaks.
#'The upper bound can be selected by increasing the
#'intercept value in the equation for the KMD plot,
#'which is y = 0.0011232 * x + b.
#'
#'This signal to noise estimation method uses the
#'equation described above to select a region of the
#'KMD plot that contains only noise and then estimates
#'the average intensity of the peaks in that region,
#'reporting this value as the noise level.
#'
#'An additional output of the function is the KMD plot
#'with the bounds of the noise estimation area highlighted
#'in red. This noise level can then be multiplied by the
#'users chosen value (3, 6, 10) in order to set the
#'signal to noise cut for formula assignment.
#'
#'The y intercept for the upper and lower bounds are set
#'to 0.05 (lower) and 0.2 (upper). The lower bound is 0.05
#'instead of 0 as mentioned previously in order to ensure
#'no analyte peaks are incorporated into the noise estimation.
#'The upper limit is set at 0.2 so that it does not
#'interact with any potentially doubly charged peaks, or the
#'"echo" of those doubly charged peaks. Both of these values
#'can be changed if desired by the user.


#'
#'@param df - dataframe of intensity and ion mass, column 1 should be mass, column 2 should be intensity
#'@param upper - sets the upper limit for the y intercept, default is 0.2
#'@param lower - sets the lower limit for the y intercept, default is 0.05
#'
#' @return list(Noise = Noise, Plot = KMD)
#'         Noise - numeric value of the noise level for the data set
#'         Plot - plot showing the KMD for the inputted data with the
#'         noise estimation area highlighted.
#'
#' @examples
#' KMDNoise(df, upper = 0.2, lower = 0.05)
#' KMDNoise(df)
#'
#' @export
#'
 #df <- Data
KMDNoise <- function(df, upper = 0.2, lower = 0.05){

  names(df)[1] <- "mass"
  names(df)[2] <- "intensity"
  df$KM <- df$mass * (14/14.01565)
  df$NM <- round(df$mass)
  df$KMD <- df$NM - df$KM
  df$int <- log(df$intensity)
  df <- df[order(df$intensity),]
  #Setting up the boundary lines for the chosen S/N area for plotting purposes
  Limits <- data.frame(Number = seq(round(min(df$mass/12)),round(max(df$mass/12)),1))
  Limits$mass <- Limits$Number*12
  Limits$KMD_low <- 0.0011232*Limits$mass + lower
  Limits$KMD_up <- 0.0011232*Limits$mass + upper

  #Filtering the actual data within the S/N area and calculating the average intensity
  SN <- df[df$KMD > 0.0011232*df$mass + lower & df$KMD < 0.0011232*df$mass + upper,]
  #SN8 <- Data8[Data8$KMD > 0.0011232*Data8$m.z +0.03 & Data8$KMD < 0.0011232*Data8$m.z + 0.1,]
  #SN$Zmod <- 0.6745*(SN$int-mean(SN$int))/(median(abs(SN$int - median(SN$int))))
  #Noise <- SN[abs(SN$Zmod) <= 3.5,]
  Noise <- mean(SN$int)
  Noise <- exp(Noise)

  df <- df[df$int > -100 & df$int < 100,]

  #Generate plot showing KMD of raw data and the boundaries for S/N cut.
  KMD <-ggplot2::ggplot() + ggplot2::geom_point(data=df,
                                               ggplot2::aes_string(x = "mass", y = "KMD",
                                                                   color = "int" ), alpha = 1/3) +
    ggplot2::geom_abline(slope = 0.0011232, intercept = lower, color = "red", size = 1)+
    ggplot2::geom_abline(slope = 0.0011232, intercept = upper, color = "red", size = 1)+
    ggplot2::scale_color_gradientn(colors = colorRamps::blue2red(50), limit = c(min(df$int), max(df$int)))+

    ggplot2::labs(x = "Ion Mass", y = "Kendrick Mass Defect",
         title = "KMD Signal to Noise Determination Plot", color = "ln(int)") + ggplot2::theme_bw() +
    ggplot2::theme(axis.title=ggplot2::element_text(size = 18, face = "bold"),
                   strip.text=ggplot2::element_text(size=18,face="bold"),
                   axis.text=ggplot2::element_text(size=17, face = "bold"),
                   legend.title=ggplot2::element_text(face="bold", size = 15),
                   legend.text=ggplot2::element_text(face="bold", size = 14),
                   strip.background = ggplot2::element_blank())
KMD
  Output <- list(Noise = Noise, KMD = KMD)
  Output
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
