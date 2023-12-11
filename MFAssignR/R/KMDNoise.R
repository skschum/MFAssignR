#' Estimates the Signal to Noise cut
#'
#' For raw mass spectral data using a method
#' based on the separation of Kendrick mass
#' defect of intense analyte peaks and low
#' intensity noise peaks.
#'
#' When the KMD is calculated for all peaks
#' in a raw mass spectrum, the more intense
#' analyte peaks form "islands" surrounded
#' by a sea of low intensity noise peaks.
#' using the mathematical limits of the KMD
#' calculation the slope of a KMD plot can be
#' calculated and used in conjunction with the
#' KMD limits of chemically feasible molecular
#' formulas to isolate regions of the plot that
#' contain only noise peaks. These noise peaks can
#' then be averaged to determine the signal to noise
#' of the spectrum.
#'
#' The upper chemical limit of the KMD plot can be
#' estimated by calculating the KMD for non-hydrogenated
#' condensed carbon formulas (C10, C20, C30, etc.).
#' This provides an x,y intercept of (0,0) when the slope
#' (0.001232) of the KMD plot is used. This provides a
#' reasonable lower bound for isolating noise peaks.
#' The upper bound can be selected by increasing the
#' intercept value in the equation for the KMD plot,
#' which is y = 0.0011232 * x + b.
#'
#' This signal to noise estimation method uses the
#' equation described above to select a region of the
#' KMD plot that contains only noise and then estimates
#' the average intensity of the peaks in that region,
#' reporting this value as the noise level.
#'
#' An additional output of the function is the KMD plot
#' with the bounds of the noise estimation area highlighted
#' in red. This noise level can then be multiplied by the
#' users chosen value (3, 6, 10) in order to set the
#' signal to noise cut for formula assignment.
#'
#' The y intercept for the upper and lower bounds are set
#' to 0.05 (lower.y) and 0.2 (upper.y). The lower bound is 0.05
#' instead of 0 as mentioned previously in order to ensure
#' no analyte peaks are incorporated into the noise estimation.
#' The upper limit is set at 0.2 so that it does not
#' interact with any potentially doubly charged peaks, or the
#' "echo" of those doubly charged peaks. Both of these values
#' can be changed if desired by the user.
#'
#' The x intercept for the upper and lower bounds are set
#' to NA, which means if the are not changed, it will default
#' to the maximum and minimum mass in the mass spectrum. If the
#' user wants to select a specific region then they can put
#' whatever limit they want on it.


#'
#' @param df - dataframe of intensity and ion mass, column 1 should be mass, column 2 should be intensity
#' @param upper.y - sets the upper limit for the y intercept, default is 0.2
#' @param lower.y - sets the lower limit for the y intercept, default is 0.05
#' @param upper.x - sets the upper limit for the x intercept, default is NA
#' @param lower.x - sets the lower limit for the x intercept, default is NA
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
# df <- Raw_Neg_ML
KMDNoise <- function(df, upper.y = 0.2, lower.y = 0.05, upper.x = NA, lower.x = NA) {
  if (is.na(upper.x)) upper.x <- max(df$mass)
  if (is.na(lower.x)) lower.x <- min(df$mass)

  # Init the dataframe with mass defect information
  # magic number is mass of CH2 which we use for the correction
  df$kendrick_mass <- df$mass * (14 / 14.01565)
  df$nominal_mass <- round(df$mass) # maybe this should actually be floor instead of round as the nominal mass is always rounded down?
  df$kendrick_mass_defect <- df$nominal_mass - df$kendrick_mass
  df$log_int <- log(df$intensity)

  df <- df[order(df$intensity), ]

  magic_number <- 0.0011232

  # Filtering the actual data within the S/N area and calculating the average intensity
  lower_bound <- magic_number * df$mass + lower.y
  upper_bound <- magic_number * df$mass + upper.y

  SN <- dplyr::filter(df, dplyr::between(kendrick_mass_defect, lower_bound, upper_bound))
  SN <- dplyr::filter(SN, dplyr::between(mass, lower.x, upper.x))

  Noise <- exp(mean(SN$log_int))

  df <- dplyr::filter(df, dplyr::between(log_int, -100, 100))

  # Generate plot showing KMD of raw data and the boundaries for S/N cut.
  KMD <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = df,
      ggplot2::aes_string(
        x = "mass",
        y = "kendrick_mass_defect",
        color = "log_int"
      ), alpha = 1 / 3
    ) +
    ggplot2::geom_abline(slope = magic_number, intercept = lower.y, color = "red", size = 1) +
    ggplot2::geom_abline(slope = magic_number, intercept = upper.y, color = "red", size = 1) +
    ggplot2::geom_vline(xintercept = lower.x, color = "red", size = 1) +
    ggplot2::geom_vline(xintercept = upper.x, color = "red", size = 1) +
    ggplot2::scale_color_gradientn(colors = colorRamps::blue2red(50), limit = c(min(df$log_int), max(df$log_int))) +
    ggplot2::labs(
      x = "Ion Mass", y = "Kendrick Mass Defect",
      title = "KMD Signal to Noise Determination Plot", color = "ln(int)"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = 18, face = "bold"),
      strip.text = ggplot2::element_text(size = 18, face = "bold"),
      axis.text = ggplot2::element_text(size = 17, face = "bold"),
      legend.title = ggplot2::element_text(face = "bold", size = 15),
      legend.text = ggplot2::element_text(face = "bold", size = 14),
      strip.background = ggplot2::element_blank()
    )
  return(list(Noise = Noise, KMD = KMD))
}

##########################################
#' Plot of Mass Spectrum with highlighted S/N cut
#'
#' SNplot which plots the mass spectrum with the S/N cut denoted by different colors
#' for the mass spectrum peaks. This is useful for a qualitative look at the effectiveness
#' of the S/N cut being used.
#'
#' @param df - dataframe of intensity and ion mass, column 1 should be intensity, column 2 should be mass
#' @param cut - numeric value of the intensity cut value being investigated
#' @param mass - numeric value setting a centerpoint to look at the mass spectrum
#' @param window.x - numeric value setting the +/- range around the mass centerpoint, default is 0.5
#' @param window.y - numeric value setting the y axis value for the plot, determined by multiplying the cut by this value
#'
#' @return S/N cut colored mass spectrum
#'
#' Zoomed Mass spectrum which shows where the cut is being applied
#'
#' @examples
#' SNplot(df, cut = 1000, mass = 300, window.x = 1, window.y = 10)
#'
#' @export



SNplot <- function(df, cut, mass, window.x = 0.5, window.y = 10) { # plots a data set displaying the SN cut around a specific mass

  df <- df[c(2, 1)]
  names(df)[2] <- "mass"
  names(df)[1] <- "Abundance"
  df$Index <- "Bad"
  df$Index <- replace(df$Index, df$Abundance > cut, "Good")
  SNplot <- ggplot2::ggplot(df, ggplot2::aes_string(x = "mass", xend = "mass", y = 0, yend = "Abundance")) +
    ggplot2::geom_segment(ggplot2::aes(color = Index), size = 0.65, alpha = 1) +
    ggplot2::geom_hline(yintercept = cut, linetype = "solid", size = 0.1) +
    ggplot2::coord_cartesian(xlim = c(mass - window.x, mass + window.x), ylim = c(0, cut * window.y))
  print(SNplot)
}
