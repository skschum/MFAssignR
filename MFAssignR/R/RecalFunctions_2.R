###############################################################
#' Generates a plot to visualize recalibrant series, recalibrates two mass lists,
#' and produces a list of the chosen recalibrants.
#'
#'
#' Recal() takes the data frame output from MFAssign functions, the outputs of
#' IsoFiltR(), and the chosen recalibrant series to generate a plot for the qualitative
#' assessment of the recalibrants, recalibrate the "Mono" and "Iso" outputs from IsoFiltR,
#' and prepare a data frame containing the chosen recalibrants.
#'
#'
#' This function can handle up to 10 homologous series, though it will work with fewer.
#' It is important for recalibrant masses to cover the entire mass range of interest,
#' and they should be among the most abundant peaks in their region of the spectrum,
#' this function helps with visualizing the recalibrants.
#'
#' This function also recalibrates up to two mass lists using the chosen recalibrants.
#' It is best to use the "Mono" and "Iso" outputs of the IsoFiltR() function.
#'
#' This recalibration method uses the first step from the paper from Kozhinov et al.
#' 2014 (Anal. Chem.) to estimate the correction term, and uses the idea of a segmented
#' or "walking" recalibration from Savory et al. 2011 (Anal. Chem) to further improve
#' the recalibration.
#'
#'
#' These parameters will be used to identify more potential recalibrants, and the
#' tallest peaks within a user defined mass range will be used as recalibrants.
#'
#' @param df data frame output from MFAssign() or MFAssignCHO()
#'
#' @param peaks data frame with two columns for recalibration, generally contains "Mono" output from IsoFiltR()
#'
#' @param isopeaks data frame with two columns, generally containing the "Iso" output from IsoFiltR()
#'
#' @param mode character "neg" or "pos" depending on the ionization mode of the data
#'
#' @param SN numeric value of the signal-to-noise cut point for making the plot, default is 0
#'
#' @param mzRange numeric value defining the mass windows used for the segmented recalibration, default is 50
#'
#' @param series(1-10) character the recalibrant series, "O7_Na_4" for example, default is NA
#'
#' @param min numeric - minimum mass range of the data being analyzed.
#'
#' @param max numeric - maximum mass range of the data being analyzed.
#'
#' @param bin numeric - sets the mass window bins to choose recalibrants from, default is 10.
#'
#' @param obs numeric - sets the number of peaks to choose as recalibrants from each bin,
#' default is 2
#'
#' @return list(MZ, Mono, Iso, RecalOut) contains the mass spectrum (MZ), recalibrated "Mono" data frame,
#'  recalibrated "Iso" data frame, and the recalibrants list (RecalOut)
#'
#' @examples
#' Recal_2(df = Data, peaks = Mono, isopeaks = Iso, mode = "neg", SN = 500, series1 = "O4_H_2", series2 = "O4_H_8", series3 = "O6_H_8")
#'
#' Recal_2(df = Data, peaks = Mono, isopeaks = Iso, mode = "pos", SN = 300, series1 = "O4_Na_2", series2 = "O4_H_8", series3 = "O6_Na_8")
#' @export


Recal_2 <- function( df, peaks, isopeaks = "none",mode, SN = 0, mzRange = 50, series1=NA, series2=NA, series3=NA, series4=NA, series5=NA,
                          series6=NA, series7=NA, series8=NA, series9=NA, series10=NA, min = 100, max = 1000,
                          bin = 10, obs = 2){
  #Preparation of the recalibrants list
  RecalList <- data.frame(series = c(series1, series2, series3, series4, series5, series6, series7, series8, series9,
                                     series10),
                          Tag = c("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10"))
  RecalList <- RecalList[!is.na(RecalList$series),]

  df$Adduct <- "H"
  df$Adduct <- replace(df$Adduct, df$M > 0, "Na")
  df$Adduct <- replace(df$Adduct, df$POE == 1, "OE")
  df$Adduct <- replace(df$Adduct, df$NOE == 1, "OE")
  df$series <- paste(df$class, df$Adduct, df$DBE, sep = "_")
  df$mode <- mode
  names(df)[1] <- "Abundance"
  names(df)[2] <- "Exp_mass"


  isopeaks <- if(isopeaks == "None") data.frame(mass = 1, Abundance = 1, Tag = "X") else isopeaks
  names(isopeaks)[1] <- "mass"
  names(isopeaks)[2] <- "Abundance"
  names(isopeaks)[3] <- "Tag"
  names(peaks)[1] <- "mass"
  names(peaks)[2] <- "Abundance"
  peaks$Tag <- "X"
  peaks$Tag2 <- "Mono"

  isopeaks$Tag2 <- "Iso"

  peaks <- if(isopeaks != "None") rbind(peaks, isopeaks) else peaks

  peaks <- peaks[c(2,1,3,4)]
  #isodummy <- data.frame
  isopeaks <- isopeaks[c(2,1,3,4)]

  #Merges recalibrant list to df in order to determine which recalibrants are in the data frame.
  RecalList2 <- merge(df, RecalList, by.x = "series", by.y = "series")


  NewRecal <- RecalList2[c("Exp_mass", "formula", "theor_mass", "C13_mass")]

  FinalRecal2 <- NewRecal

  NewRecal_mono <- FinalRecal2[c(1:3)]
  NewRecal_iso <- FinalRecal2[c(4,2,3)]
  NewRecal_iso <- NewRecal_iso[NewRecal_iso$C13_mass > 0,]
  names(NewRecal_iso)[1] <- "Exp_mass"
  NewRecal_iso$theor_mass <- NewRecal_iso$theor_mass + 1.003355
  NewRecal <- rbind(NewRecal_mono, NewRecal_iso)

  ##########################
  #Setting up the recalibrant and mass lists so they can be recalibrated by section.
  #mzRange <- 50
  peaks <- peaks[!is.na(peaks$mass),]
  Dummy <- data.frame(Values = seq(from = floor(min(peaks$mass)), to = ceiling(max(peaks$mass)), by = mzRange))
  Dummy2 <- nrow(Dummy)
  NewRecal$Range <- 1

  for(i in 1:Dummy2){
    I <- i

    NewRecal$Range <- replace(NewRecal$Range, NewRecal$Exp_mass >= min(NewRecal$Exp_mass) + (I-1) * mzRange &
                                NewRecal$Exp_mass < min(NewRecal$Exp_mass) + I*mzRange, I)

  }

  peaks$Range <- 1
  for(i in 1:Dummy2){
    I <- i

    peaks$Range <- replace(peaks$Range, peaks$mass >= min(peaks$mass) + (I-1) * mzRange &
                                peaks$mass < min(peaks$mass) + I*mzRange, I)

  }

  NewRecal$Test <- 1
  for(i in 1:Dummy2){
    I <- i

    NewRecal$Test <- replace(NewRecal$Test, NewRecal$Range == I, nrow(NewRecal[NewRecal$Range == as.numeric(I),]))


  }

  NewRecal <- NewRecal[NewRecal$Test > 2,]
  NewRecal <- NewRecal[c(1:4)]

  #This is just setting up a dataframe for final export
  Align <- NewRecal
  names(Align)[1] <- "Exp_mass"
  FinalRecal2 <- merge(FinalRecal2, Align, by.x = "Exp_mass", by.y = "Exp_mass")
  FinalRecal2 <- FinalRecal2[c(1:5)]
  names(FinalRecal2)[3] <- "formula"

  #########################
  #Calculating the weights needed for recalibration and calculating the recalibrated masses.
  names(NewRecal)[1] <- "E_mass"
  names(NewRecal)[2] <- "formula"
  names(NewRecal)[3] <- "Th_mass"
  peaks_dum <- data.frame(Abundance = -2, mass = -2, Tag = "X", Tag2 = "Y", Range = -1)
  peaks_out <- data.frame(Abundance = -2, mass = -2, Tag = "X", Tag2 = "Y", Range = -1)
  Ej_dum <- data.frame(Ejsum = -2, Range = -2)
  Ej_out <- data.frame(Ejsum = -2, Range = -2)
  Dummy3 <- data.frame(Value = 1:max(NewRecal$Range))
  Dummy3 <- as.numeric(nrow(Dummy3))

  peakshigh <- peaks[peaks$Range > max(NewRecal$Range),]
  peaks <- peaks[peaks$Range <= max(NewRecal$Range),]
  #j <- 13
  for(j in 1:(Dummy3)){
    J <- j
    while(J <= Dummy3){

      NewRecal_2 <- NewRecal[NewRecal$Range == J,]
      peaks_2 <- peaks[peaks$Range == J,]
      peaks_final <- rbind(peaks_dum, peaks_out)
      Ej_final <- rbind(Ej_dum, Ej_out)
      counter <- 1
      J <- J + Dummy3
      while(counter == 1){
        counter <- counter + 1
        NewRecal_2$num <- 0:(nrow(NewRecal_2)-1)
        NewRecal_2$weight <- factorial(nrow(NewRecal_2)-1)/
          (factorial(NewRecal_2$num)*factorial((nrow(NewRecal_2)-1)-NewRecal_2$num))

        NewRecal_2$mzweight <- NewRecal_2$weight*NewRecal_2$E_mass
        NewRecal_2$Ejweight <- NewRecal_2$weight*(NewRecal_2$E_mass-NewRecal_2$Th_mass)/NewRecal_2$Th_mass
        NewRecal_2$Ejsum <- sum(NewRecal_2$Ejweight)/2^(nrow(NewRecal_2)-1)
        NewRecal_2$masssum <- sum(NewRecal_2$mzweight)/2^(nrow(NewRecal_2)-1)
        Ejsum <- mean(NewRecal_2$Ejsum)
        peaks_out <- peaks_2
        peaks_out$mass <- peaks_out$mass/(1+Ejsum)
        peaks_out <- rbind(peaks_final, peaks_out)
        Ej_out <- data.frame(Ejsum = Ejsum, Range = mean(NewRecal_2$Range))
        Ej_out <- rbind(Ej_final, Ej_out)
      }

    }
    peaks_final <- peaks_out
    Ej_final <- Ej_out
  }


  peaks <- peaks_final[peaks_final$mass > 0,]
  Ej_final <- Ej_final[Ej_final$Range > 0,]

  Ej <- Ej_final[Ej_final$Range == max(Ej_final$Range),]

  peakshigh$mass <- peakshigh$mass/(1+Ej$Ejsum)

  peaks <- rbind(peaks, peakshigh)
  #Separating the isotope and monoisotope masses for output.
  isopeaks2 <- peaks[peaks$Tag2 != "Mono",]
  isopeaks2 <- isopeaks2[c(1,2,3)]
  names(isopeaks2)[2] <- "Iso_mass"
  names(isopeaks2)[1] <- "Iso_Abund"
  peaks <- peaks[peaks$Tag2 == "Mono",]
  peaks <- peaks[c(1,2)]

  ################################################
  #Recalibration for error comparison, calculating the average error before and after recalibration.
  FinalRecal3 <- merge(df, FinalRecal2, by.x = "Exp_mass", by.y = "Exp_mass", all = T)

  FinalRecal3$Range <- 1
  for(i in 1:Dummy2){
    I <- i

    FinalRecal3$Range <- replace(FinalRecal3$Range, FinalRecal3$Exp_mass >= min(FinalRecal3$Exp_mass) + (I-1) * mzRange &
                                   FinalRecal3$Exp_mass < min(FinalRecal3$Exp_mass) + I*mzRange, I)

  }

  FinalRecal3$Range <- replace(FinalRecal3$Range, FinalRecal3$Range > Dummy3, Dummy3)

  FinalRecal3 <- merge(FinalRecal3, Ej_final, by.x = "Range", by.y = "Range")
  FinalRecal3$R_mass <- FinalRecal3$Exp_mass/(1+FinalRecal3$Ejsum)
  FinalRecal3$R_err <- abs((FinalRecal3$R_mass-FinalRecal3$theor_mass)/FinalRecal3$theor_mass * 10^6)

  print(mean(FinalRecal3$AE_ppm))
  print(mean(FinalRecal3$R_err))

  FinalRecal3 <- FinalRecal3[!is.na(FinalRecal3$Bin),]
  plotpeak <- peaks[peaks$Abundance > SN,]

  RecalOut <- merge(Align, Ej_final, by.x = "Range", by.y = "Range")
  RecalOut$R_mass <- RecalOut$Exp_mass/(1+RecalOut$Ejsum)
  RecalOut$Recal_err <- abs((RecalOut$R_mass-RecalOut$theor_mass)/RecalOut$theor_mass * 10^6)
  RecalOut$Orig_err <- abs((RecalOut$Exp_mass-RecalOut$theor_mass)/RecalOut$theor_mass * 10^6)

  AbundM <- df[c("Exp_mass", "Abundance")]
  AbundI <- df[c("C13_mass", "C13_abund")]
  names(AbundI)[1] <- "Exp_mass"
  names(AbundI)[2] <- "Abundance"
  Abund <- rbind(AbundM, AbundI)
  Abund <- Abund[Abund$Abundance > 0,]

  RecalPlot <- merge(RecalOut, Abund, by.x = "Exp_mass", by.y = "Exp_mass")

  RecalOut <- RecalOut[c(2, 3, 4, 7, 8)]

  #Plot highlighting the recalibrant ions for qualitative assessment.
  MZ<-ggplot2::ggplot() + ggplot2::geom_segment(data=plotpeak, size=0.7,ggplot2::aes_string(x = "mass", xend = "mass", y = 0, yend = "Abundance"), alpha = 0.3, color = "grey57")+
    ggplot2::geom_segment(data=RecalPlot, size=1.2,ggplot2::aes_string(x = "Exp_mass", xend = "Exp_mass", y = 0, yend = "Abundance"), color = "blue")+

    ggplot2::theme_bw()+ggplot2::labs(x = "Ion Mass", y = "Abundance", title = "Assignment Mass Spectrum", color = "Series")+
    ggplot2::theme(axis.title=ggplot2::element_text(size = 15, face = "bold"), strip.text=ggplot2::element_text(size=15,face="bold"),
                   axis.text=ggplot2::element_text(size=15, face = "bold"), legend.title=ggplot2::element_text(face="bold", size = 15),
                   legend.text=ggplot2::element_text(face="bold", size = 15),  panel.grid.minor.x=ggplot2::element_blank(),
                   panel.grid.major.x=ggplot2::element_blank(), strip.background = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(size = 16, face = "bold"))

  #Preparing the final output of the function, the recalibrants list.

  # RecalOut <- merge(Align, Ej_final, by.x = "Range", by.y = "Range")
  # RecalOut$R_mass <- RecalOut$Exp_mass/(1+RecalOut$Ejsum)
  # RecalOut$Recal_err <- abs((RecalOut$R_mass-RecalOut$theor_mass)/RecalOut$theor_mass * 10^6)
  # RecalOut$Orig_err <- abs((RecalOut$Exp_mass-RecalOut$theor_mass)/RecalOut$theor_mass * 10^6)
  #
  # RecalOut <- RecalOut[c(2, 3, 4, 7, 8)]



  peaks <- peaks[c(2,1)]
  names(peaks)[2] <- "abundance"
  isopeaks2 <- isopeaks2[c(2,1,3)]
  names(isopeaks2)[2] <- "iso_abund"
  names(isopeaks2)[1] <- "iso_mass"
  names(isopeaks2)[3] <- "tag"
  names(RecalOut)[1] <- "exp_mass"

  Output <- list(Plot = MZ, Mono = peaks, Iso = isopeaks2, RecalList = RecalOut)
  Output
}
