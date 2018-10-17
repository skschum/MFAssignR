#' Identifies canditate series for recalibration
#'
#' MFRecalList() takes the output from MFAssign() and identifies the
#' homologous series that could be used as recalibrants.
#'
#' It returns a dataframe that contains the CH2 homologous series that
#' contain more than 3 members.
#'
#' The columns of the returned dataframe are as follows:
#'
#' Series - reports the homologous series according to class, adduct, and DBE, the
#' format is "class_adduct_DBE", for example a homologous series with class = "O6, adduct of Na+,
#' and DBE = 4 would be "O6_Na_4"
#'
#' Number Observed - reports the number of members of each homologous series.
#'
#' Series Index - represents the order of the series when ordered by length of
#' homologous series.
#'
#' Mass Range - reports the minimum and maximum mass for the compounds within a
#' homologous series.
#'
#' Tall Peak - reports the mass of the most abundant peak in each series.
#’
#' Abundance Score - reports the percentage difference between the mean abundance
#' of a homologous series and the median abundance within the mass range the
#' "Tall Peak" falls in (for example m/z 200-300). A higher score is generally better.
#'
#' Peak Score - This column compares the intensity of the tallest peak in a given series to the
#' second tallest peak in the series This comparison is calculated by log10(Max Peak Intensity/
#' Second Peak Intensity) The closer to 0 this value is the better, in general.
#'
#' Peak Distance - This column shows the number of CH2 units between the tallest and second
#' tallest peak in each series. In general, it is better for the value to be as close to 1 as possible.
#'
#' Series Score - This column compares the number of actual observations in each series to
#' the theoretical maximum number based on the CH2 homologous series. The closer to one
#' this value is, the better.
#'
#'
#' While the function does make some minor decisions as to which series to choose
#' from (more than 3 members of homologous series) the rest of the decision
#' making is left to the user.
#'
#' @param df data frame output from MFAssign() or MFAssignCHO()
#'
#' @return data frame
#'
#' @examples
#' MFRecalList(df = Data)
#' @export

MFRecalList <- function(df){
  df$number <- 1
  df$Adduct <- "H"
  df$Adduct <- replace(df$Adduct, df$M > 0, "Na")
  df$Adduct <- replace(df$Adduct, df$POE == 1, "OE")
  df$SeriesAdd <- paste(df$class, df$Adduct, sep = "_")

  df1 <- subset(aggregate(number ~ SeriesAdd + DBE, df,
                          function(x) number=sum(x, na.rm = TRUE)),
                na.action = NULL)

  longseries <- df1[order(-df1$number),]
  longseries$Index <- 1:nrow(longseries)
  Recal <- merge(df, longseries, by.x = c("SeriesAdd", "DBE"), by.y = c("SeriesAdd", "DBE"))
  Recal <- tidyr::unite(Recal, Series, SeriesAdd, DBE, sep = "_", remove = FALSE)

  Recal <- dplyr::group_by(Recal, Index)
  Recal <- dplyr::mutate(Recal, Min = min(Exp_mass), Max = max(Exp_mass), MInt = mean(Abundance),
                         Maxmass = ifelse(Abundance == max(Abundance), Exp_mass, NA), Maxint = max(Abundance),
                         Secint = sort(Abundance, TRUE)[2], Secmass = ifelse(Abundance ==sort(Abundance, TRUE)[2],Exp_mass,NA ))

  Maxmass1 <- Recal[c(1,50)]
  Maxmass1 <- Maxmass1[!is.na(Maxmass1$Maxmass),]
  Secmass1 <- Recal[c(1,53)]
  Secmass1 <- Secmass1[!is.na(Secmass1$Secmass),]

  Recal <- merge(Recal, Maxmass1, by.x = c("Series"), by.y = c("Series"))
  Recal <- merge(Recal, Secmass1, by.x = c("Series"), by.y = c("Series"))

  Recal <- Recal[Recal$number.y > 3,]
  names(Recal)[50] <- "Maxmass"
  Recal <- Recal[!is.na(Recal$Maxmass),]
  Recal$M.window <- "a"
  Recal$M.window[Recal$Maxmass> 0 & Recal$Maxmass < 200] <- "0-200"
  Recal$M.window[Recal$Maxmass> 200 & Recal$Maxmass < 300] <- "200-300"
  Recal$M.window[Recal$Maxmass> 300 & Recal$Maxmass < 400] <- "300-400"
  Recal$M.window[Recal$Maxmass> 400 & Recal$Maxmass < 500] <- "400-500"
  Recal$M.window[Recal$Maxmass> 500 & Recal$Maxmass < 600] <- "500-600"
  Recal$M.window[Recal$Maxmass> 600 & Recal$Maxmass < 700] <- "600-700"
  Recal$M.window[Recal$Maxmass> 700 & Recal$Maxmass < 800] <- "700-800"
  Recal$M.window[Recal$Maxmass> 800] <- ">800"

  Recal1 <- aggregate(MInt ~ M.window, Recal,
                      function(x) number = median(x, na.rm = TRUE))
  Recal <- merge(Recal, Recal1, by = "M.window")

  Recal$IntRel <- (Recal$MInt.x-Recal$MInt.y)/Recal$MInt.y *100
  Recal$Peakcomp <- log10(Recal$Maxint/Recal$Secint)
  Recal$Nextpeak <- abs((Recal$Maxmass.y-Recal$Secmass.y)/14)
  Recal <- Recal[c(2,46:60)]
  Recal$Min <- round(Recal$Min, 3)
  Recal$Max <- round(Recal$Max, 3)
  Recal$SerScor <- ((Recal$Max-Recal$Min)/14+1)/Recal$number.y
  Recal <- tidyr::unite(Recal, Range, Min, Max, sep = "-")
  Recal <- Recal[-c(5,7,8,9,10,11,12)]
  names(Recal)[3] <- "Series Index"
  names(Recal)[2] <- "Number Observed"
  names(Recal)[4] <- "Mass Range"
  names(Recal)[5] <- "Tall Peak"
  names(Recal)[6] <- "Abundance Score"
  names(Recal)[9] <- "Series Score"
  names(Recal)[8] <- "Peak Distance"
  names(Recal)[7] <- "Peak Score"
  Recal <- Recal[!duplicated(Recal), ]
  Recal
}

###############################################################
#' Generates a plot to check chosen recalibrant series, recalibrates two mass lists,
#' and produces a list of the chosen recalibrants.
#'
#'
#' MFRecalCheck() takes the data frame output from MFAssign, the outputs of IsoFiltR(), and
#' the chosen recalibrant series to generate a plot for the qualitative assessment of the recalibrants,
#' recalibrate the "Mono" and "Iso" outputs from IsoFiltR, and prepare a data frame containing the chosen
#' recalibrants.
#'
#'
#' This function can handle up to 10 homologous series, though it will work with fewer. It is important
#' for recalibrant masses to cover the entire mass range of interest, and they should be among the
#' most abundant peaks in their region of the spectrum, this function helps with choosing and
#' visualizing the recalibrants.
#'
#' This function also recalibrates up to two mass lists using the chosen recalibrants.
#' It is best to use the "Mono" and "Iso" outputs of the IsoFiltR() function.
#'
#'
#' These parameters will be used to identify more potential recalibrants, and the tallest peaks within
#' a user defined mass range will be used as recalibrants.
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
#' @param series(1-10) character the recalibrant series, "O7_Na_4" for example, default is NA
#'
#' @param min numeric - minimum mass range of the data being analyzed.
#'
#' @param max numeric - maximum mass range of the data being analyzed.
#'
#' @param bin numeric - sets the mass window bins to choose recalibrants from, default is 20.
#'
#' @param obs numeric - sets the number of peaks to choose as recalibrants from each bin,
#' default is 2
#'
#' @return list(MZ, Mono, Iso, RecalOut) contains the mass spectrum (MZ), recalibrated "Mono" data frame,
#'  recalibrated "Iso" data frame, and the recalibrants list (RecalOut)
#'
#' @examples
#' MFRecalCheck(df = Data, peaks = Mono, isopeaks = Iso, mode = "neg", SN = 500, series1 = "O4_H_2", series2 = "O4_H_8", series3 = "O6_H_8")
#'
#' MFRecalCheck(df = Data, peaks = Mono, isopeaks = Iso, mode = "pos", SN = 300, series1 = "O4_Na_2", series2 = "O4_H_8", series3 = "O6_Na_8")
#' @export


MFRecalCheck <- function( df, peaks, isopeaks = "None",mode, SN = 0, series1=NA, series2=NA, series3=NA, series4=NA, series5=NA,
                           series6=NA, series7=NA, series8=NA, series9=NA, series10=NA, min = 100, max = 1000,
                           bin = 20, obs = 2){
  #Preparation of the recalibrants list
  RecalList <- data.frame(series = c(series1, series2, series3, series4, series5, series6, series7, series8, series9,
                                     series10),
                          Tag = c("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10"))
  RecalList <- RecalList[!is.na(RecalList$series),]

  df$Adduct <- "H"
  df$Adduct <- replace(df$Adduct, df$M > 0, "Na")
  df$Adduct <- replace(df$Adduct, df$POE == 1, "OE")
  df$series <- paste(df$class, df$Adduct, df$DBE, sep = "_")
  df$mode <- mode


  peaks <- peaks[c(2,1)]
  isopeaks <- isopeaks[c(2,1)]

  #Merges recalibrant list to df in order to determine which recalibrants are in the data frame.
  RecalList2 <- merge(df, RecalList, by.x = "series", by.y = "series")

  #Prepares the recalibrant masses for use in recalibration steps.
  RecalList <- RecalList2[c("Abundance", "Exp_mass", "C", "H", "O", "N", "S", "P", "E",
                            "S34", "N15", "D", "Cl", "Cl37", "M", "NH4", "POE", "Z")]
  RecalList$NM <- round(RecalList$Exp_mass)

  RecalList$KM_O <- RecalList$Exp_mass * (16/15.9949146223)
  RecalList$KMD_O <- round(RecalList$NM - RecalList$KM_O, 3)
  RecalList$z_O <- round(RecalList$Exp_mass)%%16 - 16

  RecalList$KM_H2 <- RecalList$Exp_mass * (2/2.01565)
  RecalList$KMD_H2 <- round(RecalList$NM - RecalList$KM_H2, 3)
  RecalList$z_H2 <- round(RecalList$Exp_mass)%%2 - 2

  Rest <- df[c("Abundance", "Exp_mass")]
  Rest$NM <- round(Rest$Exp_mass)

  Rest$KM_O <- Rest$Exp_mass * (16/15.9949146223)
  Rest$KMD_O <- round(Rest$NM - Rest$KM_O, 3)
  Rest$z_O <- round(Rest$Exp_mass)%%16 - 16

  Rest$KM_H2 <- Rest$Exp_mass * (2/2.01565)
  Rest$KMD_H2 <- round(Rest$NM - Rest$KM_H2, 3)
  Rest$z_H2 <- round(Rest$Exp_mass)%%2 - 2

  ############
  #Picking recalibrants with series
  knownO <- RecalList[c(1:18,21,22)]

  names(knownO)[2] <- "base_mass"
  Step2 <- merge(Rest, knownO, by.x = c("KMD_O", "z_O"), by.y = c("KMD_O", "z_O"))
  Step2$O_num <- round(((Step2$Exp_mass - Step2$base_mass))/15.9949146223)
  Step2$O <- Step2$O + Step2$O_num
  Step2$Type <- "O"
  Step2$form <- paste(Step2$C, Step2$H, Step2$O, Step2$N, Step2$S, Step2$P, Step2$E, Step2$S34,
                      Step2$N15, Step2$D, Step2$Cl37, Step2$Cl37, Step2$M, Step2$NH4, Step2$POE, sep = "_")
  Step2 <- Step2[-c(10,28)]

  knownH2 <- RecalList[c(1:18,24,25)]
  names(knownH2)[2] <- "base_mass"
  Step3 <- merge(Rest, knownH2, by.x = c("KMD_H2", "z_H2"), by.y = c("KMD_H2", "z_H2"))
  Step3$H2_num <- round(((Step3$Exp_mass - Step3$base_mass))/2.01565)
  Step3$H <- Step3$H + 2*Step3$H2_num
  Step3$Type <- "H2"
  Step3$form <- paste(Step3$C, Step3$H, Step3$O, Step3$N, Step3$S, Step3$P, Step3$E, Step3$S34,
                      Step3$N15, Step3$D, Step3$Cl37, Step3$Cl37, Step3$M, Step3$NH4, Step3$POE, sep = "_")
  Step3 <- Step3[-c(10,28)]

  Out <- rbind(Step2, Step3)
  Out2 <- dplyr::distinct(Out, Exp_mass)

  NewRecal <- merge(df, Out2, by = "Exp_mass")
  NewRecal <- NewRecal[c("Abundance", "Exp_mass", "formula", "theor_mass")]
  #####################
  #Picking recalibrants by abundnce
  NewRecal$Bin<-cut(NewRecal$Exp_mass,breaks=seq(min, max, by=bin))

  FinalRecal <- dplyr::group_by(NewRecal, Bin)
  FinalRecal2 <- dplyr::top_n(FinalRecal, obs, Abundance)
  NewRecal <- FinalRecal2[c(2:4)]

  #Recalibration for all masses
  #Prepare the mass lists for recalibration
  isopeaks <- if(isopeaks == "None") data.frame(Abundance = 1, mass = 1) else isopeaks
  names(isopeaks)[1] <- "Abundance"
  names(isopeaks)[2] <- "mass"
  names(peaks)[1] <- "Abundance"
  names(peaks)[2] <- "mass"
  peaks$Tag <- "Mono"
  isopeaks$Tag <- "Iso"
  peaks <- if(isopeaks != "None") rbind(peaks, isopeaks) else peaks

  #Calculating the weights needed for recalibration and calculating the recalibrated masses.
  names(NewRecal)[1] <- "E_mass"
  names(NewRecal)[2] <- "formula"
  names(NewRecal)[3] <- "Th_mass"

  NewRecal$num <- 0:(nrow(NewRecal)-1)
  NewRecal$weight <- factorial(nrow(NewRecal)-1)/
    (factorial(NewRecal$num)*factorial((nrow(NewRecal)-1)-NewRecal$num))

  NewRecal$mzweight <- NewRecal$weight*NewRecal$E_mass
  NewRecal$Ejweight <- NewRecal$weight*(NewRecal$E_mass-NewRecal$Th_mass)/NewRecal$Th_mass
  NewRecal$Ejsum <- sum(NewRecal$Ejweight)/2^(nrow(NewRecal)-1)
  NewRecal$masssum <- sum(NewRecal$mzweight)/2^(nrow(NewRecal)-1)
  Ejsum <- mean(NewRecal$Ejsum)
  peaks$mass <- peaks$mass/(1+Ejsum)

  #Separating the isotope and monoisotope masses for output.
  isopeaks2 <- peaks[peaks$Tag != "Mono",]
  isopeaks2 <- isopeaks2[c(1,2)]
  names(isopeaks2)[2] <- "Iso_mass"
  names(isopeaks2)[1] <- "Iso_Abund"
  peaks <- peaks[peaks$Tag == "Mono",]
  peaks <- peaks[c(1,2)]

  ################################################
  #Recalibration for error comparison, calculating the average error before and after recalibration.
  FinalRecal3 <- merge(df, FinalRecal2, by.x = "Exp_mass", by.y = "Exp_mass", all = T)

  FinalRecal3$R_mass <- FinalRecal3$Exp_mass/(1+Ejsum)
  FinalRecal3$R_err <- abs((FinalRecal3$R_mass-FinalRecal3$theor_mass.x)/FinalRecal3$theor_mass.x * 10^6)

  print(mean(FinalRecal3$AE_ppm))
  print(mean(FinalRecal3$R_err))

  FinalRecal3 <- FinalRecal3[!is.na(FinalRecal3$Bin),]
  plotpeak <- peaks[peaks$Abundance > SN,]

  #Plot highlighting the recalibrant ions for qualitative assessment.
  MZ<-ggplot2::ggplot() + ggplot2::geom_segment(data=plotpeak, size=0.7,ggplot2::aes_string(x = "mass", xend = "mass", y = 0, yend = "Abundance"), alpha = 0.3, color = "grey57")+
    ggplot2::geom_segment(data=FinalRecal3, size=1.2,ggplot2::aes_string(x = "Exp_mass", xend = "Exp_mass", y = 0, yend = "Abundance.x"), color = "blue")+

    ggplot2::theme_bw()+ggplot2::labs(x = "Ion Mass", y = "Abundance", title = "Assignment Mass Spectrum", color = "Series")+
    ggplot2::theme(axis.title=ggplot2::element_text(size = 15, face = "bold"), strip.text=ggplot2::element_text(size=15,face="bold"),
                   axis.text=ggplot2::element_text(size=15, face = "bold"), legend.title=ggplot2::element_text(face="bold", size = 15),
                   legend.text=ggplot2::element_text(face="bold", size = 15),  panel.grid.minor.x=ggplot2::element_blank(),
                   panel.grid.major.x=ggplot2::element_blank(), strip.background = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(size = 16, face = "bold"))

  #Preparing the final output of the function, the recalibrants list.

  RecalOut <- FinalRecal3[c(1, 3, 23, 49, 25)]
  names(RecalOut)[2] <- "formula"
  names(RecalOut)[3] <- "theor_mass"
  names(RecalOut)[4] <- "Recal_err"
  names(RecalOut)[5] <- "Orig_err"

  peaks <- peaks[c(2,1)]
  isopeaks2 <- isopeaks2[c(2,1)]

  Output <- list(Plot = MZ, Mono = peaks, Iso = isopeaks2, RecalList = RecalOut)
  Output
}
