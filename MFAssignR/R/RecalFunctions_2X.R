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
#' @param bin numeric - sets the mass window bins to choose recalibrants from, default is 14.
#'
#' @param obs numeric - sets the number of peaks to choose as recalibrants from each bin.
#' Default is 1
#'
#' @param num numeric - sets the number of peaks on either side of defined recalibrant to choose
#' as additional recalibrants. Default is 5.
#'
#' @return list(MZ, Mono, Iso, RecalOut) contains the mass spectrum (MZ), recalibrated "Mono" data frame,
#'  recalibrated "Iso" data frame, and the recalibrants list (RecalOut)
#'
#' @examples
#' Recal_2X(df = Data, peaks = Mono, isopeaks = Iso, mode = "neg", SN = 500, series1 = "O4_H_2", series2 = "O4_H_8", series3 = "O6_H_8")
#'
#' Recal_2X(df = Data, peaks = Mono, isopeaks = Iso, mode = "pos", SN = 300, series1 = "O4_Na_2", series2 = "O4_H_8", series3 = "O6_Na_8")
#' @export

Recal_2X <- function( df, peaks, isopeaks = "none",mode, SN = 0, mzRange = 50, series1=NA, series2=NA, series3=NA, series4=NA, series5=NA,
                          series6=NA, series7=NA, series8=NA, series9=NA, series10=NA, min = 100, max = 1000,
                          bin = 14, obs = 1, num = 5){
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

  df$form <- paste(df$formula, df$Adduct) #New line to ensure no cross formula matching

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


  NewRecal <- RecalList2[c("Exp_mass", "form", "theor_mass", "C13_mass", "Abundance")]


  ##########################
  ##Full Spectrum Recal Term
  NewRecal$Bin<-cut(NewRecal$Exp_mass,breaks=seq(min, max, by=bin))

  FinalRecal <- dplyr::group_by(NewRecal, Bin)
  FullRecal <- dplyr::top_n(FinalRecal, obs, Abundance)

  #FullRecal <- dplyr::top_n(FinalRecal2, obs, Abundance)
  #low_recal <- dplyr::top_n(low, num, Abundance)
  names(FullRecal)[1] <- "E_mass"
  names(FullRecal)[3] <- "Th_mass"
  FullRecal$num <- 0:(nrow(FullRecal)-1)
  FullRecal$weight <- factorial(nrow(FullRecal)-1)/
    (factorial(FullRecal$num)*factorial((nrow(FullRecal)-1)-FullRecal$num))

  FullRecal$mzweight <- FullRecal$weight*FullRecal$E_mass
  FullRecal$Ejweight <- FullRecal$weight*(FullRecal$E_mass-FullRecal$Th_mass)/FullRecal$Th_mass
  FullRecal$Ejsum <- sum(FullRecal$Ejweight)/2^(nrow(FullRecal)-1)
  FullRecal$masssum <- sum(FullRecal$mzweight)/2^(nrow(FullRecal)-1)
  FullEjsum <- mean(FullRecal$Ejsum)

  #Setting up for second round of recalibration
  Forms <- df[c("Exp_mass", "form")]  #For switching the exp masses
  names(Forms)[1] <- "mass"
  Forms2 <- df[c("theor_mass", "form", "M", "POE", "NOE", "class", "DBE")]  #For adding the theor masses back to continue recalibration

  peaks_full <- merge(peaks, Forms, by.x = "mass", by.y= "mass", all = T)
  peaks_full <- unique(peaks_full)
  peaks_full$mass <- peaks_full$mass/(1+FullEjsum)
  df_abund <- peaks_full[!is.na(peaks_full$form),]
  df_abund <- df_abund[c("mass", "Abundance", "form", "Tag2")]

  df1 <- merge(df_abund, Forms2, by.x = "form", by.y = "form")
  df1 <- unique(df1)
  df1$Adduct <- "H"
  df1$Adduct <- replace(df1$Adduct, df1$M > 0, "Na")
  df1$Adduct <- replace(df1$Adduct, df1$POE == 1, "OE")
  df1$Adduct <- replace(df1$Adduct, df1$NOE == 1, "OE")
  df1$series <- paste(df1$class, df1$Adduct, df1$DBE, sep = "_")
  df1$mode <- mode
  names(df1)[3] <- "Abundance"
  names(df1)[2] <- "Exp_mass"

  RecalList3 <- merge(df1, RecalList, by.x = "series", by.y = "series")

  ##Breaking the spectrum into pieces for the piecewise recalibration
  pieces <- list()
  RecalList3 <- RecalList3[order(RecalList3$Exp_mass),]
  RecalListmass <- data.frame(RecalList3$Exp_mass)
  #j <- 7
  for(j in 1:nrow(RecalList3)){

    center1 <- RecalList3[j,]
    centermass1 <- center1$Exp_mass
    I <- j+1
    J <- j-1

  #Upper bound
    if(I <= nrow(RecalList3)){
      center2 <- RecalList3[I,]
    } else{
        center2 <- max(peaks_int$mass)
    }

     if(I <= nrow(RecalList3)){
       centermass2 <- center2$Exp_mass
     } else{
       centermass2 <- center2
     }

   #Lower bound
    if(J == 0){
      center3 <- 0
    } else{
      center3 <- RecalList3[J,]
    }

    if(J == 0){
      centermass3 <- 0
    } else{
      centermass3 <- center3$Exp_mass
    }


    peaks_int <- peaks_full
    peaks_int$limit1 <- (peaks_int$mass - centermass1)
    peaks_int$limit2 <- (peaks_int$mass - centermass2)
    peaks_int$limit3 <- (peaks_int$mass - centermass3)

    if(I ==2){
    piecesx <- peaks_int[abs(peaks_int$limit1) <= abs(peaks_int$limit2),]
    }

    if(I != 2 & centermass2 != max(peaks_int$mass)){

      piecesx <- peaks_int[peaks_int$mass >= centermass3 + abs(centermass1-centermass3)/2 &
                             peaks_int$mass <= centermass2 - abs(centermass1-centermass2)/2,]
    }

    if(centermass2 == max(peaks_int$mass)){   #To account for all masses beyond the highest mass recalibrant

      piecesx <- peaks_int[peaks_int$mass > RecalListmass[j-1,]+(RecalListmass[j,] - RecalListmass[j-1,])/2 ,]
    }

    piecesx$sect <- paste("S",j, sep = "")
    pieces[[j]] <- piecesx

  }

  ##Sub function to do the recalibration of the spectrum sections

  theor_form <- df1[c("form", "theor_mass", "Abundance")]


  ##Selecting recalibrants within each piece and performing recalibration
  SubRecal <- function(dfx){
    df_form <- dfx[!is.na(dfx$form),]
    df_form <- merge(df_form, theor_form, by.x = "form", by.y = "form")
    df_form <- df_form[(df_form$Abundance.x == df_form$Abundance.y),]
    df_form <- df_form[-11]
    names(df_form)[3] <- "Abundance"
    names(df_form)[2] <- "E_mass"
    names(df_form)[10] <- "Th_mass"

    recalpeak <- df_form[df_form$limit1 == 0,]
    low <- df_form[df_form$E_mass < recalpeak$E_mass,]
    high <- df_form[df_form$E_mass > recalpeak$E_mass,]

    low_recal <- dplyr::top_n(low, num, Abundance)
    high_recal <- dplyr::top_n(high, num, Abundance)

    recals <- rbind(recalpeak, low_recal, high_recal)

    recals$num <- 0:(nrow(recals)-1)
    recals$weight <- factorial(nrow(recals)-1)/
      (factorial(recals$num)*factorial((nrow(recals)-1)-recals$num))

    recals$mzweight <- recals$weight*recals$E_mass
    recals$Ejweight <- recals$weight*(recals$E_mass-recals$Th_mass)/recals$Th_mass
    recals$Ejsum <- sum(recals$Ejweight)/2^(nrow(recals)-1)
    recals$masssum <- sum(recals$mzweight)/2^(nrow(recals)-1)
    Ejsum <- mean(recals$Ejsum)

    peaks_out <- dfx
    peaks_out$mass <- peaks_out$mass/(1+Ejsum)
    peaks_out$Ejsum <- Ejsum
    Output <- list(peaks_out, recals)
  }


  #Extracting and binding all the recalibrated sections into useful dataframes
  recal_pieces <- lapply(pieces, SubRecal)


  piecemass <- do.call("rbind", recal_pieces)

  masses <- piecemass[,1]
  recal_mass <- piecemass[,2]

  piecemass <- do.call("rbind", masses)
  recal_mass <- do.call("rbind", recal_mass)
########################################################################
  ###Dealing with the export of the new recalibrated error
  Assigned <- piecemass[!is.na(piecemass$form),]
  Assigned <- merge(Assigned, theor_form, by.x = "form", by.y = "form")
  Assigned <- Assigned[Assigned$Abundance.x == Assigned$Abundance.y,]
  Assigned <- Assigned[!duplicated(Assigned$mass),]

  Exp_form <- df[c("Exp_mass", "form", "Abundance")]

  Assigned <- merge(Assigned, Exp_form, by.x = "form", by.y = "form")
  Assigned <- Assigned[Assigned$Abundance.x == Assigned$Abundance,]
  Assigned <- Assigned[!duplicated(Assigned$mass),]

  Assigned$err <- abs((Assigned$Exp_mass - Assigned$theor_mass)/ Assigned$theor_mass * 10^6)
  Assigned$R_err <- abs((Assigned$mass - Assigned$theor_mass)/ Assigned$theor_mass * 10^6)

  print(mean(Assigned$err))
  print(mean(Assigned$R_err))


  ##Setting up the data for plotting the recalibrants
  RecalPlot <- merge(recal_mass, Exp_form, by.x = "form", by.y = "form")
  RecalPlot <- RecalPlot[RecalPlot$Abundance.x == RecalPlot$Abundance.y,]
  RecalPlot <- RecalPlot[!duplicated(RecalPlot$E_mass),]

  RecalPlot <- RecalPlot[c(1,2,3,9,10,17)]
  RecalPlot$err <- (RecalPlot$Exp_mass - RecalPlot$Th_mass)/ RecalPlot$Th_mass * 10^6
  RecalPlot$R_err <- (RecalPlot$E_mass - RecalPlot$Th_mass)/ RecalPlot$Th_mass * 10^6


  #Plot highlighting the recalibrant ions for qualitative assessment.
  MZ<-ggplot2::ggplot() + ggplot2::geom_segment(data=piecemass, size=0.7,ggplot2::aes_string(x = "mass", xend = "mass", y = 0, yend = "Abundance"), alpha = 0.3, color = "grey57")+
    ggplot2::geom_segment(data=RecalPlot, size=1.2,ggplot2::aes_string(x = "E_mass", xend = "E_mass", y = 0, yend = "Abundance.x"), color = "blue")+

    ggplot2::theme_bw()+ggplot2::labs(x = "Ion Mass", y = "Abundance", title = "Assignment Mass Spectrum", color = "Series")+
    ggplot2::theme(axis.title=ggplot2::element_text(size = 15, face = "bold"), strip.text=ggplot2::element_text(size=15,face="bold"),
                   axis.text=ggplot2::element_text(size=15, face = "bold"), legend.title=ggplot2::element_text(face="bold", size = 15),
                   legend.text=ggplot2::element_text(face="bold", size = 15),  panel.grid.minor.x=ggplot2::element_blank(),
                   panel.grid.major.x=ggplot2::element_blank(), strip.background = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(size = 16, face = "bold"))

  #Preparing the final output of the function, the recalibrants list.



  peaks <- piecemass[c(1,2,3, 4)]
  names(peaks)[2] <- "abundance"
  Mono <- peaks[peaks$Tag2 != "Iso",]
  Mono <- Mono[c(1,2)]

  Iso <- peaks[peaks$Tag2 == "Iso",]

  Iso <- Iso[c(1,2,3)]
  names(Iso)[2] <- "iso_abund"
  names(Iso)[1] <- "iso_mass"
  names(Iso)[3] <- "tag"
  names(RecalPlot)[2] <- "recal_mass"
  RecalOut <- RecalPlot[-c(4,6)]
  names(RecalOut)[3] <- "abundance"
  names(RecalOut)[4] <- "theor_mass"
  names(RecalOut)[6] <- "recal_err"

  Output <- list(Plot = MZ, Mono = Mono, Iso = Iso, RecalList = RecalOut)
  Output
}
