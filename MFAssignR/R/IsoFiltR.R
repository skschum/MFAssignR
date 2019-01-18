#' Identifies and separates likely isotopic masses from monoisotopic masses
#'
#' IsoFiltR separates likely isotopic masses from monoisotopic masses in a
#' mass list. It can identify likely 13C isotopic
#' masses and put them in a separate mass list from the monoisotopic
#' masses. This should be done prior to formula assignment in order
#' to lessen the chances of incorrectly assigned formulas, where, for
#' example a 13C containing CHO formula can be assigned to a
#' monoistopic CHNO formula.
#'
#' The only necessary input is a two column data frame with the abundance in the
#' first column and the measured ion mass in the second column. This should be the
#' raw mass list output.
#'
#' At this time only the 13C isotope can be identified, although plans to improve
#' the number of isotopes to be available are in process. The output of this function
#' is a list of dataframes. Dataframe 1 contains the flagged monoisotopic masses, and
#' all masses that did not have a matching isotopic mass. The second dataframe contains
#' the masses flagged as isotopic. These two dataframes should next be run in the
#'  \code{\link{MFAssignCHO}} or \code{\link{MFAssign}} to assign formulas to the masses.
#'
#' Note that the classification of isotopic or monoisotopic from this function is not
#' definitive.
#'
#'
#' @param peaks data frame:
#' The input 2 column data frame containing abundance and peak mass
#'
#'
#' @return list(Monolist, Isolist):
#'   Monolist - monoistopic and non-matched masses,
#'   Isolist - isotopic masses
#'
#' @examples
#' IsoFiltR(peaks=df)
#' @export


IsoFiltR <- function(peaks, SN = 0, Diffrat = 0.1) {

  names(peaks)[1] <- "Exp_mass"
  names(peaks)[2] <- "Abundance"
  Data1 <- peaks[peaks$Abundance >= SN,]
  Data1<- Data1[order(Data1$Exp_mass),]

  Sect <- ceiling(nrow(Data1))/10
  Over <- round(Sect) * 0.05

  Data0 <- Data1[1:(Sect + Over),]
  Data2 <- Data1[Sect:(2*Sect+Over),]
  Data3 <- Data1[(2*Sect):(3*Sect+Over),]
  Data4 <- Data1[(3*Sect):(4*Sect+Over),]
  Data5 <- Data1[(4*Sect):(5*Sect+Over),]
  Data6 <- Data1[(5*Sect):(6*Sect+Over),]
  Data7 <- Data1[(6*Sect):(7*Sect+Over),]
  Data8 <- Data1[(7*Sect):(8*Sect+Over),]
  Data9 <- Data1[(8*Sect):(9*Sect+Over),]
  Data10 <- Data1[(9*Sect):nrow(Data1),]


  #Data frame set up
  End <- Data1
  ###############
  Sulflist <- list(Data0, Data2, Data3, Data4, Data5, Data6, Data7, Data8, Data9, Data10)

  IsoOutS_final <- data.frame(Exp_mass = -42, Abundance = -42)
  MonoOutS_final <- data.frame(Exp_mass = -42, Abundance = -42)
#i <- 3
  for(i in 1:10){
    I <- i
    Data <- Sulflist[[I]]
    Dataend <- Data
    mass <- c(unlist(Data[c(1)]))
    Data <- expand.grid(x = mass, y = mass)
    names(Data)[1] <- "Exp_mass"
    names(Data)[2] <- "Exp_mass1"

    Data$mdiff <- Data$Exp_mass - Data$Exp_mass1
    Data <- Data[Data$mdiff < 0,]
    Data <- Data[Data$mdiff > -2 & Data$mdiff < -1.990,]
    Data$S34_err <- abs(((Data$Exp_mass + 1.995797)-Data$Exp_mass1)/Data$Exp_mass1 * 10^6)
    Data <- Data[Data$S34_err <= 5,]
    Data <- Data[c(1:3)]

    Data$KM <- Data$Exp_mass * (2 / 1.995797)
    Data$KMD <- round((round(Data$Exp_mass) - Data$KM), 3)
    Data$KMr <- Data$Exp_mass * ((round((14.01565/12))/((14.01565/12))))
    Data$KMDr <- round((round(Data$KMr)-Data$KMr),3)

    Data$KM1 <- Data$Exp_mass1 * (2 / 1.995797)
    Data$KMD1 <- round((round(Data$Exp_mass1) - Data$KM1),3)
    Data$KMr1 <- Data$Exp_mass1 * ((round((14.01565/12))/((14.01565/12))))
    Data$KMDr1 <- round((round(Data$KMr1)-Data$KMr1),3)

    Data$KMDrdiff <- Data$KMDr - Data$KMDr1
    Data$KMDdiff <- Data$KMD - Data$KMD1

    Pairs <- Data[abs(Data$KMDdiff) < 0.00249 & ((Data$KMDrdiff < -0.29051 & Data$KMDrdiff > -0.29349) |
                                                       (Data$KMDrdiff < 0.70949 & Data$KMDrdiff > 0.7075)),]
    ##############################################
    MonoS <- Pairs[c(1)]
    IsoS <- Pairs[c(2)]
    names(IsoS)[1] <- "Exp_mass"

    MonoS <- merge(MonoS, End, by.x = "Exp_mass", by.y = "Exp_mass")
    #check <- unique(MonoS)
    IsoS <- merge(IsoS, End, by.x = "Exp_mass", by.y = "Exp_mass")
    #check2 <- unique(IsoS)
    names(IsoS)[1] <- "Iso_mass"
    names(IsoS)[2] <- "Iso_Abund"

    ##Final Abundance Check
    ##Be sure to check this as it is a risky point
    Abund <- cbind(MonoS, IsoS)
    Abunddummy <- data.frame(Exp_mass = -42, Abundance = -1, Iso_mass = -42, Iso_Abund = -1)
    Abund <- rbind(Abund, Abunddummy)
    Abund$ratio <- Abund$Iso_Abund / Abund$Abundance * 100

    Abund <- Abund[Abund$ratio <= 10,]
    Abund <- unique(Abund)

    MonooutS <- Abund[c(1,2)]
    IsooutS <- Abund[c(3,4)]
    names(IsooutS)[1] <- "Exp_mass"
    names(IsooutS)[2] <- "Abundance"

    IsoOutS_final <- rbind(IsoOutS_final, IsooutS)
    MonoOutS_final <- rbind(MonoOutS_final, MonooutS)

  }

  ##Test <- MonoOutS_final %>% mutate(Dups = duplicated(Exp_mass))
 # Test <- Test[Test$Dups == TRUE,]

  IsoOutS_final <- unique(IsoOutS_final)
  MonoOutS_final <- unique(MonoOutS_final)
  #########################
  #Carbon Isotoping

  Carblist <- list(Data0, Data2, Data3, Data4, Data5, Data6, Data7, Data8, Data9, Data10)

  IsoOutC1_final <- data.frame(Exp_mass = -42, Abundance = -42)
  MonoOutC_final <- data.frame(Exp_mass = -42, Abundance = -42)
  IsoOutC2_final <- data.frame(Exp_mass = -42, Abundance = -42)
  ###################################################
  #dao <- (-5/10^6)*1.0033548380+1.0033548380

  for(i in 1:10){
    I <- i
    Data <- Carblist[[I]]
    Dataend <- Data
  mass <- c(unlist(Data[c(1)]))

  Data <- expand.grid(x = mass, y = mass)
  names(Data)[1] <- "Exp_mass"
  names(Data)[2] <- "Exp_mass1"
  Data$mdiff <- Data$Exp_mass - Data$Exp_mass1
  Data <- Data[Data$mdiff < 1.0015 & Data$mdiff > -1.005, ]
  Data$C13_err <- abs(((Data$Exp_mass + 1.0033548380)-Data$Exp_mass1)/Data$Exp_mass1 * 10^6)
  Data <- Data[Data$C13_err <= 5,]
  Data <- Data[c(1:3)]

  #Variable calculation
  Data$KM <- Data$Exp_mass * (1 / 1.0033548380)
  Data$KMD <- round((round(Data$Exp_mass) - Data$KM),3)
  Data$KMr <- Data$Exp_mass * ((round((14.01565/21))/((14.01565/21))))
  Data$KMDr <- round((round(Data$KMr)-Data$KMr),3)

  Data$KM1 <- Data$Exp_mass1 * (1 / 1.0033548380)
  Data$KMD1 <- round((round(Data$Exp_mass1) - Data$KM1),3)
  Data$KMr1 <- Data$Exp_mass1 * ((round((14.01565/21))/((14.01565/21))))
  Data$KMDr1 <- round((round(Data$KMr1)-Data$KMr1),3)

  Data$KMDrdiff <- Data$KMDr - Data$KMDr1
  Data$KMDdiff <- Data$KMD - Data$KMD1
  Data$mdiff <- Data$Exp_mass - Data$Exp_mass1

  #The numbers chosen for filtering are based on the results of positive mode ESI for BB burning aerosol
  Pairs <- Data[abs(Data$KMDdiff) <= 0.00149 & ((Data$KMDrdiff < -0.494501&Data$KMDrdiff > -0.4975) |
                                                      (Data$KMDrdiff < 0.5045 & Data$KMDrdiff > 0.501501)),]

  Monopair <- Pairs[c(1,2)]
  names(Monopair)[1] <- "Mono_mass"
  names(Monopair)[2] <- "Iso_mass1"

  Isopair <- Pairs[c(1,2)]
  names(Isopair)[1] <- "Iso_mass1"
  names(Isopair)[2] <- "Iso_mass2"

  #middle <- cbind(Monopair, Isopair)

  Newpair <- merge(Monopair, Isopair, by.x = "Iso_mass1", by.y = "Iso_mass1", all = T)

  All3 <- Newpair[!is.na(Newpair$Mono_mass)& !is.na(Newpair$Iso_mass2),]

  Iso1 <- All3[c(1)]
  names(Iso1)[1] <- "Mono_mass"
  Isodummy <- data.frame(Mono_mass = -42)
  Iso1 <- rbind(Iso1, Isodummy)
  Iso1$Tag <- "Iso1"


  Finalset <- merge(Monopair, Iso1, by.x = c("Mono_mass"), by.y = c("Mono_mass"), all = T)

  Finalpair <- Finalset[is.na(Finalset$Tag),]

  Finalset2 <- merge(Finalpair, Isopair, by.x = c("Iso_mass1"), by.y = c("Iso_mass1"), all = T)

  Finalpair2 <- Finalset2[!is.na(Finalset2$Mono_mass),]

  Finalpair2 <- Finalpair2[-3]
  Finalpair2 <- Finalpair2[c(2,1,3)]

  #Adding the Abundances back in
  MonoAbund <- End
  names(MonoAbund)[1] <- "Mono_mass"
  names(MonoAbund)[2] <- "Mono_Abund"

  IsoAbund1 <- End
  names(IsoAbund1)[1] <- "Iso_mass1"
  names(IsoAbund1)[2] <- "Iso_Abund1"

  IsoAbund2 <- End
  names(IsoAbund2)[1] <- "Iso_mass2"
  names(IsoAbund2)[2] <- "Iso_Abund2"

  Align1 <- merge(Finalpair2, MonoAbund, by.x = "Mono_mass", by.y = "Mono_mass")
  Align1 <- merge(Align1, IsoAbund1, by.x = "Iso_mass1", by.y = "Iso_mass1")
  Align1 <- merge(Align1, IsoAbund2, by.x = "Iso_mass2", by.y = "Iso_mass2", all= T)
  Align1 <- Align1[!is.na(Align1$Mono_mass),]

  #Diffrat <- 0
  #Calculating Abundance Ratios for QA purposes
  #Ratios estimated with Sisweb using multiples of C

  Align1$Abund <- Align1$Iso_Abund1 / Align1$Mono_Abund
  Align100 <- Align1[Align1$Abund < 0.0973 &  Align1$Mono_mass <=100,]
  Align200 <- Align1[Align1$Abund < 0.1839 & Align1$Abund > Diffrat * 0.0433 & Align1$Mono_mass <=200&
                       Align1$Mono_mass > 100,]
  Align300 <- Align1[Align1$Abund < 0.2704 & Align1$Abund > Diffrat * 0.0757 & Align1$Mono_mass <=300&
                       Align1$Mono_mass > 200,]
  Align400 <- Align1[Align1$Abund < 0.3677 & Align1$Abund > Diffrat * 0.1082 & Align1$Mono_mass <=400&
                       Align1$Mono_mass > 300,]
  Align500 <- Align1[Align1$Abund < 0.4543 & Align1$Abund > Diffrat * 0.1406 & Align1$Mono_mass <=500&
                       Align1$Mono_mass > 400,]
  Align600 <- Align1[Align1$Abund < 0.5408 & Align1$Abund > Diffrat * 0.1731 & Align1$Mono_mass <=600&
                       Align1$Mono_mass > 500,]
  Align700 <- Align1[Align1$Abund < 0.6381 & Align1$Abund > Diffrat * 0.2055 & Align1$Mono_mass <=700&
                       Align1$Mono_mass > 600,]
  Align800 <- Align1[Align1$Abund < 0.7247 & Align1$Abund > Diffrat * 0.2379 & Align1$Mono_mass <=800&
                       Align1$Mono_mass > 700,]
  Align900 <- Align1[Align1$Abund < 0.8112 & Align1$Abund > Diffrat * 0.2704 & Align1$Mono_mass <=900&
                       Align1$Mono_mass > 800,]
  Align1000 <- Align1[Align1$Abund < 0.9085 & Align1$Abund > Diffrat * 0.3167 & Align1$Mono_mass <=1000&
                        Align1$Mono_mass > 900,]

  Align1 <- rbind(Align100, Align200, Align300, Align400, Align500, Align600, Align700,
                  Align800, Align900, Align1000)

  #Now to address the abundances of 13C_2 isotopes
  #Ratios estimated with Sisweb using multiples of CH4O

  IsoPair <- Align1[is.na(Align1$Iso_mass2),]
  IsoTri <- Align1[!is.na(Align1$Iso_mass2),]
  IsoTri$Abund <- IsoTri$Iso_Abund2 / IsoTri$Mono_Abund

  IsoTri100 <- IsoTri[IsoTri$Abund < 0.0042 &  IsoTri$Mono_mass <=100,]
  IsoTri200 <- IsoTri[IsoTri$Abund < 0.0159 & IsoTri$Abund > Diffrat * 0.00007 & IsoTri$Mono_mass <=200&
                        IsoTri$Mono_mass > 100,]
  IsoTri300 <- IsoTri[IsoTri$Abund < 0.0351 & IsoTri$Abund > Diffrat * 0.0025 & IsoTri$Mono_mass <=300&
                        IsoTri$Mono_mass > 200,]
  IsoTri400 <- IsoTri[IsoTri$Abund < 0.0656 & IsoTri$Abund > Diffrat * 0.0053 & IsoTri$Mono_mass <=400&
                        IsoTri$Mono_mass > 300,]
  IsoTri500 <- IsoTri[IsoTri$Abund < 0.1007 & IsoTri$Abund > Diffrat * 0.0091 & IsoTri$Mono_mass <=500&
                        IsoTri$Mono_mass > 400,]
  IsoTri600 <- IsoTri[IsoTri$Abund < 0.1433 & IsoTri$Abund > Diffrat * 0.014 & IsoTri$Mono_mass <=600&
                        IsoTri$Mono_mass > 500,]
  IsoTri700 <- IsoTri[IsoTri$Abund < 0.2002 & IsoTri$Abund > Diffrat * 0.02 & IsoTri$Mono_mass <=700&
                        IsoTri$Mono_mass > 600,]
  IsoTri800 <- IsoTri[IsoTri$Abund < 0.2586 & IsoTri$Abund > Diffrat * 0.027 & IsoTri$Mono_mass <=800&
                        IsoTri$Mono_mass > 700,]
  IsoTri900 <- IsoTri[IsoTri$Abund < 0.3246 & IsoTri$Abund > Diffrat * 0.0351 & IsoTri$Mono_mass <=900&
                        IsoTri$Mono_mass > 800,]
  IsoTri1000 <- IsoTri[IsoTri$Abund < 0.4078 & IsoTri$Abund > Diffrat * 0.0475 & IsoTri$Mono_mass <=1000&
                         IsoTri$Mono_mass > 900,]

  IsoTri <- rbind(IsoTri100, IsoTri200, IsoTri300, IsoTri400, IsoTri500, IsoTri600, IsoTri700,
                  IsoTri800, IsoTri900, IsoTri1000)

  FinalAlign <- rbind(IsoPair, IsoTri)
  #Now need to use the masses to determine which ones have not been accounted for.
  #FinalAlign <- Align1
  MonoC <- FinalAlign[c(3,4)]
  names(MonoC)[1] <- "Exp_mass"
  names(MonoC)[2] <- "Abundance"
  IsoC1 <- FinalAlign[c(2,5)]
  names(IsoC1)[1] <- "Exp_mass"
  names(IsoC1)[2] <- "Abundance"
  IsoC2 <- FinalAlign[c(1,6)]
  names(IsoC2)[1] <- "Exp_mass"
  names(IsoC2)[2] <- "Abundance"

  IsoOutC1_final <- rbind(IsoOutC1_final, IsoC1)
  IsoOutC2_final <- rbind(IsoOutC2_final, IsoC2)
  MonoOutC_final <- rbind(MonoOutC_final, MonoC)


}

  ##############################################
  IsoOutC1_final <- unique(IsoOutC1_final)
  IsoOutC1_final$Tag <- "C13"
  IsoOutC2_final <- unique(IsoOutC2_final)
  IsoOutC2_final$Tag <- "2C13"
  MonoOutC_final <- unique(MonoOutC_final)
  MonoOutC_final$Tag <- "C"

  MonoOutS_final$Tag <- "S"
  IsoOutS_final$Tag <- "S34"

  ###########
  #Set up the isotopes
  Isotopes <- rbind(IsoOutC1_final, IsoOutC2_final, IsoOutS_final)
  Isotopes <- Isotopes[Isotopes$Exp_mass > 0,]
  Isotopes$order <- 1:nrow(Isotopes)

  Isotopes <- Isotopes[order(Isotopes$order),]
  Isotopes$Dups1 <- duplicated(Isotopes$Exp_mass)

  Isotopes <- Isotopes[order(-Isotopes$order),]
  Isotopes$Dups2 <- duplicated(Isotopes$Exp_mass)

  Doubles <- Isotopes[(Isotopes$Dups1 == TRUE| Isotopes$Dups2 == TRUE),]
  DS34 <- Doubles[Doubles$Tag == "S34",]
  DC13 <- Doubles[Doubles$Tag != "S34",]

  NewDub <- merge(DS34, DC13, by.x = "Exp_mass", by.y = "Exp_mass")
  NewDub <- NewDub[c(1,2,3,8)]
  NewDub$Tag <- paste(NewDub$Tag.y, NewDub$Tag.x, sep = "_")

  NewDub <- NewDub[c(1,2,5)]
  names(NewDub)[2] <- "Abundance"

  Singles <- Isotopes[(Isotopes$Dups1 == FALSE& Isotopes$Dups2 == FALSE),]
  Singles <- Singles[c(1,2,3)]

  Iso_Out <- rbind(Singles, NewDub)
  ############
  #This section organizes the data and helps to ensure there are no peaks being considered as
  #both monoisotopic and isotopic peaks.
  Mono_out <- rbind(MonoOutC_final, MonoOutS_final)

  Dup_mass <- merge(Iso_Out, Mono_out, by.x = "Exp_mass", by.y = "Exp_mass")

  Dup_rem <- Dup_mass[c(1,2)]

  Iso_Sup <- Dup_mass[c(1,2,3)]
  names(Iso_Sup)[2] <- "Abundance"
  names(Iso_Sup)[3] <- "Tag"

  Iso_Out <- merge(Iso_Out, Dup_rem, by.x = "Exp_mass", by.y = "Exp_mass", all = TRUE)

  Iso_Out <- Iso_Out[is.na(Iso_Out$Abundance.x),]
  Iso_Out <- Iso_Out[c(1,2,3)]




  #The Iso_mass in this section is the final output isotope list

  Iso_mass <- unique(Iso_Out)

  #This data frame is necessary to remove the flagged isotope peaks from
  #the overall data frame.
  Iso_align <- Iso_mass[c(1,2)]

  Aligned <- merge(End, Iso_align, by.x = "Exp_mass", by.y = "Exp_mass", all = TRUE)

  Mono_final <- Aligned[is.na(Aligned$Abundance.y),]
  names(Mono_final)[1] <- "exp_mass"
  names(Mono_final)[2] <- "abundance"
  Mono_final <- Mono_final[c(1,2)]
  Mono_final <- unique(Mono_final)
  Mono_final <- Mono_final[!is.na(Mono_final$abundance),]

  Iso_final <- rbind(Iso_mass, Iso_Sup)
  Iso_final <- unique(Iso_final)
  names(Iso_final)[1] <- "exp_mass"
  names(Iso_final)[2] <- "abundance"
  names(Iso_final)[3] <- "tag"
  Iso_final <- Iso_final[!is.na(Iso_final$abundance),]

  #####################################


  Output <- list(Mono = Mono_final, Iso = Iso_final)
}


