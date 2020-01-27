#' Identifies and separates likely isotopic masses from monoisotopic masses
#'
#' IsoFiltR separates likely isotopic masses from monoisotopic masses in a
#' mass list. It can identify likely 13C and 34S isotopic
#' masses and put them in a separate mass list from the monoisotopic
#' masses. This should be done prior to formula assignment in order
#' to lessen the chances of incorrectly assigned formulas, where, for
#' example a 13C containing CHO formula can be assigned to a
#' monoistopic CHNOS formula.
#'
#' The only necessary input is a two column data frame with the abundance in the
#' first column and the measured ion mass in the second column. This should be the
#' raw mass list output.
#'
#' The output of this function is a list of dataframes. Dataframe 1 contains the flagged
#' monoisotopic masses, and all masses that did not have a matching isotopic mass. The
#' second dataframe contains the masses flagged as isotopic. These two dataframes should next be run in the
#'  \code{\link{MFAssignCHO}} or \code{\link{MFAssign}} to assign formulas to the masses.
#'
#' Note that the classification of isotopic or monoisotopic from this function is not
#' definitive.
#'
#'
#' @param peaks data frame:
#' The input 2 column data frame containing abundance and peak mass
#'
#' @param SN numeric:
#' Sets the noise cut for the data, peaks below this value will not be evaluated
#'
#' @param Carbrat numeric:
#' Sets the maximum 13C/12C ratio that is allowed for matching, default is 60
#'
#' @param Sulfrat numeric:
#' Sets the maximum 34S/32S ratio that is allowed for matching, default is 30
#'
#' @param Sulferr numeric:
#' Sets the maximum allowed error (ppm) for 34S mass matching, default is 5
#'
#' @param Carberr numeric:
#' Sets the maximum allowed error (ppm) for 13C mass matching, default is 5
#'
#'
#' @return list(Monolist, Isolist):
#'   Monolist - monoistopic and non-matched masses,
#'   Isolist - isotopic masses
#'
#' @examples
#' IsoFiltR(peaks=df)
#' @export

#peaks <- DataX[DataX$Intensity > 3*20,]

IsoFiltR <- function(peaks, SN = 0, Carbrat = 60, Sulfrat = 30, Sulferr = 5, Carberr = 5) {

  names(peaks)[1] <- "Exp_mass"
  names(peaks)[2] <- "Abundance"
  Data1 <- peaks[peaks$Abundance >= SN,]
  Data1<- Data1[order(Data1$Exp_mass),]

  Sect <- ceiling(nrow(Data1))/10
  Over <- round(Sect) * 0.15

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
  #Data <- data.frame(Exp_mass = 384.09413, Exp_mass1 = 385.09672)

  Data$mdiff <- Data$Exp_mass - Data$Exp_mass1
  Data <- Data[Data$mdiff < -1.0015 & Data$mdiff > -1.005, ]
  Data$C13_err <- abs(((Data$Exp_mass + 1.0033548380)-Data$Exp_mass1)/Data$Exp_mass1 * 10^6)
  Data <- Data[Data$C13_err <= Carberr,]
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
  #Before change on 6/25/19, the second KMDrdiff was -0.4975
  Pairs <- Data[abs(Data$KMDdiff) <= 0.00149 & ((Data$KMDrdiff < -0.494501&Data$KMDrdiff > -0.498001) |
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

  Align1 <- Align1[Align1$Iso_Abund1 < (Carbrat/100) * Align1$Mono_Abund,]  #New Adjustment 12/06/19


  IsoPair <- Align1[is.na(Align1$Iso_mass2),]
  IsoTri <- Align1[!is.na(Align1$Iso_mass2),]
  IsoTri <- IsoTri[IsoTri$Iso_Abund2 < (Carbrat/100) * IsoTri$Iso_Abund1,]  #New Adjustment 12/06/19


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

  #######Sulfur Section
  ###############
  Sulflist <- list(Data0, Data2, Data3, Data4, Data5, Data6, Data7, Data8, Data9, Data10)
  #Sulflist <- list(Data1)
  IsoOutS_final <- data.frame(Exp_mass = -42, Abundance = -42)
  MonoOutS_final <- data.frame(Exp_mass = -42, Abundance = -42)
  #i <- 1
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
    Data <- Data[Data$S34_err <= Sulferr,]
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

    Abund <- Abund[Abund$ratio <= Sulfrat,]
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
  ######End of Sulfur Section
  #########################

  ##############################################
  IsoOutC1_final <- unique(IsoOutC1_final)
  IsoOutC1_final$Tag <- "C13"
  IsoOutC2_final <- unique(IsoOutC2_final)
  IsoOutC2_final$Tag <- "2C13"
  MonoOutC_final <- unique(MonoOutC_final)
  MonoOutC_final$Tag <- "C"

  MonoOutS_final$Tag <- "S"
  IsoOutS_final$Tag <- "S34"

  ##Sulfur Check 6/24/19
  #Filtering the 34S masses that do not have a corresponding 13C
  Sulfs <- merge(MonoOutC_final, MonoOutS_final, by.x = c("Exp_mass", "Abundance"),
                    by.y = c("Exp_mass", "Abundance"), all = TRUE)
  Sulfgood <- Sulfs[!is.na(Sulfs$Tag.x) & !is.na(Sulfs$Tag.y),]
  Sulfbad <- Sulfs[is.na(Sulfs$Tag.x) & !is.na(Sulfs$Tag.y),]
  ###Comparison of "bad" 34S to C13
  Carbs1 <- merge(IsoOutC1_final, Sulfbad, by.x = c("Exp_mass", "Abundance"),
                 by.y = c("Exp_mass", "Abundance"), all = TRUE)
  IsoOutC1_final <- Carbs1[(!is.na(Carbs1$Tag) & is.na(Carbs1$Tag.y)) |
                             (!is.na(Carbs1$Tag) & !is.na(Carbs1$Tag.y)),]
  IsoOutC1_final <- IsoOutC1_final[c(1:3)]
  Sulfbad <- Carbs1[(is.na(Carbs1$Tag) & !is.na(Carbs1$Tag.y)),]
  ###Comparison of "bad"34S to 2C13
  Sulfbad <- Sulfbad[-3]
  Carbs2 <- merge(IsoOutC2_final, Sulfbad, by.x = c("Exp_mass", "Abundance"),
                  by.y = c("Exp_mass", "Abundance"), all = TRUE)
  IsoOutC2_final <- Carbs2[(!is.na(Carbs2$Tag) & is.na(Carbs2$Tag.y)) |
                             (!is.na(Carbs2$Tag) & !is.na(Carbs2$Tag.y)),]
  IsoOutC2_final <- IsoOutC2_final[c(1:3)]
  Sulfbad <- Carbs2[(is.na(Carbs2$Tag) & !is.na(Carbs2$Tag.y)),]
  ##Sulfbad are the possible monoisotopic sulfur that had a matching 34S, but no 13C, and
  ##these monoisotopic masses did not match as a possible 13C for a different mass, so they
  ##are staying as they are and are being added to the "Sulfgood" mass list.
  Sulfbad <- Sulfbad[c(1,2,5)]
  names(Sulfbad)[3] <- "Tag"

  Sulfgood <- Sulfgood[c(1,2,4)]
  names(Sulfgood)[3] <- "Tag"

  MonoOutS_final <- rbind(Sulfgood, Sulfbad)
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



