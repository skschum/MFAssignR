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


IsoFiltR <- function(peaks) {

  peaks <- peaks[c(2,1)]
  names(peaks)[2] <- "mass"
  names(peaks)[1] <- "RA"
  peakskeep <- peaks
  peakskeep <- dplyr::distinct(peakskeep, mass, .keep_all = TRUE)
  peaks$KMC <- peaks$mass* (1/1.0033548380)
  peaks$KMDC <- round(peaks$mass)-peaks$KMC
  peaks$zstar <- round(peaks$mass)%%14 - 14
  peaks$KMDTestC <- round(peaks$KMDC, 3)

  peaks$EvenOdd <- Even(peaks$mass)
  Odd <- peaks[peaks$EvenOdd == FALSE,]
  Odd$zstarAl <- Odd$zstar + 1
  Odd$zstarAl[Odd$zstarAl == 0] <- -14
  #Odd$KMDAl <- Odd$KMDTest - 0.002
  Odd <- dplyr::distinct(Odd, mass, .keep_all = TRUE)

  EvenAl <- peaks[peaks$EvenOdd == TRUE,]
  names(EvenAl)[2] <- "Iso_mass"
  names(EvenAl)[1] <- "Iso_RA"
  names(EvenAl)[5] <- "zstarAl"
  #names(EvenAl)[6] <- "KMDAl"
  #Evenkeep <- peaks[peaks$EvenOdd == TRUE,]


  Bind1 <- merge(Odd, EvenAl, by.x = c("zstarAl", "KMDTestC"), by.y = c("zstarAl", "KMDTestC"))
  Bind1$Comp_mass <- Bind1$mass + 1.0033548380
  Bind1$Err <- abs((Bind1$Comp_mass - Bind1$Iso_mass) / Bind1$Comp_mass) * 10^6
  #Selecting the good isotope masses.
  GoodIso1 <- Bind1[Bind1$Err <= 3,]
  GoodIso1$Abund <- GoodIso1$Iso_RA / GoodIso1$RA
  #OutIso1 <- GoodIso1[GoodIso1$Abund >= 0.368,] #Collect the ones that don't fit the limit.
  GoodIso100 <- GoodIso1[GoodIso1$Abund < 0.078 & GoodIso1$mass <=100,]
  GoodIso200 <- GoodIso1[GoodIso1$Abund < 0.156 & GoodIso1$mass <=200&GoodIso1$mass > 100,]
  GoodIso300 <- GoodIso1[GoodIso1$Abund < 0.234 & GoodIso1$mass <=300&GoodIso1$mass > 200,]
  GoodIso400 <- GoodIso1[GoodIso1$Abund < 0.312 & GoodIso1$mass <=400&GoodIso1$mass > 300,]
  GoodIso500 <- GoodIso1[GoodIso1$Abund < 0.39 & GoodIso1$mass <=500&GoodIso1$mass > 400,]
  GoodIso600 <- GoodIso1[GoodIso1$Abund < 0.468 & GoodIso1$mass <=600&GoodIso1$mass > 500,]
  GoodIso700 <- GoodIso1[GoodIso1$Abund < 0.557 & GoodIso1$mass <=700&GoodIso1$mass > 600,]
  GoodIso800 <- GoodIso1[GoodIso1$Abund < 0.635 & GoodIso1$mass <=800&GoodIso1$mass > 700,]
  GoodIso900 <- GoodIso1[GoodIso1$Abund < 0.713 & GoodIso1$mass <=900&GoodIso1$mass > 800,]
  GoodIso1000 <- GoodIso1[GoodIso1$Abund < 0.791 & GoodIso1$mass <=1000&GoodIso1$mass > 900,]

  GoodIso1 <- rbind(GoodIso100, GoodIso200, GoodIso300, GoodIso400, GoodIso500, GoodIso600, GoodIso700,
                    GoodIso800, GoodIso900, GoodIso1000)

  #GoodIso1 <- GoodIso1[GoodIso1$Abund < 0.368,] #Filter to ensure the Iso abundance is lower than the max allowable abundance.

  GoodMono1 <- GoodIso1[c(3,4)]
  GoodMono1 <- dplyr::distinct(GoodMono1, mass, .keep_all = TRUE)
  GoodMono1$Tag <- "Mono"
  #OutMono1 <- OutIso1[c(3,4)]
  #OutMono1 <- dplyr::distinct(OutMono1, mass, .keep_all = TRUE)

  GoodIso1 <- GoodIso1[c(9,10)]
  names(GoodIso1)[2] <- "mass"
  names(GoodIso1)[1] <- "RA"


  GoodIso1$Tag <- "Iso"
  GoodIso1 <- dplyr::distinct(GoodIso1, mass, .keep_all = TRUE)
  #OutIso1 <- dplyr::distinct(OutIso1, mass, .keep_all = TRUE)

  #Finding the non-matching peaks
  GoodPeaks <- rbind(GoodIso1, GoodMono1)

  PeakAlign <- dplyr::left_join(peakskeep, GoodPeaks, by = "mass")
  NoMatch <- PeakAlign[is.na(PeakAlign$Tag),]
  NoMatch <- NoMatch[c(1,2)]
  names(NoMatch)[1] <- "RA"
  #############################################
  NoMatch$KMC <- NoMatch$mass* (1/1.0033548380)
  NoMatch$KMDC <- round(NoMatch$mass)-NoMatch$KMC
  NoMatch$zstar <- round(NoMatch$mass)%%14 - 14
  NoMatch$KMDTestC <- round(NoMatch$KMDC, 3)

  NoMatch$EvenOdd <- Even(NoMatch$mass)
  Odd3 <- NoMatch[NoMatch$EvenOdd == FALSE,]
  names(Odd3)[2] <- "Iso_mass"
  names(Odd3)[1] <- "Iso_RA"
  names(Odd3)[5] <- "zstarAl"

  Odd3 <- dplyr::distinct(Odd3, Iso_mass, .keep_all = TRUE)

  Even3 <- NoMatch[NoMatch$EvenOdd == TRUE,]
  Even3$zstarAl <- Even3$zstar + 1
  Even3$zstarAl[Even3$zstarAl == 0] <- -14
  #EvenKeep <- peaks[peaks$EvenOdd == TRUE,]


  Bind2 <- merge(Even3, Odd3, by.x = c("zstarAl", "KMDTestC"), by.y = c("zstarAl", "KMDTestC"))
  Bind2$Comp_mass <- Bind2$mass + 1.0033548380
  Bind2$Err <- abs((Bind2$Comp_mass - Bind2$Iso_mass) / Bind2$Comp_mass) * 10^6
  #Selecting the good isotope masses.
  GoodIso2 <- Bind2[Bind2$Err <= 3,]
  GoodIso2$Abund <- GoodIso2$Iso_RA / GoodIso2$RA
  #OutIso1 <- GoodIso1[GoodIso1$Abund >= 0.368,] #Collect the ones that don't fit the limit.
  GoodIso100 <- GoodIso2[GoodIso2$Abund < 0.078 & GoodIso2$mass <=100,]
  GoodIso200 <- GoodIso2[GoodIso2$Abund < 0.156 & GoodIso2$mass <=200&GoodIso2$mass > 100,]
  GoodIso300 <- GoodIso2[GoodIso2$Abund < 0.234 & GoodIso2$mass <=300&GoodIso2$mass > 200,]
  GoodIso400 <- GoodIso2[GoodIso2$Abund < 0.312 & GoodIso2$mass <=400&GoodIso2$mass > 300,]
  GoodIso500 <- GoodIso2[GoodIso2$Abund < 0.39 & GoodIso2$mass <=500&GoodIso2$mass > 400,]
  GoodIso600 <- GoodIso2[GoodIso2$Abund < 0.468 & GoodIso2$mass <=600&GoodIso2$mass > 500,]
  GoodIso700 <- GoodIso2[GoodIso2$Abund < 0.557 & GoodIso2$mass <=700&GoodIso2$mass > 600,]
  GoodIso800 <- GoodIso2[GoodIso2$Abund < 0.635 & GoodIso2$mass <=800&GoodIso2$mass > 700,]
  GoodIso900 <- GoodIso2[GoodIso2$Abund < 0.713 & GoodIso2$mass <=900&GoodIso2$mass > 800,]
  GoodIso1000 <- GoodIso2[GoodIso2$Abund < 0.791 & GoodIso2$mass <=1000&GoodIso2$mass > 900,]

  GoodIso2 <- rbind(GoodIso100, GoodIso200, GoodIso300, GoodIso400, GoodIso500, GoodIso600, GoodIso700,
                    GoodIso800, GoodIso900, GoodIso1000) #Filter to ensure the Iso abundance is lower than the max allowable abundance.

  GoodMono2 <- GoodIso2[c(3,4)]
  GoodMono2 <- dplyr::distinct(GoodMono2, mass, .keep_all = TRUE)
  GoodMono2$Tag <- "Mono"
  #OutMono1 <- OutIso1[c(3,4)]
  #OutMono1 <- dplyr::distinct(OutMono1, mass, .keep_all = TRUE)

  GoodIso2 <- GoodIso2[c(9,10)]
  names(GoodIso2)[2] <- "mass"
  names(GoodIso2)[1] <- "RA"


  GoodIso2$Tag <- "Iso"
  GoodIso2 <- dplyr::distinct(GoodIso2, mass, .keep_all = TRUE)
  #OutIso1 <- dplyr::distinct(OutIso1, mass, .keep_all = TRUE)

  #Finding the non-matching peaks
  GoodPeaks2 <- rbind(GoodIso2, GoodMono2, GoodPeaks)

  PeakAlign2 <- dplyr::left_join(peakskeep, GoodPeaks2, by = "mass")
  NoMatch2 <- PeakAlign2[is.na(PeakAlign2$Tag),]
  NoMatch2 <- NoMatch2[c(1,2)]
  names(NoMatch2)[1] <- "RA"

  Mono <- GoodPeaks2[GoodPeaks2$Tag == "Mono",]
  Mono <- Mono[c(1,2)]
  Iso <- GoodPeaks2[GoodPeaks2$Tag == "Iso",]
  Iso <- Iso[c(1,2)]

  Monofinal <- rbind(Mono, NoMatch2)
  names(Monofinal)[1] <- "Abundance"
  names(Iso)[1] <- "Abundance"

  Monofinal <- Monofinal[c(2,1)]
  Iso <- Iso[c(2,1)]

  output <- list(Mono = Monofinal, Iso = Iso)
  output
}
