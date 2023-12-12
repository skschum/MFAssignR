#' Assigns all possible MF to each row of input data frame with CHOFIT and formula extension
#'
#' MFAssign assigns all possible molecular formulae to each
#' mass in the input file, subject to user constraints on the moles of
#' C, H, O, N, S, P, sodium (M), C13 (E), S34, N15, Deuterium (D), Cl, F, Cl37, NH4, and Z.
#'
#' There are user inputs for heteroatoms, adducts, and charge.
#' The terms for each are Nx, Sx, Px, Ex, S34x, N15x, Dx,
#' Clx, Fx, Cl37x, Mx, NH4x, and Zx. Basic QA steps are included within the
#' function. More detail about these QA steps can be seen in the
#' vignette and user manual attached to this package.
#' Additionally, an option to remove ambiguous assignments based on choosing
#' the formula with the fewest number of heteroatoms (Ohno and Ohno, 2013), and
#' the CH4 vs O replacement nominal mass series (Koch et al., 2007) is included.
#' Most of these QA parameters can be adjusted from outside the function.
#'
#' Positive mode odd-electron formula assignment is possible using the POEx
#' parameter. If set to 1 positive mode odd-electron ionized masses can be
#' assigned.
#'
#' The Ox term sets the number of loops performed by the CHOFIT core by defining the
#' maximum oxygen allowed to be looked for. Three oxygen is equal to one loop due to
#' the low mass moiety used in CHOFIT (Perdue and Green 2015). The default setting is
#' 30 oxygen, so it is limited to 10 loops. It can be increased if necessary.
#'
#' The input dataframe should have two columns, the first column should contain the
#' relative abundance or intensity for the ion mass and the second column should contain
#' the ion mass, either positive or negative mode. The function can also handle an additional
#' input dataframe containing the masses identified as istopic by the IsoFiltR function. The
#' outputs of that function can be directly put into this function as peaks and isopeaks,
#' doing so will allow likely 13C isotope masses to be match to assigned monoisotopic masses.
#' The inclusion of an isotopic mass list is not required, however.
#'
#' @param peaks data frame, Monoisotopic Masses
#' @param isopeaks data frame, Isotopic Masses. Default is "none"
#' @param ionMode string: ("neg", "pos")
#' @param lowMW numeric:
#' Sets the lower limit of molecular mass to be assigned. Default is 100.
#' @param highMW numeric:
#' Sets the upper limit of molecular mass to be assigned. Default is 1000.
#' @param POEx numeric:
#' If set to 1 and ionMode is positive, positive mode odd electron ions can be assigned.
#' Default is 0
#' @param NOEx numeric:
#' If set to 1 and ionMode is negative, negative mode odd electron ions can be assigned.
#' Default is 0
#' @param Nx numeric:
#' Sets the maximum allowable number of Nitrogen 14 to be used in assignment. Default is 0.
#' @param Sx numeric:
#' Sets the maximum allowable number of Sulfur 32 to be used in assignment. Default is 0.
#' @param Px numeric:
#' Sets the maximum allowable number of Phosphorus 31 to be used in assignment. Default is 0.
#' @param Ex numeric:
#' Sets the amount of Carbon 13 to be used in assignment. Default is 0.
#' @param S34x numeric:
#' Sets the amount of Sulfur 34 to be used in assignment. Default is 0.
#' @param N15x numeric:
#' Sets the amount of Nitrogen 15 to be used in assignment. Default is 0.
#' @param Dx numeric:
#' Sets the amount of Deuterium to be used in assignment. Default is 0.
#' @param Clx numeric:
#' Sets the amount of Chlorine to be used in assignment. Default is 0.
#' @param Fx numeric:
#' Sets the amount of Fluorine to be used in assignment. Default is 0.
#' @param Cl37x numeric:
#' Sets the amount of Chlorine 37 to be used in assignment. Default is 0.
#' @param Brx numeric:
#' Sets the amount of Bromine 79 to be used in assignment. Default is 0.
#' @param Br81x numeric:
#' Sets the amount of Bromine 81 to be used in assignment. Default is 0.
#' @param Ix numeric:
#' Sets the amount of Iodine 127 to be used in assignment. Default is 0.
#' @param Mx numeric:
#' Sets the amount of Sodium adduct to be used in assignment. Default is 0.
#' @param NH4x numeric:
#' Sets the amount of Ammonium adduct to be used in assignment. Default is 0.
#' @param Zx numeric:
#' Sets the amount of charge to be used in assignment. Default is 1.
#' @param Sval numeric:
#' Sets the valence of Sulfur. Default is 2.
#' @param Nval numeric:
#' Sets the valence of Nitrogen. Default is 3.
#' @param S34val numeric:
#' Sets the valence of Sulfur 34. Default is 2.
#' @param N15val numeric:
#' Sets the valence of Nitrogen 15. Default is 3.
#' @param Pval numeric:
#' Sets the valence of Phosphorus. Default is 5.
#' @param Ox numeric:
#' Ox sets the maximum number of oxygen looked for in the CHOFIT core, it limits the number of loops performed.
#' @param ppm_err numeric:
#' ppm_err parameter sets the error tolerance (ppm) for formula assignment. Default is 3.
#' @param iso_err numeric:
#' iso_err parameter sets the error tolerance (ppm) for polyisotope matching. Default is 3.
#' @param SN numeric:
#' SN parameter set the signal to noise cut for formula assignment. Default is 0.
#' @param O_Cmin numeric:
#' O_Cmin parameter sets the minimum allowed oxygen to carbon ratio. Default is 0.
#' @param O_Cmax numeric:
#' The O_Cmax parameter sets the upper limit for oxygen to carbon ratio. Default is 2.5.
#' @param H_Cmin numeric:
#' H_Cmin parameter sets lower limit for hydrogen to carbon ratio. Default is 0.1.
#' @param H_Cmax numeric:
#' H_Cmax parameter sets upper limit for hydrogen to carbon ratio for assigned formulas.
#' Default is 3.
#' @param DBEOmin numeric:
#' DBEOmin parameter sets lower limit for DBE minus oxygen QA parameter. Default is -13.
#' @param DBEOmax numeric:
#' DBEOmax parameter sets upper limit for DBE minus oxygen QA parameter. Default is 13.
#' @param Omin numeric:
#' Omin parameter sets lower limit for oxygen number in assigned formula. Default is 0.
#' @param max_def numeric:
#' Value for upper limit of mass defect for using floor() instead of round() for KMD. Default is 0.9
#' @param min_def numeric:
#' Value for lower limit of mass defect for using floor() instead of round() for KMD. Default is 0.5
#' @param HetCut character:
#' HetCut turns on or off the high heteroatom QA parameter. Default is "off"
#' @param NMScut character:
#' NMScut turns on or off the nominal mass series QA parameter. Default is "on".
#' @param DeNovo numeric:
#' DeNovo sets the de novo cut point for the data. Default is 300.
#' @param nLoop numeric:
#' nLoop sets the number of times the KMD and z* formula extension loops. Default is 5.
#' @param SulfCheck character:
#' Turns on or off the sulfur isotope check QA parameter. Default is "on".
#' @param Ambig character:
#' Turns on or off increased ambiguity for assignments. Default is "off".
#' @param MSMS character:
#' Turns on or off CH2 KMD prescreening before initial assignment. Default is "off".
#' @param S34_abund numeric:
#' Sets the maximum 34S/32S isotope ratio for isotope matching. Default is 30.
#' @param C13_abund numeric:
#' Sets the maximum 13C/12C isotope ratio for isotope matching. Default is 60.
#' @param N3corr character:
#' Turns on or off correction of N3OS monoisotopic assignments to 13C assignment. Default is "on"
#' @return list(Unambig = Unambig, Ambig = Ambigout, None = unassigned, MSAssign = MZ,
#'          Error = Error, MSgroups = MZgroups, VK = VK)
#'
#'   Unambig - data frame containing unambiguous assignments
#'   Ambig - data frame containing ambiguous assignments
#'   None - data frame containing unassigned masses
#'   MSAssign - ggplot of mass spectrum highlighting assigned/unassigned
#'   Error - ggplot of the Error vs. m/z
#'   MSgroups - ggplot of mass spectrum colored by molecular group
#'   VK - ggplot of van Krevelen plot, colored by molecular group
#'
#'
#' @examples
#' MFAssign_RMD(peaks = Mono_df, isopeaks = Iso_df, "neg", lowMW = 200, highMW = 700)
#' MFAssign_RMD(peaks = Mono_df, isopeaks = Iso_df, "neg", lowMW = 100, highMW = 1000, Nx = 3, Sx = 1)
#' @export


MFAssign <- function(peaks, isopeaks = "none", ionMode, lowMW=100,highMW=1000, POEx = 0, NOEx = 0, Nx=0,Sx=0, Px=0, S34x=0,
                         N15x=0, Dx=0,Ex=0, Clx=0, Cl37x=0, Fx = 0, Brx = 0, Br81x = 0, Ix = 0, Mx=0, NH4x=0, Zx=1,
                         Sval = 2, Nval = 3, S34val = 2,
                         N15val = 3, Pval = 5, Ox = 30, ppm_err = 3, iso_err = 3, SN = 0, O_Cmin = 0,
                         O_Cmax = 2.5, H_Cmin = 0.3, H_Cmax = 3, DBEOmin = -13, DBEOmax = 13, Omin = 0, max_def = 0.9,
                         min_def = 0.5, HetCut = "off",NMScut = "on", DeNovo = 300, nLoop = 5, SulfCheck = "on",
                         Ambig = "off", MSMS = "off", S34_abund = 30, C13_abund = 60, N3corr = "on") {

  if(POEx >1) print('WARNING: Positive Odd Electron (POEx) is greater than 1, are you sure that is what you want?')
  if(NOEx >1) print('WARNING: Positive Odd Electron (NOEx) is greater than 1, are you sure that is what you want?')

  if(ionMode != "pos" & ionMode != "neg") print("WARNING: ionMode should be 'pos' or 'neg' ")

  if(Nx > 5 | Sx > 5|Px >5|S34x>5|N15x >5|Dx>5|Ex>5|Clx > 5|Cl37x>5|Mx>5|NH4x>5|Fx > 5)
    print("WARNING: One or more heteroatoms are set greater than 5, this will cause the function to perform more slowly.")

  if(Ox !=30) print("WARNING: Ox is not at its default value, this will cause the core formula algorithm to perform additional
                    or fewer loops, are you sure you want it changed?")

  if(ppm_err > 3) print("WARNING: The maximum allowed error (ppm_err) is greater than 3, is this what you want?")

  # Constants
  components <- factor(c("C", "H", "O", "N", "S", "P", "Cl", "Fl","E", "S34", "N15", "D", "Cl37", "Br", "Br81",
                         "I","M", "NH4", "POE","NOE", "Z"),

                levels=c("C", "H", "O", "N", "S", "P", "Cl", "Fl",  "E", "S34", "N15", "D", "Cl37","Br", "Br81",
                          "I","M", "NH4", "POE", "NOE", "Z"))

  numComps <- length(components)
  proton = 1.00727645216
  electron =  0.000548597
  fitMode <- "ppm"
  maxErr <- ppm_err
  numDigits <- 6

  # Initialize data for analysis
  totFormulae <- 0

  #############
  #LCMS adjustment
  cols <- ncol(peaks)
  ifelse(isopeaks != "none", cols2 <- ncol(isopeaks), cols2 <- 0)
  if(cols == 3){
    if(cols == 3) {monoSave <- peaks[c(1,2,3)]}
    ifelse(cols2 == 4 & isopeaks != "none",isoSave <- isopeaks[c(1,2,3)], print("No Iso List Included"))
    ifelse(cols2 == 0 & isopeaks == "none",isoSave <- data.frame(mass = -42, abund = -42, RT = -42), isoSave <- isoSave)
    names(isoSave)[1] <- "exp_mass"
    names(isoSave)[2] <- "abundance"
    names(isoSave)[3] <- "RT"
    names(monoSave)[1] <- "exp_mass"
    names(monoSave)[2] <- "abundance"
    names(monoSave)[3] <- "RT"
  }
  #############

  peaks <- peaks[c(2,1)]

  names(peaks)[2] <- "mass"
  names(peaks)[1] <- "RA"
  rawpeaks <- peaks
  peaks <- peaks[peaks$mass >= lowMW,]
  peaks <- peaks[peaks$mass <= highMW,]

  ifelse(isopeaks != "none", isopeaks2 <- isopeaks, isopeaks2 <- data.frame(x=0,y=0,Tag=0))

  if(cols2 == 4){isopeaks2 <- isopeaks2[c(2,1,4)]}  #LC Change
  if(cols2 == 3){ isopeaks2 <- isopeaks2[c(2,1,3)]} #LC Change


  names(isopeaks2)[2] <- "Iso_mass"
  names(isopeaks2)[1] <- "Iso_RA"
  names(isopeaks2)[3] <- "Tag"
  isopeaks2 <- isopeaks2[isopeaks2$Iso_RA > SN,]

  peaksAll <- peaks

  #Splitting the peaks into those below S/N and those above
  SNcut <- peaks[peaks$RA < SN,]
  peaks <- peaks[peaks$RA >= SN,]
  peaksAll2 <- peaks
  #Splitting peaks into those below the DeNovo cut and those above
  DeNovocut <- peaks[peaks$mass > DeNovo,]
  peaks <- peaks[peaks$mass <= DeNovo,]

  Ambigcheck <- Ambig #Renames the input term so it doesn't interfere with other things

  ################################# Inital Kendrick Series implementation
  if(MSMS == "off"){
    peaks$KM <- peaks$mass* (14/14.01565)

    peaks$KMD <- ifelse(abs(floor(peaks$mass)-peaks$mass) >= min_def & abs(floor(peaks$mass)-peaks$mass) <= max_def,
    peaks$KMD <- floor(peaks$mass)-peaks$KM, peaks$KMD <- round(peaks$mass)-peaks$KM) #New 1/6/20

    peaks$zstar <- ifelse(abs(floor(peaks$mass)-peaks$mass) >= min_def & abs(floor(peaks$mass)-peaks$mass) <= max_def,
                          peaks$zstar <- floor(peaks$mass)%%14 - 14, peaks$zstar <- round(peaks$mass)%%14 - 14) #New 1/6/20

    peaks$KMDTest <- round(peaks$KMD, 3)

    Test <-  dplyr::group_by(peaks, KMDTest, zstar)
    Test <- dplyr::mutate(Test, CH2_num = round(mass - min(mass))/14)
    peaksend <- dplyr::filter(Test, CH2_num !=0 & CH2_num != (min(CH2_num[CH2_num!=min(CH2_num)])+1)&
                                CH2_num != (min(CH2_num[CH2_num!=min(CH2_num)])+3))

    names(peaksend)[1] <- "RA"
    names(peaksend)[2] <- "Exp_mass"
    peaksend <- peaksend[c(1,2)]

    peaks <- dplyr::filter(Test, CH2_num ==0 | CH2_num == (min(CH2_num[CH2_num!=min(CH2_num)])+1) |
                             CH2_num == (min(CH2_num[CH2_num!=min(CH2_num)])+3))

    peaks <- data.frame(RA = peaks[1], mass = peaks[2])
  }else{
    peaks$KM <- peaks$mass* (14/14.01565)
    peaks$KMD <- ifelse(abs(floor(peaks$mass)-peaks$mass) >= min_def & abs(floor(peaks$mass)-peaks$mass) <= max_def,
                        peaks$KMD <- floor(peaks$mass)-peaks$KM, peaks$KMD <- round(peaks$mass)-peaks$KM) #New 1/6/20

    peaks$zstar <- ifelse(abs(floor(peaks$mass)-peaks$mass) >= min_def & abs(floor(peaks$mass)-peaks$mass) <= max_def,
                          peaks$zstar <- floor(peaks$mass)%%14 - 14, peaks$zstar <- round(peaks$mass)%%14 - 14) #New 1/6/20
    peaks$KMDTest <- round(peaks$KMD, 3)

    Test <-  dplyr::group_by(peaks, KMDTest, zstar)
    Test <- dplyr::mutate(Test, CH2_num = round(mass - min(mass))/14)
    peaksend <- dplyr::filter(Test, CH2_num !=0 & CH2_num != (min(CH2_num[CH2_num!=min(CH2_num)])+1)&
                                CH2_num != (min(CH2_num[CH2_num!=min(CH2_num)])+3))

    names(peaksend)[1] <- "RA_CH2"
    names(peaksend)[2] <- "mass_CH2"
    peaksend <- peaksend[c(1,2,5,6,7)]
    #
    peaks <- dplyr::filter(Test, RA > 0)
    #
    peaks <- data.frame(RA = peaks[1], mass = peaks[2])

  }
  #################################
  Dummy <- data.frame(RA = c(-42,-42), mass = c(421.1147, 423.1293))
  peaks <- dplyr::bind_rows(Dummy, peaks)

  records <- vector('list')

  records2 <- vector('list')



  #p <-1

  pb <- txtProgressBar(min = 0, max = nrow(peaks), style = 3)

  for (p in 1:nrow(peaks)) {

    fit       <- FALSE
    formulaOK <- FALSE
    coreXEM   <- 0
    coreRNM   <- 0
    ionEM     <- peaks[p,2]
    RA        <- peaks[p,1]
    xemErr    <- 0

    moles <- vector(mode="numeric", length=numComps)
    loop <- vector(mode="numeric", length=numComps)
    #}#
    # Convert the (presumably) single-charged ion to a molecule

    if (ionMode=="neg") {
      exactEM <- ionEM + proton
    } else {
      exactEM <- ionEM - proton
    }



    # Check that exactEM is within bounds of LowMW and HighMW bounds
    if ((round(exactEM) >= lowMW) & (round(exactEM) <= highMW+1)) {

      # Start looping through components
      loop[CompFactorToInt2("C")] <- 1 #LowMoles("C")
      loop[CompFactorToInt2("H")] <- 4 #LowMoles("H")
      loop[CompFactorToInt2("O")] <- 0 #LowMoles("O")
      loop[CompFactorToInt2("Z")] <- 1 #LowMoles("Z")

      repeat {
        if (loop[CompFactorToInt2("Z")] > 2) {
          exactEM = exactEM*(loop[CompFactorToInt2("Z")] /
                               loop[CompFactorToInt2("Z")] - 1)
        }

        loop[CompFactorToInt2("NOE")] <- 0 #LowMoles("NOE")
        repeat {

          loop[CompFactorToInt2("POE")] <- 0 #LowMoles("POE")
          repeat {

            loop[CompFactorToInt2("NH4")] <- 0 #LowMoles("NH4")
            repeat {

              loop[CompFactorToInt2("M")] <- 0 #LowMoles("M")
              repeat {

                loop[CompFactorToInt2("I")] <- 0 #LowMoles("I")
                repeat {
                  loop[CompFactorToInt2("Br81")] <- 0 #LowMoles("Br81")
                  repeat {
                    loop[CompFactorToInt2("Br")] <- 0 #LowMoles("Br")
                    repeat {


                loop[CompFactorToInt2("Cl37")] <- 0 #LowMoles("Cl37")
                repeat {

                  loop[CompFactorToInt2("Fl")] <- 0 #LowMoles("Fl")
                  repeat {

                    loop[CompFactorToInt2("Cl")] <- 0 #LowMoles("Cl")
                    repeat {

                      loop[CompFactorToInt2("D")] <- 0 #LowMoles("D")
                      repeat {

                        loop[CompFactorToInt2("N15")] <- 0 #LowMoles("N15")
                        repeat {

                          loop[CompFactorToInt2("S34")] <- 0 #LowMoles("S34")
                          repeat {

                            loop[CompFactorToInt2("P")] <- 0 #LowMoles("P")
                            repeat {

                              loop[CompFactorToInt2("S")] <- 0 #LowMoles("S")
                              repeat {

                                loop[CompFactorToInt2("N")] <- 0 #LowMoles("N")
                                repeat {

                                  loop[CompFactorToInt2("E")] <- 0 #LowMoles("E")
                                  repeat {

                                    # strip the exact mass of all loop constituents to give CHO
                                    # core of formula
                                    coreXEM <- exactEM
records <- vector("list")

                                    ###Make sure this shouldn't start at E
                                    for (step in CompFactorToInt2("N"):CompFactorToInt2("NOE")) {
                                      coreXEM = coreXEM - unlist(loop[step])*EM2(CompIntToFactor2(step))

                                    }

                                    if (coreXEM >= 16.0313) {

                                      coreRNM <- round(coreXEM)
                                      formulaOK <- FALSE

                                      if (Even(coreRNM) == TRUE) {  ##Change on 4/16/20
                                        env <- environment()
                                        #if( exists(env) ) stop('Check mass list to make sure the masses are correct, or add more masses')
                                        FindCoreFormulae2_Halo(env)

                                      } #else {
                                      #   if (coreXEM >= 436.5008) {
                                      #     env2 <- environment()
                                      #     FindCoreFormulae2(env2)
                                      #
                                      #   }
                                      # }

#coreRNM <- 7


                                    } else {
                                      formulaOK <- FALSE
                                    }

                                    if(length(unlist(env$records[4])) > 0){

                                    env$records2[[length(env$records2) +1]] = env$records

                                    } else{
                                      env$records2 = env$records2
                                    }

                                    #For List option

                                    # ******##############################
                                    # #####################################
                                    loop[CompFactorToInt2("E")] <- unlist(loop[CompFactorToInt2("E")]) + 1
                                    if (loop[CompFactorToInt2("E")] > HighMoles("E", E=Ex)) {
                                      break
                                    }
                                  } # E loop

                                  loop[CompFactorToInt2("N")] <- unlist(loop[CompFactorToInt2("N")]) + 1
                                  if (loop[CompFactorToInt2("N")] > HighMoles("N",N=Nx)) {
                                    break
                                  }
                                } # N loop

                                loop[CompFactorToInt2("S")] <- unlist(loop[CompFactorToInt2("S")]) + 1
                                if (loop[CompFactorToInt2("S")] > HighMoles("S",S=Sx)) {
                                  break
                                }
                              } # S loop

                              loop[CompFactorToInt2("P")] <- unlist(loop[CompFactorToInt2("P")]) + 1
                              if (loop[CompFactorToInt2("P")] > HighMoles("P",P=Px)) {
                                break
                              }
                            } # P loop

                            loop[CompFactorToInt2("S34")] <- unlist(loop[CompFactorToInt2("S34")]) + 1
                            if (loop[CompFactorToInt2("S34")] > HighMoles("S34",S34=S34x)) {
                              break
                            }
                          } # S34 loop

                          loop[CompFactorToInt2("N15")] <- unlist(loop[CompFactorToInt2("N15")]) + 1
                          if (loop[CompFactorToInt2("N15")] > HighMoles("N15",N15=N15x)) {
                            break
                          }
                        } # N15 loop

                        loop[CompFactorToInt2("D")] <- unlist(loop[CompFactorToInt2("D")]) + 1
                        if (loop[CompFactorToInt2("D")] > HighMoles("D",D=Dx)) {
                          break
                        }
                      } # D loop

                      loop[CompFactorToInt2("Cl")] <- unlist(loop[CompFactorToInt2("Cl")]) + 1
                      if (loop[CompFactorToInt2("Cl")] > HighMoles("Cl",Cl=Clx)) {
                        break
                      }
                    } # Cl loop

                    loop[CompFactorToInt2("Fl")] <- unlist(loop[CompFactorToInt2("Fl")]) + 1
                    if (loop[CompFactorToInt2("Fl")] > HighMoles("Fl",Fl=Fx)) {
                      break
                    }
                  } # Fl loop



                  loop[CompFactorToInt2("Cl37")] <- unlist(loop[CompFactorToInt2("Cl37")]) + 1
                  if (loop[CompFactorToInt2("Cl37")] > HighMoles("Cl37",Cl37=Cl37x)) {
                    break
                  }
                }

                loop[CompFactorToInt2("Br")] <- unlist(loop[CompFactorToInt2("Br")]) + 1
                if (loop[CompFactorToInt2("Br")] > HighMoles2("Br",Br=Brx)) {
                  break
                }
                    }

                    loop[CompFactorToInt2("Br81")] <- unlist(loop[CompFactorToInt2("Br81")]) + 1
                    if (loop[CompFactorToInt2("Br81")] > HighMoles2("Br81",Br81=Br81x)) {
                      break
                    }
                  }

                  loop[CompFactorToInt2("I")] <- unlist(loop[CompFactorToInt2("I")]) + 1
                  if (loop[CompFactorToInt2("I")] > HighMoles2("I",I=Ix)) {
                    break
                  }
                }

                loop[CompFactorToInt2("M")] <- unlist(loop[CompFactorToInt2("M")]) + 1
                if (loop[CompFactorToInt2("M")] > HighMoles("M", M=Mx)) {
                  break
                }
              } # M Loop

              loop[CompFactorToInt2("NH4")] <- unlist(loop[CompFactorToInt2("NH4")]) + 1
              if (loop[CompFactorToInt2("NH4")] > HighMoles("NH4", NH4=NH4x)) {
                break
              }
            }
            loop[CompFactorToInt2("POE")] <- unlist(loop[CompFactorToInt2("POE")]) + 1
            if (loop[CompFactorToInt2("POE")] > HighMoles("POE", POE=POEx)) {
              break
            }
          } # Positive Odd Electron loop

          loop[CompFactorToInt2("NOE")] <- unlist(loop[CompFactorToInt2("NOE")]) + 1
          if (loop[CompFactorToInt2("NOE")] > HighMoles("NOE", NOE=NOEx)) {
            break
          }
        } # Negative Odd Electron loop

        if (fit) {
          loop[which(components=="Z")] <- HighMoles("Z",Z=Zx)
        }
        loop[which(components=="Z")] <- unlist(loop[which(components=="Z")]) + 1
        if ( loop[which(components=="Z")] > HighMoles("Z",Z=Zx)) {
          break
        }
      } # Z loop

    } # exactEM in lowMW and highMW bounds

    if (!fit) {
      exactEM <- exactEM/HighMoles("Z",Z=Zx)
      cem <- 0
      cnm <- 0
      xemErr <- 0

    }
    setTxtProgressBar(pb, p)
  } # for each peak in data


  recordsdf <- data.frame((do.call('rbind.data.frame', records2)))


  recordsdf$H <- recordsdf$H + recordsdf$N + recordsdf$N15 + recordsdf$P +
    2*recordsdf$POE + recordsdf$Cl + recordsdf$Cl37 + recordsdf$Fl - 2*recordsdf$NOE + recordsdf$Br +
    recordsdf$Br81 + recordsdf$I

  records1 <- recordsdf[recordsdf$C > 1 & recordsdf$O >=0 & recordsdf$H > 0 & recordsdf$RA >=0,]

  records1$mode <- ionMode

  df1 <- records1[records1$mode == "pos" & records1$M > 0 & records1$POE == 0 & records1$NH4 == 0,]
  df1$Neutral_mass <- df1$Exp_mass - df1$M * 22.989221

  df1x <- records1[records1$mode == "pos" & records1$M == 0 & records1$POE == 0 & records1$NH4 == 1,]
  df1x$Neutral_mass <- df1x$Exp_mass - df1x$NH4 * 18.033823

  df2 <- records1[records1$mode == "pos" & records1$M == 0 & records1$POE == 0& records1$NH4 == 0,]
  df2$Neutral_mass <- df2$Exp_mass - 1.00727645216

  df3 <- records1[records1$mode == "neg" & records1$NOE == 0,]
  df3$Neutral_mass <- df3$Exp_mass + 1.00727645216

  df4 <- records1[records1$mode == "neg" & records1$NOE == 1,]
  df4$Neutral_mass <- df4$Exp_mass - 0.000548597

  df5 <- records1[records1$mode == "pos" & records1$M == 0 & records1$POE == 1& records1$NH4 == 0,]
  df5$Neutral_mass <- df5$Exp_mass + 0.000548597

  records1 <- rbind(df1, df1x, df2, df3, df4, df5)

  records1 <- records1[-29]   ##CHECK THIS 6/2/20

  records1 <- dplyr::mutate(records1, O_C = O/(C+E), H_C =H/(C+E),

                            theor_mass1 = EM2("C") * C + EM2("H") * H + EM2("O") * O + N * EM2("N14") + S * EM2("S") + P * EM2("P31") +
                              Cl * EM2("Cl35") + Fl * EM2("Fl19") + E * EM2("E") + S34 * EM2("S34") + Cl37 * EM2("Cl37m") + N15 * EM2("N15H") +
                              D * EM2("D") + M * EM2("M") + NH4 * EM2("NH4+") + NOE * electron - POE * electron + Br * EM2("Br79") +
                              Br81 * EM2("Br81m") + I * EM2("I127"),

                            theor_mass = EM2("C") * C + EM2("H") * H + EM2("O") * O + N * EM2("N14") + S * EM2("S") + P * EM2("P31") +
                              Cl * EM2("Cl35") + Fl * EM2("Fl19") + E * EM2("E") + S34 * EM2("S34") + Cl37 * EM2("Cl37m") + N15 * EM2("N15H") +
                              D * EM2("D") + Br * EM2("Br79") +
                              Br81 * EM2("Br81m") + I * EM2("I127"), #

                            C = C + E,

                            DBE = C - 0.5 * (H + Cl + Cl37 + Fl + Br + Br81 + I) + 0.5 * (N + N15 + P) + 1,

                            err_ppm = ((Neutral_mass - theor_mass) / Neutral_mass * 10^6),

                            AE_ppm = round(abs((Neutral_mass - theor_mass) / Neutral_mass * 10^6),2),

                            #NM = floor(Exp_mass),

                            KM = Exp_mass * (14 / 14.01565), #KMD = NM - KM,

                            max_LA = theor_mass1 / 13, actual_LA = ((C - E) + N + S + O + E + S34 + P + Cl +Fl+Cl37+N15
                                                                    +Br + Br81 + I) ,

                            rule_13= round(actual_LA/max_LA,1),

                            Senior1 = H + P + N + Cl + Cl37 + N15+Fl + Br + Br81 + I  ,

                            STest = S + S34, ClTest = Cl + Cl37, NTest = N + N15, DBEO = DBE-O, BrTest = Br + Br81,

                            max_H = C * 2 + 2, H_test = round(H / max_H, 1),

                            Senior2 = Pval + Nval + N15val + Valence("H")  +
                              Valence("Cl") + Valence("Cl37")+Valence("Fl") + Valence2("Br") + Valence2("Br81") +
                              Valence2("I"),

                            Senior3Atom = C + H + O + N + S + P + N15 + Cl + Cl37 + S34 + Br + Br81 + I,

                            Senior3Val = C*Valence("C") + H*Valence("H") + O*Valence("O") + N*Nval +
                              S*Sval + P*Pval + S34*S34val + Fl*Valence("Fl")+
                              N15*N15val + Cl*Valence("Cl") + Cl37*Valence("Cl37") + Br *Valence2("Br")+
                              Br81 *Valence2("Br81") + I * Valence2("I")

  )

  records1$NM <- ifelse(abs(floor(records1$Exp_mass)-records1$Exp_mass) >= min_def & abs(floor(records1$Exp_mass)-records1$Exp_mass) <= max_def,
                        records1$NM <- floor(records1$Exp_mass), records1$NM <- round(records1$Exp_mass)) #New 1/6/20

  records1$KMD <- ifelse(abs(floor(records1$Exp_mass)-records1$Exp_mass) >= min_def & abs(floor(records1$Exp_mass)-records1$Exp_mass) <= max_def,
                         records1$KMD <- floor(records1$Exp_mass)-records1$KM, records1$KMD <- round(records1$Exp_mass)-records1$KM) #New 1/6/20



  #recordssave <- records1


  records1 <- dplyr::filter(records1, C>0& H>0&O>=Omin& H >= D)
  records1 <- unique(records1)
  #check <- records1 %>% filter(AE_ppm <= 3 & S34 == 1)
  records1 <- dplyr::filter(records1, O_C < O_Cmax & H_C <= H_Cmax & H_C > H_Cmin & O_C >= O_Cmin &

                              DBEO >= DBEOmin & DBEO <= DBEOmax &

                              #!((ClTest) > HighMoles("Cl",Cl=Clx)|(ClTest) > HighMoles("Cl37",Cl37=Cl37x)) &
                              #!((STest) > HighMoles("S",S=Sx)|(STest) > HighMoles("S34",S34=S34x)) &
                              #!((NTest) > HighMoles("N",N=Nx)|(NTest) > HighMoles("N15",N15=N15x)) &
                              (STest <= Sx + S34x) &
                              (ClTest <= Clx + Cl37x) &
                              (NTest <= Nx + N15x) &
                              (BrTest <= Brx + Br81x) &

                              H_test <= 1 & rule_13 <= 1 &

                              AE_ppm <= ppm_err &

                              Even(Senior1)==TRUE & DBE >= 0 & DBE <= round(0.9 * (C + N)) &

                              O <= 2 * C + 3 * (N+N15) + 4 * P + 4 * (S+S34)&

                              O >= Omin&

                              #RA >= 0 &

                              Senior2 >= 2*Valence("C") &

                              Senior3Val >= (2*Senior3Atom - 1)
  )

  #check <- records1 %>% distinct(Exp_mass, .keep_all = TRUE)
  #recordssave <- records1
  #NOEx working perfect to this point.
  #Sulf <- records1[records1$S == 1,]
  ##################
  ###S34 isotope check QA
  check <-0
  ifelse(isopeaks != "none" , check <- 1, print(""))
  #ifelse(cols2 == 4 & isopeaks != "none",isoSave <- isopeaks[c(1,2,3)], print("No Iso List Included"))

  if(check == 1 & SulfCheck == "on"){
    #The next line was changed for the new isotoping
    SIso <- isopeaks2[isopeaks2$Tag == "S34"|isopeaks2$Tag == "C13_S34"|isopeaks2$Tag == "2C13_S34",]
    SIso <- unlist(SIso[2])
    recordsS <- records1[records1$S > 0,]
    recordsS <- unlist(recordsS[3])
    rest <- records1[records1$S == 0,]

    Sulf <- expand.grid(recordsS, SIso)
    names(Sulf)[1] <- "Exp_mass"
    names(Sulf)[2] <- "Iso_mass"

    Sulf$S_mass <- Sulf$Exp_mass + 1.995797   ##New as of 6/19/2019
    Sulf$err <- ((Sulf$S_mass - Sulf$Iso_mass)/Sulf$S_mass) * 10^6   ##New as of 6/19/2019

    Sulf <- Sulf[abs(Sulf$err) < iso_err,]   ##New as of 6/19/2019

    # Sulf$mdiff <- Sulf$Exp_mass - Sulf$Iso_mass
    # Sulf <- Sulf[Sulf$mdiff > -2 & Sulf$mdiff < -1.98,]
    #
    #
    # Sulf$KM2 <- Sulf$Iso_mass * (2 / 1.995797)
    # Sulf$KMD2 <- round((round(Sulf$Iso_mass) - Sulf$KM2),3)
    # Sulf$KM <- Sulf$Exp_mass * (2 / 1.995797)
    # Sulf$KMD <- round((round(Sulf$Exp_mass) - Sulf$KM),3)
    #
    # Sulf$KMDdiff <- Sulf$KMD - Sulf$KMD2
    #
    # Sulf <- Sulf[abs(Sulf$KMDdiff) < 0.002,]

    Sulf <- Sulf[1]
    Sulfdata <- records1[records1$S > 0,]
    Sulfur <- merge(Sulf, Sulfdata, by.x = "Exp_mass", by.y = "Exp_mass")


    restalign <- rest[c(3,1)]
    restalign <- merge(Sulfur, restalign, by.x = "Exp_mass", by.y = "Exp_mass", all = T)
    restmass <- restalign[is.na(restalign$RA.x),]
    restmass <- restmass[c(1,53)]  ##CHECK THIS POINT 6/2/20
    restfinal <- merge(rest, restmass, by.x = "Exp_mass", by.y = "Exp_mass", all = T)
    restfinal <- restfinal[!is.na(restfinal$RA.y),]
    restfinal <- restfinal[-53]  ##CHECK THIS POINT 6/2/20

    records1 <- rbind(restfinal, Sulfur)
    records1 <- records1[c(2,3,1, 4:52)]  ##CHECK THIS POINT 6/2/20
  }
  ###################################################
  #records1 <- recordssave
  ###Series analysis
  ##Determining Ambiguity
  records1$Test <- paste(records1$Exp_mass, records1$RA, sep = "_")
  records1$num <- 1:nrow(records1)
  records1$dups <- duplicated(records1$Test)
  records1 <-records1[order(-records1$num),]
  records1$dups2 <- duplicated(records1$Test)
  records1 <- records1[-c(53)]  #removes Test  ##CHECK THIS POINT 6/2/20


  Unambig <- records1[records1$dups == FALSE & records1$dups2 == FALSE,]
  Unamatch <- Unambig[c(1,3, 4)]

  Ambig <- records1[records1$dups == TRUE | records1$dups2 == TRUE,]
  Ambig <- Ambig[c(1,3)]
  names(peaksAll2)[2] <- "Exp_mass"
  Ambig <- rbind(Ambig, peaksAll2)

  Ambig <- merge(Ambig, Unamatch, by.x = c("RA", "Exp_mass"), by.y = c("RA", "Exp_mass"), all = T)
  Ambig <- Ambig[is.na(Ambig$C),]
  Ambig <- Ambig[-c(3)]
  Ambig <- unique(Ambig)
  ##Unambiguous formulas prep

  Unambig <- Unambig[c("RA", "Exp_mass", "C", "H", "O", "N", "S", "P", "E", "S34", "N15", "D",
                       "Cl", "Fl", "Cl37", "Br", "Br81", "I", "M", "NH4", "POE", "NOE", "Z")]
 ####
   Unambig$NM <- ifelse(abs(floor(Unambig$Exp_mass)-Unambig$Exp_mass) >= min_def & abs(floor(Unambig$Exp_mass)-Unambig$Exp_mass) <= max_def,
                      Unambig$NM <- floor(Unambig$Exp_mass), Unambig$NM <- round(Unambig$Exp_mass)) #New 1/6/20
 ####

   peaks$zstar <- ifelse(abs(floor(peaks$mass)-peaks$mass) >= min_def & abs(floor(peaks$mass)-peaks$mass) <= max_def,
                         peaks$zstar <- floor(peaks$mass)%%14 - 14, peaks$zstar <- round(peaks$mass)%%14 - 14) #New 1/6/20


  Unambig$NM <- floor(Unambig$Exp_mass)

  Unambig$KM_CH2 <- Unambig$Exp_mass * (14/14.01565)
  Unambig$KMD_CH2 <- round(Unambig$NM - Unambig$KM_CH2, 3)
  ###
  Unambig$z_CH2 <- ifelse(abs(floor(Unambig$Exp_mass)-Unambig$Exp_mass) >= min_def & abs(floor(Unambig$Exp_mass)-Unambig$Exp_mass) <= max_def,
                        Unambig$z_CH2 <- floor(Unambig$Exp_mass)%%14 - 14, Unambig$z_CH2 <- round(Unambig$Exp_mass)%%14 - 14) #New 1/6/20
  ###
  Unambig$KM_O <- Unambig$Exp_mass * (16/15.9949146223)
  Unambig$KMD_O <- round(Unambig$NM - Unambig$KM_O, 3)
  ###
  Unambig$z_O <- ifelse(abs(floor(Unambig$Exp_mass)-Unambig$Exp_mass) >= min_def & abs(floor(Unambig$Exp_mass)-Unambig$Exp_mass) <= max_def,
                          Unambig$z_O <- floor(Unambig$Exp_mass)%%16 - 16, Unambig$z_O <- round(Unambig$Exp_mass)%%16 - 16) #New 1/6/20
  ###

  Unambig$KM_H2 <- Unambig$Exp_mass * (2/2.01565)
  Unambig$KMD_H2 <- round(Unambig$NM - Unambig$KM_H2, 3)
  ###
  Unambig$z_H2 <- ifelse(abs(floor(Unambig$Exp_mass)-Unambig$Exp_mass) >= min_def & abs(floor(Unambig$Exp_mass)-Unambig$Exp_mass) <= max_def,
                          Unambig$z_H2 <- floor(Unambig$Exp_mass)%%2 - 2, Unambig$z_H2 <- round(Unambig$Exp_mass)%%2 - 2) #New 1/6/20
  ###

  Unambig$KM_H2O <- Unambig$Exp_mass * (18/18.01056468)
  Unambig$KMD_H2O <- round(Unambig$NM - Unambig$KM_H2O, 3)
  ###
  Unambig$z_H2O <- ifelse(abs(floor(Unambig$Exp_mass)-Unambig$Exp_mass) >= min_def & abs(floor(Unambig$Exp_mass)-Unambig$Exp_mass) <= max_def,
                          Unambig$z_H2O <- floor(Unambig$Exp_mass)%%18 - 18, Unambig$z_H2O <- round(Unambig$Exp_mass)%%18 - 18) #New 1/6/20
  ###

  Unambig$KM_CH2O <- Unambig$Exp_mass * (30/30.01056468)
  Unambig$KMD_CH2O <- round(Unambig$NM - Unambig$KM_CH2O, 3)
  ###
  Unambig$z_CH2O <- ifelse(abs(floor(Unambig$Exp_mass)-Unambig$Exp_mass) >= min_def & abs(floor(Unambig$Exp_mass)-Unambig$Exp_mass) <= max_def,
                          Unambig$z_CH2O <- floor(Unambig$Exp_mass)%%30 - 30, Unambig$z_CH2O <- round(Unambig$Exp_mass)%%30 - 30) #New 1/6/20
  ###
  #Good to this point

  ##Ambiguous
  Ambig <- Ambig[c(1,2)]
  ####
  Ambig$NM <- ifelse(abs(floor(Ambig$Exp_mass)-Ambig$Exp_mass) >= min_def & abs(floor(Ambig$Exp_mass)-Ambig$Exp_mass) <= max_def,
                       Ambig$NM <- floor(Ambig$Exp_mass), Ambig$NM <- round(Ambig$Exp_mass)) #New 1/6/20
  ####

  Ambig$KM_CH2 <- Ambig$Exp_mass * (14/14.01565)
  Ambig$KMD_CH2 <- round(Ambig$NM - Ambig$KM_CH2, 3)
  ###
  Ambig$z_CH2 <- ifelse(abs(floor(Ambig$Exp_mass)-Ambig$Exp_mass) >= min_def & abs(floor(Ambig$Exp_mass)-Ambig$Exp_mass) <= max_def,
                          Ambig$z_CH2 <- floor(Ambig$Exp_mass)%%14 - 14, Ambig$z_CH2 <- round(Ambig$Exp_mass)%%14 - 14) #New 1/6/20
  ###

  Ambig$KM_O <- Ambig$Exp_mass * (16/15.9949146223)
  Ambig$KMD_O <- round(Ambig$NM - Ambig$KM_O, 3)
  ###
  Ambig$z_O <- ifelse(abs(floor(Ambig$Exp_mass)-Ambig$Exp_mass) >= min_def & abs(floor(Ambig$Exp_mass)-Ambig$Exp_mass) <= max_def,
                        Ambig$z_O <- floor(Ambig$Exp_mass)%%16 - 16, Ambig$z_O <- round(Ambig$Exp_mass)%%16 - 16) #New 1/6/20
  ###

  Ambig$KM_H2 <- Ambig$Exp_mass * (2/2.01565)
  Ambig$KMD_H2 <- round(Ambig$NM - Ambig$KM_H2, 3)
  ###
  Ambig$z_H2 <- ifelse(abs(floor(Ambig$Exp_mass)-Ambig$Exp_mass) >= min_def & abs(floor(Ambig$Exp_mass)-Ambig$Exp_mass) <= max_def,
                         Ambig$z_H2 <- floor(Ambig$Exp_mass)%%2 - 2, Ambig$z_H2 <- round(Ambig$Exp_mass)%%2 - 2) #New 1/6/20
  ###

  Ambig$KM_H2O <- Ambig$Exp_mass * (18/18.01056468)
  Ambig$KMD_H2O <- round(Ambig$NM - Ambig$KM_H2O, 3)
  ###
  Ambig$z_H2O <- ifelse(abs(floor(Ambig$Exp_mass)-Ambig$Exp_mass) >= min_def & abs(floor(Ambig$Exp_mass)-Ambig$Exp_mass) <= max_def,
                          Ambig$z_H2O <- floor(Ambig$Exp_mass)%%18 - 18, Ambig$z_H2O <- round(Ambig$Exp_mass)%%18 - 18) #New 1/6/20
  ###

  Ambig$KM_CH2O <- Ambig$Exp_mass * (30/30.01056468)
  Ambig$KMD_CH2O <- round(Ambig$NM - Ambig$KM_CH2O, 3)
  ###
  Ambig$z_CH2O <- ifelse(abs(floor(Ambig$Exp_mass)-Ambig$Exp_mass) >= min_def & abs(floor(Ambig$Exp_mass)-Ambig$Exp_mass) <= max_def,
                           Ambig$z_CH2O <- floor(Ambig$Exp_mass)%%30 - 30, Ambig$z_CH2O <- round(Ambig$Exp_mass)%%30 - 30) #New 1/6/20
  ###
  #Good to this point
  ###Looping series assignment
  knowndummy <- data.frame(KMD_CH2 = -42, KMD_O = -42, KMD_CH2O = -42, KMD_H2O = -42, KMD_H2 = -42, z_CH2 = -42,
                           z_O = -42,
                           z_CH2O = -42, z_H2O = -42, z_H2 = -42, RA = -42, Exp_mass = -42, C = 4, H = 4, O = 0,
                           N = 0, S = 0,
                           S34 = 0, P = 0, N15 = 0, Cl = 0, Fl = 0, Cl37 = 0, Br = 0, Br81 = 0, I = 0,
                           M = 0, NH4 = 0, Z= 0, POE = 0,
                           NOE = 0,D = 0, E = 0, KM_O = -42,
                           KM_CH2 = -42, KM_CH2O = -42, KM_H2O= -42, KM_H2 = -42, NM = 0)
  known <- Unambig
  known <- rbind(known, knowndummy)

  unknowndummy <- data.frame(KMD_CH2 = -42, KMD_O = -42, KMD_CH2O = -42, KMD_H2O = -42, KMD_H2 = -42, z_CH2 = -42, z_O = -42,
                             z_CH2O = -42, z_H2O = -42, z_H2 = -42, RA = -42, Exp_mass = -42, KM_O = -42,
                             KM_CH2 = -42, KM_CH2O = -42, KM_H2O= -42, KM_H2 = -42, NM = 0)
  unknown <- Ambig
  unknown <- rbind(unknown, unknowndummy)
  Ambigreturn <- unknowndummy

  DummyOut <- data.frame(KMD_CH2 = -42, KMD_O = -42, RA.x = -42,KMD_CH2O = -42,  KMD_H2O = -42, KMD_H2 = -42, z_CH2 = -42, z_O = -42,
                         z_CH2O = -42, z_H2O = -42, z_H2 = -42,  Exp_mass = -42, C = 4, H = 4, O = 0, N = 0, S = 0,
                         S34 = 0, P = 0, N15 = 0, Cl = 0, Fl = 0, Cl37 = 0,Br = 0, Br81 = 0, I = 0,
                         M = 0, NH4 = 0, Z= 0, POE = 0, NOE = 0, D = 0, E = 0, KM_O = -42,
                         KM_CH2 = -42, KM_CH2O = -42, KM_H2O= -42, KM_H2 = -42, NM = 0, RA.y = -42, base_mass = -42,
                         Type = "X", form = "E")

  pb2 <- txtProgressBar(min = 0, max = nLoop, style = 3)

  ###############
  Allmasses <- peaksAll[!is.na(peaksAll$mass),]
  nloop1 <- ceiling((max(Allmasses$mass)-DeNovo)/200)
  #Needs to be saved outside loop so it stays intact
  Unambigsave <- Unambig
  Ambigsave <- Ambig
  Ambigreturn <- Ambigsave[Ambigsave$Exp_mass < DeNovo,]
  #j = 0
  for(j in 0:(nloop1)){
    loop <- j
    masstrim <- DeNovo + 200*loop
    loop1 <- j-1
    if(loop1 < 0) {loop1 <- 0}
    masstrim1 <- DeNovo + 200*loop1
    Ambig2 <- Ambigsave[Ambigsave$Exp_mass < masstrim,] #Resets a comparison df each time

    Ambig4 <- merge(Ambigreturn, Ambig2, by.x = c("Exp_mass", "RA"), by.y = c("Exp_mass", "RA"), all = T) #LCMS

    Common <- Ambig4[!is.na(Ambig4$NM.x)&!is.na(Ambig4$NM.y),]  #LCMS
    Common <- Common[c(1:18)]      ###CHECK TO MAKE SURE THIS IS GOOD 6/2/20
    names(Common) <- gsub(".x","",names(Common),fixed = TRUE)
    New <- Ambig4[is.na(Ambig4$NM.x) & Ambig4$Exp_mass >= masstrim1,]  #LCMS
    New <- New[c(1,2, 19:ncol(New))] #LCMS
    names(New) <- gsub(".y","",names(New),fixed = TRUE)

    Ambig <- rbind(Common, New)
    #Ambig <- unique(Ambig)

    #Seems good to this point
    #i <- 1
    for(i in 1:nLoop){

      known <- Unambig
      known <- rbind(known, knowndummy)

      unknown <- Ambig
      unknown <- rbind(unknown, unknowndummy)
      x <- 1
      repeat{
        x = x+1
        knownCH2 <- known[c("RA", "Exp_mass", "KMD_CH2", "z_CH2", "C", "H", "O", "N", "S", "P", "E",
                            "S34", "N15", "D", "Cl", "Fl", "Cl37", "Br", "Br81", "I",
                            "M", "NH4", "POE", "NOE", "Z")]
        names(knownCH2)[2] <- "base_mass"
        Step1 <- merge(unknown, knownCH2, by.x = c("KMD_CH2", "z_CH2"), by.y = c("KMD_CH2", "z_CH2"))
        Step1$CH2_num <- round(((Step1$Exp_mass - Step1$base_mass))/14.01565)
        Step1$C <- Step1$C + Step1$CH2_num
        Step1$H <- Step1$H + 2 * Step1$CH2_num
        Step1$Type <- "CH2"
        Step1$form <- paste(Step1$C, Step1$H, Step1$O, Step1$N, Step1$S, Step1$P, Step1$E, Step1$S34,
                            Step1$N15, Step1$D, Step1$Cl, Step1$Fl, Step1$Cl37, Step1$Br, Step1$Br81,
                            Step1$I, Step1$M, Step1$NH4,
                            Step1$POE, Step1$NOE, sep = "_")
        Step1 <- Step1[-c(42)]


        knownO <- known[c("RA", "Exp_mass", "KMD_O", "z_O", "C", "H", "O", "N", "S", "P", "E",
                          "S34", "N15", "D", "Cl", "Fl", "Cl37","Br", "Br81", "I", "M", "NH4", "POE", "NOE", "Z")]
        names(knownO)[2] <- "base_mass"
        Step2 <- merge(unknown, knownO, by.x = c("KMD_O", "z_O"), by.y = c("KMD_O", "z_O"))
        Step2$O_num <- round(((Step2$Exp_mass - Step2$base_mass))/15.9949146223)
        Step2$O <- Step2$O + Step2$O_num
        Step2$Type <- "O"
        Step2$form <- paste(Step2$C, Step2$H, Step2$O, Step2$N, Step2$S, Step2$P, Step2$E, Step2$S34,
                            Step2$N15, Step2$D, Step2$Cl, Step2$Fl, Step2$Cl37,Step2$Br, Step2$Br81,
                            Step2$I, Step2$M, Step2$NH4,
                            Step2$POE, Step2$NOE, sep = "_")
        Step2 <- Step2[-c(42)]

        knownH2 <- known[c("RA", "Exp_mass", "KMD_H2", "z_H2", "C", "H", "O", "N", "S", "P", "E",
                           "S34", "N15", "D", "Cl", "Fl", "Cl37","Br", "Br81", "I", "M", "NH4", "POE", "NOE", "Z")]
        names(knownH2)[2] <- "base_mass"
        Step3 <- merge(unknown, knownH2, by.x = c("KMD_H2", "z_H2"), by.y = c("KMD_H2", "z_H2"))
        Step3$H2_num <- round(((Step3$Exp_mass - Step3$base_mass))/2.01565)
        Step3$H <- Step3$H + 2*Step3$H2_num
        Step3$Type <- "H2"
        Step3$form <- paste(Step3$C, Step3$H, Step3$O, Step3$N, Step3$S, Step3$P, Step3$E, Step3$S34,
                            Step3$N15, Step3$D, Step3$Cl, Step3$Fl, Step3$Cl37, Step3$Br, Step3$Br81,
                            Step3$I,Step3$M, Step3$NH4,
                            Step3$POE, Step3$NOE, sep = "_")
        Step3 <- Step3[-c(42)]

        knownH2O <- known[c("RA", "Exp_mass", "KMD_H2O", "z_H2O", "C", "H", "O", "N", "S", "P", "E",
                            "S34", "N15", "D", "Cl", "Fl", "Cl37","Br", "Br81", "I", "M", "NH4", "POE", "NOE", "Z")]
        names(knownH2O)[2] <- "base_mass"
        Step4 <- merge(unknown, knownH2O, by.x = c("KMD_H2O", "z_H2O"), by.y = c("KMD_H2O", "z_H2O"))
        Step4$H2O_num <- round(((Step4$Exp_mass - Step4$base_mass))/18.01056468)
        Step4$H <- Step4$H + 2*Step4$H2O_num
        Step4$O <- Step4$O + Step4$H2O_num
        Step4$Type <- "H2O"
        Step4$form <- paste(Step4$C, Step4$H, Step4$O, Step4$N, Step4$S, Step4$P, Step4$E, Step4$S34,
                            Step4$N15, Step4$D, Step4$Cl, Step4$Fl, Step4$Cl37,Step4$Br, Step4$Br81,
                            Step4$I, Step4$M, Step4$NH4,
                            Step4$POE, Step4$NOE, sep = "_")
        Step4 <- Step4[-c(42)]

        knownCH2O <- known[c("RA", "Exp_mass", "KMD_CH2O", "z_CH2O", "C", "H", "O", "N", "S", "P", "E",
                             "S34", "N15", "D", "Cl", "Fl", "Cl37","Br", "Br81", "I", "M", "NH4", "POE", "NOE", "Z")]
        names(knownCH2O)[2] <- "base_mass"
        Step5 <- merge(unknown, knownCH2O, by.x = c("KMD_CH2O", "z_CH2O"), by.y = c("KMD_CH2O", "z_CH2O"))
        Step5$CH2O_num <- round(((Step5$Exp_mass - Step5$base_mass))/30.01056468)
        Step5$H <- Step5$H + 2*Step5$CH2O_num
        Step5$O <- Step5$O + Step5$CH2O_num
        Step5$C <- Step5$C + Step5$CH2O_num
        Step5$Type <- "CH2O"
        Step5$form <- paste(Step5$C, Step5$H, Step5$O, Step5$N, Step5$S, Step5$P, Step5$E, Step5$S34,
                            Step5$N15, Step5$D, Step5$Cl, Step5$Fl, Step5$Cl37, Step5$Br, Step5$Br81,
                            Step5$I,Step5$M, Step5$NH4,
                            Step5$POE, Step5$NOE, sep = "_")
        Step5 <- Step5[-c(42)]

        Out <- rbind(Step1, Step2, Step3, Step4, Step5)
        Out <- Out[(Out$C >= 2 & Out$H >= 4 & Out$O >= 0),]

        Out$H_C <- Out$H/Out$C   #Quick internal QA to limit bad assignments
        Out$O_C <- Out$O/Out$C
        Out <- Out[Out$H_C <= H_Cmax & Out$H_C >= H_Cmin &
                     Out$O_C <= O_Cmax & Out$O_C >= O_Cmin,]
        Out <- Out[!names(Out) %in% c("H_C", "O_C")]

        Out <- rbind(Out, DummyOut)

        Out_form <- dplyr::group_by(Out, Exp_mass, RA.x, form)  #LCMS

        Out_form$number <- 1
        Out_form <- dplyr::summarize_at(Out_form, "number", sum, na.rm = TRUE)

        #Turns on or off ambiguity based on user input
        if(Ambigcheck == "off") {
          Out_form <- dplyr::filter(Out_form, number == max(number))
        }

        Out_form <- unique(Out_form)
        Out2<- merge(Out, Out_form, by.x = c("Exp_mass", "RA.x", "form"), by.y = c("Exp_mass", "RA.x", "form")) #LCMS
        Out3 <- dplyr::distinct(Out2, form, Exp_mass, RA.x, .keep_all = TRUE) #LCMS
        Out3 <- Out3[!names(Out3) %in% c("number")]


        Next <- Out3[!names(Out3) %in% c("base_mass", "Type", "form", "RA.y")]
        colnames(Next)[colnames(Next) == "RA.x"] <- "RA"
        Next <- rbind(known, Next)
        Next <- unique(Next)
        Unambig <- Next

        masses <- Unambig[c("Exp_mass", "RA", "C")] #LCMS
        names(masses)[3] <- "Var"
        Ambig <- merge(unknown, masses, by.x = c("Exp_mass", "RA"), by.y = c("Exp_mass", "RA"), all = T) #LCMS
        Ambig <- Ambig[is.na(Ambig$Var),]
        Ambig <- Ambig[-19]    ##CHECK THIS 6/2/20
        Ambigreturn <- unique(Ambig)

        if(x == 2){
          break
        }
      }
      setTxtProgressBar(pb2, i)
      Unambig
    }
  }
  records1 <- Unambig[c(1:23)]  #CHECK THIS 6/2/20

  records1$mode <- ionMode



  df1 <- records1[records1$mode == "pos" & records1$M > 0 & records1$POE == 0 & records1$NH4 == 0,]
  df1$Neutral_mass <- df1$Exp_mass - df1$M * 22.989221

  df1x <- records1[records1$mode == "pos" & records1$M == 0 & records1$POE == 0 & records1$NH4 == 1,]
  df1x$Neutral_mass <- df1x$Exp_mass - df1x$NH4 * 18.033823

  df2 <- records1[records1$mode == "pos" & records1$M == 0 & records1$POE == 0& records1$NH4 == 0,]
  df2$Neutral_mass <- df2$Exp_mass - 1.00727645216

  df3 <- records1[records1$mode == "neg" & records1$NOE == 0,]
  df3$Neutral_mass <- df3$Exp_mass + 1.00727645216

  df4 <- records1[records1$mode == "neg" & records1$NOE == 1,]
  df4$Neutral_mass <- df4$Exp_mass - 0.000548597

  df5 <- records1[records1$mode == "pos" & records1$M == 0 & records1$POE == 1& records1$NH4 == 0,]
  df5$Neutral_mass <- df5$Exp_mass + 0.000548597

  records1 <- rbind(df1, df1x, df2, df3, df4, df5)

  records1 <- records1[-c(24)]  ##CHECK THIS 6/2/20

  #recordssave <- records1

  ###Standard QA steps, second round

  records1 <- dplyr::mutate(records1, O_C = O/(C+E), H_C =H/(C+E),

                            #Neutral_mass = Neutral_mass + POE * (2.0156500638/2)- NOE * (2.0156500638/2),


                            theor_mass1 = EM2("C") * C + EM2("H") * H + EM2("O") * O + N * EM2("N14") + S * EM2("S") + P * EM2("P31") +
                              Cl * EM2("Cl35") + Fl * EM2("Fl19") + E * EM2("E") + S34 * EM2("S34") + Cl37 * EM2("Cl37m") + N15 * EM2("N15H") +
                              D * EM2("D") + M * EM2("M") + NH4 * EM2("NH4+") + NOE * electron - POE * electron +
                              Br * EM2("Br79") + Br81 * EM2("Br81m") + I * EM2("I127"),

                            theor_mass = EM2("C") * C + EM2("H") * H + EM2("O") * O + N * EM2("N14") +
                              S * EM2("S") + P * EM2("P31") + Fl * EM2("Fl19") +
                              Cl * EM2("Cl35") +  E * EM2("E2") + S34 * EM2("S34") + Cl37 * EM2("Cl37m") +
                              N15 * EM2("N15H") +
                              D * EM2("D")+
                              Br * EM2("Br79") + Br81 * EM2("Br81m") + I * EM2("I127"),

                            #C = C + E, #It is added back so that formulas are more accurate.

                            DBE = C - 0.5 * (H + Cl + Cl37 +Fl + Br + Br81 + I) + 0.5 * (N +N15+ P) + 1,

                            err_ppm = ((Neutral_mass - theor_mass) / Neutral_mass * 10^6),

                            AE_ppm = round(abs((Neutral_mass - theor_mass) / Neutral_mass * 10^6),2),

                            #NM = floor(Exp_mass),

                            KM = Exp_mass * (14 / 14.01565), #KMD = NM - KM,

                            max_LA = theor_mass1 / 13, actual_LA = ((C - E) + N + S + Fl + O + E + S34 + P +
                                                                      Cl +Cl37+N15 + Br + Br81 + I) ,

                            rule_13=round(actual_LA/max_LA,1),

                            Senior1 = H + P + N + Cl + Fl + Cl37 + N15 + Br + Br81 + I  ,

                            STest = S + S34, ClTest = Cl + Cl37, NTest = N + N15, DBEO = DBE-O, BrTest = Br + Br81,

                            max_H = C * 2 + 2, H_test = round(H / max_H,1),

                            Senior2 = Pval + Nval + N15val + Valence("H")  +
                              Valence("Cl") + Valence("Cl37")+Valence("Fl") + Valence2("Br") + Valence2("Br81") +
                              Valence2("I"),

                            Senior3Atom = C + H + O + N + S + P + N15 + Cl + Cl37 + S34 + Br + Br81 + I,

                            Senior3Val = C*Valence("C") + H*Valence("H") + O*Valence("O") + N*Nval +
                              S*Sval + P*Pval + S34*S34val + Fl*Valence("Fl")+
                              N15*N15val + Cl*Valence("Cl") + Cl37*Valence("Cl37") + Br * Valence2("Br") +
                              Br81 * Valence2("Br81") + I * Valence2("I")
  )

  records1$NM <- ifelse(abs(floor(records1$Exp_mass)-records1$Exp_mass) >= min_def & abs(floor(records1$Exp_mass)-records1$Exp_mass) <= max_def,
                        records1$NM <- floor(records1$Exp_mass), records1$NM <- round(records1$Exp_mass)) #New 1/6/20

  records1$KMD <- ifelse(abs(floor(records1$Exp_mass)-records1$Exp_mass) >= min_def & abs(floor(records1$Exp_mass)-records1$Exp_mass) <= max_def,
                         records1$KMD <- floor(records1$Exp_mass)-records1$KM, records1$KMD <- round(records1$Exp_mass)-records1$KM) #New 1/6/20


  #recordssave <- records1
  #records1 <- recordssave

  records1 <- dplyr::filter(records1, C>0, H>0,O>=Omin, H >= D)
  records1 <- unique(records1)
  records1 <- dplyr::filter(records1, O_C < O_Cmax & H_C <= H_Cmax & H_C > H_Cmin & O_C >= O_Cmin &

                              DBEO >= DBEOmin & DBEO <= DBEOmax &

                              (STest <= Sx + S34x) &
                              (ClTest <= Clx + Cl37x) &
                              (NTest <= Nx + N15x) &
                              (BrTest <= Brx + Br81x) &

                              H_test <= 1 & rule_13 <= 1 &

                              AE_ppm <= ppm_err &

                              Even(Senior1)==TRUE & DBE >= 0 & DBE <= round(0.9 * (C + N)) &

                              O <= 2 * C + 3 * (N+N15) + 4 * P + 4 * (S+S34)&

                              O >= Omin&

                              RA >= 0 &

                              Senior2 >= 2*Valence("C") &

                              Senior3Val >= (2*Senior3Atom - 1)
  )



  records1 <- records1[!names(records1) %in% c("Senior1", "Senior2", "Senior3Val", "Senior3Atom")]

  #records1 <- records1 %>% filter(AE_ppm <= 2)
  #check <- records1 %>% distinct(Exp_mass, .keep_all = TRUE) %>% filter(AE_ppm <= 2)
  ###Formula generation
  records1 <-
    dplyr::mutate(records1, Cform = ifelse(C == 0 , "",
                                            ifelse(C == 1 , "C", paste("C",C, sep = ""))),
                  Hform = ifelse(H == 0 , "",
                                 ifelse(H == 1 , "H", paste("H",H, sep = ""))),
                  Nform = ifelse(NTest == 0 , "",
                                 ifelse(NTest == 1 , "N", paste("N",NTest, sep = ""))),
                  Oform = ifelse(O == 0 , "",
                                 ifelse(O == 1 , "O", paste("O",O, sep = ""))),
                  Sform = ifelse(STest == 0 , "",
                                 ifelse(STest == 1 , "S", paste("S",STest, sep = ""))),
                  Pform = ifelse(P == 0 , "",
                                 ifelse(P == 1 , "P", paste("P",P, sep = ""))),
                  Clform = ifelse(ClTest == 0 , "",
                                  ifelse(ClTest == 1 , "Cl", paste("Cl",ClTest, sep = ""))),
                  Flform = ifelse(Fl == 0 , "",
                                  ifelse(Fl == 1 , "F", paste("F",Fl, sep = ""))),
                  Brform = ifelse(BrTest == 0 , "",
                                  ifelse(BrTest == 1 , "Br", paste("Br",BrTest, sep = ""))),
                  Iform = ifelse(I == 0 , "",
                                 ifelse(I == 1 , "I", paste("I",I, sep = ""))))

  records1 <- tidyr::unite(records1, class, Nform, Oform, Sform, Pform, Clform, Flform, Brform, Iform,
                            sep = "", remove = FALSE)

  records1 <- tidyr::unite(records1, formula, Cform, Hform, Nform, Oform, Sform, Pform, Clform,
                            Flform, Brform, Iform, sep = "")

  records1 <-
    dplyr::mutate(records1, Cform = ifelse(C == 0 , "", "C"),
                  Hform = ifelse(H == 0 , "", "H"),

                  Nform = ifelse(NTest == 0 , "", "N"),

                  Oform = ifelse(O == 0 , "", "O"),

                  Sform = ifelse(STest == 0 , "","S"),
                  Pform = ifelse(P == 0 , "", "P"),
                  Clform = ifelse(ClTest == 0 , "", "Cl"),
                  Flform = ifelse(Fl == 0 , "", "F"),
                  Brform = ifelse(BrTest == 0 , "", "Br"),
                  Iform = ifelse(I == 0 , "", "I"))

  records1 <- tidyr::unite(records1, group, Cform, Hform, Nform, Oform, Sform, Pform,
                            Clform, Flform, Brform, Iform, sep = "")

  ###Supplemental Specialized QA Steps
  records1<-dplyr::mutate(records1, HA = NTest + STest + P + ClTest + E + Fl + I + BrTest)

  records1<-dplyr::group_by(records1, Exp_mass, RA) #LCMS

  ifelse(HetCut == "on", records1<-dplyr::filter(records1, HA == (min(HA))), records1<- records1)

  records1<-dplyr::group_by(records1, Exp_mass, RA) #LCMS

  records1 <- dplyr::distinct(records1, formula, RA, .keep_all = TRUE)  ###Change Here 05/14/19

  records1 <- dplyr::ungroup(records1)

  records1 <- records1[!names(records1) %in% c("HA")]



  #records3 <- dplyr::rename(records1, mass = Exp_mass)

  cut <- (peaksAll2)
  cut <- unique(cut)
  cut$Test <- paste(cut$Exp_mass, cut$RA, sep = "_")   #New 05/14/19
  records1$Test <- paste(records1$Exp_mass, records1$RA, sep = "_") #New 05/14/19

  unassigned <- dplyr::left_join(cut, records1, by = "Test")  #New 05/14/19
  unassigned <- unassigned[is.na(unassigned$formula),]
  unassigned <- unassigned[c("RA.x", "Exp_mass.x")]
  names(unassigned)[1] <- "RA"
  names(unassigned)[2] <- "mass"
  unassigned <- unique(unassigned)  #Good to this point
  records1 <- records1[-c(48)] #Removes Test  ##CHECK THIS 6/2/20

  records1$mode <- ionMode


  df1 <- records1[records1$mode == "pos",]
  df1$theor_mass1 <- df1$theor_mass1 + proton

  df2 <- records1[records1$mode == "neg",]
  df2$theor_mass1 <- df2$theor_mass1 - proton

  records1 <- rbind(df1, df2)
  records1 <- records1[-c(48)]  #removes mode  ##CHECK THIS 6/2/20



  records1 <- records1[c(1,2,45:47,3:24,27,30:34, 29, 25:26,39, 35, 41:42,43, 44)]  ##CHECK THIS 6/2/20

  #recordssave <- records1
  #records1 <- recordssave
  #Columns still good at this point
  ######################################################################


  ##Aligning Isotope masses back into the mass spectrum
  ##Align single C13 masses
  records1$C13_mass <- records1$Exp_mass + 1.0033548380
  err <- iso_err*10^-6
  #The following line was changed for isotoping
  C13Iso <- isopeaks2[isopeaks2$Tag == "C13"|isopeaks2$Tag == "C13_S34",]
  names(C13Iso)[2] <- "C13_mass"
  names(C13Iso)[1] <- "C13_Abund"
  records1$C13_mass <- sapply(records1$C13_mass, function(x){
    # First check if the element lies within tolerance limits of any element in df2
    ifelse(min(abs(C13Iso$C13_mass - x), na.rm=TRUE) < err * x,
           # If yes, replace that element in df1 with the matching element in df2
           C13Iso[which.min(abs(C13Iso$C13_mass - x)), "C13_mass"], 0)
  })
  ## New as of 12/13/19
  records1$C13_mass_2 <- records1$theor_mass1 + 1.0033548380
  err <- iso_err*10^-6
  #The following line was changed for isotoping
  C13Iso <- isopeaks2[isopeaks2$Tag == "C13"|isopeaks2$Tag == "C13_S34",]
  names(C13Iso)[2] <- "C13_mass_2"
  names(C13Iso)[1] <- "C13_Abund"
  records1$C13_mass_2 <- sapply(records1$C13_mass_2, function(x){
    # First check if the element lies within tolerance limits of any element in df2
    ifelse(min(abs(C13Iso$C13_mass_2 - x), na.rm=TRUE) < err * x,
           # If yes, replace that element in df1 with the matching element in df2
           C13Iso[which.min(abs(C13Iso$C13_mass_2 - x)), "C13_mass_2"], 0)
  })

  No_iso <- records1[records1$C13_mass == 0 & records1$C13_mass_2 == 0,]
  No_iso <- No_iso[-44]  ##CHECK
  Both_iso <- records1[records1$C13_mass != 0 & records1$C13_mass_2 != 0,]
  Both_iso <- Both_iso[-44] ##CHECK
  Match_iso <- records1[records1$C13_mass != 0 & records1$C13_mass_2 == 0,]
  Match_iso <- Match_iso[-44] ##CHECK
  Assign_iso <- records1[records1$C13_mass == 0 & records1$C13_mass_2 != 0,]
  Assign_iso <- Assign_iso[-43] ##CHECK
  names(Assign_iso)[43] <- "C13_mass"  ##CHECK

  records1 <- rbind(No_iso, Both_iso, Match_iso, Assign_iso)
  C13Iso <- isopeaks2[isopeaks2$Tag == "C13"|isopeaks2$Tag == "C13_S34",]
  names(C13Iso)[2] <- "C13_mass"
  names(C13Iso)[1] <- "C13_Abund"

  ##
  #Ccheck <- records1[records1$C13_mass > 0,]
  #check1 <- records1 %>% mutate(C13 = Exp_mass + 1.0033548380, err = abs((C13_mass- C13)/C13_mass*10^6))
  #check2 <- check1 %>% filter(err < 3)
  #Ccheck <- records1[records1$C13_mass > 0,]



  records1 <- dplyr::left_join(records1, C13Iso, by = "C13_mass")
  records1 <- records1[-45]  ##CHECK
  #########
  #Align double C13 masses
  records1$C13_mass2 <- records1$Exp_mass + 2.006709676
  err <- iso_err*10^-6
  #The following line was changed for isotoping
  C13Iso2 <- isopeaks2[isopeaks2$Tag == "2C13"| isopeaks2$Tag == "2C13_S34",]
  names(C13Iso2)[2] <- "C13_mass2"
  names(C13Iso2)[1] <- "C13_Abund2"
  records1$C13_mass2 <- sapply(records1$C13_mass2, function(x){
    # First check if the element lies within tolerance limits of any element in df2
    ifelse(min(abs(C13Iso2$C13_mass2 - x), na.rm=TRUE) < err * x,
           # If yes, replace that element in df1 with the matching element in df2
           C13Iso2[which.min(abs(C13Iso2$C13_mass2 - x)), "C13_mass2"], 0)
  })

  ## New as of 12/13/19
  records1$C13_mass2_2 <- records1$theor_mass1 + 2* 1.0033548380
  err <- iso_err*10^-6

  #The following line was changed for isotoping
  C13Iso2 <- isopeaks2[isopeaks2$Tag == "2C13"| isopeaks2$Tag == "2C13_S34",]
  names(C13Iso2)[2] <- "C13_mass2_2"
  names(C13Iso2)[1] <- "C13_Abund2_2"
  records1$C13_mass2_2 <- sapply(records1$C13_mass2_2, function(x){
    # First check if the element lies within tolerance limits of any element in df2
    ifelse(min(abs(C13Iso2$C13_mass2_2 - x), na.rm=TRUE) < err * x,
           # If yes, replace that element in df1 with the matching element in df2
           C13Iso2[which.min(abs(C13Iso2$C13_mass2_2 - x)), "C13_mass2_2"], 0)
  })


  No_iso <- records1[records1$C13_mass2 == 0 & records1$C13_mass2_2 == 0,]
  No_iso <- No_iso[-46]  ##CHECK
  Both_iso <- records1[records1$C13_mass2 != 0 & records1$C13_mass2_2 != 0,]
  Both_iso <- Both_iso[-46]  ##CHECK
  Match_iso <- records1[records1$C13_mass2 != 0 & records1$C13_mass2_2 == 0,]
  Match_iso <- Match_iso[-46] ##CHECK
  Assign_iso <- records1[records1$C13_mass2 == 0 & records1$C13_mass2_2 != 0,]
  Assign_iso <- Assign_iso[-45]  ##CHECK
  names(Assign_iso)[45] <- "C13_mass2"   ##CHECK

  records1 <- rbind(No_iso, Both_iso, Match_iso, Assign_iso)
  C13Iso2 <- isopeaks2[isopeaks2$Tag == "2C13"|isopeaks2$Tag == "2C13_S34",]
  names(C13Iso2)[2] <- "C13_mass2"
  names(C13Iso2)[1] <- "C13_Abund2"

  ##


  #Ccheck <- records1[records1$C13_mass > 0,]
  records1 <- dplyr::left_join(records1, C13Iso2, by = "C13_mass2")
  records1$C13_mass2 <- ifelse(records1$C13_mass == 0, 0, records1$C13_mass2)
  records1$C13_Abund2 <- ifelse(records1$C13_mass2 == 0, 0, records1$C13_Abund2)
  # This is removing the peaks because it cannot find the single C13
  #records1 <- records1[!(records1$C13_mass == 0 & records1$C13_mass2 > 0),]
  records1 <- records1[-47]  ##CHECK
  #recordssave <- records1

  ###########
  #Align S34 masses  #This is a potential error zone if Sx >0 and no sulfur assigned, or if Sx >0 and
  #no non-sulfur are assigned. Standard conditions should not have a problem.
  if(Sx > 0){
    recordsdummy <- records1[1,]
    recordsdummy[is.na(recordsdummy)] <- 0
    recordsdummy$group[!is.na(recordsdummy$group)] <- "Dummy"
    recordsSulf <- records1[records1$S >0,]
    recordsSulf <- rbind(recordsSulf, recordsdummy)
    recordsrest <- records1[records1$S ==0,]
    recordsrest <- rbind(recordsrest, recordsdummy)
    recordsrest$S34_mass <- 0
    recordsrest$S34_Abund <- 0
    recordsSulf$S34_mass <- recordsSulf$Exp_mass + 1.995797
    err <- iso_err*10^-6
    recordsSulf <- recordsSulf[!is.na(recordsSulf$S),]

    #The following line was changed for the isotoping
    S34Iso <- isopeaks2[isopeaks2$Tag == "S34"|isopeaks2$Tag == "C13_S34"|isopeaks2$Tag == "2C13_S34",]
    names(S34Iso)[2] <- "S34_mass"
    names(S34Iso)[1] <- "S34_Abund"
    recordsSulf$S34_mass <- sapply(recordsSulf$S34_mass, function(x){
      # First check if the element lies within tolerance limits of any element in df2
      ifelse(min(abs(S34Iso$S34_mass - x), na.rm=TRUE) < err * x,
             # If yes, replace that element in df1 with the matching element in df2
             S34Iso[which.min(abs(S34Iso$S34_mass - x)), "S34_mass"], 0)
    })

    recordsSulf <- recordsSulf[!is.na(recordsSulf$C),]
    ## New as of 12/13/19
    recordsSulf$S34_mass_2 <- recordsSulf$theor_mass1 + 1.995797
    err <- iso_err*10^-6

    #The following line was changed for isotoping
    S34Iso <- isopeaks2[isopeaks2$Tag == "S34"|isopeaks2$Tag == "C13_S34"|isopeaks2$Tag == "2C13_S34",]
    names(S34Iso)[2] <- "S34_mass_2"
    names(S34Iso)[1] <- "S34_Abund_2"
    recordsSulf$S34_mass_2 <- sapply(recordsSulf$S34_mass_2, function(x){
      # First check if the element lies within tolerance limits of any element in df2
      ifelse(min(abs(S34Iso$S34_mass_2 - x), na.rm=TRUE) < err * x,
             # If yes, replace that element in df1 with the matching element in df2
             S34Iso[which.min(abs(S34Iso$S34_mass_2 - x)), "S34_mass_2"], 0)
    })


    No_iso <- recordsSulf[recordsSulf$S34_mass == 0 & recordsSulf$S34_mass_2 == 0,]
    No_iso <- No_iso[-48]  ##CHECK
    Both_iso <- recordsSulf[recordsSulf$S34_mass != 0 & recordsSulf$S34_mass_2 != 0,]
    Both_iso <- Both_iso[-48]  ##CHECK
    Match_iso <- recordsSulf[recordsSulf$S34_mass != 0 & recordsSulf$S34_mass_2 == 0,]
    Match_iso <- Match_iso[-48]  ##CHECK
    Assign_iso <- recordsSulf[recordsSulf$S34_mass == 0 & recordsSulf$S34_mass_2 != 0,]
    Assign_iso <- Assign_iso[-47]  ##CHECK
    names(Assign_iso)[47] <- "S34_mass"  ##CHECK

    recordsSulf <- rbind(No_iso, Both_iso, Match_iso, Assign_iso)
    S34Iso <- isopeaks2[isopeaks2$Tag == "S34"|isopeaks2$Tag == "C13_S34"|isopeaks2$Tag == "2C13_S34",]
    names(S34Iso)[2] <- "S34_mass"
    names(S34Iso)[1] <- "S34_Abund"
    ##





    recordsSulf <- dplyr::left_join(recordsSulf, S34Iso, by = "S34_mass")
    #recordsSulf <- merge(recordsSulf, S34Iso, by.x = "S34_mass", by.y = "S34_mass")
    recordsSulf <- recordsSulf[-49]  ##CHECK
    records1 <- rbind(recordsrest, recordsSulf)
    records1 <- records1[records1$group !="Dummy",]
  } else{records1$S34_mass <- 0;
  records1$S34_Abund <- 0}



  #Checking to see if any polyisotope masses match an assigned monoisotope mass
  Iso_check <- isopeaks2
  names(Iso_check)[2] <- "Exp_mass"
  records1$Iso_RA <- records1$RA
  Iso_check2 <- merge(records1, Iso_check, by.x = c("Exp_mass", "Iso_RA"), by.y = c("Exp_mass", "Iso_RA")) #LCMS
  Mono_check2 <- merge(records1, Iso_check, by.x = c("Exp_mass", "Iso_RA"), by.y = c("Exp_mass", "Iso_RA"), all = T) #LCMS
  Mono_check2 <- Mono_check2[!is.na(Mono_check2$C)& is.na(Mono_check2$Tag),]

  I1 <- Mono_check2[c("C13_mass", "C13_Abund")]  #LCMS
  names(I1)[1] <- "Exp_mass"    #LCMS
  names(I1)[2] <- "Iso_RA"      #LCMS
  I2 <- Mono_check2[c("C13_mass2", "C13_Abund2")] #LCMS
  names(I2)[1] <- "Exp_mass"    #LCMS
  names(I2)[2] <- "Iso_RA"    #LCMS
  I3 <- Mono_check2[c("S34_mass", "S34_Abund")]   #LCMS
  names(I3)[1] <- "Exp_mass"    #LCMS
  names(I3)[2] <- "Iso_RA"     #LCMS

  IM <- rbind(I1, I2, I3)
  IM$Tag2 <- "Iso"
  IM <- unique(IM)

  Iso_check3 <- merge(Iso_check2, IM, by.x = c("Exp_mass", "Iso_RA"), by.y = c("Exp_mass", "Iso_RA"), all = T) #LCMS
  Iso_check3 <- Iso_check3[!is.na(Iso_check3$C),]
  MonoG1 <- Iso_check3[Iso_check3$C13_mass > 0 | Iso_check3$C13_mass2 > 0 | Iso_check3$S34_mass > 0,]
  MonoG2 <- Iso_check3[is.na(Iso_check3$Tag2) & (Iso_check3$C13_mass == 0 &
                                                   Iso_check3$C13_mass2 == 0 & Iso_check3$S34_mass == 0),]
  MonoGF <- rbind(MonoG1, MonoG2)
  MonoGF <- MonoGF[c(1,3:50)]   ##LCMS
  MonoRest <- Mono_check2[c(1,3:50)]  ##LCMS
  records1 <- rbind(MonoRest, MonoGF)
  records1 <- records1[c(2,1,3:49)]  ##CHECK

  Iso_save <- records1[c("Exp_mass", "RA", "C13_mass", "C13_Abund", "C13_mass2", "C13_Abund2", "S34_mass", "S34_Abund")] #LCMS
  ###########################################################################################################
  #Additional formula extension for the isotope masses not matched to monoisotopic mass

  #check <- records1[records1$S34_mass > 0,]
  #Preparation of unmatched "isotope" mass list
  #Identifying the peaks that match first
  C13 <- records1[c("C13_mass", "C13_Abund")]   #LCMS
  C13_2 <- records1[c("C13_mass2", "C13_Abund2")]  #LCMS
  S34 <- records1[records1$S34_mass > 0 & records1$C13_mass > 0,] #Chosen as most likely to be true 34S
  S34 <- S34[c("S34_mass", "S34_Abund")]  #LCMS
  C13 <- unique(C13)
  C13_2 <- unique(C13_2)
  S34 <- unique(S34)
  names(C13)[1] <- "Iso_mass"
  names(C13)[2] <- "Iso_RA"   #LCMS
  names(C13_2)[1] <- "Iso_mass"
  names(C13_2)[2] <- "Iso_RA"   #LCMS
  names(S34)[1] <- "Iso_mass"
  names(S34)[2] <- "Iso_RA"   #LCMS
  Iso_match <- rbind(C13, C13_2, S34)
  Iso_match <- unique(Iso_match)
  Iso_match$Tag2 <- "Match"
  Iso_match <- Iso_match[Iso_match$Iso_mass > 0,] #Removes rows with 0 in them

  #Finding the "isotope" masses that did not get assigned
  Iso_align <- merge(isopeaks2, Iso_match, by.x = c("Iso_mass", "Iso_RA"), by.y = c("Iso_mass", "Iso_RA"), all = T) #LCMS
  Iso_nomatch <- Iso_align[is.na(Iso_align$Tag2),]
  Iso_nomatch <- Iso_nomatch[c(1,2)]
  names(Iso_nomatch)[1] <- "Exp_mass"
  names(Iso_nomatch)[2] <- "RA"


  ###################
  #New 12/10/19
  Unassigned <- records1[c(1,2,3)] #LCMS
  Unassigned <- merge(peaksAll2, Unassigned, by.x = c("Exp_mass", "RA"), by.y = c("Exp_mass", "RA"), all = T) #LCMS
  Unassigned <- Unassigned[is.na(Unassigned$formula),]  #LCMS
  Unassigned <- Unassigned[c(1,2)]
  #names(Unassigned)[2] <- "RA" #LCMS

  Iso_nomatch <- rbind(Iso_nomatch, Unassigned)
  Iso_nomatch <- unique(Iso_nomatch)
  ###################
  #Kendrick Series Preparation
  Iso_nomatch <- Iso_nomatch[c(1,2)]
  ####
  Iso_nomatch$NM <- ifelse(abs(floor(Iso_nomatch$Exp_mass)-Iso_nomatch$Exp_mass) >= min_def & abs(floor(Iso_nomatch$Exp_mass)-Iso_nomatch$Exp_mass) <= max_def,
                       Iso_nomatch$NM <- floor(Iso_nomatch$Exp_mass), Iso_nomatch$NM <- round(Iso_nomatch$Exp_mass)) #New 1/6/20
  ####

  Iso_nomatch$KM_CH2 <- Iso_nomatch$Exp_mass * (14/14.01565)
  Iso_nomatch$KMD_CH2 <- round(Iso_nomatch$NM - Iso_nomatch$KM_CH2, 3)
  ###
  Iso_nomatch$z_CH2 <- ifelse(abs(floor(Iso_nomatch$Exp_mass)-Iso_nomatch$Exp_mass) >= min_def & abs(floor(Iso_nomatch$Exp_mass)-Iso_nomatch$Exp_mass) <= max_def,
                          Iso_nomatch$z_CH2 <- floor(Iso_nomatch$Exp_mass)%%14 - 14, Iso_nomatch$z_CH2 <- round(Iso_nomatch$Exp_mass)%%14 - 14) #New 1/6/20
  ###

  Iso_nomatch$KM_O <- Iso_nomatch$Exp_mass * (16/15.9949146223)
  Iso_nomatch$KMD_O <- round(Iso_nomatch$NM - Iso_nomatch$KM_O, 3)
  ###
  Iso_nomatch$z_O <- ifelse(abs(floor(Iso_nomatch$Exp_mass)-Iso_nomatch$Exp_mass) >= min_def & abs(floor(Iso_nomatch$Exp_mass)-Iso_nomatch$Exp_mass) <= max_def,
                              Iso_nomatch$z_O <- floor(Iso_nomatch$Exp_mass)%%16 - 16, Iso_nomatch$z_O <- round(Iso_nomatch$Exp_mass)%%16 - 16) #New 1/6/20
  ###

  Iso_nomatch$KM_H2 <- Iso_nomatch$Exp_mass * (2/2.01565)
  Iso_nomatch$KMD_H2 <- round(Iso_nomatch$NM - Iso_nomatch$KM_H2, 3)
  ###
  Iso_nomatch$z_H2 <- ifelse(abs(floor(Iso_nomatch$Exp_mass)-Iso_nomatch$Exp_mass) >= min_def & abs(floor(Iso_nomatch$Exp_mass)-Iso_nomatch$Exp_mass) <= max_def,
                              Iso_nomatch$z_H2 <- floor(Iso_nomatch$Exp_mass)%%2 - 2, Iso_nomatch$z_H2 <- round(Iso_nomatch$Exp_mass)%%2 - 2) #New 1/6/20
  ###

  Iso_nomatch$KM_H2O <- Iso_nomatch$Exp_mass * (18/18.01056468)
  Iso_nomatch$KMD_H2O <- round(Iso_nomatch$NM - Iso_nomatch$KM_H2O, 3)
  ###
  Iso_nomatch$z_H2O <- ifelse(abs(floor(Iso_nomatch$Exp_mass)-Iso_nomatch$Exp_mass) >= min_def & abs(floor(Iso_nomatch$Exp_mass)-Iso_nomatch$Exp_mass) <= max_def,
                             Iso_nomatch$z_H2O <- floor(Iso_nomatch$Exp_mass)%%18 - 18, Iso_nomatch$z_H2O <- round(Iso_nomatch$Exp_mass)%%18 - 18) #New 1/6/20
  ###

  Iso_nomatch$KM_CH2O <- Iso_nomatch$Exp_mass * (30/30.01056468)
  Iso_nomatch$KMD_CH2O <- round(Iso_nomatch$NM - Iso_nomatch$KM_CH2O, 3)
  ###
  Iso_nomatch$z_CH2O <- ifelse(abs(floor(Iso_nomatch$Exp_mass)-Iso_nomatch$Exp_mass) >= min_def & abs(floor(Iso_nomatch$Exp_mass)-Iso_nomatch$Exp_mass) <= max_def,
                              Iso_nomatch$z_CH2O <- floor(Iso_nomatch$Exp_mass)%%30 - 30, Iso_nomatch$z_CH2O <- round(Iso_nomatch$Exp_mass)%%30 - 30) #New 1/6/20
  ###


  recordsx <- records1[c("RA", "Exp_mass", "C", "H", "O", "N", "S", "P", "E", "S34", "N15", "D",
                         "Cl", "Fl", "Cl37", "Br", "Br81", "I", "M", "NH4", "POE", "NOE", "Z")]
  ####
  recordsx$NM <- ifelse(abs(floor(recordsx$Exp_mass)-recordsx$Exp_mass) >= min_def & abs(floor(recordsx$Exp_mass)-recordsx$Exp_mass) <= max_def,
                           recordsx$NM <- floor(recordsx$Exp_mass), recordsx$NM <- round(recordsx$Exp_mass)) #New 1/6/20
  ####

  recordsx$KM_CH2 <- recordsx$Exp_mass * (14/14.01565)
  recordsx$KMD_CH2 <- round(recordsx$NM - recordsx$KM_CH2, 3)
  ###
  recordsx$z_CH2 <- ifelse(abs(floor(recordsx$Exp_mass)-recordsx$Exp_mass) >= min_def & abs(floor(recordsx$Exp_mass)-recordsx$Exp_mass) <= max_def,
                              recordsx$z_CH2 <- floor(recordsx$Exp_mass)%%14 - 14, recordsx$z_CH2 <- round(recordsx$Exp_mass)%%14 - 14) #New 1/6/20
  ###

  recordsx$KM_O <- recordsx$Exp_mass * (16/15.9949146223)
  recordsx$KMD_O <- round(recordsx$NM - recordsx$KM_O, 3)
  ###
  recordsx$z_O <- ifelse(abs(floor(recordsx$Exp_mass)-recordsx$Exp_mass) >= min_def & abs(floor(recordsx$Exp_mass)-recordsx$Exp_mass) <= max_def,
                            recordsx$z_O <- floor(recordsx$Exp_mass)%%16 - 16, recordsx$z_O <- round(recordsx$Exp_mass)%%16 - 16) #New 1/6/20
  ###

  recordsx$KM_H2 <- recordsx$Exp_mass * (2/2.01565)
  recordsx$KMD_H2 <- round(recordsx$NM - recordsx$KM_H2, 3)
  ###
  recordsx$z_H2 <- ifelse(abs(floor(recordsx$Exp_mass)-recordsx$Exp_mass) >= min_def & abs(floor(recordsx$Exp_mass)-recordsx$Exp_mass) <= max_def,
                             recordsx$z_H2 <- floor(recordsx$Exp_mass)%%2 - 2, recordsx$z_H2 <- round(recordsx$Exp_mass)%%2 - 2) #New 1/6/20
  ###

  recordsx$KM_H2O <- recordsx$Exp_mass * (18/18.01056468)
  recordsx$KMD_H2O <- round(recordsx$NM - recordsx$KM_H2O, 3)
  ###
  recordsx$z_H2O <- ifelse(abs(floor(recordsx$Exp_mass)-recordsx$Exp_mass) >= min_def & abs(floor(recordsx$Exp_mass)-recordsx$Exp_mass) <= max_def,
                              recordsx$z_H2O <- floor(recordsx$Exp_mass)%%18 - 18, recordsx$z_H2O <- round(recordsx$Exp_mass)%%18 - 18) #New 1/6/20
  ###

  recordsx$KM_CH2O <- recordsx$Exp_mass * (30/30.01056468)
  recordsx$KMD_CH2O <- round(recordsx$NM - recordsx$KM_CH2O, 3)
  ###
  recordsx$z_CH2O <- ifelse(abs(floor(recordsx$Exp_mass)-recordsx$Exp_mass) >= min_def & abs(floor(recordsx$Exp_mass)-recordsx$Exp_mass) <= max_def,
                               recordsx$z_CH2O <- floor(recordsx$Exp_mass)%%30 - 30, recordsx$z_CH2O <- round(recordsx$Exp_mass)%%30 - 30) #New 1/6/20
  ###
  Unambig <- recordsx
  ##Formula Extension for isotope masses

  pb2 <- txtProgressBar(min = 0, max = nLoop, style = 3)

  ###############
  Allmasses <- peaksAll[!is.na(peaksAll$mass),]
  nloop1 <- ceiling((max(Allmasses$mass)-DeNovo)/200)
  #Needs to be saved outside loop so it stays intact
  Unambigsave <- recordsx
  Ambigsave <- Iso_nomatch
  Ambigreturn <- Ambigsave[Ambigsave$Exp_mass < DeNovo,]
  #j = 0
  for(j in 0:(nloop1)){
    loop <- j
    masstrim <- DeNovo + 200*loop
    loop1 <- j-1
    if(loop1 < 0) {loop1 <- 0}
    masstrim1 <- DeNovo + 200*loop1
    Ambig2 <- Ambigsave[Ambigsave$Exp_mass < masstrim,] #Resets a comparison df each time

    Ambig4 <- merge(Ambigreturn, Ambig2, by.x = c("Exp_mass", "RA"), by.y = c("Exp_mass", "RA"), all = T) #LCMS

    Common <- Ambig4[!is.na(Ambig4$NM.x)&!is.na(Ambig4$NM.y),]  #LCMS
    Common <- Common[c(1:18)]   ##CHECK
    names(Common) <- gsub(".x","",names(Common),fixed = TRUE)
    New <- Ambig4[is.na(Ambig4$NM.x) & Ambig4$Exp_mass >= masstrim1,]   #LCMS
    New <- New[c(1,2,19:ncol(New))]  #LCMS
    names(New) <- gsub(".y","",names(New),fixed = TRUE)

    Ambig <- rbind(Common, New)
    #Ambig <- unique(Ambig)
    #i <- 1
    #Seems good to this point
    for(i in 1:nLoop){

      known <- Unambig
      known <- rbind(known, knowndummy)

      unknown <- Ambig
      unknown <- rbind(unknown, unknowndummy)
      x <- 1
      repeat{
        x = x+1
        knownCH2 <- known[c("RA", "Exp_mass", "KMD_CH2", "z_CH2", "C", "H", "O", "N", "S", "P", "E",
                            "S34", "N15", "D", "Cl", "Fl", "Cl37", "Br", "Br81", "I", "M", "NH4", "POE", "NOE", "Z")]
        names(knownCH2)[2] <- "base_mass"
        Step1 <- merge(unknown, knownCH2, by.x = c("KMD_CH2", "z_CH2"), by.y = c("KMD_CH2", "z_CH2"))
        Step1$CH2_num <- round(((Step1$Exp_mass - Step1$base_mass))/14.01565)
        Step1$C <- Step1$C + Step1$CH2_num
        Step1$H <- Step1$H + 2 * Step1$CH2_num
        Step1$Type <- "CH2"
        Step1$form <- paste(Step1$C, Step1$H, Step1$O, Step1$N, Step1$S, Step1$P, Step1$E, Step1$S34,
                            Step1$N15, Step1$D, Step1$Cl, Step1$Fl, Step1$Cl37, Step1$Br, Step1$Br81,
                            Step1$I, Step1$M, Step1$NH4,
                            Step1$POE, Step1$NOE, sep = "_")
        Step1 <- Step1[-c(42)]  ##CHECK THIS


        knownO <- known[c("RA", "Exp_mass", "KMD_O", "z_O", "C", "H", "O", "N", "S", "P", "E",
                          "S34", "N15", "D", "Cl", "Fl", "Cl37", "Br", "Br81", "I", "M", "NH4", "POE", "NOE", "Z")]
        names(knownO)[2] <- "base_mass"
        Step2 <- merge(unknown, knownO, by.x = c("KMD_O", "z_O"), by.y = c("KMD_O", "z_O"))
        Step2$O_num <- round(((Step2$Exp_mass - Step2$base_mass))/15.9949146223)
        Step2$O <- Step2$O + Step2$O_num
        Step2$Type <- "O"
        Step2$form <- paste(Step2$C, Step2$H, Step2$O, Step2$N, Step2$S, Step2$P, Step2$E, Step2$S34,
                            Step2$N15, Step2$D, Step2$Cl, Step2$Fl, Step2$Cl37,Step2$Br, Step2$Br81,
                            Step2$I, Step2$M, Step2$NH4,
                            Step2$POE, Step2$NOE, sep = "_")
        Step2 <- Step2[-c(42)] ##CHECK THIS

        knownH2 <- known[c("RA", "Exp_mass", "KMD_H2", "z_H2", "C", "H", "O", "N", "S", "P", "E",
                           "S34", "N15", "D", "Cl", "Fl", "Cl37", "Br", "Br81", "I", "M", "NH4", "POE", "NOE", "Z")]
        names(knownH2)[2] <- "base_mass"
        Step3 <- merge(unknown, knownH2, by.x = c("KMD_H2", "z_H2"), by.y = c("KMD_H2", "z_H2"))
        Step3$H2_num <- round(((Step3$Exp_mass - Step3$base_mass))/2.01565)
        Step3$H <- Step3$H + 2*Step3$H2_num
        Step3$Type <- "H2"
        Step3$form <- paste(Step3$C, Step3$H, Step3$O, Step3$N, Step3$S, Step3$P, Step3$E, Step3$S34,
                            Step3$N15, Step3$D, Step3$Cl, Step3$Fl, Step3$Cl37,Step3$Br, Step3$Br81,
                            Step3$I, Step3$M, Step3$NH4,
                            Step3$POE, Step3$NOE, sep = "_")
        Step3 <- Step3[-c(42)]  ##CHECK THIS

        knownH2O <- known[c("RA", "Exp_mass", "KMD_H2O", "z_H2O", "C", "H", "O", "N", "S", "P", "E",
                            "S34", "N15", "D", "Cl", "Fl", "Cl37", "Br", "Br81", "I", "M", "NH4", "POE", "NOE", "Z")]
        names(knownH2O)[2] <- "base_mass"
        Step4 <- merge(unknown, knownH2O, by.x = c("KMD_H2O", "z_H2O"), by.y = c("KMD_H2O", "z_H2O"))
        Step4$H2O_num <- round(((Step4$Exp_mass - Step4$base_mass))/18.01056468)
        Step4$H <- Step4$H + 2*Step4$H2O_num
        Step4$O <- Step4$O + Step4$H2O_num
        Step4$Type <- "H2O"
        Step4$form <- paste(Step4$C, Step4$H, Step4$O, Step4$N, Step4$S, Step4$P, Step4$E, Step4$S34,
                            Step4$N15, Step4$D, Step4$Cl, Step4$Fl, Step4$Cl37, Step4$Br, Step4$Br81,
                            Step4$I,Step4$M, Step4$NH4,
                            Step4$POE, Step4$NOE, sep = "_")
        Step4 <- Step4[-c(42)]  ##CHECK THIS

        knownCH2O <- known[c("RA", "Exp_mass", "KMD_CH2O", "z_CH2O", "C", "H", "O", "N", "S", "P", "E",
                             "S34", "N15", "D", "Cl", "Fl", "Cl37","Br", "Br81", "I", "M", "NH4", "POE", "NOE", "Z")]
        names(knownCH2O)[2] <- "base_mass"
        Step5 <- merge(unknown, knownCH2O, by.x = c("KMD_CH2O", "z_CH2O"), by.y = c("KMD_CH2O", "z_CH2O"))
        Step5$CH2O_num <- round(((Step5$Exp_mass - Step5$base_mass))/30.01056468)
        Step5$H <- Step5$H + 2*Step5$CH2O_num
        Step5$O <- Step5$O + Step5$CH2O_num
        Step5$C <- Step5$C + Step5$CH2O_num
        Step5$Type <- "CH2O"
        Step5$form <- paste(Step5$C, Step5$H, Step5$O, Step5$N, Step5$S, Step5$P, Step5$E, Step5$S34,
                            Step5$N15, Step5$D, Step5$Cl, Step5$Fl, Step5$Cl37, Step5$Br, Step5$Br81,
                            Step5$I,Step5$M, Step5$NH4,
                            Step5$POE, Step5$NOE, sep = "_")
        Step5 <- Step5[-c(42)] ##CHECK THIS

        Out <- rbind(Step1, Step2, Step3, Step4, Step5)
        Out <- Out[(Out$C >= 2 & Out$H >= 4 & Out$O >= 0),]

        Out$H_C <- Out$H/Out$C   #Quick internal QA to limit bad assignments
        Out$O_C <- Out$O/Out$C
        Out <- Out[Out$H_C <= H_Cmax & Out$H_C >= H_Cmin &
                     Out$O_C <= O_Cmax & Out$O_C >= O_Cmin,]
        Out <- Out[!names(Out) %in% c("H_C", "O_C")]


        Out <- rbind(Out, DummyOut)

        Out_form <- dplyr::group_by(Out, Exp_mass, RA.x, form) #LCMS

        Out_form$number <- 1
        Out_form <- dplyr::summarize_at(Out_form, "number", sum, na.rm = TRUE)

        #Turns on or off ambiguity based on user input
        if(Ambigcheck == "off") {
          Out_form <- dplyr::filter(Out_form, number == max(number))
        }

        Out_form <- unique(Out_form)
        Out2<- merge(Out, Out_form, by.x = c("Exp_mass", "form", "RA.x"), by.y = c("Exp_mass", "form", "RA.x")) #LCMS
        Out3 <- dplyr::distinct(Out2, form, Exp_mass, RA.x, .keep_all = TRUE)  #LCMS
        Out3 <- Out3[!names(Out3) %in% c("number")]


        Next <- Out3[!names(Out3) %in% c("base_mass", "Type", "form", "RA.y")]
        colnames(Next)[colnames(Next) == "RA.x"] <- "RA"

        Next <- rbind(known, Next)
        Next <- unique(Next)
        Unambig <- Next

        masses <- Unambig[c("Exp_mass", "RA", "C")]  #LCMS
        names(masses)[3] <- "Var"   #LCMS
        Ambig <- merge(unknown, masses, by.x = c("Exp_mass", "RA"), by.y = c("Exp_mass", "RA"),all = T)   #LCMS
        Ambig <- Ambig[is.na(Ambig$Var),]
        Ambig <- Ambig[-19]  ##CHECK THIS
        Ambigreturn <- unique(Ambig)

        if(x == 2){
          break
        }
      }
      setTxtProgressBar(pb2, i)
      Unambig
    }
  }
  records1X <- Unambig[c(1:23)]  ##CHECK THIS

  records1X <- unique(records1X)

#check <- records1X %>% distinct(Exp_mass, .keep_all = TRUE)
  records1X$mode <- ionMode


  df1 <- records1X[records1X$mode == "pos" & records1X$M > 0 & records1X$POE == 0 & records1X$NH4 == 0,]
  df1$Neutral_mass <- df1$Exp_mass - df1$M * 22.989221

  df1x <- records1X[records1X$mode == "pos" & records1X$M == 0 & records1X$POE == 0 & records1X$NH4 == 1,]
  df1x$Neutral_mass <- df1x$Exp_mass - df1x$NH4 * 18.033823

  df2 <- records1X[records1X$mode == "pos" & records1X$M == 0 & records1X$POE == 0& records1X$NH4 == 0,]
  df2$Neutral_mass <- df2$Exp_mass - 1.00727645216

  df3 <- records1X[records1X$mode == "neg" & records1X$NOE == 0,]
  df3$Neutral_mass <- df3$Exp_mass + 1.00727645216

  df4 <- records1X[records1X$mode == "neg" & records1X$NOE == 1,]
  df4$Neutral_mass <- df4$Exp_mass - 0.000548597

  df5 <- records1X[records1X$mode == "pos" & records1X$M == 0 & records1X$POE == 1& records1X$NH4 == 0,]
  df5$Neutral_mass <- df5$Exp_mass + 0.000548597

  records1X <- rbind(df1, df1x, df2, df3, df4, df5)


  records1X <- records1X[-c(24)]  ##CHECK THIS

  #NOEx seems to be good to this point.
  ###Standard QA steps, second round
  #if(records1X$E > 0) {records1X$C <- records1X$C + records1X$E} # This is to fix the theor. masses
  records1X <- dplyr::mutate(records1X, O_C = O/(C+E), H_C =H/(C+E),

                             #Neutral_mass = Neutral_mass + POE * (2.0156500638/2)- NOE * (2.0156500638/2),


                             theor_mass1 = EM2("C") * C + EM2("H") * H + EM2("O") * O + N * EM2("N14") +
                               S * EM2("S") + P * EM2("P31") +
                               Cl * EM2("Cl35") + Fl * EM2("Fl19") + E * EM2("E2") + S34 * EM2("S34") +
                               Cl37 * EM2("Cl37m") + N15 * EM2("N15H") + Br * EM2("Br79") + Br81 * EM2("Br81m") +
                               I * EM2("I127") +
                               D * EM2("D") + M * EM2("M") + NH4 * EM2("NH4") - POE * electron + NOE*electron,

                             theor_mass = EM2("C") * C + EM2("H") * H + EM2("O") * O + N * EM2("N14") +
                               S * EM2("S") + P * EM2("P31") + Fl * EM2("Fl19") +Br * EM2("Br79") + Br81 * EM2("Br81m") +
                               I * EM2("I127") +
                               Cl * EM2("Cl35") +  E * EM2("E2") + S34 * EM2("S34") + Cl37 * EM2("Cl37m") +
                               N15 * EM2("N15H") +
                               D * EM2("D"),

                             #C = C + E, #It is added back so that formulas are more accurate.

                             DBE = C - 0.5 * (H + Cl + Cl37 +Fl +Br + Br81 + I) + 0.5 * (N +N15+ P) + 1,

                             err_ppm = ((Neutral_mass - theor_mass) / Neutral_mass * 10^6),

                             AE_ppm = round(abs((Neutral_mass - theor_mass) / Neutral_mass * 10^6),2))
  ######
  records1X$NM <- ifelse(abs(floor(records1X$Exp_mass)-records1X$Exp_mass) >= min_def & abs(floor(records1X$Exp_mass)-records1X$Exp_mass) <= max_def,
                         records1X$NM <- floor(records1X$Exp_mass), records1X$NM <- round(records1X$Exp_mass)) #New 1/6/20

  records1X$KM = records1X$Exp_mass * (14 / 14.01565)

  records1X$KMD <- ifelse(abs(floor(records1X$Exp_mass)-records1X$Exp_mass) >= min_def & abs(floor(records1X$Exp_mass)-records1X$Exp_mass) <= max_def,
                          records1X$KMD <- floor(records1X$Exp_mass)-records1X$KM, records1X$KMD <- round(records1X$Exp_mass)-records1X$KM) #New 1/6/20
  ##########

  records1X <- dplyr::mutate(records1X, max_LA = theor_mass1 / 13, actual_LA = ((C - E) + N + S + Fl + Br + Br81 + I + O + E + S34 + P + Cl +Cl37+N15) ,

                             rule_13= round(actual_LA/max_LA,1),

                             Senior1 = H + P + N + Cl + Fl + Cl37 + N15 + Br + Br81 + I  ,

                             STest = S + S34, ClTest = Cl + Cl37, NTest = N + N15, DBEO = DBE-O, BrTest = Br + Br81,

                             max_H = C * 2 + 2, H_test = round(H / max_H,1),

                             Senior2 = Pval + Nval + N15val + Valence("H")  +
                               Valence("Cl") + Valence("Cl37")+Valence("Fl") + Valence2("Br") + Valence2("Br81") +
                               Valence2("I"),

                             Senior3Atom = C + H + O + N + S + P + N15 + Cl + Cl37 + S34,

                             Senior3Val = C*Valence("C") + H*Valence("H") + O*Valence("O") + N*Nval +
                               S*Sval + P*Pval + S34*S34val + Fl*Valence("Fl")+
                               N15*N15val + Cl*Valence("Cl") + Cl37*Valence("Cl37") + Br * Valence2("Br") +
                               Br81 * Valence2("Br81") + I * Valence2("I")
  )


  #recordssave <- records1X
  #records1X <- recordssave

  records1X <- dplyr::filter(records1X, C>0, H>0,O>=Omin, H >= D)
  records1X <- unique(records1X)
  records1X <- dplyr::filter(records1X, O_C < O_Cmax & H_C <= H_Cmax & H_C > H_Cmin & O_C >= O_Cmin &

                               DBEO >= DBEOmin & DBEO <= DBEOmax &

                               !((ClTest) > (HighMoles("Cl",Cl=Clx) + HighMoles("Cl37",Cl37=Cl37x))) &
                               !((STest) > (HighMoles("S",S=Sx)+HighMoles("S34", S34 = S34x))) &
                               !((NTest) > (HighMoles("N",N=Nx) +HighMoles("N15",N15=N15x))) &
                               BrTest <= Brx + Br81x &

                               H_test <= 1 & rule_13 <= 1 &

                               AE_ppm <= ppm_err &

                               Even(Senior1)==TRUE & DBE >= 0 & DBE <= round(0.9 * (C + N)) &

                               O <= 2 * C + 3 * (N+N15) + 4 * P + 4 * (S+S34)&

                               O >= Omin&

                               RA > 0 &

                               Senior2 >= 2*Valence("C") &

                               Senior3Val >= (2*Senior3Atom - 1)
  )



  records1X <- records1X[!names(records1X) %in% c("Senior1", "Senior2", "Senior3Val", "Senior3Atom")]

  #record1sX <- records1X %>% filter(AE_ppm <= 2)
  #check <- record1sX %>% distinct(Exp_mass, .keep_all = TRUE)

  ###Formula generation
  records1X <-
    dplyr::mutate(records1X, Cform = ifelse(C == 0 , "",
                                            ifelse(C == 1 , "C", paste("C",C, sep = ""))),
                  Hform = ifelse(H == 0 , "",
                                 ifelse(H == 1 , "H", paste("H",H, sep = ""))),
                  Nform = ifelse(NTest == 0 , "",
                                 ifelse(NTest == 1 , "N", paste("N",NTest, sep = ""))),
                  Oform = ifelse(O == 0 , "",
                                 ifelse(O == 1 , "O", paste("O",O, sep = ""))),
                  Sform = ifelse(STest == 0 , "",
                                 ifelse(STest == 1 , "S", paste("S",STest, sep = ""))),
                  Pform = ifelse(P == 0 , "",
                                 ifelse(P == 1 , "P", paste("P",P, sep = ""))),
                  Clform = ifelse(ClTest == 0 , "",
                                  ifelse(ClTest == 1 , "Cl", paste("Cl",ClTest, sep = ""))),
                  Flform = ifelse(Fl == 0 , "",
                                  ifelse(Fl == 1 , "F", paste("F",Fl, sep = ""))),
                 Brform = ifelse(BrTest == 0 , "",
                                  ifelse(BrTest == 1 , "Br", paste("Br",BrTest, sep = ""))),
                 Iform = ifelse(I == 0 , "",
                                  ifelse(I == 1 , "I", paste("I",I, sep = ""))))

  records1X <- tidyr::unite(records1X, class, Nform, Oform, Sform, Pform, Clform, Flform, Brform, Iform,
                            sep = "", remove = FALSE)

  records1X <- tidyr::unite(records1X, formula, Cform, Hform, Nform, Oform, Sform, Pform, Clform,
                            Flform, Brform, Iform, sep = "")

  records1X <-
    dplyr::mutate(records1X, Cform = ifelse(C == 0 , "", "C"),
                  Hform = ifelse(H == 0 , "", "H"),

                  Nform = ifelse(NTest == 0 , "", "N"),

                  Oform = ifelse(O == 0 , "", "O"),

                  Sform = ifelse(STest == 0 , "","S"),
                  Pform = ifelse(P == 0 , "", "P"),
                  Clform = ifelse(ClTest == 0 , "", "Cl"),
                  Flform = ifelse(Fl == 0 , "", "F"),
                  Brform = ifelse(BrTest == 0 , "", "Br"),
                  Iform = ifelse(I == 0 , "", "I"))

  records1X <- tidyr::unite(records1X, group, Cform, Hform, Nform, Oform, Sform, Pform,
                            Clform, Flform, Brform, Iform, sep = "")

  ###Supplemental Specialized QA Steps
  records1X<-dplyr::mutate(records1X, HA = NTest + STest + P + ClTest + E + Fl + BrTest + I)

  records1X<-dplyr::group_by(records1X, Exp_mass, RA) #LCMS

  ifelse(HetCut == "on", records1X<-dplyr::filter(records1X, HA == (min(HA))), records1X<- records1X)

  records1X <- dplyr::group_by(records1X, Exp_mass, RA) #LCMS

  records1X <- dplyr::distinct(records1X, formula, RA, .keep_all = TRUE)  #Fix on 05/14/19

  records1X <- dplyr::ungroup(records1X)

  records1X <- records1X[!names(records1X) %in% c("HA")]



  #records3 <- dplyr::rename(records1X, mass = Exp_mass)

  cut <- (peaksAll2)
  cut <- unique(cut)
  cut$Test <- paste(cut$Exp_mass, cut$RA, sep = "_")   #New 05/14/19
  records1X$Test <- paste(records1X$Exp_mass, records1X$RA, sep = "_") #New 05/14/19
  unassigned <- dplyr::left_join(cut, records1X, by = "Test")  #New 05/14/19
  unassigned <- unassigned[is.na(unassigned$formula),]
  unassigned <- unassigned[c("RA.x", "Exp_mass.x")]
  names(unassigned)[1] <- "RA"
  names(unassigned)[2] <- "mass"
  unassigned <- unique(unassigned)  #Good to this point
  records1X <- records1X[-c(48)]  #removes Test   ###CHECK THIS

  records1X$mode <- ionMode


  df1 <- records1X[records1X$mode == "pos" & records1X$POE == 0,]
  df1$theor_mass1 <- df1$theor_mass1 + proton

  df2 <- records1X[records1X$mode == "pos" & records1X$POE == 1,]
  df2$theor_mass1 <- df2$theor_mass1

  df3 <- records1X[records1X$mode == "neg" & records1X$NOE == 0,]
  df3$theor_mass1 <- df3$theor_mass1 - proton

  df4 <- records1X[records1X$mode == "neg" & records1X$NOE == 1,]
  df4$theor_mass1 <- df4$theor_mass1

  records1X <- rbind(df1, df2, df3, df4)

  records1X <- records1X[-c(48)]  #removes mode     ###CHECK THIS

  records1X <- records1X[c(1,2,45:47,3:24,27,30:34, 29, 25:26,41, 35:37,43, 44)]  ###CHECK THIS

  records1 <- records1X

  rec_mass <- records1[c(2,1)]
  rec_mass$tag2 <- "Mono"
  names(rec_mass)[1] <- "Iso_mass"
  names(rec_mass)[2] <- "Iso_RA"  #LCMS

  isopeaks3 <- isopeaks2  ##New 12/11/19
  peaksAll3 <- peaksAll2
  names(peaksAll3)[2] <- "Iso_mass"
  names(peaksAll3)[1] <- "Iso_RA"
  peaksAll3$Tag <- "Mono2"
  isopeaks3 <- rbind(isopeaks3, peaksAll3)
  isopeaks3 <- isopeaks3[!duplicated(isopeaks3$Iso_mass),]

  final_iso <- merge(isopeaks3, rec_mass, by.x = c("Iso_mass", "Iso_RA"), by.y = c("Iso_mass", "Iso_RA"), all = T)##LCMS
  final_iso <- final_iso[is.na(final_iso$tag2),]
  final_iso <- final_iso[c(1:3)]

  C13 <- final_iso[final_iso$Tag == "C13" | final_iso$Tag == "C13_S34" | final_iso$Tag == "Mono2",]
  DC13 <- final_iso[final_iso$Tag == "2C13" | final_iso$Tag == "2C13_S34"| final_iso$Tag == "Mono2",]
  S34 <- final_iso[final_iso$Tag == "S34" | final_iso$Tag == "C13_S34"| final_iso$Tag == "2C13_S34"|
                     final_iso$Tag == "Mono2",]

  #recordssave <- records1
  #records1 <- recordssave
  #########C13 Isotope Alignment
 # iso_err = 1
  records1$Iso_mass <- records1$Exp_mass + 1.0033548380
  err <- iso_err * 10^-6
  records1$Iso_mass <- sapply(records1$Iso_mass, function(x){
    # First check if the element lies within tolerance limits of any element in df2
    ifelse(min(abs(C13$Iso_mass - x), na.rm=TRUE) < err * x,
           # If yes, replace that element in df1 with the matching element in df2
           C13[which.min(abs(C13$Iso_mass - x)), "Iso_mass"], 0)
  })

  #####New 12/13/15
  C13_2 <- final_iso[final_iso$Tag == "C13" | final_iso$Tag == "C13_S34" | final_iso$Tag == "Mono2",]
  names(C13_2)[1] <- "Iso_mass_2"
  #names(C13_2)[2] <- "Iso_RA_2"  #LCMS

  records1$Iso_mass_2 <- records1$theor_mass1 + 1.0033548380
  err <- iso_err * 10^-6
  records1$Iso_mass_2 <- sapply(records1$Iso_mass_2, function(x){
    # First check if the element lies within tolerance limits of any element in df2
    ifelse(min(abs(C13_2$Iso_mass_2 - x), na.rm=TRUE) < err * x,
           # If yes, replace that element in df1 with the matching element in df2
           C13_2[which.min(abs(C13_2$Iso_mass_2 - x)), "Iso_mass_2"], 0)
  })


  No_iso <- records1[records1$Iso_mass == 0 & records1$Iso_mass_2 == 0,]
  No_iso <- No_iso[-44]  ###CHECK THIS
  Both_iso <- records1[records1$Iso_mass != 0 & records1$Iso_mass_2 != 0,]
  Both_iso <- Both_iso[-44]  ###CHECK THIS
  Match_iso <- records1[records1$Iso_mass != 0 & records1$Iso_mass_2 == 0,]
  Match_iso <- Match_iso[-44]  ###CHECK THIS
  Assign_iso <- records1[records1$Iso_mass == 0 & records1$Iso_mass_2 != 0,]
  Assign_iso <- Assign_iso[-43]  ###CHECK THIS
  names(Assign_iso)[43] <- "Iso_mass"  ###CHECK THIS

  records1 <- rbind(No_iso, Both_iso, Match_iso, Assign_iso)



  ####
  records1 <- dplyr::left_join(records1, C13, by = "Iso_mass")
  #chceck <- records1 %>% filter(Iso_mass > 0)

  names(records1)[43] <- "C13_mass"   ###CHECK THIS
  names(records1)[44] <- "C13_Abund"  ###CHECK THIS
  C13_next <- records1[!is.na(records1$Tag),]
  C13_next <- C13_next[C13_next$Tag == "Mono2",]
  C13_next <- C13_next[c(43:45)]   ###LCMS
  names(C13_next)[1] <- "Iso_mass"   #Highlighting Mono peaks to not consider further
  names(C13_next)[2] <- "Iso_RA"  #LCMS
  records1 <- records1[-c(45)]   ###CHECK THIS

  #check <- records1 %>% filter(C13_mass > 0)
  ##########
  ## 2 C13 Isotope Assignment

  DC13 <- merge(DC13, C13_next, by.x = c("Iso_mass", "Iso_RA"), by.y = c("Iso_mass", "Iso_RA"), all = T)  #LCMS
  DC13 <- DC13[is.na(DC13$Tag.y),]
  DC13 <- DC13[c(1:3)]
  names(DC13)[3] <- "Tag"

  records1$Iso_mass <- records1$Exp_mass + 2 * 1.0033548380
  err <- iso_err * 10^-6
  records1$Iso_mass <- sapply(records1$Iso_mass, function(x){
    # First check if the element lies within tolerance limits of any element in df2
    ifelse(min(abs(DC13$Iso_mass - x), na.rm=TRUE) < err * x,
           # If yes, replace that element in df1 with the matching element in df2
           DC13[which.min(abs(DC13$Iso_mass - x)), "Iso_mass"], 0)
  })

  #####New 12/13/15
  DC13_2 <- DC13
  names(DC13_2)[1] <- "Iso_mass_2"
  names(DC13_2)[2] <- "Iso_RA_2"
  records1$Iso_mass_2 <- records1$theor_mass1 + 2* 1.0033548380
  err <- iso_err * 10^-6
  records1$Iso_mass_2 <- sapply(records1$Iso_mass_2, function(x){
    # First check if the element lies within tolerance limits of any element in df2
    ifelse(min(abs(DC13_2$Iso_mass_2 - x), na.rm=TRUE) < err * x,
           # If yes, replace that element in df1 with the matching element in df2
           DC13_2[which.min(abs(DC13_2$Iso_mass_2 - x)), "Iso_mass_2"], 0)
  })

  No_C13dummy <- records1[1,]
  #No_C13dummy <- data.frame(as.numeric(No_C13dummy))
  #No_C13dummy[!is.na(No_C13dummy)] <- 0
  No_C13dummy$group <- "Dummy"   #New  8/21/20

  No_iso <- records1[records1$Iso_mass == 0 & records1$Iso_mass_2 == 0,]
  No_iso <- No_iso[-46]  ###CHECK THIS
  Both_iso <- records1[records1$Iso_mass != 0 & records1$Iso_mass_2 != 0 & records1$C13_mass != 0,]
  Both_iso <- Both_iso[-46]  ###CHECK THIS
  Match_iso <- records1[records1$Iso_mass != 0 & records1$Iso_mass_2 == 0& records1$C13_mass != 0,]
  Match_iso <- Match_iso[-46]  ###CHECK THIS
  Assign_iso <- records1[records1$Iso_mass == 0 & records1$Iso_mass_2 != 0& records1$C13_mass != 0,]
  Assign_iso <- Assign_iso[-45]  ###CHECK THIS
  names(Assign_iso)[45] <- "Iso_mass"   ###CHECK THIS
  No_C13 <- records1[(records1$Iso_mass != 0 | records1$Iso_mass_2 != 0)& records1$C13_mass == 0,]
  No_C13 <- rbind(No_C13dummy, No_C13)
  No_C13$Iso_mass <- 0
  No_C13 <- No_C13[-46]   ###CHECK THIS
  records1 <- rbind(No_iso, Both_iso, Match_iso, Assign_iso, No_C13)



  ####
#########################Problem starting here.
  records1 <- dplyr::left_join(records1, DC13, by = "Iso_mass")
  #chceck <- records1 %>% filter(Iso_mass > 0)

  names(records1)[45] <- "C13_mass2"  ###CHECK THIS
  names(records1)[46] <- "C13_Abund2"  ###CHECK THIS

  records1$C13_mass2[records1$C13_mass == 0] <- 0
  records1$C13_Abund2[records1$C13_mass == 0] <- 0


  C13_next2 <- records1[!is.na(records1$Tag),]
  C13_next2 <- C13_next2[C13_next2$Tag == "Mono2",]
  C13_next2 <- C13_next2[c(45:47)]   ###LCMS
  C13_next2 <- C13_next2[C13_next2$C13_mass2 > 0,]
  names(C13_next2)[1] <- "Iso_mass"
  names(C13_next2)[2] <- "Iso_RA"  ###LCMS
  C13_next2 <- rbind(C13_next, C13_next2)

  records1 <- records1[-c(47)]   ###CHECK THIS
  ##########
  #THIS IS THE PROBLEM SECTION

  ## S34 Isotope Assignment
  if(Sx >0 & sum(records1$S) > 0){


    S34 <- merge(S34, C13_next2, by.x = c("Iso_mass", "Iso_RA"), by.y = c("Iso_mass", "Iso_RA"), all = T)  #LCMS
    S34 <- S34[is.na(S34$Tag.y),]
    S34 <- S34[c(1:3)]
    names(S34)[3] <- "Tag"



    records1K <- records1[records1$S ==0,]
    records1S <- records1[records1$S > 0,]
    records1S$Iso_mass <- records1S$Exp_mass + 1.995797
    err <- iso_err * 10^-6
    records1S$Iso_mass <- sapply(records1S$Iso_mass, function(x){
      # First check if the element lies within tolerance limits of any element in df2
      ifelse(min(abs(S34$Iso_mass - x), na.rm=TRUE) < err * x,
             # If yes, replace that element in df1 with the matching element in df2
             S34[which.min(abs(S34$Iso_mass - x)), "Iso_mass"], 0)
    })

    ###New 12/13/19
    S34_2 <- S34
    names(S34_2)[1] <- "Iso_mass_2"
    names(S34_2)[2] <- "Iso_RA"  #LCMS

    records1S$Iso_mass_2 <- records1S$theor_mass1 + 1.995797

    records1S$Iso_mass_2 <- sapply(records1S$Iso_mass_2, function(x){
      # First check if the element lies within tolerance limits of any element in df2
      ifelse(min(abs(S34_2$Iso_mass_2 - x), na.rm=TRUE) < err * x,
             # If yes, replace that element in df1 with the matching element in df2
             S34_2[which.min(abs(S34_2$Iso_mass_2 - x)), "Iso_mass_2"], 0)
    })

    records1S <- merge(records1S, S34_2, by.x = c("Iso_mass_2"), by.y = c("Iso_mass_2"), all = T) #LCMS
    records1S <- records1S[!is.na(records1S$DBEO),] #LCMS

    No_iso <- records1S[records1S$Iso_mass == 0 & records1S$Iso_mass_2 == 0,]
    No_iso <- No_iso[-1]    ###CHECK THIS  #LCMS
    Both_iso <- records1S[records1S$Iso_mass != 0 & records1S$Iso_mass_2 != 0,]
    Both_iso <- Both_iso[-1]    ###CHECK THIS  #LCMS
    Match_iso <- records1S[records1S$Iso_mass != 0 & records1S$Iso_mass_2 == 0,]
    Match_iso <- Match_iso[-1]   ###CHECK THIS  #LCMS
    Assign_iso <- records1S[records1S$Iso_mass == 0 & records1S$Iso_mass_2 != 0,]
    Assign_iso <- Assign_iso[-48]   ###CHECK THIS   #LCMS
    names(Assign_iso)[1] <- "Iso_mass"    ###CHECK THIS  #LCMS

    records1S <- rbind(No_iso, Both_iso, Match_iso, Assign_iso)


    S34dummy <- data.frame(Iso_mass = 0, Iso_RA = 0, Tag = "X")
    S34 <- rbind(S34,S34dummy)

    #records1Sdummy <- data.frame(Iso_mass = 1, Iso_RA = -42)  #LCMS
    records1Sdummy <- records1[1,]   #LCMS
    records1Sdummy$Iso_mass <- -42   #LCMS
    records1Sdummy$Iso_RA <- -42   #LCMS
    records1Sdummy$Tag <- "X"   #LCMS

    records1S <-rbind(records1S, records1Sdummy)
    records1S <- dplyr::left_join(records1S, S34, by = c("Iso_mass", "Iso_RA"))  #LCMS
    names(records1S)[47] <- "S34_mass"     ###CHECK THIS
    names(records1S)[48] <- "S34_Abund"    ###CHECK THIS
    records1S <- records1S[-c(49,50)]      ###CHECK THIS
    records1K$S34_mass <- 0
    records1K$S34_Abund <- 0
    records1 <- rbind(records1K, records1S)
    records1 <- records1[records1$S34_mass >= 0,]

  }else{
    records1$S34_mass <- 0;
    records1$S34_Abund <- 0
  }


  #New section to make sure the isotope abundance does not exceed the monoisotope abundance
#recordssave <- records1

  records1$C13_mass[records1$C13_Abund > (C13_abund/100) * records1$RA] <- 0
  records1$C13_Abund[records1$C13_Abund > (C13_abund/100) * records1$RA] <- 0

  records1$C13_mass2[records1$C13_Abund2 > (C13_abund/100) * records1$C13_Abund] <- 0
  records1$C13_Abund2[records1$C13_Abund2 > (C13_abund/100) * records1$C13_Abund] <- 0

  records1$S34_mass[records1$S34_Abund > (S34_abund/100) * records1$RA] <- 0
  records1$S34_Abund[records1$S34_Abund > (S34_abund/100) * records1$RA] <- 0

  records1 <- records1[!is.na(records1$C),]
  #check <- records1 %>% filter(S34_mass > 0 )
  ##########################################
#check <- records1 %>% filter(C13_mass > 0)
  ##############################################################
  records1$Test <- paste(records1$Exp_mass, records1$RA, sep = "_") #Added 05/14/19
  records1$Dups <- duplicated(records1$Test) #Changed 05/14/19
  records1$num <- 1:nrow(records1)
  records1 <- records1[order(-records1$num),]  #6/20/19
  records1$Dup2 <- duplicated(records1$Test)   #6/20/19
  records1 <- records1[records1$Exp_mass >= 16,]
  records1 <- records1[records1$group != "Dummy",]    #8/21/19

  Unambig <- records1[records1$Dups == FALSE & records1$Dup2 == FALSE,]
  Ambig <- records1[records1$Dups == TRUE | records1$Dup2 == TRUE,]
  Ambigout <- Ambig[-c(49:52)]   ###CHECK THIS
  Unambig <- Unambig[-c(49:52)]  ###CHECK THIS

  Ambigout2 <- data.frame(Exp_mass = 1)
  Ambigout <- dplyr::bind_rows(Ambigout, Ambigout2)
  Ambigout <- unique(Ambigout)
  Ambigout$Tag <- "Ambiguous"
  Unambig$Tag <- "Unambiguous"
  #Everything is good to this point, columns too, 05/14/19

  #############
  ##Nominal Mass Series QA Step
  if(NMScut == "on"){
    NewKMD <- dplyr::group_by(Ambigout, NM)
    NewKMD <- dplyr:: mutate(NewKMD, mDa = format(round((Exp_mass-min(Exp_mass))/0.0363855,2), nsmall = 2))

    NewKMD <- tidyr::separate(NewKMD, mDa, c("Whole", "Dec"), sep = -2)

    NewKMD$AddForm <- paste(NewKMD$formula, NewKMD$M, NewKMD$POE, NewKMD$NOE, sep = "_")
    NewKMD <- unique(NewKMD)

    Ambigout2 <- dplyr::group_by(NewKMD,NM, Dec, DBEO)
    Ambigout2 <- dplyr::mutate(Ambigout2, number = 1)

    Ambigout2 <- dplyr::summarize_at(Ambigout2, "number", sum, na.rm = TRUE)
    Ambigout2 <- dplyr::filter(Ambigout2, number == max(number))

    Ambigout3 <- merge(NewKMD, Ambigout2, by.x = c("NM", "DBEO", "Dec"), by.y = c("NM", "DBEO", "Dec"))

    Ambigout3$num <- 1:nrow(Ambigout3)
    Ambigout3$dups <- duplicated(Ambigout3$Exp_mass)
    Ambigout3 <-Ambigout3[order(-Ambigout3$num),]
    Ambigout3$dups2 <- duplicated(Ambigout3$Exp_mass)

    Unambig2 <- Ambigout3[Ambigout3$dups == FALSE & Ambigout3$dups2 == FALSE,]
    Ambig <- Ambigout3[Ambigout3$dups == TRUE | Ambigout3$dups2 == TRUE,]
    Ambig <- unique(Ambig)

    Ambigout <- Ambig[-c(3,51:56)]   ###CHECK THIS
    Ambigout <- Ambigout[c(3:32,1,33:37, 2, 38:49)]    ###CHECK THIS
    Unambigout <- Unambig2[-c(3,51:56)]     ###CHECK THIS
    Unambigout <- Unambigout[c(3:32,1,33:37, 2, 38:49)]    ###CHECK THIS
    Unambig <- rbind(Unambig, Unambigout)
  }



  #Columns still good
  # Unambig$theor_mass1 <- Unambig$theor_mass1 - Unambig$POE * 2.0156500638 +
  #   Unambig$NOE * 2.0156500638 + Unambig$NOE * electron - Unambig$POE * electron
  # #Everything is good to this point

  Unambig <- Unambig[Unambig$Exp_mass >= 16,]
  Unambig[is.na(Unambig)] <- 0
  Ambigout <- Ambigout[Ambigout$Exp_mass >= 16,]
  Ambigout[is.na(Ambigout)] <- 0
  records1[is.na(records1)] <- 0


  #check1 <- Unambig %>% filter(C13_mass > 0)

  ######################
  #N3OS to 13C conversion
  if(N3corr == "on"){

  Unambigkeep <- Unambig
  N3OS <- Unambig[(Unambig$group == "CHNOS"|Unambig$group == "CHNS") & Unambig$N==3,]#&Unambig$O_C < 0.5,]
  N3OSkeep <- Unambig[(Unambig$group == "CHNOS"|Unambig$group == "CHNS")& Unambig$N==3,]#&Unambig$O_C < 0.5,]
  Rest <- Unambig[!(((Unambig$group == "CHNOS"|Unambig$group == "CHNS"))& Unambig$N==3),]#&Unambig$O_C < 0.5),]
  N3OS$Cy <- N3OS$C - 2
  N3OS$Hy <- N3OS$H + 1
  N3OS$Ny <- N3OS$N - 3
  N3OS$Sy <- N3OS$S - 1
  N3OS$Oy <- N3OS$O + 6

  N3OS <-
    dplyr::mutate(N3OS, Cform = ifelse(Cy == 0 , "",
                                       ifelse(Cy == 1 , "C", paste("C",Cy, sep = ""))),
                  Hform = ifelse(Hy == 0 , "",
                                 ifelse(Hy == 1 , "H", paste("H",Hy, sep = ""))),
                  Nform = ifelse(Ny == 0 , "",
                                 ifelse(Ny == 1 , "N", paste("N",Ny, sep = ""))),
                  Oform = ifelse(Oy == 0 , "",
                                 ifelse(Oy == 1 , "O", paste("O",Oy, sep = ""))),
                  Sform = ifelse(Sy == 0 , "",
                                 ifelse(Sy == 1 , "S", paste("S",Sy, sep = "")))
    )

  N3OS <- tidyr::unite(N3OS, formula, Cform, Hform, Nform, Oform, Sform, sep = "")

  N3OSAlign <- N3OS[c("formula", "Exp_mass", "RA")] #Fix at end
  names(N3OSAlign)[2] <- "C13_mass1"
  names(N3OSAlign)[3] <- "C13_Abund1"

  Aligned <- merge(Unambig, N3OSAlign, by.x = "formula", by.y = "formula")
  Aligned <- Aligned[Aligned$RA > Aligned$C13_Abund1,]

  GoodN3OS <- merge(N3OSAlign, Aligned, by.x = c("C13_mass1", "C13_Abund1"),
                    by.y = c("C13_mass1", "C13_Abund1"), all = T) #LCMS
  GoodN3OS <- GoodN3OS[is.na(GoodN3OS$formula.y),]

  GoodN3OS <- GoodN3OS[c(1,2)]  #LCMS
  names(GoodN3OS)[1] <- "Exp_mass"
  names(GoodN3OS)[2] <- "RA"   #LCMS
  GoodN3OS <- merge(N3OSkeep, GoodN3OS, by.x = c("Exp_mass", "RA"), by.y = c("Exp_mass", "RA"))  #LCMS

  New13C <- Aligned[-c(43,44)]  ###CHECK THIS
  names(New13C)[48] <- "C13_mass"   ###CHECK THIS
  names(New13C)[49] <- "C13_Abund"  ###CHECK THIS

  Final13C <- rbind(New13C, Rest)
  Final13C <- Final13C[!duplicated(Final13C[c("formula", "Exp_mass", "RA")]),]

  Unambig <- rbind(GoodN3OS, Final13C)
  #This was fixed after new isotoping
  Unambig <- Unambig[c(2,3, 1,4:49)]   ###CHECK THIS

  ####Fix on 12/12/19, updated 12/13/19
#Save <- Unambig
  Unambig_C13 <- Unambig
  Unambig_C13$C13 <-  Unambig_C13$Exp_mass + 1.0033548380

  Unambig_C13$C13_2 <-  Unambig_C13$theor_mass1 + 1.0033548380

  Unambig_C13$err <-  round(abs((Unambig_C13$C13_mass- Unambig_C13$C13)/Unambig_C13$C13_mass*10^6),2)
  Unambig_C13$err_2 <-  round(abs((Unambig_C13$C13_mass- Unambig_C13$C13_2)/Unambig_C13$C13_mass*10^6),2)

  Unambig_rest <- Unambig_C13[Unambig_C13$C13_mass == 0,]
  Unambig_C13 <- Unambig_C13[Unambig_C13$C13_mass != 0,]


  Unambig_C13$C13_mass[Unambig_C13$err > iso_err & Unambig_C13$err_2 > iso_err] <- 0
  Unambig_C13$C13_Abund[Unambig_C13$err > iso_err& Unambig_C13$err_2 > iso_err] <- 0

  Unambig2 <- rbind(Unambig_rest, Unambig_C13)
  Unambig2 <- Unambig2[c(1:49)]    ###CHECK THIS
  N3OS_d <- Unambigkeep[(Unambigkeep$group == "CHNOS"|Unambigkeep$group == "CHNS") & Unambigkeep$N==3,]


  N3OS_check <- N3OS_d[c(1,3)]
  N3OS_merge <- merge(Unambig2, N3OS_check, by.x = "formula", by.y = "formula", all = T)
  N3OS_shift <- N3OS_merge[is.na(N3OS_merge$C),]
  N3OS_shift <- N3OS_shift[1]
  N3OS_rem <- merge(N3OS_d, N3OS_shift, by.x = "formula", by.y = "formula")
  N3OS_rem <- N3OS_rem[c(3,2,1)]  #LCMS
  N3OS_rem2 <- merge(Unambig2, N3OS_rem, by.x = c("C13_mass", "C13_Abund"),
                     by.y = c("Exp_mass", "RA"), all = T)  #LCMS
  N3OS_fin <- N3OS_rem2[!is.na(N3OS_rem2$formula.y) & is.na(N3OS_rem2$C),]
  N3OS_fin <- N3OS_fin[50]   ###CHECK THIS
  N3OS_fin2 <- merge(N3OS_d, N3OS_fin, by.x = "formula", by.y = "formula.y")

  #N3OS_merge <- merge(Unambigkeep, N3OS_check, by.x = "Exp_mass", by.y = "Exp_mass", all = T)
  #N3OS_good <- N3OS_merge[is.na(N3OS_merge$RA.y),]
  #N3OS_good <- N3OS_good[c(1:46)]
  #names(N3OS_good)[2] <- "RA"

  Unambig <- rbind(Unambig2, N3OS_fin2)
}
#########################



  ###########
  #Checking iso/mono matches
  #Likely Obsolete
  # Iso_check <- isopeaks2
  # names(Iso_check)[2] <- "Exp_mass"
  # Iso_check2_X <- merge(Unambig, Iso_check, by.x = "Exp_mass", by.y = "Exp_mass")


  ##############
  ##Plot data preparation
  Unambig <- Unambig[!is.na(Unambig$C)&Unambig$Exp_mass > 2,]
  Ambigout <- Ambigout[Ambigout$Exp_mass > 2,]

  #check1 <- Unambig %>% mutate(C13 = Exp_mass + 1.0033548380, err = abs((C13_mass- C13)/C13_mass*10^6))
  #check2 <- check1 %>% filter(err < 3)%>% mutate(C13_t = theor_mass1 + 1.0033548380, err_t = abs((C13_mass- C13_t)/C13_mass*10^6))



  #Final Unassigned Peaks
  P1 <- Unambig[c("Exp_mass", "RA")]
  P1$Test <- paste(P1$Exp_mass, P1$RA, sep = "_")
  names(P1)[1] <- "mass"
  names(P1)[2] <- "RA"

  P2 <- Unambig[c("C13_mass", "C13_Abund")]
  P2$Test <- paste(P2$C13_mass, P2$C13_Abund, sep = "_")
  names(P2)[1] <- "mass"
  names(P2)[2] <- "RA"

  P3 <- Unambig[c("C13_mass2", "C13_Abund2")]
  P3$Test <- paste(P3$C13_mass2, P3$C13_Abund2, sep = "_")
  names(P3)[1] <- "mass"
  names(P3)[2] <- "RA"

  P4 <- Unambig[c("S34_mass", "S34_Abund")]
  P4$Test <- paste(P4$S34_mass, P4$S34_Abund, sep = "_")
  names(P4)[1] <- "mass"
  names(P4)[2] <- "RA"

  P5 <- Ambigout[c("Exp_mass", "RA")]
  P5$Test <- paste(P5$Exp_mass, P5$RA, sep = "_")
  names(P5)[1] <- "mass"
  names(P5)[2] <- "RA"

  P6 <- Ambigout[c("C13_mass", "C13_Abund")]
  P6$Test <- paste(P6$C13_mass, P6$C13_Abund, sep = "_")
  names(P6)[1] <- "mass"
  names(P6)[2] <- "RA"

  P7 <- Ambigout[c("C13_mass2", "C13_Abund2")]
  P7$Test <- paste(P7$C13_mass2, P7$C13_Abund2, sep = "_")
  names(P7)[1] <- "mass"
  names(P7)[2] <- "RA"

  P8 <- Ambigout[c("S34_mass", "S34_Abund")]
  P8$Test <- paste(P8$S34_mass, P8$S34_Abund, sep = "_")
  names(P8)[1] <- "mass"
  names(P8)[2] <- "RA"

  AI <- isopeaks2[c(1,2)]
  names(AI)[2] <- "mass"
  names(AI)[1] <- "RA"
  AI$Test <- paste(AI$mass, AI$RA, sep = "_")
  AM <- peaksAll2
  names(AM)[2] <- "mass"
  names(AM)[1] <- "RA"
  AM$Test <- paste(AM$mass, AM$RA, sep = "_")

  AP <- rbind(AM, AI)
  AP <- unique(AP)

  GP <- rbind(P1, P2, P3, P4, P5, P6, P7, P8)
  GP <- unique(GP)
  GP$Tag <- "Good"

  UP <- merge(AP, GP, by.x = "Test", by.y = "Test", all = T)
  unassigned <- UP[is.na(UP$Tag),]
  unassigned <- unassigned[c(2,3)]
  names(unassigned)[1] <- "RA"
  names(unassigned)[2] <- "mass"

  ##############################################

  Ambigout <- Ambigout[Ambigout$group != "Dummy",]
  Unambig <- Unambig[Unambig$Tag != "Ambiguous",]
  PD <- rbind(Ambigout, Unambig)
  PDG <- subset( PD, group == "CHO"|group == "CHNO"|group == "CHOS"|group == "CHNOS"|group == "CH"|group == "CHN")
  PDB <- subset( PD, group != "CHO"&group != "CHNO"&group != "CHOS"&group != "CHNOS"&group != "CH"&group != "CHN")
  PDBdummy <- data.frame(Exp_mass = 1)
  PDB <- dplyr::bind_rows(PDB, PDBdummy)
  PDG$Tag2 <- PDG$group
  PDB$Tag2 <- "Other"
  PD <- rbind(PDG, PDB)
  PD <- PD[!is.na(PD$Tag),]
  records1 <- records1[records1$Exp_mass > 0,]
  records1 <- dplyr::distinct(records1, Test, .keep_all = TRUE)
  unassigned <- unassigned[unassigned$mass > 0,]


  ##########Palettes for plots
  group_colors <- data.frame(group = c("CHO", "CHNO", "CHOS", "CHNOS", "CH", "CHN", "Other"),
                             color = c("green", "blue", "red", "purple", "gold", "cyan", "grey67"))



  form_group <- data.frame(group = unique(PD$Tag2))

  form_palette <- merge(form_group, group_colors, by.x = "group", by.y = "group")
  form_palette <- setNames(form_palette$color, form_palette$group)
  #print(form_palette)
  ###############

  MZ<-ggplot2::ggplot() + ggplot2::geom_segment(data=records1, size=0.7,ggplot2::aes_string(x = "Exp_mass", xend = "Exp_mass", y = 0, yend = "RA"), color = "green")+
    ggplot2::geom_segment(data=records1, size=0.7,ggplot2::aes_string(x = "C13_mass", xend = "C13_mass", y = 0, yend = "C13_Abund"), color = "blue")+
    ggplot2::geom_segment(data=records1, size=0.7,ggplot2::aes_string(x = "C13_mass2", xend = "C13_mass2", y = 0, yend = "C13_Abund2"), color = "blue")+
    ggplot2::geom_segment(data=records1, size=0.7,ggplot2::aes_string(x = "S34_mass", xend = "S34_mass", y = 0, yend = "S34_Abund"), color = "blue")+
    ggplot2::geom_segment(data=unassigned, size=0.7,ggplot2::aes_string(x = "mass", xend = "mass", y = 0, yend = "RA"), color = "red")+
    ggplot2::coord_cartesian(xlim = c(min(rawpeaks$mass), max(rawpeaks$mass)))+
    ggplot2::theme_bw()+ggplot2::labs(x = "Ion Mass", y = "Abundance", title = "Assignment Mass Spectrum", color = "DBE")+
    ggplot2::theme(axis.title=ggplot2::element_text(size = 15, face = "bold"), strip.text=ggplot2::element_text(size=15,face="bold"),
                   axis.text=ggplot2::element_text(size=15, face = "bold"), legend.title=ggplot2::element_text(face="bold", size = 15),
                   legend.text=ggplot2::element_text(face="bold", size = 15),  panel.grid.minor.x=ggplot2::element_blank(),
                   panel.grid.major.x=ggplot2::element_blank(), strip.background = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(size = 16, face = "bold"))

  Error <- ggplot2::ggplot() + ggplot2::geom_point(data=Unambig, ggplot2::aes_string(x = "Exp_mass", y = "AE_ppm", color = "Tag"), alpha = 1/3) +
    ggplot2::geom_point(data=Ambigout, ggplot2::aes_string(x = "Exp_mass", y = "AE_ppm", color = "Tag"), alpha = 1/3) +
    ggplot2::coord_cartesian(xlim = c(min(records1$Exp_mass), max(records1$Exp_mass)), ylim = c(min(records1$AE_ppm), max(records1$AE_ppm))) +
    ggplot2::scale_colour_manual(name = "Ambiguity", values = c(Unambiguous = "blue", Ambiguous = "red")) +
    ggplot2::labs(x = "Ion Mass", y = "Absolute Error (ppm)", color = "Ambiguity", title = "Error Plot") + ggplot2::theme_bw() +
    ggplot2::theme(axis.title=ggplot2::element_text(size = 15, face = "bold"), strip.text=ggplot2::element_text(size=15,face="bold"),
                   axis.text=ggplot2::element_text(size=15, face = "bold"), legend.title=ggplot2::element_text(face="bold", size = 11),
                   legend.text=ggplot2::element_text(face="bold", size = 10),  panel.grid.minor.x=ggplot2::element_blank(),
                   panel.grid.major.x=ggplot2::element_blank(), strip.background = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(size = 16, face = "bold"))


  MZgroups<-ggplot2::ggplot() + ggplot2::geom_segment(data=PD, size=0.7,ggplot2::aes_string(x = "Exp_mass", xend = "Exp_mass", y = 0, yend = "RA", color = "Tag2"))+
    ggplot2::facet_wrap(~Tag, ncol = 1, scales = 'free_y')+
    ggplot2::scale_colour_manual(name = "Groups", values = form_palette) +
    ggplot2::theme_bw()+ggplot2::labs(x = "Ion Mass", y = "Abundance", title = "Assignment Mass Spectrum", color = "DBE")+
    ggplot2::theme(axis.title=ggplot2::element_text(size = 15, face = "bold"), strip.text=ggplot2::element_text(size=15,face="bold"),
                   axis.text=ggplot2::element_text(size=15, face = "bold"), legend.title=ggplot2::element_text(face="bold", size = 12),
                   legend.text=ggplot2::element_text(face="bold", size = 12),  panel.grid.minor.x=ggplot2::element_blank(),
                   panel.grid.major.x=ggplot2::element_blank(), strip.background = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(size = 16, face = "bold"))

  VK <- ggplot2::ggplot() + ggplot2::geom_point(data=PD, ggplot2::aes_string(x = "O_C", y = "H_C", color = "Tag2"), alpha = 1/3) +
    ggplot2::facet_wrap(~Tag, ncol = 2)+
    #ggplot2::coord_cartesian(xlim = c(min(PD$O_C), max(PD$O_C), ylim = c(min(PD$H_C), max(PD$H_C)))) +
    ggplot2::scale_colour_manual(name = "Groups", values = form_palette) +
    ggplot2::labs(x = "Oxygen-to-Carbon Ratio", y = "Hydrogen-to-Carbon Ratio", color = "Groups", title = "van Krevelen Plot") + ggplot2::theme_bw() +
    ggplot2::theme(axis.title=ggplot2::element_text(size = 15, face = "bold"), strip.text=ggplot2::element_text(size=15,face="bold"),
                   axis.text=ggplot2::element_text(size=15, face = "bold"), legend.title=ggplot2::element_text(face="bold", size = 11),
                   legend.text=ggplot2::element_text(face="bold", size = 10),  panel.grid.minor.x=ggplot2::element_blank(),
                   panel.grid.major.x=ggplot2::element_blank(), strip.background = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(size = 16, face = "bold"))

  colnames(Unambig)[colnames(Unambig)=="RA"] <- "abundance"
  colnames(Unambig)[colnames(Unambig)=="Exp_mass"] <- "exp_mass"
  colnames(Unambig)[colnames(Unambig)=="Neutral_mass"] <- "neutral_mass"
  colnames(Unambig)[colnames(Unambig)=="C13_Abund"] <- "C13_abund"
  colnames(Unambig)[colnames(Unambig)=="C13_Abund2"] <- "C13_abund2"
  colnames(Unambig)[colnames(Unambig)=="S34_Abund"] <- "S34_abund"
  colnames(Unambig)[colnames(Unambig)=="Tag"] <- "tag"
  CU12 <- Unambig[c("abundance", "exp_mass", "formula")]
  CURest <- Unambig[c(4:49)]  ###CHECK THIS
  Unambig <- cbind(CU12, CURest)
  names(Unambig)[28] <- "theor_mass"   ###CHECK THIS
  Unambig[is.na(Unambig)] <- 0

  colnames(Ambigout)[colnames(Ambigout)=="RA"] <- "abundance"
  colnames(Ambigout)[colnames(Ambigout)=="Exp_mass"] <- "exp_mass"
  colnames(Ambigout)[colnames(Ambigout)=="Neutral_mass"] <- "neutral_mass"
  colnames(Ambigout)[colnames(Ambigout)=="Tag"] <- "tag"
  colnames(Ambigout)[colnames(Ambigout)=="C13_Abund"] <- "C13_abund"
  colnames(Ambigout)[colnames(Ambigout)=="C13_Abund2"] <- "C13_abund2"
  colnames(Ambigout)[colnames(Ambigout)=="S34_Abund"] <- "S34_abund"
  CA12 <- Ambigout[c("abundance", "exp_mass", "formula")]
  CARest <- Ambigout[c(4:49)]   ###CHECK THIS
  Ambigout <- cbind(CA12, CARest)
  names(Ambigout)[28] <- "theor_mass"   ###CHECK THIS
  Unambig[is.na(Unambig)] <- 0

  names(unassigned)[1] <- "abundance"
  names(unassigned)[2] <- "exp_mass"
  unassigned <- unassigned[c(2,1)]
  unassigned <- unassigned[unassigned$abundance > SN,]

  #LCMS Changes
  if(cols == 3){
    peaksSave <- rbind(monoSave, isoSave)
    peaksSave <- unique(peaksSave)
    peaksC13 <- peaksSave
    peaks2C13 <- peaksSave
    peaksS34 <- peaksSave
    peaksUna <- peaksSave
    names(peaksC13)[1] <- "C13_mass"
    names(peaksC13)[2] <- "C13_abund"
    names(peaksC13)[3] <- "C13_RT"
    names(peaks2C13)[1] <- "C13_mass2"
    names(peaks2C13)[2] <- "C13_abund2"
    names(peaks2C13)[3] <- "C13_RT2"
    names(peaksS34)[1] <- "S34_mass"
    names(peaksS34)[2] <- "S34_abund"
    names(peaksS34)[3] <- "S34_RT"



    Unambig <- merge(Unambig, peaksSave, by.x = c("exp_mass", "abundance"), by.y = c("exp_mass", "abundance"))
    Unambig <- merge(Unambig, peaksC13, by.x = c("C13_mass", "C13_abund"), by.y = c("C13_mass", "C13_abund"), all = T)
    Unambig <- Unambig[!is.na(Unambig$exp_mass),]
    Unambig <- merge(Unambig, peaks2C13, by.x = c("C13_mass2", "C13_abund2"), by.y = c("C13_mass2", "C13_abund2"), all = T)
    Unambig <- Unambig[!is.na(Unambig$exp_mass),]
    Unambig <- merge(Unambig, peaksS34, by.x = c("S34_mass", "S34_abund"), by.y = c("S34_mass", "S34_abund"), all = T)
    Unambig <- Unambig[!is.na(Unambig$exp_mass),]

    Ambigout <- merge(Ambigout, peaksSave, by.x = c("exp_mass", "abundance"), by.y = c("exp_mass", "abundance"))
    Ambigout <- merge(Ambigout, peaksC13, by.x = c("C13_mass", "C13_abund"), by.y = c("C13_mass", "C13_abund"), all = T)
    Ambigout <- Ambigout[!is.na(Ambigout$exp_mass),]
    Ambigout <- merge(Ambigout, peaks2C13, by.x = c("C13_mass2", "C13_abund2"), by.y = c("C13_mass2", "C13_abund2"), all = T)
    Ambigout <- Ambigout[!is.na(Ambigout$exp_mass),]
    Ambigout <- merge(Ambigout, peaksS34, by.x = c("S34_mass", "S34_abund"), by.y = c("S34_mass", "S34_abund"), all = T)
    Ambigout <- Ambigout[!is.na(Ambigout$exp_mass),]

    unassigned <- merge(unassigned, peaksSave, by.x = c("exp_mass", "abundance"), by.y = c("exp_mass", "abundance"))

    Unambig <- Unambig[c(7,8,50,9:48,5,6,51,3,4,52,1,2,53,49)]  #Will need fixing
    Unambig[is.na(Unambig)] <- 0
    Ambigout <- Ambigout[c(7,8,50,9:48,5,6,51,3,4,52,1,2,53,49)]
    Ambigout[is.na(Ambigout)] <- 0
  }
  ##########

  #.rs.restartR()

  #Ambigout <- Ambigout[Ambigout$group != "Dummy",]
  ##Final Output list
  output <- list(Unambig = Unambig, Ambig = Ambigout, None = unassigned, MSAssign = MZ,
                 Error = Error, MSgroups = MZgroups, VK = VK)

  output

}






