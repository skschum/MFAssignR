#' Assigns all possible MF to each row of input data frame with CHOFIT algorithm
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
#' assigned. This is only really needed for APCI and APPI ionization modes.
#' This option is currently unavailable, but will be addressed in a future release.
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
#'  @param Fx numeric:
#' Sets the amount of Fluorine to be used in assignment. Default is 0.
#' @param Cl37x numeric:
#' Sets the amount of Chlorine 37 to be used in assignment. Default is 0.
#' @param Mx numeric:
#' Sets the amount of Sodium adduct to be used in assignment. Default is 0.
#' @param NH4x numeric:
#' Sets the amount of Ammonium adduct to be used in assignment. Default is 0.
#' @param Zx numeric:
#' Sets the amount of charge to be used in assignment. Default is 1.
#' @param Ox numeric:
#' Ox sets the maximum number of oxygen looked for in the CHOFIT core, it limits the number of loops performed
#' @param ppm_err numeric:
#' ppm_err parameter sets the error tolerance (ppm) for formula assignment. Default is 3.
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
#' @param HetCut character:
#' HetCut turns on or off the high heteroatom QA parameter. Default is "off"
#' @param NMScut character:
#' NMScut turns on or off the nominal mass series QA parameter. Default is “on”.
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
#' MFAssignAll_MSMS(peaks = Mono_df, isopeaks = Iso_df, "neg", lowMW = 200, highMW = 700)
#' MFAssignALl_MSMS(peaks = Mono_df, isopeaks = Iso_df, "neg", lowMW = 100, highMW = 1000, Nx = 3, Sx = 1)
#' @export
MFAssignAll_MSMS <- function(peaks, isopeaks = "none", ionMode, lowMW=100,highMW=1000, POEx = 0, NOEx = 0, Nx=0,Sx=0, Px=0, S34x=0,
                                  N15x=0, Dx=0,Ex=0, Clx=0, Fx = 0, Cl37x=0, Mx=0, NH4x=0, Zx=1, Ox = 30, ppm_err = 3, SN = 0, O_Cmin = 0,
                                  O_Cmax = 2.5, H_Cmin = 0.3, H_Cmax = 3, DBEOmin = -13, DBEOmax = 13, Omin = 0, HetCut = "off",
                                  NMScut = "on") {

  if(POEx >1) print('WARNING: Positive Odd Electron (POEx) is greater than 1, are you sure that is what you want?')
  if(NOEx >1) print('WARNING: Positive Odd Electron (NOEx) is greater than 1, are you sure that is what you want?')

  if(ionMode != "pos" & ionMode != "neg") print("WARNING: ionMode should be 'pos' or 'neg' ")

  #if(ionMode != "neg") print("WARNING: ionMode should be 'neg'")

  if(Nx > 5 | Sx > 5|Px >5|S34x>5|N15x >5|Dx>5|Ex>5|Clx > 5|Cl37x>5|Mx>5|NH4x>5|Fx > 5)
    print("WARNING: One or more heteroatoms are set greater than 5, this will cause the function to perform more slowly.")

  if(Ox !=30) print("WARNING: Ox is not at its default value, this will cause the core formula algorithm to perform additional
                      or fewer loops, are you sure you want it changed?")

  if(ppm_err > 3) print("WARNING: The maximum allowed error (ppm_err) is greater than 3, is this what you want?")

  # Constants
  components <- factor(c("C", "H", "O", "N", "S", "P", "Cl", "Fl", "E", "S34", "N15", "D", "Cl37",
                         "M", "NH4", "POE", "NOE", "Z"),
                       levels=c("C", "H", "O", "N", "S", "P", "Cl", "Fl",  "E", "S34", "N15", "D", "Cl37",
                                "M", "NH4", "POE", "NOE", "Z"))
  numComps <- length(components)
  proton = 1.00727645216
  electron =  0.000548597
  fitMode <- "ppm"
  maxErr <- ppm_err
  numDigits <- 6

  # Initialize data for analysis
  totFormulae <- 0

  peaks <- peaks[c(2,1)]

  names(peaks)[2] <- "mass"
  names(peaks)[1] <- "RA"

  isopeaks2 <- if(isopeaks != "none") isopeaks else data.frame(x=0,y=0, Tag = "X")# ,isopeaks = data.frame(x=0, y=0), isopeaks = data.frame(isopeaks))

  isopeaks2 <- isopeaks2[c(2,1,3)]

  names(isopeaks2)[2] <- "Iso_mass"
  names(isopeaks2)[1] <- "Iso_RA"
  names(isopeaks2)[3] <- "Tag"

  peaksAll <- peaks


  SNcut <- peaks[peaks$RA < SN,]
  peaks <- peaks[peaks$RA >= SN,]
  peaksAll2 <- peaks

  ################################# Inital Kendrick Series implementation
  peaks$KM <- peaks$mass* (14/14.01565)
  peaks$KMD <- round(peaks$mass)-peaks$KM
  peaks$zstar <- round(peaks$mass)%%14 - 14
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
  # Test <- dplyr::ungroup(Test)
  #################################
  Dummy <- data.frame(RA = c(-42,-42), mass = c(421.1147, 423.1293))
  peaks <- dplyr::bind_rows(Dummy, peaks)

  records <- vector('list')

  records2 <- vector('list')

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


    # Convert the (presumably) single-charged ion to a molecule

    if (ionMode=="neg") {
      exactEM <- ionEM + proton
    } else {
      exactEM <- ionEM - proton
    }



    # Check that exactEM is within bounds of LowMW and HighMW bounds
    if ((round(exactEM) >= lowMW) & (round(exactEM) <= highMW+1)) {

      # Start looping through components
      loop[CompFactorToInt("C")] <- 1 #LowMoles("C")
      loop[CompFactorToInt("H")] <- 4 #LowMoles("H")
      loop[CompFactorToInt("O")] <- 0 #LowMoles("O")
      loop[CompFactorToInt("Z")] <- 1 #LowMoles("Z")

      repeat {
        if (loop[CompFactorToInt("Z")] > 2) {
          exactEM = exactEM*(loop[CompFactorToInt("Z")] /
                               loop[CompFactorToInt("Z")] - 1)
        }

        loop[CompFactorToInt("NOE")] <- 0 #LowMoles("NOE")
        repeat {
        loop[CompFactorToInt("POE")] <- 0 #LowMoles("POE")
        repeat {

          loop[CompFactorToInt("NH4")] <- 0 #LowMoles("NH4")
          repeat {

            loop[CompFactorToInt("M")] <- 0 #LowMoles("M")
            repeat {

              loop[CompFactorToInt("Cl37")] <- 0 #LowMoles("Cl37")
              repeat {

                loop[CompFactorToInt("Fl")] <- 0 #LowMoles("Fl")
                repeat {

                loop[CompFactorToInt("Cl")] <- 0 #LowMoles("Cl")
                repeat {

                  loop[CompFactorToInt("D")] <- 0 #LowMoles("D")
                  repeat {

                    loop[CompFactorToInt("N15")] <- 0 #LowMoles("N15")
                    repeat {

                      loop[CompFactorToInt("S34")] <- 0 #LowMoles("S34")
                      repeat {

                        loop[CompFactorToInt("P")] <- 0 #LowMoles("P")
                        repeat {

                          loop[CompFactorToInt("S")] <- 0 #LowMoles("S")
                          repeat {

                            loop[CompFactorToInt("N")] <- 0 #LowMoles("N")
                            repeat {

                              loop[CompFactorToInt("E")] <- 0 #LowMoles("E")
                              repeat {

                                # strip the exact mass of all loop constituents to give CHO
                                # core of formula
                                coreXEM <- exactEM

                                ###Make sure this shouldn't start at E
                                for (step in CompFactorToInt("N"):CompFactorToInt("NOE")) {
                                  coreXEM = coreXEM - unlist(loop[step])*EM(CompIntToFactor(step))
                                }

                                if (coreXEM >= 16.0313) {

                                  coreRNM <- round(coreXEM)
                                  formulaOK <- FALSE

                                  if (Even(coreRNM)) {
                                    env <- environment()
                                    #if( exists(env) ) stop('Check mass list to make sure the masses are correct, or add more masses')
                                    FindCoreFormulae(env)

                                  } else {
                                    if (coreXEM >= 436.5008) {
                                      env2 <- environment()
                                      FindCoreFormulae(env2)

                                    }
                                  }




                                } else {
                                  formulaOK <- FALSE
                                }

                                env$records2[[length(env$records2) +1]] = env$records
                                #For List option

                                # ******##############################
                                # #####################################
                                loop[CompFactorToInt("E")] <- unlist(loop[CompFactorToInt("E")]) + 1
                                if (loop[CompFactorToInt("E")] > HighMoles("E", E=Ex)) {
                                  break
                                }
                              } # E loop

                              loop[CompFactorToInt("N")] <- unlist(loop[CompFactorToInt("N")]) + 1
                              if (loop[CompFactorToInt("N")] > HighMoles("N",N=Nx)) {
                                break
                              }
                            } # N loop

                            loop[CompFactorToInt("S")] <- unlist(loop[CompFactorToInt("S")]) + 1
                            if (loop[CompFactorToInt("S")] > HighMoles("S",S=Sx)) {
                              break
                            }
                          } # S loop

                          loop[CompFactorToInt("P")] <- unlist(loop[CompFactorToInt("P")]) + 1
                          if (loop[CompFactorToInt("P")] > HighMoles("P",P=Px)) {
                            break
                          }
                        } # P loop

                        loop[CompFactorToInt("S34")] <- unlist(loop[CompFactorToInt("S34")]) + 1
                        if (loop[CompFactorToInt("S34")] > HighMoles("S34",S34=S34x)) {
                          break
                        }
                      } # S34 loop

                      loop[CompFactorToInt("N15")] <- unlist(loop[CompFactorToInt("N15")]) + 1
                      if (loop[CompFactorToInt("N15")] > HighMoles("N15",N15=N15x)) {
                        break
                      }
                    } # N15 loop

                    loop[CompFactorToInt("D")] <- unlist(loop[CompFactorToInt("D")]) + 1
                    if (loop[CompFactorToInt("D")] > HighMoles("D",D=Dx)) {
                      break
                    }
                  } # D loop

                  loop[CompFactorToInt("Cl")] <- unlist(loop[CompFactorToInt("Cl")]) + 1
                  if (loop[CompFactorToInt("Cl")] > HighMoles("Cl",Cl=Clx)) {
                    break
                  }
                } # Cl loop

                loop[CompFactorToInt("Fl")] <- unlist(loop[CompFactorToInt("Fl")]) + 1
                if (loop[CompFactorToInt("Fl")] > HighMoles("Fl",Fl=Fx)) {
                  break
                }
                } # Fl loop

                loop[CompFactorToInt("Cl37")] <- unlist(loop[CompFactorToInt("Cl37")]) + 1
                if (loop[CompFactorToInt("Cl37")] > HighMoles("Cl37",Cl37=Cl37x)) {
                  break
                }
              }

              loop[CompFactorToInt("M")] <- unlist(loop[CompFactorToInt("M")]) + 1
              if (loop[CompFactorToInt("M")] > HighMoles("M", M=Mx)) {
                break
              }
            } # M Loop

            loop[CompFactorToInt("NH4")] <- unlist(loop[CompFactorToInt("NH4")]) + 1
            if (loop[CompFactorToInt("NH4")] > HighMoles("NH4", NH4=NH4x)) {
              break
            }
          }
          loop[CompFactorToInt("POE")] <- unlist(loop[CompFactorToInt("POE")]) + 1
          if (loop[CompFactorToInt("POE")] > HighMoles("POE", POE=POEx)) {
            break
          }
        } # Positive Odd Electron loop

        loop[CompFactorToInt("NOE")] <- unlist(loop[CompFactorToInt("NOE")]) + 1
        if (loop[CompFactorToInt("NOE")] > HighMoles("NOE", NOE=NOEx)) {
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
#.rs.restartR()

    recordsdf <- data.frame((do.call('rbind', records2)))

    recordsdf <- data.frame(RA = unlist(recordsdf$RA), coreNM = unlist(recordsdf$coreNM),

                            Exp_mass = unlist(recordsdf$Exp_mass),
                            C = unlist(recordsdf$C), H = unlist(recordsdf$H),
                            O = unlist(recordsdf$O), N = unlist(recordsdf$N),
                            S = unlist(recordsdf$S), P = unlist(recordsdf$P),
                            E = unlist(recordsdf$E), S34 = unlist(recordsdf$S34),
                            N15 = unlist(recordsdf$N15),
                            D = unlist(recordsdf$D), Cl = unlist(recordsdf$Cl),
                            Fl = unlist(recordsdf$Fl), Cl37 = unlist(recordsdf$Cl37),
                            M = unlist(recordsdf$M), NH4 = unlist(recordsdf$NH4),POE = unlist(recordsdf$POE),
                            NOE = unlist(recordsdf$NOE),
                            Z = unlist(recordsdf$Z), Neutral_mass = unlist(recordsdf$Neutral_mass),
                            CHO_mass = unlist(recordsdf$CHO_mass), CHO_Err = unlist(recordsdf$CHO_Err),
                            Ratio = unlist(recordsdf$Ratio))

    records1 <- dplyr::mutate(env$recordsdf, C = C+1*Ratio, H = H+4*Ratio+N+N15+P+2*POE+Cl+Cl37
                              + Fl - 2*NOE, O = O-1*Ratio)

    records1$KM <- records1$Exp_mass* (14/14.01565)
    records1$KMD <- round(records1$Exp_mass)-records1$KM
    records1$zstar <- round(records1$Exp_mass)%%14 - 14
    records1$KMDTest <- round(records1$KMD, 3)
    records1Isol <- records1 #New

########################
    ###############################################
    ###S34 isotope check QA
    if(isopeaks != "none"){
      SIso <- isopeaks2[isopeaks2$Tag == "S34"|isopeaks2$Tag == "C13_S34"|isopeaks2$Tag == "2C13_S34",]
      SIso <- unlist(SIso[2])
      recordsS <- records1[records1$S > 0,]
      recordsS <- unlist(recordsS[3])
      rest <- records1[records1$S == 0,]

      Sulf <- expand.grid(recordsS, SIso)
      names(Sulf)[1] <- "Exp_mass"
      names(Sulf)[2] <- "Iso_mass"

      Sulf$mdiff <- Sulf$Exp_mass - Sulf$Iso_mass
      Sulf <- Sulf[Sulf$mdiff > -2 & Sulf$mdiff < -1.98,]


      Sulf$KM2 <- Sulf$Iso_mass * (2 / 1.995797)
      Sulf$KMD2 <- round((round(Sulf$Iso_mass) - Sulf$KM2),3)
      Sulf$KM <- Sulf$Exp_mass * (2 / 1.995797)
      Sulf$KMD <- round((round(Sulf$Exp_mass) - Sulf$KM),3)

      Sulf$KMDdiff <- Sulf$KMD - Sulf$KMD2

      Sulf <- Sulf[abs(Sulf$KMDdiff) < 0.002,]
      Sulf <- Sulf[1]  #Pulls out the Exp_mass
      Sulfdata <- records1[records1$S > 0,]
      Sulfur <- merge(Sulf, Sulfdata, by.x = "Exp_mass", by.y = "Exp_mass")


      restalign <- rest[c(1,3)]  #Select the Exp_mass and RA
      restalign <- merge(Sulfur, restalign, by.x = "Exp_mass", by.y = "Exp_mass", all = T)
      restmass <- restalign[is.na(restalign$RA.x),]
      restmass <- restmass[c(1,30)]
      restfinal <- merge(rest, restmass, by.x = "Exp_mass", by.y = "Exp_mass", all = T)
      restfinal <- restfinal[!is.na(restfinal$RA.y),]
      restfinal <- restfinal[-30]

      records1 <- rbind(restfinal, Sulfur)
      records1 <- records1[c(2,3,1, 4:29)]
      records1 <- unique(records1)
    }
    ########################################

    records1 <- merge(records1, peaksend, by.x = c("zstar", "KMDTest"), by.y = c("zstar", "KMDTest"))
    records1$CH2_num <- round((records1$mass_CH2 - records1$Exp_mass)/14)

    records1 <- dplyr::mutate(records1, Cr = C+1*CH2_num, Hr = H+2*CH2_num, Or = O, Nr = N, Sr = S,
                              Pr = P, Dr = D,
                              S34r = S34, N15r = N15, Er = E, Mr = M, NH4r = NH4, Zr = Z,
                              POEr = POE, NOEr = NOE,
                              Clr = Cl, Flr = Fl, Cl37r = Cl37)
    Isolated <- records1Isol[!records1Isol$Exp_mass %in% records1$Exp_mass,] #New



    CH2 <- records1[c(1,2,30:50)]
    CH2 <- dplyr::rename(CH2, C= Cr, H = Hr, O=Or, N=Nr, D = Dr, S=Sr, P= Pr, S34 = S34r, N15 = N15r, E = Er, M=Mr, NH4 = NH4r,
                         Z=Zr, POE = POEr, NOE = NOEr, Cl = Clr, Fl = Flr, Cl37 = Cl37r, Exp_mass = mass_CH2, RA = RA_CH2)

    if (ionMode=="neg") {
      CH2$Neutral_mass <- CH2$Exp_mass + proton
    } else {
      CH2$Neutral_mass <- CH2$Exp_mass - proton
    }

    records1 <- records1[c(1:29)]
    records1 <- dplyr::bind_rows(records1, CH2)
    records1 <- dplyr::bind_rows(records1, Isolated)



    records1 <- dplyr::mutate(records1, O_C = O/(C+E), H_C =H/(C+E),

                              Neutral_mass = Neutral_mass + POE * 2.0156500638 - NOE * 2.0156500638,

           theor_mass1 = EM("C") * C + EM("H") * H + EM("O") * O + N * EM("N14") + S * EM("S") +
             P * EM("P31") +
           Cl * EM("Cl35") + Fl * EM("Fl19")+ E * EM("E") + S34 * EM("S34") + Cl37 * EM("Cl37m") +
             N15 * EM("N15H") +
           D * EM("D") + M * EM("M") + NH4 * EM("NH4") +POE * EM("POE") + NOE * EM("NOE"),

           theor_mass = EM("C") * C + EM("H") * H + EM("O") * O + N * EM("N14") + S * EM("S") +
             P * EM("P31") +
             Cl * EM("Cl35") + Fl * EM("Fl19") + E * EM("E") + S34 * EM("S34") + Cl37 * EM("Cl37m") +
             N15 * EM("N15H") +
             D * EM("D"),

           C = C + E,

           DBE = C - 0.5 * (H + Cl + Cl37 +Fl) + 0.5 * (N + P) + 1,

           err_ppm = ((Neutral_mass - theor_mass1) / Neutral_mass * 10^6),

           AE_ppm = abs((Neutral_mass - theor_mass1) / Neutral_mass * 10^6),

           NM = round(Exp_mass),

           KM = Exp_mass * (14 / 14.01565), KMD = round(Exp_mass) - KM,

           max_LA = theor_mass1 / 13, actual_LA = ((C - E) + N + S + O + E + S34 + P + Cl +Cl37+N15+ Fl) ,

           rule_13=actual_LA/max_LA,

           Senior1 = H + P + N + Cl + Cl37 + N15+Fl  ,

           STest = S + S34, ClTest = Cl + Cl37, NTest = N + N15, DBEO = DBE-O,

           max_H = C * 2 + 2, H_test = H / max_H,

           Senior2 = Valence("P") + Valence("N") + Valence("N15") + Valence("H")  + Valence("Cl")
           + Valence("Cl37")+ Valence("Fl"),
           Senior3Atom = C + H + O + N + S + P + N15 + E + Cl + Cl37 + S34 + Fl,
           Senior3Val = C*Valence("C") + H*Valence("H") + O*Valence("O") + N*Valence("N") +
             S*Valence("S") + P*Valence("P") + S34*Valence("S34") +
             N15*Valence("N15") + Cl*Valence("Cl") + Cl37*Valence("Cl37")+ Fl*Valence("Fl")
           )



#recordssave <- records1

     records1 <- dplyr::filter(records1, C>0, H>0,O>=Omin, H >= D)
     records1 <- unique(records1)
   records1 <- dplyr::filter(records1, O_C < O_Cmax & H_C <= H_Cmax & H_C > H_Cmin & O_C >= O_Cmin &

                               DBEO >= DBEOmin & DBEO <= DBEOmax &

                                !((ClTest) > HighMoles("Cl",Cl=Clx)) & !((STest) > HighMoles("S",S=Sx)) &

                               !((NTest) > HighMoles("N",N=Nx)) &

                               H_test <= 1 & rule_13 <= 1 &

                                AE_ppm <= ppm_err &

                               Even(Senior1)==TRUE & DBE >= 0 & DBE <= 0.9 * (C + N) &

                               O <= 2 * C + 3 * (N+N15) + 4 * P + 4 * (S+S34)&

                               #O>= P*4&

                               RA > 0 &

                               Senior2 >= 2*Valence("C") &

                               Senior3Val >= (2*Senior3Atom - 1)
                            )


   ###################################################

  records1 <- records1[-c(42,49:51)] #Removes the Senior columns

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
                                  ifelse(Fl == 1 , "F", paste("F",Fl, sep = ""))))

  records1 <- tidyr::unite(records1, class, Nform, Oform, Sform, Pform, Clform, Flform, sep = "", remove = FALSE)

  records1 <- tidyr::unite(records1, formula, Cform, Hform, Nform, Oform, Sform, Pform, Clform, Flform, sep = "")

  records1 <-
    dplyr::mutate(records1, Cform = ifelse(C == 0 , "", "C"),
                  Hform = ifelse(H == 0 , "", "H"),

                  Nform = ifelse(NTest == 0 , "", "N"),

                  Oform = ifelse(O == 0 , "", "O"),

                  Sform = ifelse(STest == 0 , "","S"),
                  Pform = ifelse(P == 0 , "", "P"),
                  Clform = ifelse(ClTest == 0 , "", "Cl"),
                  Flform = ifelse(Fl == 0 , "", "F"))

  records1 <- tidyr::unite(records1, group, Cform, Hform, Nform, Oform, Sform, Pform, Clform, Flform, sep = "")


  records1<-dplyr::mutate(records1, HA = NTest + STest + P + ClTest + E)

  records1<-dplyr::group_by(records1, Exp_mass)

  ifelse(HetCut == "On", records1<-dplyr::filter(records1, HA == (min(HA))), records1<- records1)

  records1 <- dplyr::group_by(records1, Exp_mass)

  records1 <- dplyr::distinct(records1, formula, .keep_all = TRUE)

  records1 <- dplyr::ungroup(records1)

  records1 <- dplyr::select(records1, -c(coreNM, CHO_mass, CHO_Err, Ratio, HA))
###################################################################


  records3 <- dplyr::rename(records1, mass = Exp_mass)

  cut <- dplyr::bind_rows(peaksAll, SNcut)
  unassigned <- dplyr::left_join(cut, records3, by = "mass")
  unassigned <- unassigned[is.na(unassigned$formula),]  #dplyr::filter(unassigned, is.na(formula))
  unassigned <- dplyr::select(unassigned, RA.x, mass)
  unassigned <- dplyr::rename(unassigned, RA = RA.x)
  unassigned <- unassigned[unassigned$RA >=0 & unassigned$RA < SN,] #dplyr::filter(unassigned, RA >= 0)
  unassignedH2O <- unassigned[unassigned$RA >= SN,]


  #################################
  #H2O series

  H2O_una <- unassignedH2O
  names(H2O_una)[2] <- "mass_H2O"
  names(H2O_una)[1] <- "RA_H2O"

  H2O_una$KMH2O <- H2O_una$mass_H2O* (18/18.0105)
  H2O_una$KMDH2O <- round(H2O_una$mass_H2O)-H2O_una$KMH2O
  H2O_una$zstarH2O <- round(H2O_una$mass_H2O)%%18 - 18
  H2O_una$KMDTestH2O <- round(H2O_una$KMDH2O, 3)

  records1$KMH2O <- records1$Exp_mass* (18/18.0105)
  records1$KMDH2O <- round(records1$Exp_mass)-records1$KMH2O
  records1$zstarH2O <- round(records1$Exp_mass)%%18 - 18
  records1$KMDTestH2O <- round(records1$KMDH2O, 3)

  records13 <- merge(records1, H2O_una, by.x = c("zstarH2O", "KMDTestH2O"), by.y = c("zstarH2O", "KMDTestH2O"))
  IsolatedH2O <- records1[!records1$Exp_mass %in% records13$Exp_mass,]
  records13$H2O_num <- round((records13$mass_H2O-records13$Exp_mass)/18)




  records13 <- dplyr::mutate(records13, Cr = C, Hr = H+2*H2O_num, Or = O + H2O_num, Nr = N, Sr = S, Pr = P, Dr = D,
                            S34r = S34, N15r = N15, Er = E, Mr = M, NH4r = NH4, Zr = Z, POEr = POE,
                            Clr = Cl, Cl37r = Cl37, Flr = Fl, NOEr = NOE)

  H2O <- records13[c(1:4,25:73)]
  H2O <- dplyr::rename(H2O, C= Cr, H = Hr, O=Or, N=Nr, D = Dr, S=Sr, P= Pr, S34 = S34r, N15 = N15r, E = Er, M=Mr, NH4 = NH4r,
                       Z=Zr, POE = POEr, Cl = Clr, Cl37 = Cl37r, NOE = NOEr, Fl = Flr, Exp_mass = mass_H2O, RA = RA_H2O)

  if (ionMode=="neg") {
    H2O$Neutral_mass <- H2O$Exp_mass + proton
  } else {
    H2O$Neutral_mass <- H2O$Exp_mass - proton
  }
  records13 <- dplyr::bind_rows(records13, H2O)
  records13 <- dplyr::bind_rows(records13, IsolatedH2O)

  records14 <- records13[c(5:25, 1:4)]
  ##################################
  records1 <- records14
  records1 <- dplyr::mutate(records1, O_C = O/(C), H_C =H/(C),

                            Neutral_mass = Neutral_mass + POE * (-2.0156500638/2) + NOE * (2.0156500638/2),

                            theor_mass1 = EM("C") * C + EM("H") * H + EM("O") * O + N * EM("N14") + S * EM("S") + P * EM("P31") +
                              Cl * EM("Cl35") + Fl*EM("Fl19") +  E * EM("E2") + S34 * EM("S34") + Cl37 * EM("Cl37m") + N15 * EM("N15H") +
                              D * EM("D") + M * EM("M") + NH4 * EM("NH4"), #+POE * EM("POE"),

                            theor_mass = EM("C") * C + EM("H") * H + EM("O") * O + N * EM("N14") + S * EM("S") + P * EM("P31") +
                              Cl * EM("Cl35") + Fl*EM("Fl19") +  E * EM("E2") + S34 * EM("S34") + Cl37 * EM("Cl37m") + N15 * EM("N15H") +
                              D * EM("D"),

                            DBE = C - 0.5 * (H + Cl + Cl37 + Fl) + 0.5 * (N +N15+ P) + 1,

                            #C = C + E,

                            err_ppm = ((Neutral_mass - theor_mass1) / Neutral_mass * 10^6),

                            AE_ppm = abs((Neutral_mass - theor_mass1) / Neutral_mass * 10^6),

                            NM = round(Exp_mass),

                            KM = Exp_mass * (14 / 14.01565), KMD = round(Exp_mass) - KM,

                            max_LA = theor_mass1 / 13, actual_LA = ((C - E) + N + S + O + E +
                                                                      S34 + P + Cl +Cl37+N15+Fl) ,

                            rule_13=actual_LA/max_LA,

                            Senior1 = H + P + N + Cl + Cl37 + N15+ Fl  ,

                            STest = S + S34, ClTest = Cl + Cl37, NTest = N + N15, DBEO = DBE-O,

                            max_H = C * 2 + 2, H_test = H / max_H,

                            Senior2 = Valence("P") + Valence("N") + Valence("N15") + Valence("H")  +
                              Valence("Cl") + Valence("Cl37") + Valence("Fl"),
                            Senior3Atom = C + H + O + N + S + P + N15 + E + Cl + Cl37 + S34 + Fl,
                            Senior3Val = C*Valence("C") + H*Valence("H") + O*Valence("O") + N*Valence("N") +
                              S*Valence("S") + P*Valence("P") + S34*Valence("S34") +
                              N15*Valence("N15") + Cl*Valence("Cl") + Cl37*Valence("Cl37") + Fl*Valence("Fl")
  )






  records1 <- dplyr::filter(records1, O_C < O_Cmax & H_C <= H_Cmax & H_C > H_Cmin & O_C >= O_Cmin &

                              DBEO >= DBEOmin & DBEO <= DBEOmax &

                              !((ClTest) > HighMoles("Cl",Cl=Clx)) & !((STest) > HighMoles("S",S=Sx)) &

                              !((NTest) > HighMoles("N",N=Nx)) &

                              H_test <= 1 & rule_13 <= 1 &

                              AE_ppm <= ppm_err &

                              Even(Senior1)==TRUE & DBE >= 0 & DBE <= 0.9*(C+N) &

                              O <= 2 * C + 3 * (N+N15) + 4 * P + 4 * (S+S34)&

                              #O>= P*4&

                              RA > 0 &

                              Senior2 >= 2*Valence("C") &

                              Senior3Val >= (2*Senior3Atom - 1)

  )



  records1 <- dplyr::select(records1, -Senior1, -Senior2, -Senior3Val, -Senior3Atom)

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
                                  ifelse(Fl == 1 , "F", paste("F",Fl, sep = ""))))

  records1 <- tidyr::unite(records1, class, Nform, Oform, Sform, Pform, Clform, Flform, sep = "", remove = FALSE)

  records1 <- tidyr::unite(records1, formula, Cform, Hform, Nform, Oform, Sform, Pform, Clform, Flform, sep = "")

  records1 <-
    dplyr::mutate(records1, Cform = ifelse(C == 0 , "", "C"),
                  Hform = ifelse(H == 0 , "", "H"),

                  Nform = ifelse(NTest == 0 , "", "N"),

                  Oform = ifelse(O == 0 , "", "O"),

                  Sform = ifelse(STest == 0 , "","S"),
                  Pform = ifelse(P == 0 , "", "P"),
                  Clform = ifelse(ClTest == 0 , "", "Cl"),
                  Flform = ifelse(Fl == 0 , "", "F"))

  records1 <- tidyr::unite(records1, group, Cform, Hform, Nform, Oform, Sform, Pform, Clform, Flform, sep = "")


  records1<-dplyr::mutate(records1, HA = NTest + STest + P + ClTest + E)

  records1<-dplyr::group_by(records1, Exp_mass)

  ifelse(HetCut == "on", records1<-dplyr::filter(records1, HA == (min(HA))), records1<- records1)

  records1 <- dplyr::group_by(records1, Exp_mass)

  records1 <- dplyr::distinct(records1, formula, .keep_all = TRUE)

  records1 <- dplyr::ungroup(records1)

  records1 <- dplyr::select(records1, -c( HA))

  records1 <- records1[c(1,2,45,46,47,3:20,28,29,31:35,30, 26,27, 42, 36:38,43:44)]

  #records3 <- dplyr::rename(records1, mass = Exp_mass)

  cut <- dplyr::bind_rows(peaksAll)
  cut <- unique(cut)
  names(cut)[2] <- "Exp_mass"
  unassigned <- dplyr::left_join(cut, records1, by = "Exp_mass")
  unassigned <- unassigned[is.na(unassigned$formula),]
  unassigned <- unassigned[c("RA.x", "Exp_mass")]
  names(unassigned)[1] <- "RA"
  names(unassigned)[2] <- "mass"
  unassigned <- unique(unassigned)  #Good to this point
###################################################

  ##Aligning Isotope masses back into the mass spectrum
  ##Align single C13 masses
  records1$C13_mass <- records1$Exp_mass + 1.0033548380
  err <- ppm_err*10^-6
  C13Iso <- isopeaks2[isopeaks2$Tag == "C13"|isopeaks2$Tag == "C13_S34",]
  names(C13Iso)[2] <- "C13_mass"
  names(C13Iso)[1] <- "C13_Abund"
  records1$C13_mass <- sapply(records1$C13_mass, function(x){
    # First check if the element lies within tolerance limits of any element in df2
    ifelse(min(abs(C13Iso$C13_mass - x), na.rm=TRUE) < err * x,
           # If yes, replace that element in df1 with the matching element in df2
           C13Iso[which.min(abs(C13Iso$C13_mass - x)), "C13_mass"], 0)
  })


  records1 <- dplyr::left_join(records1, C13Iso, by = "C13_mass")
  records1 <- records1[-42] #Removes Tag column
  #########
  #Align double C13 masses
  records1$C13_mass2 <- records1$Exp_mass + 2.006709676
  err <- ppm_err*10^-6
  C13Iso2 <- isopeaks2[isopeaks2$Tag == "2C13" | isopeaks2$Tag == "2C13_S34",]
  names(C13Iso2)[2] <- "C13_mass2"
  names(C13Iso2)[1] <- "C13_Abund2"
  records1$C13_mass2 <- sapply(records1$C13_mass2, function(x){
    # First check if the element lies within tolerance limits of any element in df2
    ifelse(min(abs(C13Iso2$C13_mass2 - x), na.rm=TRUE) < err * x,
           # If yes, replace that element in df1 with the matching element in df2
           C13Iso2[which.min(abs(C13Iso2$C13_mass2 - x)), "C13_mass2"], 0)
  })

  #Ccheck <- records1[records1$C13_mass > 0,]
  records1 <- dplyr::left_join(records1, C13Iso2, by = "C13_mass2")
  records1 <- records1[!(records1$C13_mass == 0 & records1$C13_mass2 > 0),]
  records1 <- records1[-44] #Removes Tag column
  ###########
  #Align S34 masses  #This is a potential error zone if Sx >0 and no sulfur assigned, or if Sx >0 and
  #no non-sulfur are assigned. Standard conditions should not have a problem.
  if(Sx > 0){
    recordsdummy <- records1[1,]
    recordsdummy[!is.na(recordsdummy)] <- NA
    recordsSulf <- records1[records1$S >0,]
    recordsSulf <- rbind(recordsSulf, recordsdummy)
    recordsrest <- records1[records1$S ==0,]
    recordsrest <- rbind(recordsrest, recordsdummy)
    recordsrest$S34_mass <- 0
    recordsrest$S34_Abund <- 0
    recordsSulf$S34_mass <- recordsSulf$Exp_mass + 1.995797
    err <- ppm_err*10^-6

    S34Iso <- isopeaks2[isopeaks2$Tag == "S34" | isopeaks2$Tag == "C13_S34" | isopeaks2$Tag == "2C13_S34",]
    names(S34Iso)[2] <- "S34_mass"
    names(S34Iso)[1] <- "S34_Abund"
    recordsSulf$S34_mass <- sapply(recordsSulf$S34_mass, function(x){
      # First check if the element lies within tolerance limits of any element in df2
      ifelse(min(abs(S34Iso$S34_mass - x), na.rm=TRUE) < err * x,
             # If yes, replace that element in df1 with the matching element in df2
             S34Iso[which.min(abs(S34Iso$S34_mass - x)), "S34_mass"], 0)
    })


    recordsSulf <- dplyr::left_join(recordsSulf, S34Iso, by = "S34_mass")
    recordsSulf <- recordsSulf[-46]
    records1 <- rbind(recordsrest, recordsSulf)
  } else{records1$S34_mass <- 0;
  records1$S34_Abund <- 0}

  ###########
  #Checking to see if any polyisotope masses match an assigned monoisotope mass
  Iso_check <- isopeaks2
  names(Iso_check)[2] <- "Exp_mass"
  Iso_check2 <- merge(records1, Iso_check, by.x = "Exp_mass", by.y = "Exp_mass")
  Mono_check2 <- merge(records1, Iso_check, by.x = "Exp_mass", by.y = "Exp_mass", all = T)
  Mono_check2 <- Mono_check2[!is.na(Mono_check2$C)& is.na(Mono_check2$Tag),]

  I1 <- Mono_check2[c("C13_mass")]
  names(I1)[1] <- "Exp_mass"
  I2 <- Mono_check2[c("C13_mass2")]
  names(I2)[1] <- "Exp_mass"
  I3 <- Mono_check2[c("S34_mass")]
  names(I3)[1] <- "Exp_mass"

  IM <- rbind(I1, I2, I3)
  IM$Tag2 <- "Iso"
  IM <- unique(IM)

  Iso_check3 <- merge(Iso_check2, IM, by.x = "Exp_mass", by.y = "Exp_mass", all = T)
  Iso_check3 <- Iso_check3[!is.na(Iso_check3$C),]
  MonoG1 <- Iso_check3[Iso_check3$C13_mass > 0 | Iso_check3$C13_mass2 > 0 | Iso_check3$S34_mass > 0,]
  MonoG2 <- Iso_check3[is.na(Iso_check3$Tag2) & (Iso_check3$C13_mass == 0 &
                                                   Iso_check3$C13_mass2 == 0 & Iso_check3$S34_mass == 0),]
  MonoGF <- rbind(MonoG1, MonoG2)
  MonoGF <- MonoGF[c(1:45)]
  MonoRest <- Mono_check2[c(1:45)]
  records1 <- rbind(MonoRest, MonoGF)
  records1 <- records1[c(2,1,3:45)]
  ######################################################################
  #ID unmatched isotope masses.
  C13 <- records1[c(41,40)]
  names(C13)[1] <- "RA"
  names(C13)[2] <- "Iso_mass"
  C13_2 <- records1[c(43,42)]
  names(C13_2)[1] <- "RA"
  names(C13_2)[2] <- "Iso_mass"
  S34 <- records1[c(45,44)]
  names(S34)[1] <- "RA"
  names(S34)[2] <- "Iso_mass"

  Isomass <- rbind(C13, C13_2, S34)
  Isomass <- Isomass[Isomass$Iso_mass > 0,]

  unassignedIso <- merge(isopeaks2, Isomass, by.x = "Iso_mass", by.y = "Iso_mass", all = T)
  unassignedIso <- unassignedIso[is.na(unassignedIso$RA),]

  unassignedIso <- unassignedIso[c(2,1)]
  names(unassignedIso)[2] <- "mass"
  names(unassignedIso)[1] <- "RA"

  ######################################
#Old Isotope alignment
  # records1$Iso_mass <- records1$Exp_mass + 1.0033548380
  # err <- ppm_err*10^-6
  # records1$Iso_mass <- sapply(records1$Iso_mass, function(x){
  #   # First check if the element lies within tolerance limits of any element in df2
  #   ifelse(min(abs(isopeaks2$Iso_mass - x), na.rm=TRUE) < err * x,
  #          # If yes, replace that element in df1 with the matching element in df2
  #          isopeaks2[which.min(abs(isopeaks2$Iso_mass - x)), "Iso_mass"], 0)
  # })
  # #sapply(records1, class)
  # records1 <- dplyr::left_join(records1, isopeaks2, by = "Iso_mass")
  # unassignedIso <- dplyr::left_join(isopeaks2, records1, by = "Iso_mass")
  # unassignedIso <- unassignedIso[is.na(unassignedIso$Exp_mass),]
  # unassignedIso <- unassignedIso[c(1,2)]
  #
  # names(unassignedIso)[2] <- "mass"
  # names(unassignedIso)[1] <- "RA"
   unassigned <- rbind(unassigned, unassignedIso)
#########################################################
  records1$Dups <- duplicated(records1$Exp_mass)
  records2 <- records1[records1$Dups == FALSE,]
  Ambig <- records1[records1$Dups == TRUE,]
  Ambigcheck <- Ambig[c(2)]
  Ambigdummy <- data.frame(Exp_mass = 1)
  Ambigcheck <- rbind(Ambigcheck, Ambigdummy)
  Ambigcheck$Tag <- "Ambig"

  Ambigfinal <- dplyr::left_join(records2, Ambigcheck, by = "Exp_mass")
  Unambig <- Ambigfinal[is.na(Ambigfinal$Tag),]
  Unambig <- Unambig[-c(46,47)]
  Ambigfinal <- Ambigfinal[!is.na(Ambigfinal$Tag),]
  Ambigfinal <- Ambigfinal[-47]
  Ambigout <- rbind(Ambig, Ambigfinal)
  Ambigout <- Ambigout[-46]
  Ambigout2 <- data.frame(Exp_mass = 1)
  Ambigout <- dplyr::bind_rows(Ambigout, Ambigout2)
  Ambigout <- unique(Ambigout)
  Ambigout$Tag <- "Ambiguous"
  Unambig$Tag <- "Unambiguous"

  #############
  if(NMScut == "on"){
  NewKMD <- dplyr::group_by(Ambigout, NM)
  NewKMD <- dplyr:: mutate(NewKMD, mDa = format(round((Exp_mass-min(Exp_mass))/0.0363855,2), nsmall = 2))

  NewKMD <- tidyr::separate(NewKMD, mDa, c("Whole", "Dec"), sep = -2)

  NewKMD$AddForm <- paste(NewKMD$formula, NewKMD$M, sep = "_")
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

  Ambigout <- Ambig[-c(3,48:53)]
  Ambigout <- Ambigout[c(3:27,1,28:32, 2, 33:46)]
  Unambigout <- Unambig2[-c(3,48:53)]
  Unambigout <- Unambigout[c(3:27,1,28:32, 2, 33:46)]
  Unambig <- rbind(Unambig, Unambigout)
  }


  #######################
  #Mass correction Unambig
  Unambig$mode <- ionMode

  df1 <- Unambig[Unambig$mode == "pos" & Unambig$M == 0,]
   df1$theor_mass <- df1$theor_mass + proton

   df2 <- Unambig[Unambig$mode == "neg",]
   df2$theor_mass <- df2$theor_mass - proton

   #Unambig6$mode <- ionMode
   #Unambig <- Unambig6

    df3 <- Unambig[Unambig$mode == "pos" & Unambig$M > 0,]
    df3$theor_mass <- df3$theor_mass +  df3$M*22.989221

   Unambig <- rbind(df1, df2, df3)
   #Unambig <- records1[-c(41)]


  df1 <- Unambig[Unambig$mode == "pos" & Unambig$M > 0,]
  df1$Neutral_mass <- df1$Exp_mass - df1$M * 22.989221

  df2 <- Unambig[Unambig$mode == "pos" & Unambig$M == 0,]
  df2$Neutral_mass <- df2$Exp_mass - 1.00727645216

  df3 <- Unambig[Unambig$mode == "neg",]
  df3$Neutral_mass <- df3$Exp_mass + 1.00727645216

  Unambig <- rbind(df1, df2, df3)
  Unambig <- Unambig[-c(47)]

  #######################
  #Mass correction Ambigout
  Ambigdummy2 <- data.frame(Exp_mass = 21)
  Ambigout <- dplyr::bind_rows(Ambigout, Ambigdummy2)
  Ambigout$mode <- ionMode

  df1 <- Ambigout[Ambigout$mode == "pos" & Ambigout$M == 0,]
  df1$theor_mass <- df1$theor_mass + proton

  df2 <- Ambigout[Ambigout$mode == "neg",]
  df2$theor_mass <- df2$theor_mass - proton

   df3 <- Ambigout[Ambigout$mode == "pos" & Ambigout$M > 0,]
   df3$theor_mass <- df3$theor_mass +  df3$M*22.989221

  Ambigout <- rbind(df1, df2)
  #Ambigout <- records1[-c(41)]


  df1 <- Ambigout[Ambigout$mode == "pos" & Ambigout$M > 0,]
  df1$Neutral_mass <- df1$Exp_mass - df1$M * 22.989221

  df2 <- Ambigout[Ambigout$mode == "pos" & Ambigout$M == 0,]
  df2$Neutral_mass <- df2$Exp_mass - 1.00727645216

  df3 <- Ambigout[Ambigout$mode == "neg",]
  df3$Neutral_mass <- df3$Exp_mass + 1.00727645216

  Ambigout <- rbind(df1, df2, df3)
  Ambigout <- Ambigout[-c(47)]

  Unambig <- Unambig[!is.na(Unambig$N15),]
  Unambig$theor_mass <- Unambig$theor_mass - Unambig$POE * (2.0156500638/2) +
     Unambig$NOE * (2.0156500638/2) - Unambig$NOE * electron + Unambig$POE * electron

  Unambig <- Unambig[c(1:23, 47, 25:46)]
  Ambigout <- Ambigout[c(1:23, 47, 25:46)]
  ###########

  Unambig <- Unambig[!is.na(Unambig$C)&Unambig$Exp_mass > 2,]
  Ambigout <- Ambigout[Ambigout$Exp_mass > 2,]

  #Final Unassigned Peaks
  P1 <- Unambig[c("Exp_mass")]
  names(P1)[1] <- "mass"
  P2 <- Unambig[c("C13_mass")]
  names(P2)[1] <- "mass"
  P3 <- Unambig[c("C13_mass2")]
  names(P3)[1] <- "mass"
  P4 <- Unambig[c("S34_mass")]
  names(P4)[1] <- "mass"
  P5 <- Ambigout[c("Exp_mass")]
  names(P5)[1] <- "mass"
  P6 <- Ambigout[c("C13_mass")]
  names(P6)[1] <- "mass"
  P7 <- Ambigout[c("C13_mass2")]
  names(P7)[1] <- "mass"
  P8 <- Ambigout[c("S34_mass")]
  names(P8)[1] <- "mass"

  AI <- isopeaks2[c(1,2)]
  names(AI)[2] <- "mass"
  names(AI)[1] <- "RA"
  AM <- peaksAll2
  names(AM)[2] <- "mass"
  names(AM)[1] <- "RA"

  AP <- rbind(AM, AI)
  AP <- unique(AP)

  GP <- rbind(P1, P2, P3, P4, P5, P6, P7, P8)
  GP <- unique(GP)
  GP$Tag <- "Good"

  UP <- merge(AP, GP, by.x = "mass", by.y = "mass", all = T)
  unassigned <- UP[is.na(UP$Tag),]
  unassigned <- unassigned[c(1,2)]

  ##############################################
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
  unassigned <- unassigned[unassigned$mass > 0,]

  Ambigout <- Ambigout[!is.na(Ambigout$RA),]
  Unambig <- Unambig[!is.na(Unambig$RA),]
  records1 <- records1[!is.na(records1$RA),]



  MZ<-ggplot2::ggplot() + ggplot2::geom_segment(data=records1, size=0.7,ggplot2::aes_string(x = "Exp_mass", xend = "Exp_mass", y = 0, yend = "RA"), color = "green")+
    ggplot2::geom_segment(data=records1, size=0.7,ggplot2::aes_string(x = "C13_mass", xend = "C13_mass", y = 0, yend = "C13_Abund"), color = "blue")+
    ggplot2::geom_segment(data=records1, size=0.7,ggplot2::aes_string(x = "C13_mass2", xend = "C13_mass2", y = 0, yend = "C13_Abund2"), color = "blue")+
    ggplot2::geom_segment(data=records1, size=0.7,ggplot2::aes_string(x = "S34_mass", xend = "S34_mass", y = 0, yend = "S34_Abund"), color = "blue")+
    ggplot2::geom_segment(data=unassigned, size=0.7,ggplot2::aes_string(x = "mass", xend = "mass", y = 0, yend = "RA"), color = "red")+
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
    ggplot2::scale_colour_manual(name = "Groups", values = c(CHO = "green", CHNO = "blue", CHOS = "red", CHNOS = "purple",
                                                             CH = "gold", CHN = "cyan", Other = "grey67")) +
    ggplot2::theme_bw()+ggplot2::labs(x = "Ion Mass", y = "Abundance", title = "Assignment Mass Spectrum", color = "DBE")+
    ggplot2::theme(axis.title=ggplot2::element_text(size = 15, face = "bold"), strip.text=ggplot2::element_text(size=15,face="bold"),
                   axis.text=ggplot2::element_text(size=15, face = "bold"), legend.title=ggplot2::element_text(face="bold", size = 12),
                   legend.text=ggplot2::element_text(face="bold", size = 12),  panel.grid.minor.x=ggplot2::element_blank(),
                   panel.grid.major.x=ggplot2::element_blank(), strip.background = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(size = 16, face = "bold"))

  VK <- ggplot2::ggplot() + ggplot2::geom_point(data=PD, ggplot2::aes_string(x = "O_C", y = "H_C", color = "Tag2"), alpha = 1/3) +
    ggplot2::facet_wrap(~Tag, ncol = 2)+
    #ggplot2::coord_cartesian(xlim = c(min(PD$O_C), max(PD$O_C), ylim = c(min(PD$H_C), max(PD$H_C)))) +
    ggplot2::scale_colour_manual(name = "Groups", values = c(CHO = "green", CHNO = "blue", CHOS = "red", CHNOS = "purple",
                                                             CH = "gold", CHN = "cyan", Other = "grey67")) +
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
  CURest <- Unambig[c(4:46)]
  Unambig <- cbind(CU12, CURest)
  names(Unambig)[25] <- "theor_mass"
  Unambig[is.na(Unambig)] <- 0

  colnames(Ambigout)[colnames(Ambigout)=="RA"] <- "abundance"
  colnames(Ambigout)[colnames(Ambigout)=="Exp_mass"] <- "exp_mass"
  colnames(Ambigout)[colnames(Ambigout)=="Neutral_mass"] <- "neutral_mass"
  colnames(Ambigout)[colnames(Ambigout)=="Tag"] <- "tag"
  colnames(Ambigout)[colnames(Ambigout)=="C13_Abund"] <- "C13_abund"
  colnames(Ambigout)[colnames(Ambigout)=="C13_Abund2"] <- "C13_abund2"
  colnames(Ambigout)[colnames(Ambigout)=="S34_Abund"] <- "S34_abund"
  CA12 <- Ambigout[c("abundance", "exp_mass", "formula")]
  CARest <- Ambigout[c(4:46)]
  Ambigout <- cbind(CA12, CARest)
  names(Ambigout)[25] <- "theor_mass"
  Unambig[is.na(Unambig)] <- 0

  names(unassigned)[2] <- "abundance"
  unassigned <- unassigned[unassigned$abundance > SN,]
  .rs.restartR()


  output <- list(Unambig = Unambig, Ambig = Ambigout, None = unassigned, MSAssign = MZ,
                 Error = Error, MSgroups = MZgroups, VK = VK)

  output


}



