#' Assigns all possible CHO molecular formulae to each row of input data frame
#'
#' MFAssignCHO() assigns all possible molecular formulae to each
#' mass in the input file, subject to user constraints on the moles of
#' C, H, O
#'
#' This function operates in the same way as MFAssign() but only assigns
#' carbon, hydrogen, and oxygen formulas, heteroatoms are not included.
#' Basic QA steps are included within the function,more detail about
#' these QA steps can be seen in the vignette and user manual attached to this package.
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
#' @param Ex numeric:
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
#' HetCut turns on or off the high heteroatom QA parameter. Default is "off".
#' @param NMScut character:
#' NMScut turns on or off the nominal mass series QA parameter. Default is “on”.
#' @param DeNovo numeric:
#' DeNovo sets the de novo cut point for the data. Default is 1000.
#' @param nLoop numeric:
#' nLoops sets the number of times the KMD and z* series assignment loops. Default is 5.
#' @param Ambig character:
#' Turns on or off increased ambiguity for assignments. Default is "off".
#' @param MSMS character:
#' Turns on or off CH2 KMD prescreening before initial assignment. Default is "off".
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
#' @examples
#' MFAssignCHO(peaks = Mono_df, isopeaks = Iso_df, "neg", lowMW = 200, highMW = 700)
#' MFAssignCHO(peaks = Mono_df, isopeaks = Iso_df, "neg", lowMW = 100, highMW = 1000)
#' @export

MFAssignCHO_RMD <- function(peaks, isopeaks = "none", ionMode, lowMW=100,highMW=1000, POEx = 0, NOEx = 0, Ex = 0, Mx = 0, NH4x = 0, Zx=1, Ox = 30,
                          ppm_err = 3, SN = 0, O_Cmin = 0, O_Cmax = 2.5, H_Cmin = 0.3, H_Cmax = 3, DBEOmin = -13, DBEOmax = 13,
                          Omin = 0, HetCut = "off", NMScut = "on", DeNovo = 1000, nLoop = 5, Ambig = "off", MSMS = "off") {

  if(POEx >1) print('WARNING: Positive Odd Electron (POEx) is greater than 1, are you sure that is what you want?')
  if(NOEx >1) print('WARNING: Positive Odd Electron (NOEx) is greater than 1, are you sure that is what you want?')

  if(ionMode != "pos" & ionMode != "neg") print("WARNING: ionMode should be 'pos' or 'neg' ")


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

  isopeaks2 <- if(isopeaks != "none") isopeaks else data.frame(x=0,y=0,Tag = 0)

  isopeaks2 <- isopeaks2[c(2,1,3)]

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
    peaks$KMD <- round(peaks$mass)-peaks$KM
    peaks$zstar <- round(peaks$mass)%%14 - 14
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

  }

  #################################
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

        loop[CompFactorToInt("POE")] <- 0 #LowMoles("POE")
        repeat {

        loop[CompFactorToInt("POE")] <- 0 #LowMoles("POE")
        repeat {

          loop[CompFactorToInt("NH4")] <- 0 #LowMoles("NH4")
          repeat {

            loop[CompFactorToInt("M")] <- 0 #LowMoles("M")
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


    recordsdf <- data.frame((do.call('rbind', records2)))

    recordsdf <- data.frame(RA = unlist(recordsdf$RA), coreNM = unlist(recordsdf$coreNM),

                            Exp_mass = unlist(recordsdf$Exp_mass),
                            C = unlist(recordsdf$C), H = unlist(recordsdf$H),
                            O = unlist(recordsdf$O), N = unlist(recordsdf$N),
                            S = unlist(recordsdf$S), P = unlist(recordsdf$P),
                            E = unlist(recordsdf$E), S34 = unlist(recordsdf$S34),
                            N15 = unlist(recordsdf$N15),
                            D = unlist(recordsdf$D), Cl = unlist(recordsdf$Cl),
                            Fl = unlist(recordsdf$Fl),
                             Cl37 = unlist(recordsdf$Cl37),
                            M = unlist(recordsdf$M), NH4 = unlist(recordsdf$NH4),POE = unlist(recordsdf$POE),
                            NOE = unlist(recordsdf$NOE),
                            Z = unlist(recordsdf$Z), Neutral_mass = unlist(recordsdf$Neutral_mass),
                            CHO_mass = unlist(recordsdf$CHO_mass), CHO_Err = unlist(recordsdf$CHO_Err),
                            Ratio = unlist(recordsdf$Ratio))

    records1 <- dplyr::mutate(env$recordsdf, C = C+1*Ratio, H = H+4*Ratio+N+P+2*POE+Cl+Cl37-2*NOE, O = O-1*Ratio)

    records1 <- dplyr::mutate(records1, O_C = O/(C+E), H_C =H/(C+E),

                              Neutral_mass = Neutral_mass + POE * 2.0156500638 - NOE * 2.0156500638,


                              theor_mass1 = EM("C") * C + EM("H") * H + EM("O") * O + N * EM("N14") + S * EM("S") + P * EM("P31") +
                                Cl * EM("Cl35") + Fl * EM("Fl19") + E * EM("E") + S34 * EM("S34") + Cl37 * EM("Cl37m") + N15 * EM("N15H") +
                                D * EM("D") + M * EM("M") + NH4 * EM("NH4") +POE * EM("POE")+NOE * EM("NOE"),

                              theor_mass = EM("C") * C + EM("H") * H + EM("O") * O + N * EM("N14") + S * EM("S") + P * EM("P31") +
                                Cl * EM("Cl35") + Fl * EM("Fl19") + E * EM("E") + S34 * EM("S34") + Cl37 * EM("Cl37m") + N15 * EM("N15H") +
                                D * EM("D"),

                              C = C + E,

                              DBE = C - 0.5 * (H + Cl + Cl37) + 0.5 * (N + N15 + P) + 1,

                              err_ppm = ((Neutral_mass - theor_mass1) / Neutral_mass * 10^6),

                              AE_ppm = abs((Neutral_mass - theor_mass1) / Neutral_mass * 10^6),

                              NM = round(Exp_mass),

                              KM = Exp_mass * (14 / 14.01565), KMD = round(Exp_mass) - KM,

                              max_LA = theor_mass1 / 13, actual_LA = ((C - E) + N + S + O + E + S34 + P + Cl +Cl37+N15) ,

                              rule_13=actual_LA/max_LA,

                              Senior1 = H ,#+ P + N + Cl + Cl37 + N15  ,

                              DBEO = DBE-O,

                              max_H = C * 2 + 2, H_test = H / max_H,

                              Senior2 = H*Valence("H") + O*Valence("O") + C*Valence("C"),

                              Senior3Atom = C + H + O,

                              Senior3Val = C*Valence("C") + H*Valence("H") + O*Valence("O")
    )

    records1 <- dplyr::filter(records1, C>0, H>0,O>=Omin, H >= D)
    records1 <- unique(records1)
    records1 <- dplyr::filter(records1, O_C < O_Cmax & H_C <= H_Cmax & H_C > H_Cmin & O_C >= O_Cmin &

                                DBEO >= DBEOmin & DBEO <= DBEOmax &

                                H_test <= 1 & rule_13 <= 1 &

                                AE_ppm <= ppm_err &

                                Even(Senior1)==TRUE & DBE >= 0 & DBE <= 0.9 * (C + N) &

                                O <= 2 * C + 3 * (N+N15) + 4 * P + 4 * (S+S34)&

                                #O>= P*4&

                                RA > 0 &

                                Senior2 >= 2*Valence("C") &

                                Senior3Val >= (2*Senior3Atom - 1)
    )



    ###Series analysis
    ##Determining Ambiguity
    records1$num <- 1:nrow(records1)
    records1$dups <- duplicated(records1$Exp_mass)
    records1 <-records1[order(-records1$num),]
    records1$dups2 <- duplicated(records1$Exp_mass)

    Unambig <- records1[records1$dups == FALSE & records1$dups2 == FALSE,]
    Unamatch <- Unambig[c(1,3, 4)]

    Ambig <- records1[records1$dups == TRUE | records1$dups2 == TRUE,]
    Ambig <- Ambig[c(1,3)]
    names(peaksAll2)[2] <- "Exp_mass"
    Ambig <- rbind(Ambig, peaksAll2)

    Ambig <- merge(Ambig, Unamatch, by.x = c("RA", "Exp_mass"), by.y = c("RA", "Exp_mass"), all = T)
    Ambig <- Ambig[is.na(Ambig$C),]
    Ambig <- Ambig[-c(3)]
    #Good to this point
    ##Unambiguous formulas prep

    Unambig <- Unambig[c("RA", "Exp_mass", "C", "H", "O", "N", "S", "P", "E", "S34", "N15", "D",
                         "Cl", "Fl", "Cl37", "M", "NH4", "POE", "NOE", "Z")]
    Unambig$NM <- round(Unambig$Exp_mass)

    Unambig$KM_CH2 <- Unambig$Exp_mass * (14/14.01565)
    Unambig$KMD_CH2 <- round(Unambig$NM - Unambig$KM_CH2, 3)
    Unambig$z_CH2 <- round(Unambig$Exp_mass)%%14 - 14

    Unambig$KM_O <- Unambig$Exp_mass * (16/15.9949146223)
    Unambig$KMD_O <- round(Unambig$NM - Unambig$KM_O, 3)
    Unambig$z_O <- round(Unambig$Exp_mass)%%16 - 16

    Unambig$KM_H2 <- Unambig$Exp_mass * (2/2.01565)
    Unambig$KMD_H2 <- round(Unambig$NM - Unambig$KM_H2, 3)
    Unambig$z_H2 <- round(Unambig$Exp_mass)%%2 - 2

    Unambig$KM_H2O <- Unambig$Exp_mass * (18/18.01056468)
    Unambig$KMD_H2O <- round(Unambig$NM - Unambig$KM_H2O, 3)
    Unambig$z_H2O <- round(Unambig$Exp_mass)%%18 - 18

    Unambig$KM_CH2O <- Unambig$Exp_mass * (30/30.01056468)
    Unambig$KMD_CH2O <- round(Unambig$NM - Unambig$KM_CH2O, 3)
    Unambig$z_CH2O <- round(Unambig$Exp_mass)%%30 - 30
    #Good to this point

    ##Ambiguous
    #Ambig <- Ambig[c(1,2)]
    Ambig$NM <- round(Ambig$Exp_mass)

    Ambig$KM_CH2 <- Ambig$Exp_mass * (14/14.01565)
    Ambig$KMD_CH2 <- round(Ambig$NM - Ambig$KM_CH2, 3)
    Ambig$z_CH2 <- round(Ambig$Exp_mass)%%14 - 14

    Ambig$KM_O <- Ambig$Exp_mass * (16/15.9949146223)
    Ambig$KMD_O <- round(Ambig$NM - Ambig$KM_O, 3)
    Ambig$z_O <- round(Ambig$Exp_mass)%%16 - 16

    Ambig$KM_H2 <- Ambig$Exp_mass * (2/2.01565)
    Ambig$KMD_H2 <- round(Ambig$NM - Ambig$KM_H2, 3)
    Ambig$z_H2 <- round(Ambig$Exp_mass)%%2 - 2

    Ambig$KM_H2O <- Ambig$Exp_mass * (18/18.01056468)
    Ambig$KMD_H2O <- round(Ambig$NM - Ambig$KM_H2O, 3)
    Ambig$z_H2O <- round(Ambig$Exp_mass)%%18 - 18

    Ambig$KM_CH2O <- Ambig$Exp_mass * (30/30.01056468)
    Ambig$KMD_CH2O <- round(Ambig$NM - Ambig$KM_CH2O, 3)
    Ambig$z_CH2O <- round(Ambig$Exp_mass)%%30 - 30
    #Good to this point
    ###Looping series assignment
    known <- Unambig
    unknown <- Ambig

    knowndummy <- data.frame(KMD_CH2 = -42, KMD_O = -42, KMD_CH2O = -42, KMD_H2O = -42, KMD_H2 = -42, z_CH2 = -42,
                             z_O = -42,
                             z_CH2O = -42, z_H2O = -42, z_H2 = -42, RA = -42, Exp_mass = -42, C = 4, H = 4, O = 0,
                             N = 0, S = 0,
                             S34 = 0, P = 0, N15 = 0, Cl = 0, Fl = 0, Cl37 = 0, M = 0, NH4 = 0, Z= 0, POE = 0,
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
                           S34 = 0, P = 0, N15 = 0, Cl = 0, Fl = 0, Cl37 = 0, M = 0, NH4 = 0, Z= 0, POE = 0, NOE = 0, D = 0, E = 0, KM_O = -42,
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
    #j = 1
    for(j in 0:(nloop1)){
      loop <- j
      masstrim <- DeNovo + 200*loop
      loop1 <- j-1
      if(loop1 < 0) {loop1 <- 0}
      masstrim1 <- DeNovo + 200*loop1
      Ambig2 <- Ambigsave[Ambigsave$Exp_mass < masstrim,] #Resets a comparison df each time

      Ambig4 <- merge(Ambigreturn, Ambig2, by.x = "Exp_mass", by.y = "Exp_mass", all = T)

      Common <- Ambig4[!is.na(Ambig4$RA.x)&!is.na(Ambig4$RA.y),]
      Common <- Common[c(1:18)]
      names(Common) <- gsub(".x","",names(Common),fixed = TRUE)
      New <- Ambig4[is.na(Ambig4$RA.x) & Ambig4$Exp_mass >= masstrim1,]
      New <- New[c(1,19:ncol(New))]
      names(New) <- gsub(".y","",names(New),fixed = TRUE)

      Ambig <- rbind(Common, New)
      #Ambig <- unique(Ambig)

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
                              "S34", "N15", "D", "Cl", "Fl", "Cl37", "M", "NH4", "POE", "NOE", "Z")]
          names(knownCH2)[2] <- "base_mass"
          Step1 <- merge(unknown, knownCH2, by.x = c("KMD_CH2", "z_CH2"), by.y = c("KMD_CH2", "z_CH2"))
          Step1$CH2_num <- round(((Step1$Exp_mass - Step1$base_mass))/14.01565)
          Step1$C <- Step1$C + Step1$CH2_num
          Step1$H <- Step1$H + 2 * Step1$CH2_num
          Step1$Type <- "CH2"
          Step1$form <- paste(Step1$C, Step1$H, Step1$O, Step1$N, Step1$S, Step1$P, Step1$E, Step1$S34,
                              Step1$N15, Step1$D, Step1$Cl, Step1$Fl, Step1$Cl37, Step1$M, Step1$NH4,
                              Step1$POE, Step1$NOE, sep = "_")
          Step1 <- Step1[-c(39)]


          knownO <- known[c("RA", "Exp_mass", "KMD_O", "z_O", "C", "H", "O", "N", "S", "P", "E",
                            "S34", "N15", "D", "Cl", "Fl", "Cl37", "M", "NH4", "POE", "NOE", "Z")]
          names(knownO)[2] <- "base_mass"
          Step2 <- merge(unknown, knownO, by.x = c("KMD_O", "z_O"), by.y = c("KMD_O", "z_O"))
          Step2$O_num <- round(((Step2$Exp_mass - Step2$base_mass))/15.9949146223)
          Step2$O <- Step2$O + Step2$O_num
          Step2$Type <- "O"
          Step2$form <- paste(Step2$C, Step2$H, Step2$O, Step2$N, Step2$S, Step2$P, Step2$E, Step2$S34,
                              Step2$N15, Step2$D, Step2$Cl, Step2$Fl, Step2$Cl37, Step2$M, Step2$NH4,
                              Step2$POE, Step2$NOE, sep = "_")
          Step2 <- Step2[-c(39)]

          knownH2 <- known[c("RA", "Exp_mass", "KMD_H2", "z_H2", "C", "H", "O", "N", "S", "P", "E",
                             "S34", "N15", "D", "Cl", "Fl", "Cl37", "M", "NH4", "POE", "NOE", "Z")]
          names(knownH2)[2] <- "base_mass"
          Step3 <- merge(unknown, knownH2, by.x = c("KMD_H2", "z_H2"), by.y = c("KMD_H2", "z_H2"))
          Step3$H2_num <- round(((Step3$Exp_mass - Step3$base_mass))/2.01565)
          Step3$H <- Step3$H + 2*Step3$H2_num
          Step3$Type <- "H2"
          Step3$form <- paste(Step3$C, Step3$H, Step3$O, Step3$N, Step3$S, Step3$P, Step3$E, Step3$S34,
                              Step3$N15, Step3$D, Step3$Cl, Step3$Fl, Step3$Cl37, Step3$M, Step3$NH4,
                              Step3$POE, Step3$NOE, sep = "_")
          Step3 <- Step3[-c(39)]

          knownH2O <- known[c("RA", "Exp_mass", "KMD_H2O", "z_H2O", "C", "H", "O", "N", "S", "P", "E",
                              "S34", "N15", "D", "Cl", "Fl", "Cl37", "M", "NH4", "POE", "NOE", "Z")]
          names(knownH2O)[2] <- "base_mass"
          Step4 <- merge(unknown, knownH2O, by.x = c("KMD_H2O", "z_H2O"), by.y = c("KMD_H2O", "z_H2O"))
          Step4$H2O_num <- round(((Step4$Exp_mass - Step4$base_mass))/18.01056468)
          Step4$H <- Step4$H + 2*Step4$H2O_num
          Step4$O <- Step4$O + Step4$H2O_num
          Step4$Type <- "H2O"
          Step4$form <- paste(Step4$C, Step4$H, Step4$O, Step4$N, Step4$S, Step4$P, Step4$E, Step4$S34,
                              Step4$N15, Step4$D, Step4$Cl, Step4$Fl, Step4$Cl37, Step4$M, Step4$NH4,
                              Step4$POE, Step4$NOE, sep = "_")
          Step4 <- Step4[-c(39)]

          knownCH2O <- known[c("RA", "Exp_mass", "KMD_CH2O", "z_CH2O", "C", "H", "O", "N", "S", "P", "E",
                               "S34", "N15", "D", "Cl", "Fl", "Cl37", "M", "NH4", "POE", "NOE", "Z")]
          names(knownCH2O)[2] <- "base_mass"
          Step5 <- merge(unknown, knownCH2O, by.x = c("KMD_CH2O", "z_CH2O"), by.y = c("KMD_CH2O", "z_CH2O"))
          Step5$CH2O_num <- round(((Step5$Exp_mass - Step5$base_mass))/30.01056468)
          Step5$H <- Step5$H + 2*Step5$CH2O_num
          Step5$O <- Step5$O + Step5$CH2O_num
          Step5$C <- Step5$C + Step5$CH2O_num
          Step5$Type <- "CH2O"
          Step5$form <- paste(Step5$C, Step5$H, Step5$O, Step5$N, Step5$S, Step5$P, Step5$E, Step5$S34,
                              Step5$N15, Step5$D, Step5$Cl, Step5$Fl, Step5$Cl37, Step5$M, Step5$NH4,
                              Step5$POE, Step5$NOE, sep = "_")
          Step5 <- Step5[-c(39)]

          Out <- rbind(Step1, Step2, Step3, Step4, Step5)
          Out <- Out[(Out$C >= 2 & Out$H >= 4 & Out$O >= 0),]
          Out$H_C <- Out$H/Out$C   #Quick internal QA to limit bad assignments
          Out$O_C <- Out$O/Out$C
          Out <- Out[Out$H_C < H_Cmax & Out$H_C > H_Cmin &
                       Out$O_C < O_Cmax & Out$O_C > O_Cmin,]
          Out <- Out[!names(Out) %in% c("H_C", "O_C")]


          Out <- rbind(Out, DummyOut)

          Out_form <- dplyr::group_by(Out, Exp_mass, form)

          Out_form$number <- 1
          Out_form <- dplyr::summarize_at(Out_form, "number", sum, na.rm = TRUE)

          #Turns on or off ambiguity based on user input
          if(Ambigcheck == "off") {
            Out_form <- dplyr::filter(Out_form, number == max(number))
          }

          Out_form <- unique(Out_form)
          Out2<- merge(Out, Out_form, by.x = c("Exp_mass", "form"), by.y = c("Exp_mass", "form"))
          Out3 <- dplyr::distinct(Out2, form, Exp_mass, .keep_all = TRUE)
          Out3 <- Out3[!names(Out3) %in% c("number")]


          Next <- Out3[!names(Out3) %in% c("base_mass", "Type", "form", "RA.y")]
          colnames(Next)[colnames(Next) == "RA.x"] <- "RA"
          Next <- rbind(known, Next)
          Next <- unique(Next)
          Unambig <- Next

          masses <- Unambig[c("Exp_mass", "RA")]
          names(masses)[2] <- "Var"
          Ambig <- merge(unknown, masses, by = "Exp_mass", all = T)
          Ambig <- Ambig[is.na(Ambig$Var),]
          Ambig <- Ambig[-19]
          Ambigreturn <- unique(Ambig)

          if(x == 2){
            break
          }
        }
        setTxtProgressBar(pb2, i)
        Unambig
      }
    }
    records1 <- Unambig[c(1:20)]

    records1$mode <- ionMode


    df1 <- records1[records1$mode == "pos" & records1$M > 0,]
    df1$Neutral_mass <- df1$Exp_mass - df1$M * 22.989221

    df2 <- records1[records1$mode == "pos" & records1$M == 0,]
    df2$Neutral_mass <- df2$Exp_mass - 1.00727645216

    df3 <- records1[records1$mode == "neg",]
    df3$Neutral_mass <- df3$Exp_mass + 1.00727645216

    records1 <- rbind(df1, df2, df3)
    records1 <- records1[-c(21)]

    ###Standard QA steps
    records1 <- dplyr::mutate(records1, O_C = O/(C+E), H_C =H/(C+E),

                              Neutral_mass = Neutral_mass + POE * (2.0156500638/2)- NOE * (2.0156500638/2),

                              theor_mass1 = EM("C") * C + EM("H") * H + EM("O") * O + N * EM("N14") +
                                S * EM("S") + P * EM("P31") +
                                Cl * EM("Cl35") + Fl * EM("Fl19") + E * EM("E") + S34 * EM("S34") +
                                Cl37 * EM("Cl37m") + N15 * EM("N15H") +
                                D * EM("D") + M * EM("M") + NH4 * EM("NH4") +POE * EM("POE") + NOE*EM("NOE"),

                              theor_mass = EM("C") * C + EM("H") * H + EM("O") * O + N * EM("N14") +
                                S * EM("S") + P * EM("P31") + Fl * EM("Fl19") +
                                Cl * EM("Cl35") +  E * EM("E") + S34 * EM("S34") + Cl37 * EM("Cl37m") +
                                N15 * EM("N15H") +
                                D * EM("D"),


                              C = C + E,

                              DBE = C - 0.5 * (H + Cl + Cl37 + Fl) + 0.5 * (N+ N15 + P) + 1,

                              err_ppm = ((Neutral_mass - theor_mass) / Neutral_mass * 10^6),

                              AE_ppm = abs((Neutral_mass - theor_mass) / Neutral_mass * 10^6),

                              NM = round(Exp_mass),

                              KM = Exp_mass * (14 / 14.01565), KMD = round(Exp_mass) - KM,

                              max_LA = theor_mass1 / 13, actual_LA = ((C - E) + N + S + O + E + S34 + P + Cl +Cl37+N15) ,

                              rule_13=actual_LA/max_LA,

                              Senior1 = H ,#+ P + N + Cl + Cl37 + N15  ,

                              DBEO = DBE-O,

                              max_H = C * 2 + 2, H_test = H / max_H,

                              Senior2 = H*Valence("H") + O*Valence("O") + C*Valence("C"),
                              Senior3Atom = C + H + O,
                              Senior3Val = C*Valence("C") + H*Valence("H") + O*Valence("O")
    )

    #recordssave <- records1

    records1 <- dplyr::filter(records1, C>0, H>0,O>=Omin, H >= D)
    records1 <- unique(records1)
    records1 <- dplyr::filter(records1, O_C < O_Cmax & H_C <= H_Cmax & H_C > H_Cmin & O_C >= O_Cmin &

                                DBEO >= DBEOmin & DBEO <= DBEOmax &

                                H_test <= 1 & rule_13 <= 1 &

                                AE_ppm <= ppm_err &

                                Even(Senior1)==TRUE & DBE >= 0 & DBE <= 0.9 * (C + N) &

                                O <= 2 * C + 3 * (N+N15) + 4 * P + 4 * (S+S34)&

                                #O>= P*4&

                                RA > 0 &

                                Senior2 >= 2*Valence("C") &

                                Senior3Val >= (2*Senior3Atom - 1)
    )



    records1 <- records1[!names(records1) %in% c("Senior1", "Senior2", "Senior3Val", "Senior3Atom")]

    ###Formula generation
    records1 <-
      dplyr::mutate(records1, Cform = ifelse(C == 0 , "",
                                             ifelse(C == 1 , "C", paste("C",C, sep = ""))),
                    Hform = ifelse(H == 0 , "",
                                   ifelse(H == 1 , "H", paste("H",H, sep = ""))),
                    Nform = ifelse(N == 0 , "",
                                   ifelse(N == 1 , "N", paste("N",N, sep = ""))),
                    Oform = ifelse(O == 0 , "",
                                   ifelse(O == 1 , "O", paste("O",O, sep = ""))),
                    Sform = ifelse(S == 0 , "",
                                   ifelse(S == 1 , "S", paste("S",S, sep = ""))),
                    Pform = ifelse(P == 0 , "",
                                   ifelse(P == 1 , "P", paste("P",O, sep = ""))),
                    Clform = ifelse(Cl == 0 , "",
                                    ifelse(Cl == 1 , "Cl", paste("Cl",Cl, sep = ""))))

    records1 <- tidyr::unite(records1, class, Nform, Oform, Sform, Pform, Clform, sep = "", remove = FALSE)

    records1 <- tidyr::unite(records1, formula, Cform, Hform, Nform, Oform, Sform, Pform, Clform, sep = "")

    records1 <-
      dplyr::mutate(records1, Cform = ifelse(C == 0 , "", "C"),
                    Hform = ifelse(H == 0 , "", "H"),

                    Nform = ifelse(N == 0 , "", "N"),

                    Oform = ifelse(O == 0 , "", "O"),

                    Sform = ifelse(S == 0 , "","S"),
                    Pform = ifelse(P == 0 , "", "P"),
                    Clform = ifelse(Cl == 0 , "", "Cl"))

    records1 <- tidyr::unite(records1, group, Cform, Hform, Nform, Oform, Sform, Pform, Clform, sep = "")

    ###Supplemental Specialized QA Steps
    records1<-dplyr::mutate(records1, HA =  E)

    records1<-dplyr::group_by(records1, Exp_mass)

    ifelse(HetCut == "on", records1<-dplyr::filter(records1, HA == (min(HA))), records1<- records1)

    records1 <- dplyr::group_by(records1, Exp_mass)

    records1 <- dplyr::distinct(records1, formula, .keep_all = TRUE)

    records1 <- dplyr::ungroup(records1)

    records1 <- records1[!names(records1) %in% c("HA")]





    cut <- dplyr::bind_rows(peaksAll2)
    cut <- unique(cut)
    unassigned <- dplyr::left_join(cut, records1, by = "Exp_mass")
    unassigned <- unassigned[is.na(unassigned$formula),]
    unassigned <- unassigned[c("RA.x", "Exp_mass")]
    names(unassigned)[1] <- "RA"
    names(unassigned)[2] <- "mass"
    unassigned <- unique(unassigned)  #Good to this point

    records1$mode <- ionMode


    df1 <- records1[records1$mode == "pos",]
    df1$theor_mass1 <- df1$theor_mass1 + proton

    df2 <- records1[records1$mode == "neg",]
    df2$theor_mass1 <- df2$theor_mass1 - proton

    records1 <- rbind(df1, df2)
    records1 <- records1[-c(41)]  #removes mode


    records1 <- records1[c(1,2,38:40,3:21,24, 27:31, 26, 22:23,35, 32:34,36,37)]

    #######################

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

    #Ccheck <- records1[records1$C13_mass > 0,]
    records1 <- dplyr::left_join(records1, C13Iso, by = "C13_mass")
    records1 <- records1[-42]
    #########
    #Align double C13 masses
    records1$C13_mass2 <- records1$Exp_mass + 2.006709676
    err <- ppm_err*10^-6
    C13Iso2 <- isopeaks2[isopeaks2$Tag == "2C13"|isopeaks2$Tag == "2C13_S34",]
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
    records1 <- records1[-44]
    ###########

    records1$S34_mass <- 0
    records1$S34_Abund <- 0
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
    Iso_save <- records1[c("Exp_mass", "C13_mass", "C13_Abund", "C13_mass2", "C13_Abund2", "S34_mass", "S34_Abund")]
    ###########################################################################################################
    #Additional formula extension for the isotope masses not matched to monoisotopic mass

    #Preparation of unmatched "isotope" mass list
    C13 <- records1[c("C13_mass")]
    C13_2 <- records1[c("C13_mass2")]
    S34 <- records1[c("S34_mass")]
    C13 <- unique(C13)
    C13_2 <- unique(C13_2)
    S34 <- unique(S34)
    names(C13)[1] <- "Iso_mass"
    names(C13_2)[1] <- "Iso_mass"
    names(S34)[1] <- "Iso_mass"
    Iso_match <- rbind(C13, C13_2, S34)
    Iso_match <- unique(Iso_match)
    Iso_match$Tag2 <- "Match"
    Iso_align <- merge(isopeaks2, Iso_match, by.x = "Iso_mass", by.y = "Iso_mass", all = T)
    Iso_nomatch <- Iso_align[is.na(Iso_align$Tag2),]
    Iso_nomatch <- Iso_nomatch[c(1,2)]
    names(Iso_nomatch)[1] <- "Exp_mass"
    names(Iso_nomatch)[2] <- "RA"

    #Kendrick Series Preparation
    Iso_nomatch <- Iso_nomatch[c(1,2)]
    Iso_nomatch$NM <- round(Iso_nomatch$Exp_mass)

    Iso_nomatch$KM_CH2 <- Iso_nomatch$Exp_mass * (14/14.01565)
    Iso_nomatch$KMD_CH2 <- round(Iso_nomatch$NM - Iso_nomatch$KM_CH2, 3)
    Iso_nomatch$z_CH2 <- round(Iso_nomatch$Exp_mass)%%14 - 14

    Iso_nomatch$KM_O <- Iso_nomatch$Exp_mass * (16/15.9949146223)
    Iso_nomatch$KMD_O <- round(Iso_nomatch$NM - Iso_nomatch$KM_O, 3)
    Iso_nomatch$z_O <- round(Iso_nomatch$Exp_mass)%%16 - 16

    Iso_nomatch$KM_H2 <- Iso_nomatch$Exp_mass * (2/2.01565)
    Iso_nomatch$KMD_H2 <- round(Iso_nomatch$NM - Iso_nomatch$KM_H2, 3)
    Iso_nomatch$z_H2 <- round(Iso_nomatch$Exp_mass)%%2 - 2

    Iso_nomatch$KM_H2O <- Iso_nomatch$Exp_mass * (18/18.01056468)
    Iso_nomatch$KMD_H2O <- round(Iso_nomatch$NM - Iso_nomatch$KM_H2O, 3)
    Iso_nomatch$z_H2O <- round(Iso_nomatch$Exp_mass)%%18 - 18

    Iso_nomatch$KM_CH2O <- Iso_nomatch$Exp_mass * (30/30.01056468)
    Iso_nomatch$KMD_CH2O <- round(Iso_nomatch$NM - Iso_nomatch$KM_CH2O, 3)
    Iso_nomatch$z_CH2O <- round(Iso_nomatch$Exp_mass)%%30 - 30


    recordsx <- records1[c("RA", "Exp_mass", "C", "H", "O", "N", "S", "P", "E", "S34", "N15", "D",
                           "Cl", "Fl", "Cl37", "M", "NH4", "POE", "NOE", "Z")]
    recordsx$NM <- round(recordsx$Exp_mass)

    recordsx$KM_CH2 <- recordsx$Exp_mass * (14/14.01565)
    recordsx$KMD_CH2 <- round(recordsx$NM - recordsx$KM_CH2, 3)
    recordsx$z_CH2 <- round(recordsx$Exp_mass)%%14 - 14

    recordsx$KM_O <- recordsx$Exp_mass * (16/15.9949146223)
    recordsx$KMD_O <- round(recordsx$NM - recordsx$KM_O, 3)
    recordsx$z_O <- round(recordsx$Exp_mass)%%16 - 16

    recordsx$KM_H2 <- recordsx$Exp_mass * (2/2.01565)
    recordsx$KMD_H2 <- round(recordsx$NM - recordsx$KM_H2, 3)
    recordsx$z_H2 <- round(recordsx$Exp_mass)%%2 - 2

    recordsx$KM_H2O <- recordsx$Exp_mass * (18/18.01056468)
    recordsx$KMD_H2O <- round(recordsx$NM - recordsx$KM_H2O, 3)
    recordsx$z_H2O <- round(recordsx$Exp_mass)%%18 - 18

    recordsx$KM_CH2O <- recordsx$Exp_mass * (30/30.01056468)
    recordsx$KMD_CH2O <- round(recordsx$NM - recordsx$KM_CH2O, 3)
    recordsx$z_CH2O <- round(recordsx$Exp_mass)%%30 - 30

    ##Formula Extension for isotope masses

    pb2 <- txtProgressBar(min = 0, max = nLoop, style = 3)

    ###############
    Allmasses <- peaksAll[!is.na(peaksAll$mass),]
    nloop1 <- ceiling((max(Allmasses$mass + 1)-DeNovo)/200)
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

      Ambig4 <- merge(Ambigreturn, Ambig2, by.x = "Exp_mass", by.y = "Exp_mass", all = T)

      Common <- Ambig4[!is.na(Ambig4$RA.x)&!is.na(Ambig4$RA.y),]
      Common <- Common[c(1:18)]
      names(Common) <- gsub(".x","",names(Common),fixed = TRUE)
      New <- Ambig4[is.na(Ambig4$RA.x) & Ambig4$Exp_mass >= masstrim1,]
      New <- New[c(1,19:ncol(New))]
      names(New) <- gsub(".y","",names(New),fixed = TRUE)

      Ambig <- rbind(Common, New)
      #Ambig <- unique(Ambig)

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
                              "S34", "N15", "D", "Cl", "Fl", "Cl37", "M", "NH4", "POE", "NOE", "Z")]
          names(knownCH2)[2] <- "base_mass"
          Step1 <- merge(unknown, knownCH2, by.x = c("KMD_CH2", "z_CH2"), by.y = c("KMD_CH2", "z_CH2"))
          Step1$CH2_num <- round(((Step1$Exp_mass - Step1$base_mass))/14.01565)
          Step1$C <- Step1$C + Step1$CH2_num
          Step1$H <- Step1$H + 2 * Step1$CH2_num
          Step1$Type <- "CH2"
          Step1$form <- paste(Step1$C, Step1$H, Step1$O, Step1$N, Step1$S, Step1$P, Step1$E, Step1$S34,
                              Step1$N15, Step1$D, Step1$Cl, Step1$Fl, Step1$Cl37, Step1$M, Step1$NH4,
                              Step1$POE, Step1$NOE, sep = "_")
          Step1 <- Step1[-c(39)]


          knownO <- known[c("RA", "Exp_mass", "KMD_O", "z_O", "C", "H", "O", "N", "S", "P", "E",
                            "S34", "N15", "D", "Cl", "Fl", "Cl37", "M", "NH4", "POE", "NOE", "Z")]
          names(knownO)[2] <- "base_mass"
          Step2 <- merge(unknown, knownO, by.x = c("KMD_O", "z_O"), by.y = c("KMD_O", "z_O"))
          Step2$O_num <- round(((Step2$Exp_mass - Step2$base_mass))/15.9949146223)
          Step2$O <- Step2$O + Step2$O_num
          Step2$Type <- "O"
          Step2$form <- paste(Step2$C, Step2$H, Step2$O, Step2$N, Step2$S, Step2$P, Step2$E, Step2$S34,
                              Step2$N15, Step2$D, Step2$Cl, Step2$Fl, Step2$Cl37, Step2$M, Step2$NH4,
                              Step2$POE, Step2$NOE, sep = "_")
          Step2 <- Step2[-c(39)]

          knownH2 <- known[c("RA", "Exp_mass", "KMD_H2", "z_H2", "C", "H", "O", "N", "S", "P", "E",
                             "S34", "N15", "D", "Cl", "Fl", "Cl37", "M", "NH4", "POE", "NOE", "Z")]
          names(knownH2)[2] <- "base_mass"
          Step3 <- merge(unknown, knownH2, by.x = c("KMD_H2", "z_H2"), by.y = c("KMD_H2", "z_H2"))
          Step3$H2_num <- round(((Step3$Exp_mass - Step3$base_mass))/2.01565)
          Step3$H <- Step3$H + 2*Step3$H2_num
          Step3$Type <- "H2"
          Step3$form <- paste(Step3$C, Step3$H, Step3$O, Step3$N, Step3$S, Step3$P, Step3$E, Step3$S34,
                              Step3$N15, Step3$D, Step3$Cl, Step3$Fl, Step3$Cl37, Step3$M, Step3$NH4,
                              Step3$POE, Step3$NOE, sep = "_")
          Step3 <- Step3[-c(39)]

          knownH2O <- known[c("RA", "Exp_mass", "KMD_H2O", "z_H2O", "C", "H", "O", "N", "S", "P", "E",
                              "S34", "N15", "D", "Cl", "Fl", "Cl37", "M", "NH4", "POE", "NOE", "Z")]
          names(knownH2O)[2] <- "base_mass"
          Step4 <- merge(unknown, knownH2O, by.x = c("KMD_H2O", "z_H2O"), by.y = c("KMD_H2O", "z_H2O"))
          Step4$H2O_num <- round(((Step4$Exp_mass - Step4$base_mass))/18.01056468)
          Step4$H <- Step4$H + 2*Step4$H2O_num
          Step4$O <- Step4$O + Step4$H2O_num
          Step4$Type <- "H2O"
          Step4$form <- paste(Step4$C, Step4$H, Step4$O, Step4$N, Step4$S, Step4$P, Step4$E, Step4$S34,
                              Step4$N15, Step4$D, Step4$Cl, Step4$Fl, Step4$Cl37, Step4$M, Step4$NH4,
                              Step4$POE, Step4$NOE, sep = "_")
          Step4 <- Step4[-c(39)]

          knownCH2O <- known[c("RA", "Exp_mass", "KMD_CH2O", "z_CH2O", "C", "H", "O", "N", "S", "P", "E",
                               "S34", "N15", "D", "Cl", "Fl", "Cl37", "M", "NH4", "POE", "NOE", "Z")]
          names(knownCH2O)[2] <- "base_mass"
          Step5 <- merge(unknown, knownCH2O, by.x = c("KMD_CH2O", "z_CH2O"), by.y = c("KMD_CH2O", "z_CH2O"))
          Step5$CH2O_num <- round(((Step5$Exp_mass - Step5$base_mass))/30.01056468)
          Step5$H <- Step5$H + 2*Step5$CH2O_num
          Step5$O <- Step5$O + Step5$CH2O_num
          Step5$C <- Step5$C + Step5$CH2O_num
          Step5$Type <- "CH2O"
          Step5$form <- paste(Step5$C, Step5$H, Step5$O, Step5$N, Step5$S, Step5$P, Step5$E, Step5$S34,
                              Step5$N15, Step5$D, Step5$Cl, Step5$Fl, Step5$Cl37, Step5$M, Step5$NH4,
                              Step5$POE, Step5$NOE, sep = "_")
          Step5 <- Step5[-c(39)]

          Out <- rbind(Step1, Step2, Step3, Step4, Step5)
          Out <- Out[(Out$C >= 2 & Out$H >= 4 & Out$O >= 0),]
          Out$H_C <- Out$H/Out$C   #Quick internal QA to limit bad assignments
          Out$O_C <- Out$O/Out$C
          Out <- Out[Out$H_C < H_Cmax & Out$H_C > H_Cmin &
                       Out$O_C < O_Cmax & Out$O_C > O_Cmin,]
          Out <- Out[!names(Out) %in% c("H_C", "O_C")]


          Out <- rbind(Out, DummyOut)

          Out_form <- dplyr::group_by(Out, Exp_mass, form)

          Out_form$number <- 1
          Out_form <- dplyr::summarize_at(Out_form, "number", sum, na.rm = TRUE)

          #Turns on or off ambiguity based on user input
          if(Ambigcheck == "off") {
            Out_form <- dplyr::filter(Out_form, number == max(number))
          }

          Out_form <- unique(Out_form)
          Out2<- merge(Out, Out_form, by.x = c("Exp_mass", "form"), by.y = c("Exp_mass", "form"))
          Out3 <- dplyr::distinct(Out2, form, Exp_mass, .keep_all = TRUE)
          Out3 <- Out3[!names(Out3) %in% c("number")]


          Next <- Out3[!names(Out3) %in% c("base_mass", "Type", "form", "RA.y")]
          colnames(Next)[colnames(Next) == "RA.x"] <- "RA"
          Next <- rbind(known, Next)
          Next <- unique(Next)
          Unambig <- Next

          masses <- Unambig[c("Exp_mass", "RA")]
          names(masses)[2] <- "Var"
          Ambig <- merge(unknown, masses, by = "Exp_mass", all = T)
          Ambig <- Ambig[is.na(Ambig$Var),]
          Ambig <- Ambig[-19]
          Ambigreturn <- unique(Ambig)

          if(x == 2){
            break
          }
        }
        setTxtProgressBar(pb2, i)
        Unambig
      }
    }
    records1X <- Unambig[c(1:20)]

    records1X <- unique(records1X)


    records1X$mode <- ionMode


    df1 <- records1X[records1X$mode == "pos" & records1X$M > 0,]
    df1$Neutral_mass <- df1$Exp_mass - df1$M * 22.989221

    df2 <- records1X[records1X$mode == "pos" & records1X$M == 0,]
    df2$Neutral_mass <- df2$Exp_mass - 1.00727645216

    df3 <- records1X[records1X$mode == "neg",]
    df3$Neutral_mass <- df3$Exp_mass + 1.00727645216

    records1X <- rbind(df1, df2, df3)
    records1X <- records1X[-c(21)]

    #NOEx seems to be good to this point.
    ###Standard QA steps, second round
    #if(records1X$E > 0) {records1X$C <- records1X$C + records1X$E} # This is to fix the theor. masses
    records1X <- dplyr::mutate(records1X, O_C = O/(C+E), H_C =H/(C+E),

                               Neutral_mass = Neutral_mass + POE * (2.0156500638/2)- NOE * (2.0156500638/2),


                               theor_mass1 = EM("C") * C + EM("H") * H + EM("O") * O + N * EM("N14") +
                                 S * EM("S") + P * EM("P31") +
                                 Cl * EM("Cl35") + Fl * EM("Fl19") + E * EM("E2") + S34 * EM("S34") +
                                 Cl37 * EM("Cl37m") + N15 * EM("N15H") +
                                 D * EM("D") + M * EM("M") + NH4 * EM("NH4") +POE * EM("POE") + NOE*EM("NOE"),

                               theor_mass = EM("C") * C + EM("H") * H + EM("O") * O + N * EM("N14") +
                                 S * EM("S") + P * EM("P31") + Fl * EM("Fl19") +
                                 Cl * EM("Cl35") +  E * EM("E2") + S34 * EM("S34") + Cl37 * EM("Cl37m") +
                                 N15 * EM("N15H") +
                                 D * EM("D"),

                               #C = C + E, #It is added back so that formulas are more accurate.

                               DBE = C - 0.5 * (H + Cl + Cl37 +Fl) + 0.5 * (N +N15+ P) + 1,

                               err_ppm = ((Neutral_mass - theor_mass) / Neutral_mass * 10^6),

                               AE_ppm = abs((Neutral_mass - theor_mass) / Neutral_mass * 10^6),

                               NM = round(Exp_mass),

                               KM = Exp_mass * (14 / 14.01565), KMD = round(Exp_mass) - KM,

                               max_LA = theor_mass1 / 13, actual_LA = ((C - E) + N + S + Fl + O + E + S34 + P + Cl +Cl37+N15) ,

                               rule_13=actual_LA/max_LA,

                               Senior1 = H + P + N + Cl + Fl + Cl37 + N15  ,

                               STest = S + S34, ClTest = Cl + Cl37, NTest = N + N15, DBEO = DBE-O,

                               max_H = C * 2 + 2, H_test = H / max_H,

                               Senior2 = Valence("P") + Valence("N") + Valence("N15") + Valence("H")  +
                                 Valence("Cl") + Valence("Cl37") + Valence("Fl"),
                               Senior3Atom = C + H + O + N + S + P + N15 + E + Cl + Fl + Cl37 + S34,
                               Senior3Val = C*Valence("C") + H*Valence("H") + O*Valence("O") + N*Valence("N")
                               + S*Valence("S") + P*Valence("P") + S34*Valence("S34") + Fl*Valence("Fl") +
                                 N15*Valence("N15") + Cl*Valence("Cl") + Cl37*Valence("Cl37")
    )



    #recordssave <- records1X
    #records1X <- recordssave

    records1X <- dplyr::filter(records1X, C>0, H>0,O>=Omin, H >= D)
    records1X <- unique(records1X)
    records1X <- dplyr::filter(records1X, O_C < O_Cmax & H_C <= H_Cmax & H_C > H_Cmin & O_C >= O_Cmin &

                                 DBEO >= DBEOmin & DBEO <= DBEOmax &



                                 H_test <= 1 & rule_13 <= 1 &

                                 AE_ppm <= ppm_err &

                                 Even(Senior1)==TRUE & DBE >= 0 & DBE <= 0.9 * (C + N) &

                                 O <= 2 * C + 3 * (N+N15) + 4 * P + 4 * (S+S34)&

                                 #O>= P*4&

                                 RA > 0 &

                                 Senior2 >= 2*Valence("C") &

                                 Senior3Val >= (2*Senior3Atom - 1)
    )



    records1X <- records1X[!names(records1X) %in% c("Senior1", "Senior2", "Senior3Val", "Senior3Atom")]

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
                                    ifelse(Fl == 1 , "F", paste("F",Fl, sep = ""))))

    records1X <- tidyr::unite(records1X, class, Nform, Oform, Sform, Pform, Clform, Flform,
                              sep = "", remove = FALSE)

    records1X <- tidyr::unite(records1X, formula, Cform, Hform, Nform, Oform, Sform, Pform, Clform,
                              Flform, sep = "")

    records1X <-
      dplyr::mutate(records1X, Cform = ifelse(C == 0 , "", "C"),
                    Hform = ifelse(H == 0 , "", "H"),

                    Nform = ifelse(NTest == 0 , "", "N"),

                    Oform = ifelse(O == 0 , "", "O"),

                    Sform = ifelse(STest == 0 , "","S"),
                    Pform = ifelse(P == 0 , "", "P"),
                    Clform = ifelse(ClTest == 0 , "", "Cl"),
                    Flform = ifelse(Fl == 0 , "", "F"))

    records1X <- tidyr::unite(records1X, group, Cform, Hform, Nform, Oform, Sform, Pform,
                              Clform, Flform, sep = "")

    ###Supplemental Specialized QA Steps
    records1X<-dplyr::mutate(records1X, HA = NTest + STest + P + ClTest + E + Fl)

    records1X<-dplyr::group_by(records1X, Exp_mass)

    ifelse(HetCut == "on", records1X<-dplyr::filter(records1X, HA == (min(HA))), records1X<- records1X)

    records1X <- dplyr::group_by(records1X, Exp_mass)

    records1X <- dplyr::distinct(records1X, formula, .keep_all = TRUE)

    records1X <- dplyr::ungroup(records1X)

    records1X <- records1X[!names(records1X) %in% c("HA")]



    #records3 <- dplyr::rename(records1X, mass = Exp_mass)

    cut <- dplyr::bind_rows(peaksAll2)
    cut <- unique(cut)
    unassigned <- dplyr::left_join(cut, records1X, by = "Exp_mass")
    unassigned <- unassigned[is.na(unassigned$formula),]
    unassigned <- unassigned[c("RA.x", "Exp_mass")]
    names(unassigned)[1] <- "RA"
    names(unassigned)[2] <- "mass"
    unassigned <- unique(unassigned)  #Good to this point

    records1X$mode <- ionMode


    df1 <- records1X[records1X$mode == "pos",]
    df1$theor_mass1 <- df1$theor_mass1 + proton

    df2 <- records1X[records1X$mode == "neg",]
    df2$theor_mass1 <- df2$theor_mass1 - proton

    records1X <- rbind(df1, df2)
    records1X <- records1X[-c(44)]  #removes mode

    records1X <- records1X[c(1,2,41:43,3:21,24,27:31, 26, 22:23,38, 32:34,39, 40)]

    records1 <- records1X
    records1 <- merge(records1, Iso_save, by.x = "Exp_mass", by.y = "Exp_mass", all = T)
    records1 <- records1[c(2,1, 3:ncol(records1))]
    ###THIS IS THE END OF THE SECONDARY ASSIGNMENT FOR THE UNMATCHED ISOTOPE MASSES

    ##############################################################

    records1$Dups <- duplicated(records1$Exp_mass)
    records2 <- records1[records1$Dups == FALSE,]
    Ambig <- records1[records1$Dups == TRUE,]
    Ambigcheck <- Ambig[c(2)]
    Ambigdummy <- data.frame(Exp_mass = 1)
    Ambigcheck <- rbind(Ambigcheck, Ambigdummy) #This one is fine
    Ambigcheck$Tag <- "Ambig"

    Ambigfinal <- dplyr::left_join(records2, Ambigcheck, by = "Exp_mass")
    Unambig <- Ambigfinal[is.na(Ambigfinal$Tag),]
    Unambig <- Unambig[-c(46,47)]
    Ambigfinal <- Ambigfinal[!is.na(Ambigfinal$Tag),]
    Ambigfinal <- Ambigfinal[-47]
    Ambigout <- rbind(Ambig, Ambigfinal) #This one is fine
    Ambigout <- Ambigout[-46]
    Ambigout2 <- data.frame(Exp_mass = 1)
    Ambigout <- dplyr::bind_rows(Ambigout, Ambigout2)
    Ambigout <- unique(Ambigout)
    Ambigout$Tag <- "Ambiguous"
    Unambig$Tag <- "Unambiguous"
    #Everything is good to this point

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

      Ambigout <- Ambig[-c(3,48:53)]
      Ambigout <- Ambigout[c(3:29,1,30:34, 2, 35:46)]
      Unambigout <- Unambig2[-c(3,48:53)]
      Unambigout <- Unambigout[c(3:29,1,30:34, 2, 35:46)]
      Unambig <- rbind(Unambig, Unambigout)
    }

    Unambig <- Unambig[!is.na(Unambig$N15),]
    Unambig$theor_mass1 <- Unambig$theor_mass1 - Unambig$POE * 2.0156500638+
      Unambig$NOE * 2.0156500638 + Unambig$NOE * electron - Unambig$POE * electron
    #Everything is good to this point

    ###########
    ##Plot data preparation
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


    #.rs.restartR()

    ##Final Output list
    output <- list(Unambig = Unambig, Ambig = Ambigout, None = unassigned, MSAssign = MZ,
                   Error = Error, MSgroups = MZgroups, VK = VK)

    output


}



