#' Assigns all possible combinations of C, H, and O for a given MW
#'
#' FindCoreFormulae assigns all possible combinations of C, H, and O for a
#' given MW and subject to constraints in the \code{\link{ValidFormula}} function.
#' This is a subfunction of the MFAssign function and will not work
#' independently of that function.
#'
#' @param env environment from \code{\link{MFAssign}} function
#'
#' @examples
#' env <- environment()
#' FindCoreFormulae(env)
#'
#' @export
FindCoreFormulae <- function(env) {
  # constants
  #coreXEM <- 162.127970
  #coreRNM <- 162
  #RECORDS <- vector("list", 1000000)
  moietyMass <- 0.0363855087200022;
  maxAttempts <- env$Ox/3;

  env$formulaOK <- FALSE

  if (Even(env$coreRNM)) {
    #try <- 4
    for(try in 1:maxAttempts) {
    # Find the hydrocarbon having this NM and maximum number of moles of C.

    C0 <- floor(round(env$coreRNM) / NM("C"))
    env$loop[CompFactorToInt("H")] <- env$coreRNM - C0*NM("C")
    env$loop[CompFactorToInt("C")] <- floor(round(env$coreRNM) / NM("C")) - 4* (try-1)
    env$loop[CompFactorToInt("O")] <- 0 + 3*(try-1)
    env$coreCEM <- 0

    for (step in CompFactorToInt("C"):CompFactorToInt("O")) {
      env$coreCEM <- env$coreCEM + unlist(env$loop[step])*EM(CompIntToFactor(step))
    }
    #switch(env$fitMode,
           #"ppm" = { env$xemErr <-
            # round(1e6*(env$coreCEM-env$coreXEM)/env$coreCEM, env$numDigits) }#,
           #"mDa" = { env$xemErr <-
             #round(1e3*(env$coreCEM-env$coreXEM), env$numDigits) }
           #)

    env$Ratio <- round((env$coreXEM - env$coreCEM)/moietyMass *10)/10
    if(abs(env$Ratio - round(env$Ratio))<=0.10001){  #Changed 5/21/19 This has an impact on what is assigned
      env$Ratio <- round(env$Ratio, 0)

      env$xemErr <- round(1e6*(env$coreCEM-(env$coreXEM - env$Ratio*moietyMass))/env$coreCEM, env$numDigits)

      env$records = list(RA = (env$RA), coreNM = env$coreRNM, Exp_mass = env$ionEM,
                                  C = env$loop[CompFactorToInt("C")],
                                  H = env$loop[CompFactorToInt("H")], O = env$loop[CompFactorToInt("O")],
                                  N= env$loop[CompFactorToInt("N")], S=env$loop[CompFactorToInt("S")],
                                  P = env$loop[CompFactorToInt("P")], E = env$loop[CompFactorToInt("E")],
                                  S34 = env$loop[CompFactorToInt("S34")], N15 = env$loop[CompFactorToInt("N15")],
                                  D = env$loop[CompFactorToInt("D")], Cl = env$loop[CompFactorToInt("Cl")],
                                  Fl = env$loop[CompFactorToInt("Fl")],
                                  Cl37 = env$loop[CompFactorToInt("Cl37")], M = env$loop[CompFactorToInt("M")],
                                  NH4 = env$loop[CompFactorToInt("NH4")], POE = env$loop[CompFactorToInt("POE")],
                                  NOE = env$loop[CompFactorToInt("NOE")],
                                  Z = env$loop[CompFactorToInt("Z")],
                                  Neutral_mass = env$exactEM, CHO_mass = env$coreCEM,
                                  CHO_Err = env$xemErr, Ratio = env$Ratio)
      break  #This causes the function to pick the first matching hit.
    }

    #env$RECORDS[[step]] <- env$records
    }


    #env$formulaOK <- ( (abs(env$xemErr) <= env$maxErr*10))
    #ValidFormula(env$loop) &
  }

}


