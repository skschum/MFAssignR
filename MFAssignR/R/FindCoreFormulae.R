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
  moietyMass <- 0.0363855087200022;
  maxAttempts <- env$Ox/3;

  env$formulaOK <- FALSE

  if (Even(env$coreRNM)) {
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
    switch(env$fitMode,
           "ppm" = { env$xemErr <-
             round(1e6*(env$coreCEM-env$coreXEM)/env$exactEM, env$numDigits) },
           "mDa" = { env$xemErr <-
             round(1e3*(env$coreCEM-env$coreXEM), env$numDigits) })

    env$Ratio <- round((env$coreXEM - env$coreCEM)/moietyMass *10)/10
    if(env$Ratio - round(env$Ratio)==0){

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
                                  Neutral_mass = env$exactEM, CHO_mass = env$coreCEM, CHO_Err = env$xemErr, Ratio = env$Ratio)}}


    env$formulaOK <- (ValidFormula(env$loop) & (abs(env$xemErr) <= env$maxErr))

  }

}


