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
#'
#'
#'
#Ox <- 30
#maxAttempts  <- 10
FindCoreFormulae2_Halo <- function(env) {
  # constants
  #coreXEM <- 232.052663
  #coreRNM <- 232
  #RECORDS <- vector("list", 1000000)
  moietyMass <- 0.0363855087200022;
  maxAttempts <- env$Ox/3;

  env$formulaOK <- FALSE

  #env$coreRNM <- 102
  #env$coreXEM <- 102.0322
  #maxAttempts  <- 10
  if (Even(env$coreRNM) == TRUE) {
  C0 <- floor(env$coreRNM / 12)
  H0 <- env$coreRNM - C0 * 12
  O0 <- 0
  Output <- data.frame(C = 0, H = 0, O = 1, XEM = 0, CEM = 0, Ratio = 0, Ratio2 = 2.3)
  for(try in 0:(maxAttempts-1)) {

   # try <- 1
    C1 <- C0 - 4 * try
    H1 <- H0
    O1 <- O0 + 3 * try
    coreCEM <- C1* 12 + H1 * 1.0078250319 + O1 * 15.9949146223
    Ratio <- round((env$coreXEM - coreCEM)/moietyMass, 2)
    Ratio2 <- round(Ratio)
    Out <- dplyr::bind_cols(C = C1, H = H1, O = O1, XEM = env$coreXEM, CEM = coreCEM, Ratio = Ratio, Ratio2 = Ratio2)
    ifelse(abs(Out$Ratio2 - Out$Ratio) <= 0.05, Output <- Out, Output <- Output)
    if(abs(Out$Ratio2 - Out$Ratio) <= 0.05){
      #Output <- Out
      break
    }

  }



  #Output <- Output[abs(Output$Ratio2 - Output$Ratio) <= 0.05,]
  Output$C <- Output$C + Output$Ratio2
  Output$H <- Output$H + 4 * Output$Ratio2
  Output$O <- Output$O + (-1) * Output$Ratio2
  Output$CEM <- Output$C* 12 + Output$H * 1.0078250319 + Output$O * 15.9949146223
  Output$xemErr <- ((Output$CEM - Output$XEM)/Output$CEM * 10^6)
  Output <- Output[abs(Output$xemErr) <= env$ppm_err & Output$C > 0 & Output$O >= 0,]

  env$records = list(RA = (env$RA), coreNM = env$coreRNM, Exp_mass = env$ionEM,
                     C = Output$C,
                     H = Output$H, O = Output$O,
                     N= env$loop[CompFactorToInt2("N")], S=env$loop[CompFactorToInt2("S")],
                     P = env$loop[CompFactorToInt2("P")], E = env$loop[CompFactorToInt2("E")],
                     S34 = env$loop[CompFactorToInt2("S34")], N15 = env$loop[CompFactorToInt2("N15")],
                     D = env$loop[CompFactorToInt2("D")], Cl = env$loop[CompFactorToInt2("Cl")],
                     Fl = env$loop[CompFactorToInt2("Fl")],
                     Cl37 = env$loop[CompFactorToInt2("Cl37")], Br = env$loop[CompFactorToInt2("Br")],
                     Br81 = env$loop[CompFactorToInt2("Br81")], I = env$loop[CompFactorToInt2("I")],
                     M = env$loop[CompFactorToInt2("M")],
                     NH4 = env$loop[CompFactorToInt2("NH4")], POE = env$loop[CompFactorToInt2("POE")],
                     NOE = env$loop[CompFactorToInt2("NOE")],
                     Z = env$loop[CompFactorToInt2("Z")],
                     Neutral_mass = Output$XEM, CHO_mass = Output$CEM,
                     CHO_Err = Output$xemErr, Ratio = Output$Ratio2)

#
#}
}
}

