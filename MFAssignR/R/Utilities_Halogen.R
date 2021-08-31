# -----------------------------------------------------------------------------
# Functions for replacing components with given values:
#   NM - nominal mass
#   EM - exact mass
#   MinMoles
#   MaxMoles
#   Valence
# -----------------------------------------------------------------------------


#' Component factor to integer
#'
#' Given a component as a factor returns the correponding integer representation.
#' For example, "C" -> 1, "H" -> 2, ..., "Z" -> 9
#'
#' This is an internal fuction that will not work outside the \code{\link{MFAssign}} function
#'  environment.
#'
#' @param x component; factor (\code{"C", "H", "O", "N", "S", "P", "E", "S34", "N15", "D", "Cl", "Cl37", "M", "NH4", "Z"})
#'
#' @return numeric; integer representation for passed in factor component
#'
#' @examples
#' CompFactorToInt("C")
#'
#' CompFactorToInt('C')
#'
#' @export
#'
CompFactorToInt2 <- function(x) {
  # CompFactorToInt - Component factor to integer
  # input:  factor {"C", "H", "O", "N", "S", "P", "E", "S34", "N15", "D", "Cl", "Cl37", "M", "NH4", "Z"}
  # output: numeric
  y <- 0
  switch(as.character(x),
         "C"   = {y <- 1},
         "H"   = {y <- 2},
         "O"   = {y <- 3},
         "N"   = {y <- 4},
         "S"   = {y <- 5},
         "P"   = {y <- 6},
         "E"   = {y <- 7},
         "S34" = {y <- 8},
         "N15" = {y <- 9},
         "D"   = {y <- 10},
         "Cl"  = {y <- 11},
         "Fl"  = {y <- 12},
         "Cl37"= {y <- 13},
         "Br"  = {y <- 14},
         "Br81"= {y <- 15},
         "I"   = {y <- 16},
         "M"   = {y <- 17},
         "NH4" = {y <- 18},
         "POE" = {y <- 19},
         "NOE" = {y <- 20},
         "Z"   = {y <- 21},
         { stop("CompFactorToInt called on undefined component ") })
  return(y)
}


#' Component integer to factor
#'
#' Given a component as an integer returns the correponding factor representation.
#' For example, 1 -> "C", 2 -> "H", ..., 9 -> "Z"
#'
#' This is an internal function that will not work outside the \code{\link{MFAssign}} function
#'  environment.
#'
#' @param x component; numeric; integer representation for component
#'
#'
#' @return factor (\code{"C", "H", "O", "N", "S", "P", "E", "S34", "N15", "D", "Cl", "Fl",
#'  "Cl37", "M", "NH4", "POE", "NOE", "Z"})
#'
#' @examples
#' CompIntToFactor(1)
#'
#' CompIntToFactor(2)
#'
#' @export
#'
CompIntToFactor2 <- function(x) {
  # CompIntToFactor - Component integer to factor
  # input: numeric
  # output: factor {"C", "H", "O", "N", "S", "P", "E", "S34", "N15", "D", "Cl", "Fl",
  # "Cl37", "M", "NH4", "POE", "NOE", "Z"}
  y <- NA
  switch(x,
         {y <- factor("C",
                      levels=c("C", "H", "O", "N", "S", "P", "E", "S34", "N15", "D", "Cl", "Fl",
                               "Cl37", "Br", "Br81", "I", "M", "NH4","POE", "NOE", "Z"))},
         {y <- factor("H",
                      levels=c("C", "H", "O", "N", "S", "P", "E", "S34", "N15", "D", "Cl", "Fl",
                               "Cl37", "Br", "Br81", "I", "M", "NH4","POE", "NOE", "Z"))},
         {y <- factor("O",
                      levels=c("C", "H", "O", "N", "S", "P", "E", "S34", "N15", "D", "Cl", "Fl",
                               "Cl37", "Br", "Br81", "I", "M", "NH4","POE", "NOE", "Z"))},
         {y <- factor("N",
                      levels=c("C", "H", "O", "N", "S", "P", "E", "S34", "N15", "D", "Cl", "Fl",
                               "Cl37", "Br", "Br81", "I", "M", "NH4","POE", "NOE", "Z"))},
         {y <- factor("S",
                      levels=c("C", "H", "O", "N", "S", "P", "E", "S34", "N15", "D", "Cl", "Fl",
                               "Cl37", "Br", "Br81", "I", "M", "NH4","POE", "NOE", "Z"))},
         {y <- factor("P",
                      levels=c("C", "H", "O", "N", "S", "P", "E", "S34", "N15", "D", "Cl", "Fl",
                               "Cl37", "Br", "Br81", "I", "M", "NH4","POE", "NOE", "Z"))},
         {y <- factor("E",
                      levels=c("C", "H", "O", "N", "S", "P", "E", "S34", "N15", "D", "Cl", "Fl",
                               "Cl37", "Br", "Br81", "I", "M", "NH4","POE", "NOE", "Z"))},
         {y <- factor("S34",
                      levels=c("C", "H", "O", "N", "S", "P", "E", "S34", "N15", "D", "Cl", "Fl",
                               "Cl37", "Br", "Br81", "I", "M", "NH4","POE", "NOE", "Z"))},
         {y <- factor("N15",
                      levels=c("C", "H", "O", "N", "S", "P", "E", "S34", "N15", "D", "Cl", "Fl",
                               "Cl37", "Br", "Br81", "I", "M", "NH4","POE", "NOE", "Z"))},
         {y <- factor("D",
                      levels=c("C", "H", "O", "N", "S", "P", "E", "S34", "N15", "D", "Cl", "Fl",
                               "Cl37", "Br", "Br81", "I", "M", "NH4","POE", "NOE", "Z"))},
         {y <- factor("Cl",
                      levels=c("C", "H", "O", "N", "S", "P", "E", "S34", "N15", "D", "Cl", "Fl",
                               "Cl37", "Br", "Br81", "I", "M", "NH4","POE", "NOE", "Z"))},
         {y <- factor("Fl",
                      levels=c("C", "H", "O", "N", "S", "P", "E", "S34", "N15", "D", "Cl", "Fl",
                               "Cl37", "Br", "Br81", "I", "M", "NH4","POE", "NOE", "Z"))},

         {y <- factor("Cl37",
                      levels=c("C", "H", "O", "N", "S", "P", "E", "S34", "N15", "D", "Cl", "Fl",
                               "Cl37", "Br", "Br81", "I", "M", "NH4","POE", "NOE", "Z"))},
         {y <- factor("Br",
                      levels=c("C", "H", "O", "N", "S", "P", "E", "S34", "N15", "D", "Cl", "Fl",
                               "Cl37", "Br", "Br81", "I", "M", "NH4","POE", "NOE", "Z"))},
         {y <- factor("Br81",
                      levels=c("C", "H", "O", "N", "S", "P", "E", "S34", "N15", "D", "Cl", "Fl",
                               "Cl37", "Br", "Br81", "I", "M", "NH4","POE", "NOE", "Z"))},
         {y <- factor("I",
                      levels=c("C", "H", "O", "N", "S", "P", "E", "S34", "N15", "D", "Cl", "Fl",
                               "Cl37", "Br", "Br81", "I", "M", "NH4","POE", "NOE", "Z"))},

         {y <- factor("M",
                      levels=c("C", "H", "O", "N", "S", "P", "E", "S34", "N15", "D", "Cl", "Fl",
                               "Cl37", "Br", "Br81", "I", "M", "NH4","POE", "NOE", "Z"))},
         {y <- factor("NH4",
                      levels=c("C", "H", "O", "N", "S", "P", "E", "S34", "N15", "D", "Cl", "Fl",
                               "Cl37", "Br", "Br81", "I", "M", "NH4","POE", "NOE", "Z"))},
         {y <- factor("POE",
                      levels=c("C", "H", "O", "N", "S", "P", "E", "S34", "N15", "D", "Cl", "Fl",
                               "Cl37", "Br", "Br81", "I", "M", "NH4","POE", "NOE", "Z"))},
         {y <- factor("NOE",
                      levels=c("C", "H", "O", "N", "S", "P", "E", "S34", "N15", "D", "Cl", "Fl",
                               "Cl37", "Br", "Br81", "I", "M", "NH4","POE", "NOE", "Z"))},
         {y <- factor("Z",
                      levels=c("C", "H", "O", "N", "S", "P", "E", "S34", "N15", "D", "Cl", "Fl",
                               "Cl37", "Br", "Br81", "I", "M", "NH4","POE", "NOE", "Z"))},
         { stop("CompIntToFactor called on undefined component ") })
  return(y)
}


#' Nominal masses of components
#'
#' Returns the nominal masses of molecular components, using zero for charge
#'
#' This is an internal fuction that will not work outside the \code{\link{MFAssign}} function
#'  environment.
#'
#' @param x component; factor (\code{"C", "H", "O", "N", "S", "P", "E", "S34", "N15", "D", "Cl", "Cl37", "M", "NH4", "Z"})
#'
#' @return numeric; value of nominal mass
#'
#' @examples
#' EM("C")
#'
#' EM('C')
#'
#' @export
#'
NM2 <- function(x) {
  # NM - nominal masses of components, using zero for charge
  # input:  factor {"C", "H", "O", "N", "S", "P", "E", "S34", "N15", "D", "Cl", "Cl37", "M", "NH4", "Z"}
  # output: numeric
  y <- 0
  switch(as.character(x),
         "C"   = {y <- 12},
         "H"   = {y <- 1},
         "O"   = {y <- 16},
         "N"   = {y <- 15},
         "S"   = {y <- 32},
         "P"   = {y <- 32},
         "E"   = {y <- 1},
         "S34" = {y <- 34},
         "N15" = {y <- 16},
         "D"   = {y <- 1},
         "Cl"  = {y <- 35},
         "Fl"  = {y <- 19},
         "POE" = {y <- 1},
         "NOE" = {y <- (-1)},
         "Cl37"= {y <- 37},
         "Br"  = {y <- 79},
         "Br81"  = {y <- 81},
         "I"  = {y <- 127},
         "M"   = {y <- 22},
         "NH4" = {y <- 17},
         "Z"   = {y <- 0},
         { stop("NM called on undefined component ") })
  return(y)
}


#' Expected masses of components
#'
#' Returns the nominal masses of molecular components, using the mass of an
#'   electron for charge
#'
#'   This is an internal fuction that will not work outside the \code{\link{MFAssign}} function
#'  environment.
#'
#' @param x component; factor (\code{"C", "H", "O", "N", "S", "P", "E", "S34", "N15", "D", "Cl", "Cl37", "M", "NH4", "Z"})
#'
#' @return numeric; value of expected mass
#'
#' @examples
#' EM("C")
#'
#' EM('C')
#'
#' @export
#'
EM2 <- function(x) {
  # EM - exact masses of component atoms and molecules, using the mass of an electron for Z
  # input:  factor {"C", "H", "D", "O", "N14", "N", "S", "P31", "P", "Na", "M", "E", "Z", "Cl", "Cl35", "Fl", "NH4", "NH4-", "Cl37","Cl37m "S34", "C13", "N15", "N15H"}
  # output: numeric
  y <- 0
  switch(as.character(x),
         "C"  = {y <- 12.0000000000},
         "H"  = {y <- 1.0078250319},
         #"D"  = {y <- 2.014102},
         "D"  = {y <- 1.006277},
         "O"  = {y <- 15.9949146223},
         "N14"  = {y <- 14.0030740074},
         "N" = {y <- 15.0108990393},
         "S"  = {y <- 31.9720707300},
         "P31"  = {y <- 30.97376149},
         "P" = {y <- 31.9815865219},
         "Na" = {y <- 22.989770},
         "M"  = {y <-21.9819446281},
         "E"  = {y <- 13.0033548380},
         "E2"  = {y <- 1.0033548380},
         "Z"  = {y <- 0.0005485799},
         "Cl35" = {y <- 34.968853},
         "Cl" = {y <- 35.976678},
         "Fl18" = {y <- 17.990578},
         "Fl19" = {y <- 18.998403},
         "Fl" = {y <- 20.006228},
         "POE" = {y <- 1.007825049},
         "NOE" = {y <- (-1.007825049)},
         "NH4+"  = {y <- 18.033823},
         "NH4" = {y <- 17.026550},
         "Cl37" = {y <- 37.973728},
         "Cl37m" = {y <- 36.965903},
         "S34" = {y <- 33.967868},
         "C13" = {y <- 1.0033548380},
         "N15" = {y <- 16.007934},
         "N15H" = {y <- 15.000109},
         "Br" = {y <- 79.926161},
         "Br79" = {y <- 78.918336},
         "Br81" = {y <- 81.924115},
         "Br81m" = {y <- 80.916290},
         "I"    = {y <- 127.91302},
         "I127" = {y <- 126.904477},
         { stop("EM called on undefined component ") })
  return(y)
}

#mass <- 18.998403+1.0078250319
#'
#' #' Minimum limit of moles (program limit)
#' #'
#' #' Returns the minimum limit of moles used by the
#' #'   \code{\link{MFAssign}} function.
#' #'
#' #'   This is an internal fuction that will not work outside the \code{\link{MFAssign}} function
#' #'  environment.
#' #'
#' #' @param x component; factor (\code{"C", "H", "O", "N", "S", "P", "E", "S34", "N15", "D", "Cl", "Cl37", "M", "NH4", "Z"})
#' #'
#' #' @return numeric; value of minimum limit of moles
#' #'
#' #' @examples
#' #' MinMoles("C")
#' #'
#' #' MinMoles('C')
#' #'
#' #' @export
#' #'
#' MinMoles2 <- function(x) {
#'   # MinMoles - program defined minimum moles of components,
#'   # input:  factor {"C", "H", "O", "N", "S", "P", "E", "S34", "N15", "D", "Cl", "Cl37", "M", "NH4", "Z"}
#'   # output: numeric
#'   y <- 0
#'   switch(as.character(x),
#'          "C"  = {y <- 1},
#'          "H"  = {y <-  2},
#'          "O"  = {y <- 0},
#'          "N"  = {y <- 0},
#'          "S"  = {y <- 0},
#'          "P"  = {y <- 0},
#'          "M"  = {y <- 0},
#'          "E"  = {y <- 0},
#'          "Z"  = {y <- 1},
#'          "Cl" = {y <- 0},
#'          "Fl" = {y <- 0},
#'          "POE" = {y <- 0},
#'          "NOE" = {y <- 0},
#'          "Cl37" = {y <- 0},
#'          "S34" = {y <- 0},
#'          "C13" = {y <- 0},
#'          "Br" = {y <- 0},
#'          "Br81" = {y <- 0},
#'          "I"  = {y <- 0},
#'          { stop("MinMoles called on undefined component ") })
#'   return(y)
#' }
#'
#' #' Maximum limit of moles (program limit)
#' #'
#' #' Returns the maximum limit of moles used by the
#' #'   \code{\link{MFAssign}} function.
#' #'
#' #'   This is an internal fuction that will not work outside the \code{\link{MFAssign}} function
#' #'  environment.
#' #'
#' #' @param x component; factor (\code{"C", "H", "O", "N", "S", "P", "M", "E", "Z"})
#' #'
#' #' @return numeric; value of maximum limit of moles
#' #'
#' #' @examples
#' #' MaxMoles("C")
#' #'
#' #' MaxMoles('C')
#' #'
#' #' @export
#' #'
#' MaxMoles2 <- function(x) {
#'   # MaxMoles - program maximum limit of moles
#'   # input:  factor {"C", "H", "O", "N", "S", "P", "E", "S34", "N15", "D", "Cl", "Cl37", "M", "NH4", "Z"}
#'   # output: numeric
#'   y <- 0
#'   switch(as.character(x),
#'          "C"  = {y <- 166},
#'          "H"  = {y <- 284},
#'          "O"  = {y <-  72},
#'          "N"  = {y <-   24},
#'          "S"  = {y <-   8},
#'          "P"  = {y <-   5},
#'          "M"  = {y <-   2},
#'          "E"  = {y <-   3},
#'          "Z"  = {y <-   5},
#'          "Cl" = {y <-   5},
#'          "Fl" = {y <-   5},
#'          "POE" = {y <-   1},
#'          "NOE" = {y <-   1},
#'          "Cl37" = {y <- 3},
#'          "S34" = {y <-  3},
#'          { stop("MaxMoles called on undefined component ") })
#'   return(y)
#' }
#'


#' Minimum limit of moles (user limit)
#'
#' Returns the minimum limit of moles used by the
#'   \code{\link{MFAssign}} function, specified by the user.
#'
#'   This is an internal fuction that will not work outside the \code{\link{MFAssign}} function
#'  environment.
#'
#' @param x component; factor (\code{"C", "H", "O", "N", "S", "P", "E", "S34", "N15", "D", "Cl", "Cl37", "M", "NH4", "Z"})
#'
#' @return numeric; value of minimum limit of moles
#'
#' @examples
#' LowMoles("C")
#'
#' LowMoles('C')
#'
#' @export
#'
LowMoles2 <- function(x) {
  # LowMoles - user suplied limit on minimum moles
  # input:  factor {"C", "H", "O", "N", "S", "P", "E", "S34", "N15", "D", "Cl", "Cl37", "M", "NH4", "Z"}
  # output: numeric
  y <- 0
  switch(as.character(x),
         "C"  = {y <- 1},
         "H"  = {y <- 2},
         "O"  = {y <- 0},
         "N"  = {y <- 0},
         "S"  = {y <- 0},
         "P"  = {y <- 0},
         "M"  = {y <- 0},
         "E"  = {y <- 0},
         "Z"  = {y <- 1},
         "Cl" = {y <- 0},
         "Fl" = {y <- 0},
         "POE" = {y <- 0},
         "NOE" = {y <- 0},
         "N15" = {y <- 0},
         "D" = {y <- 0},
         "NH4" = {y <- 0},
         "Cl37" = {y <- 0},
         "S34" = {y <- 0},
         "Br" = {y <- 0},
         "Br81" = {y <- 0},
          "I"  = {y <- 0},
         { stop("LowMoles called on undefined component ") })
  return(y)
}




#' Maximum limit of moles (user limit)
#'
#' Returns the maximum limit of moles used by the
#'   \code{\link{MFAssign}} function, specified by the user.
#'
#'  This is an internal fuction that will not work outside the \code{\link{MFAssign}} function
#'  environment.
#'
#' @param x component; factor (\code{"C", "H", "O", "N", "S", "P", "E", "S34", "N15", "D", "Cl", "Cl37", "M", "NH4", "Z"})
#'
#' @return numeric; value of maximum limit of moles
#'
#' @examples
#' HighMoles("C")
#'
#' HighMoles('C')
#'
#' @export
#'
HighMoles2 <- function(x, N=0, S=0, P=0, Cl=0, Fl = 0, POE = 0, NOE = 0, E=0, S34=0, Cl37=0,
                      N15=0, D=0, Br = 0, Br81 = 0, I = 0, M=0, NH4=0, Z=1) {
  # max_moles - user suplied limit on maximum moles
  # input:  factor {"C", "H", "O", "N", "S", "P", "E", "S34", "N15", "D", "Cl", "Cl37", "M", "NH4", "Z"}
  # output: numeric
  y <- 0
  switch(as.character(x),
         "C"  = {y <- 100},
         "H"  = {y <- 202},
         "O"  = {y <- 60},
         "N"  = {y <- N},
         "S"  = {y <- S},
         "P"  = {y <- P},
         "M"  = {y <- M},
         "E"  = {y <- E},
         "Z"  = {y <- Z},
         "Cl" = {y <- Cl},
         "Fl" = {y <- Fl},
         "POE" = {y <- POE},
         "NOE" = {y <- NOE},
         "N15" = {y <- N15},
         "D"  = {y <- D},
         "NH4" = {y <- NH4},
         "Cl37" = {y <- Cl37},
         "S34" = {y <- S34},
         "Br" = {y <- Br},
         "Br81" = {y <- Br81},
         "I" = {y <- I},
         { stop("HighMoles called on undefined component ") })
  return(y)
}

#' Valence of component
#'
#' Returns the valence of the specified component; used by the
#'   \code{\link{MFAssign}} function.
#'  This is an internal fuction that will not work outside the \code{\link{MFAssign}} function
#'  environment.
#'
#'
#' @param x component; factor (\code{"C", "H", "O", "N", "S", "P", "E", "S34", "N15", "D", "Cl", "Cl37", "M", "NH4", "Z"})
#'
#' @return numeric; value of valence
#'
#' @examples
#' Valence("C")
#'
#' Valence('C')
#'
#' @export
#'
Valence2 <- function(x) {
  # Valence - valences of the components
  # input:  factor {"C", "H", "O", "N", "S", "P", "E", "S34", "N15", "D", "Cl", "Cl37", "M", "NH4", "Z"}
  # output: numeric
  y <- 0
  switch(as.character(x),
         "C"  = {y <- 4},
         "H"  = {y <- 1},
         "O"  = {y <- 2},
         "N"  = {y <- 3},
         "S"  = {y <- 6},
         "P"  = {y <- 5},
         "M"  = {y <- 0},
         "E"  = {y <- 0},
         "Z"  = {y <- 0},
         "Cl" = {y <- 1},
         "Fl" = {y <- 1},
         "POE" = {y <- 0},
         "NOE" = {y <- 0},
         "N15" = {y <- 3},
         "D" = {y <- 1},
         "NH4" = {y <- 0},
         "Cl37" = {y <- 1},
         "S34" = {y <- 6},
         "Br" = {y <- 1},
         "Br81" = {y <- 1},
         "I"  = {y <- 1},
         { stop("Valence called on undefined component ") })
  return(y)
}


# ----------------------------------------------------------------------------
# Utility funcitons
# ----------------------------------------------------------------------------

#' Even integer
#'
#' Returns Boolean if number is Even
#'
#' @param x numeric
#'
#' @return Boolean
#'
#' @examples
#' Even(123)
#'
#' Even(12)
#'
#' @export
#'
Even <- function(x) {
  # even - function to test whether an integer is even or odd
  (round(x)%%2)==0
}

#' Determine if molecular formula is valid
#'
#' Evaluates if a molecular formula is valid, obeying senior rules and other
#' constraints
#'
#'  This is an internal fuction that will not work outside the \code{\link{MFAssign}} function
#'  environment.
#'
#' @param moles vector of integers (length 9) - indicating number of each
#'   elements
#'
#' @return boolean
#'
#' @examples
#' ValidFormula(moles)
#'
#' ValidFormula(c(2,3,1,0,0,0,0,0,1))
#'
#' @export
#'
ValidFormula2 <- function(moles) {
  #The function Valid evaluates a molecular formula to ensure that its
  #composition obeys the Senior Rules and meets other compositional
  #constraints that X >= Low[X] for X=CHONSPME and O <= (C+2+3*N+4*S+4*P).

  ok <- TRUE

  # Minimum constraints on CHONSPME formulae
  for (i in CompFactorToInt("C"):CompFactorToInt("NH4")) {
    ok <- (ok & moles[i] >= LowMoles(CompIntToFactor(i)))
  }

  # Senior Rule #1 for CHONSPME formulae
  # The sum of atoms having odd valences must be an even number. The use of
  # components results in only H having an odd valence.
  sum <- 0
  for (i in CompFactorToInt("C"):CompFactorToInt("NH4")) {
    if (!Even(Valence(CompIntToFactor(i)))) {
      sum <- sum + unlist(moles[i])
    }
    ok <- (ok & Even(sum))
  }

  # Senior Rule #2 for CHONSPME formulae
  # The sum of the valences must equal or exceed two times the maximum
  # valence.
  sum <- 0
  for (i in CompFactorToInt("C"):CompFactorToInt("NH4")) {
    if ((i==CompFactorToInt("N")) | (i==CompFactorToInt("P"))| (i==CompFactorToInt("Cl"))|
        (i==CompFactorToInt("N15")) | (i==CompFactorToInt("Cl37"))| (i==CompFactorToInt("Fl"))) {
      sum <- sum + unlist(moles[i])*Valence(CompIntToFactor(i)) + unlist(moles[i])*2
    } else {
      sum <- sum + unlist(moles[i])*Valence(CompIntToFactor(i))
    }
  }
  ok <- (ok & (sum >= 2*Valence("C")))

  # Senior Rule #3 for CHONSPME formulae
  # The difference between the maximum and minimum number of bonds that
  # can exist in a molecular formula is known as unsaturation (U) or double
  # bond equivalents (DBE), and U (or DBE) must be >= 0.

  sum <- 2
  for (i in CompFactorToInt("C"):CompFactorToInt("P")) {
    sum <- sum + unlist(moles[i])*(Valence(CompIntToFactor(i))-2)
  }
  ok <- (ok & (sum >= 0))

  # Other compositional constraints
  ok <- (ok & (moles[CompFactorToInt("O")] <= 2+unlist(moles[CompFactorToInt("O")])
               + 3*unlist(moles[CompFactorToInt("N")]) + 4*unlist(moles[CompFactorToInt("P")])
               + 4*unlist(moles[CompFactorToInt("S")])))
  return(ok)
}

#' Example Data, masslist
#' CHNOS_ML_Ex is a list of 2121 negative ion masses and their corresponding relative
#' abundances. It contains monoisotopic CHO, CHNO, and CHOS molecules and serves as an
#' example of data

