
#' Observation error specification for Geoduck
#'
#' Controls bias and imprecision of observed catches and indices
#'
#' @param Hist Object of class 'Hist' made by runMSE(OM, Hist=T)
#' @param Catch_CV Positive real number - the coeficient of variation of catch observations
#' @param Catch_bias Positive real number - the mean bias of the catch observations. 1 is unbiased, 1.1 is a 10 percent positive bias on average
#' @param Index_CV Positive real number - the coeficient of variation of relative or absolute abundance observations
#' @param Index_bias Positive real number - the mean bias of the relatie or absolute abundance observations. 1 is unbiased, 1.1 is a 10 percent positive bias on average
#' @author T. Carruthers
#' @export
Obs_GD = function(Hist, Catch_CV = 0.01, Catch_bias = 1,  Index_CV=0.05, Index_bias=1){
  Hist@SampPars$Obs$Cobs_y[] = trlnorm(prod(dim(Hist@SampPars$Obs$Cobs_y)), Catch_bias, Catch_CV)
  Hist@SampPars$Obs$AddIerr[] = trlnorm(prod(dim(Hist@SampPars$Obs$AddIerr)), Index_bias, Index_CV)
  Hist
}
