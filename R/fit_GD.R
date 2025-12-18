
#' Condition an operating model for Geoduck
#'
#' Uses the RCM model of OpenMSE to fit an operating model to data
#'
#' @param In A list object of class ('In') that includes the operating model parameters (slot 2, class OM) and data (slot 3, class RCMinput) for conditioning
#' @param sims Integer or vector of integers - the number of simulations for the operating model (e.g. 96) or the specific vector of simulations for the operating model (e.g. 13:48)
#' @param max_F Positive real number - the maximum instantaneous mortality rate in any historical time step
#' @param comp_like The likelihood function used for composition data (age / length) c("multinomial", "lognormal", "mvlogistic", "dirmult1", "dirmult2") see ?RCM
#' @param resample Logical - should parameters be drawn from the variance covariance matrix of the MLE fit (T) or just the MLE parameter values (F)
#' @param parallel Logical - should model fitting use multiple cores
#' @param silent Logical - should all RCM messages be repressed (ie not read out onto the console)?
#' @examples
#' cond.GD("In.GD.7")
#' @author T. Carruthers
#' @seealso \link{RCM}
#' @export
cond.GD = function(RCMinput, sims = 12, max_F = 0.22, comp_like = "lognormal", resample = T, parallel=F, silent=T){
  cores=1
  if(parallel){
    setup()
    cores=parallel::detectCores()/2
  }
  if(length(sims) == 1) simy = 1:sims
  if(length(sims) > 1) simy = sims
  OM = SubCpars(RCMinput$OM, simy)
  RCMfit = RCM(OM, RCMinput[[3]], s_selectivity=spec_sel(RCMinput),
               max_F = max_F, mean_fit = T, condition = "catch", cores = cores,
               comp_like=comp_like, priors = list(M = c(0.05,0.15)), LWT=list(IAA=1),
               drop_nonconv=T, drop_highF=T, resample = resample, silent=silent);
  RCMfit@OM@Name = paste0("Geoduck Stat Area ",RCMinput[[1]])
  RCMfit@Misc$keepsims = Geoduck_OMfilter(RCMfit) # removes simulations with completely unrealistic recruitment deviations and age_selectivity
  RCMfit@OM = SubCpars(RCMfit@OM, RCMfit@Misc$keepsims)
  RCMfit

}


#  maxRecDev = 20; vrefs1 = list(index=1,age=6,yr=1,maxV=0.5); vrefs2 = list(index=1,age=80,yr=1,minV=0.95)
Geoduck_OMfilter = function(RCMfit,
                            maxRecDev = 20,
                            vrefs1 = list(index=1,age=6,yr=1,maxV=0.5),
                            vrefs2 = list(index=1,age=80,yr=1,minV=0.95)){

  nsub = RCMfit@data@Misc$nsub # number of subareas with age data (n age surveys)
  isCPUE = RCMfit@data@Misc$isCPUE # is there CPUE data?
  Bsuvpos = nsub+isCPUE+1 # position of absolute biomass survey

  keep1 = apply(RCMfit@OM@cpars$Perr_y,1,max)<maxRecDev
  #keep2 = sapply(RCMfitobj@Misc,function(x)x$ivul[vrefs1$yr,vrefs1$age,vrefs1$index]<vrefs1$maxV)

  # Biomass must be within 95% interval of current biomass estimate in index[last index]
  Bmus = RCMfit@data@Index
  Bcvs =  RCMfit@data@I_sd

  mu = Bmus[nrow(Bmus),Bsuvpos]
  cv = Bcvs[nrow(Bmus),Bsuvpos]
  trunc = qlnorm(c(0.025,0.975),log(mu),cv)
  OM = RCMfit@OM
  Best = t(sapply(RCMfit@report,function(x)x$B))
  Bfinal = Best[,OM@nyears]
  B1 = Best[,1]

  # B changes can't be more than 2/3x or 3x
  keep3 = Bfinal > trunc[1] & Bfinal < trunc[2] & Bfinal > (1/3 * B1) & Bfinal < (3 * B1)
  keep = keep1 & keep3
  #RCMfit@OM=SubCpars(RCMfit@OM,sims=keep)
  #RCMfitobj@Misc$keepsims = keep
  (1:length(keep))[keep]

}

spec_sel=function(RCMinput){
  nBsuv = RCMinput[[3]]@Misc$nBsuv # number of biomass survey time series
  nsub = RCMinput[[3]]@Misc$nsub # number of subareas with age data (n age surveys)
  nBbio = RCMinput[[3]]@Misc$nBbio
  isCPUE = RCMinput[[3]]@Misc$isCPUE # is there CPUE data?
  if(nsub==0)s_selectivity = c(rep(1,isCPUE),rep("B",nBsuv+nBbio))
  if(nsub>0)s_selectivity = c(rep("dome_age",nsub-1),"logistic_age",rep(1,isCPUE),rep("B",nBsuv+nBbio))
  s_selectivity
}
