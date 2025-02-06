
#' Management Procedure for Geoduck
#'
#' A modifiable management procedure for Manila Clam that allows for adjustments to size limits, current effort and TAC control via indices
#'
#' @param x Positive integer - the simulation number
#' @param Data Object of class 'Data'
#' @param reps positive integer - not applicable - the number of stochastic draws of advice from which to sample a percentile
#' @param Min.size Positive real number - minimum size limit mm
#' @param Max.size Positive real number - maximum size limit (NaN is no limit)
#' @param CEff.Mult Positive real number (imperfect fraction) the multiplier on current fishing effort (fishing pressure)
#' @param C_I.targ Positive real number (imperfect fraction) TAC is calculated TAC(t+1) = I(t) x C(2022)/I(2022) x C_I.targ when TAC.calc = "Ratio"
#' @param I.targ Positive real number (imperfect fraction). TAC is reduced / increased to reach I.targ (a fraction of current index) when TAC.calc = "Target"
#' @param IS.targ Real number (imperfect fraction) TAC is reduced / increased to reach index target slope (IS.targ) when TAC.calc = "Slope"
#' @param IS.yrs Positive integer - the number of years to evalute index slope for the IS rule.
#' @param IS.fac Positive real number - sensitivity of the slope - TAC change rule. Lower values are less sensitive. A value of 1 makes changes in proportion to index changes.
#' @param TAC.calc Character string. Can be NaN, "Ratio", "Target", "Slope" where TACs are either not constrained by data or set either by index ratio (C_I.targ), index target (I.targ) or index slope target (IS.targ)
#' @param maxTAC Positive real number (imperfect fraction) - annual catches cannot exceed current catches muliplied by this factor
#' @param minTAC Positive real number (imperfect fraction) - annual catches cannot be lower than current catches muliplied by this factor
#' @param TACdec Positive real number (fraction) - downward TAC changes cannot exceed this fraction (e.g. 0.2 is maximum decline of 20 per cent among management cycles)
#' @param TACinc Positive real number (fraction) - upward TAC changes cannot exceed this fraction (e.g. 0.1 is maximum increase of 10 per cent among management cycles)
#' @param I.enp Positive real number (fraction) the parameter controlling effective number of parameters for the polynomial smoother on the indices. npar = ny * I.enp so higher values mean more smoother parameters and less smoothing
#' @param I_freq Vector of positive integers - how frequently do you collect the indices (0 = do not collect, 1 = every year, 2 = every other year, 3 = once every three years, ...). Default is "auto" - the absolute abundance index is located and updated every year
#' @param rotation Vector proyears long of catch on/off (1 and 0, respectively). If a single value is specified (1, 2,..) the pattern alternates on once every year (1), every other year (2) etc. The default is "auto" where the MP tried to detect the most recent pattern
#' @param calib_yrs Positive integer - how many of the recent years are used to calibrate Index to catch ratio for TAC based MPs
#' @param HCR_CP_B A positive numeric vector two positions long of biomass control points (x axis) relative to recent index / curI_2_target for a hockey stick harvest control rule c(0,0) essentially has no control rule if HCR_CP_TAC = 1
#' @param HCR_CP_TAC A positive numeric vector two positions long that is the fraction of the recommended TAC taken below control point 1 and above control point 2 (y axis adjustment of the harvest control rule)
#' @param curI_2_target A positive real number that is the level of the current index (recent historical year) relative to the target biomass level (e.g. BMSY) 2 implies recent indices are at twice target levels (underexploited)
#' @param DR Fraction - the discard rate
#' @param Fdisc The fraction of discarded individuals that die
#' @examples
#' MP.GD(Simulation_1) # apply to a generic simulated dataset from openMSE
#' @author T. Carruthers
#' @export
MP.GD = function(x, Data, reps=1, Min.size = NaN, Max.size = NaN, CEff.Mult = NaN, C_I.targ = 1.0,
                 I.targ = 0.5, IS.targ = 0, IS.yrs = 6, IS.fac = 1, TAC.calc = "Ratio", maxTAC = 5.0, minTAC = 0.1,
                 TACdec = 0.2, TACinc = 0.2, I.enp = 0.25, I_freq = "auto", rotation ="auto",calib_yrs = 3,
                 HCR_CP_B = c(0, 0), HCR_CP_TAC = c(0,1), curI_2_target = 2,
                 DR = 0, Fdisc = 0.5){

  # x = readRDS("C:/temp/x_GD.rds"); Data = readRDS("C:/temp/Data_GD.rds"); reps=1; Min.size = NaN; Max.size = NaN; CEff.Mult = NaN; C_I.targ = 1.0;  I.targ = 0.5; IS.targ = 0; IS.yrs = 6; IS.fac = 1; TAC.calc = "Ratio"; maxTAC = 5.0; minTAC = 0.1;  TACdec = 1; TACinc = 1E5; I.enp = 0.25; I_freq = "auto"; rotation = "auto"; calib_yrs = 3;  HCR_CP_B = c(0, 0); HCR_CP_TAC = c(0,1); curI_2_target = 2;  DR = 0; Fdisc = 0.5

  dependencies = "Data@Cat, Data@AddInd"
  ny = length(Data@Year)
  if(ny==54){saveRDS(x,"C:/temp/x_GD.rds"); saveRDS(Data,"C:/temp/Data_GD.rds");stop()}

  nI = dim(Data@AddInd)[2]
  if(I_freq == "auto"){
    I_freq = rep(0,nI)
    LHind = match(Data@LHYear, Data@Year)
    absind = which.max(Data@AddInd[x,,LHind])
    I_freq[absind] = 1
  }


  ystart <- which(!is.na(Data@Cat[x, ]))[1]
  yind <- ystart:length(Data@Cat[x, ])
  LHYr = Data@LHYear
  LHYrInd = match(LHYr, Data@Year)
  proyears = dim(Data@Misc$ReferencePoints$ByYear$BMSY)[2] - LHYrInd
  CurYr = max(Data@Year)

  Year <- Data@Year[yind]
  C_hist <- Data@Cat[x, yind]

  Crel <- C_hist  / mean(C_hist)
  on_off <- Crel > 0.01
  laston = max((1:LHYrInd)[on_off[1:LHYrInd]])
  fracon = mean(on_off[LHYrInd-9:0])
  if(rotation == "auto")  rotation = ceiling(1/fracon)-1

  lastMPrec = max((1:ny)[on_off])
  MPrec = Data@Cat[x,lastMPrec] # if not available assume last historical catch observation

  if(length(rotation)==1){ # user has not specified a frequency for every proyears
    makeNA = rep(c(rep(FALSE,rotation),TRUE),proyears*2)
    lagon= LHYrInd-laston
    rotation = as.integer(makeNA[lagon+1:proyears])
  }

  I_hist <- Data@AddInd[x,,yind]  # all indices
  I_hist[I_hist < 0] <- NA        # zeros are NAs
  nkeep = sum(I_freq!=0)
  I_keep = I_obs_freq(I_hist, I_freq, LHYr, CurYr, Year) # This is the code that filters out future years where observations are not available
  I_smth = array(NA,dim(I_keep))

  caliby = LHYrInd-(calib_yrs-1):0
  calibmuI = apply(I_keep[,caliby,drop=F],1,mean,na.rm=T)
  calibmuC = mean(C_hist[caliby][on_off[caliby]],na.rm=T)
  C_I = calibmuC / calibmuI # historical catch per index (average over calib_yrs)

  for(j in 1:nkeep) I_smth[j,] = smoothy(I_keep[j,], enp_mult = I.enp) #, plot=x==1)

  Rec = new('Rec')

  if(!is.na(CEff.Mult)){  # if Effort is imposed (defaults to 1)
    Rec@Effort = CEff.Mult
  }

  if(!is.na(TAC.calc)){ # if TACs are imposed

    if(TAC.calc == "Ratio"){
      TACbyI = I_smth[,ny] * C_I                                   # basic TAC is smoothed index multiplied by current Catch / Index
      TACtemp = mean(TACbyI,na.rm=T) * C_I.targ                    # modified by tuning parameter
      if(is.na(TACtemp)|is.null(TACtemp))TACtemp = Data@MPrec[x]   #
      trial_TAC = TACfilter(TACtemp)
      ref = mean(calibmuI, na.rm=T) / curI_2_target
      est = mean(I_smth[,ny], na.rm=T)
      HCRadj_TAC = doHCR(trial_TAC, est = est, ref = ref, CP = HCR_CP_B, CPy = HCR_CP_TAC)
      mod = HCRadj_TAC/MPrec                              # implied TAC change

    }

    if(TAC.calc == "Target"){
      est = mean(I_smth[,ny], na.rm=T)
      mod = est /  (calibmuI * I.targ)
    }

    if(TAC.calc == "Slope"){

      Is0 = I_keep/calibmuI
      Is1 = Is0[length(Is0)-((IS.yrs-1):0)]
      Is2 = Is1[!is.na(Is1)]
      ns = length(Is2)
      slp = lm(y~x,data=data.frame(x=1:ns,y=Is2))$coefficients[[2]]
      mod = exp((slp-IS.targ)*IS.fac)
    }
    #
    Rec = doRec(Rec, MPrec, mod, TACdelta = c(TACdec, TACinc), TACrng = calibmuC*c(minTAC, maxTAC))

  } # end of TAC calcs

  # Minimum size limits
  if(!is.na(Min.size)){
    Rec@LR5 = Min.size*0.80
    Rec@LFR = Min.size *1.25
  }
  # Maximum size limits
  if(!is.na(Max.size)) Rec@HS = Max.size

  # Rotation
  py = ny-LHYrInd+1
  Rec@TAC = Rec@TAC * rotation[py]
  Rec@Effort = Rec@Effort*rotation[py]
  Rec

}
class(MP.GD) = c('MP', 'iMP')
