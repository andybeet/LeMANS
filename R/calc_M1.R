####################################################################################
# calculates other sources of natural mortality other than predation
# assumption: follows a beta distribution and that the small er the size class the greater the mortality
# this is based on the ratio of the sizeclass midpoint to the largest sizeclass midpoint (0,1) variable
calc_M1 <- function(nSize,nSpecies,lBound,mBound,alphaM1,betaM1,cM1,scLinfMat,scLinf,phiMin) {

  ## need to work out why scaled by cM1
  x <- outer((mBound / max(parameterValues$Linf)),rep(1,nSpecies))
  M1 <- stats::dbeta(x,alphaM1,betaM1)*cM1

  # largest size class has wrong midpoint since Linf for a species may be < the midpoint of the bin
  # the midpoint used in midway between lowerbounfd and Linf
  for (isp in 1:nSpecies) {
    midPoint <- lBound[scLinf[isp]]  +  (parameterValues$Linf[isp]-lBound[scLinf[isp]])/2
    x <- midPoint/max(parameterValues$Linf)
    M1[scLinf[isp],isp] <- stats::dbeta(x,alphaM1,betaM1)*cM1
  }
  # scaled to time step
  M1 <- M1*scLinfMat*phiMin

  return(M1)
}
