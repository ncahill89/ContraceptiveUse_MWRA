WriteModel <- function(# Writes BUGS model
  ### Writes BUGS model to output.dir (stored in mcmc.meta)
  mcmc.meta ##<< Object
  ){


#------ MODEL ----------
cat("
#--------------------------------------------------------------
# Model for contraceptive use
# Leontine Alkema, 2011
#--------------------------------------------------------------

model{

tau.tot.st <- tau.tot*(1-pow(rho.tot,2))
rho.tot ~ dunif(0,rho.max) 
tau.tot <- pow(sigma.tot, -2)
sigma.tot ~ dunif(0.01, sigma.ar.max) 

tau.rat.st <- tau.rat*(1-pow(rho.rat,2))
rho.rat ~  dunif(0,rho.max) 
tau.rat <- pow(sigma.rat, -2)
sigma.rat ~ dunif(0.01, sigma.ar.max) 

tau.unmet.st <- tau.unmet*(1-pow(rho.unmet,2))
rho.unmet ~ dunif(0,rho.max.unmet) 
tau.unmet <- pow(sigma.ar.unmet, -2)
sigma.ar.unmet ~ dunif(0.01, sigma.ar.max.unmet) 

for (c in 1:C){
  # gett.ci[c,i] gives year index for obs i in country c
  # geti.j # index of the obs. year
  # (for country c, if the obs years are sorted in increasing order), 
  # that obs j refers to
  # gett.ci[c,i] # t is the index in (start.c, end.c) that obs i refers to, 
  # where i=1 refers to first obs year
  theta.ci[c,1] ~ dnorm(0, tau.unmet.st)
  eps.ci[c, 1] ~ dnorm(0, tau.tot.st)
  eta.ci[c,1] ~ dnorm(0, tau.rat.st)
}
# AR-loop for countries with more than 1 obs
for (z in 1:n.countriesmorethan1obs){
  for (i in 2:N.unique.c[getc.z[z]]){
    thetahat.ci[getc.z[z],i] <- theta.ci[getc.z[z],i-1]*
      pow(rho.unmet, gett.ci[getc.z[z],i] - gett.ci[getc.z[z],i-1])
    epshat.ci[getc.z[z],i] <- eps.ci[getc.z[z],i-1]*
      pow(rho.tot, gett.ci[getc.z[z],i] - gett.ci[getc.z[z],i-1])
    etahat.ci[getc.z[z],i] <- eta.ci[getc.z[z],i-1]*
      pow(rho.rat, gett.ci[getc.z[z],i] - gett.ci[getc.z[z],i-1])
    tautheta.ci[getc.z[z],i] <- tau.unmet.st/
      (1-pow(rho.unmet, 2*(gett.ci[getc.z[z],i] - gett.ci[getc.z[z],i-1])))
    taueps.ci[getc.z[z],i] <- tau.tot.st/
      (1-pow(rho.tot, 2*(gett.ci[getc.z[z],i] - gett.ci[getc.z[z],i-1])))
    taueta.ci[getc.z[z],i] <- tau.rat.st/
      (1-pow(rho.rat, 2*(gett.ci[getc.z[z],i] - gett.ci[getc.z[z],i-1])))
    theta.ci[getc.z[z],i] ~ dnorm(thetahat.ci[getc.z[z],i], tautheta.ci[getc.z[z],i])
    eps.ci[getc.z[z],i] ~ dnorm(epshat.ci[getc.z[z],i], taueps.ci[getc.z[z],i])
    eta.ci[getc.z[z],i] ~ dnorm(etahat.ci[getc.z[z],i], taueta.ci[getc.z[z],i])
  }
}# end AR stuff
#### til here it can be moved into ar loop if par's are fixed as well

#--------------------------------------------------------------
for (c in 1:C){
  for (i in 1:N.unique.c[c]){
    mu.ci[c, i] <- omega.c[c]*(gett.ci[c,i] - T.c[c]) 
    pstar.ci[c, i] <- pmax.c[c]/(1+exp(-mu.ci[c, i]))
    # p.ci defined later, depends on whether AR is added or not
    Rmu.ci[c, i] <- Romega.c[c]*(gett.ci[c,i] - RT.c[c]) 
    Rstar.ci[c, i] <- Rmax.c[c]/(1+exp(-Rmu.ci[c, i]))
    # R.ci defined later, depends on whether AR is added or not
  } # end unique-obs year loop

  # Asymptotes: 1-level
  pmax.c[c] <- (exp(logitpmax.c[c])+0.5)/(1+exp(logitpmax.c[c]))
  logitpmax.c[c] ~ dnorm(lp.world, tau.lpc)
  Rmax.c[c] <- (exp(logitRmax.c[c])+0.5)/(1+exp(logitRmax.c[c]))
  logitRmax.c[c] ~ dnorm(lr.world, tau.lrc)

  # Omegas: 3-level
  #omega.c[c] <- 1/(1+exp(-logitomega.c[c]))
  #logitomega.c[c] ~ dnorm(w.subreg[subreg.c[c]], tau.wc)T(-4.5, 0) # in Bugs: use I()
  #Romega.c[c] <- 1/(1+exp(-logitRomega.c[c]))
  #logitRomega.c[c] ~ dnorm(Rw.subreg[subreg.c[c]], tau.Rwc)T(-4.5, 0) # in Bugs: use I()
  omega.c[c] <- (0.5*exp(logitomega.c[c])+0.01)/(1+exp(logitomega.c[c]))
  logitomega.c[c] ~ dnorm(w.subreg[subreg.c[c]], tau.wc)
  Romega.c[c] <- (0.5*exp(logitRomega.c[c])+0.01)/(1+exp(logitRomega.c[c]))
  logitRomega.c[c] ~ dnorm(Rw.subreg[subreg.c[c]], tau.Rwc)

  # Midyear ratio modern/total: 3-level
  RT.c[c] ~ dnorm(RT.subreg[subreg.c[c]], tau.RTc)T(1800,)

  # unmet
  unmet.intercept.c[c] ~ dnorm(unmet.subreg[subreg.c[c]], tau.unmetc)
} # end country loop

# dummy # change JR, 20131105
pmax.c[C+1] <- 0
omega.c[C+1] <- 0
RT.c[C+1] <- 0
Rmax.c[C+1] <- 0
Romega.c[C+1] <- 0
T.c[C+1] <- 0
unmet.intercept.c[C+1] <- 0

# Midyear total: rich countries
for (i in 1:n.rich){
   T.c[crich.index[i]] ~ dnorm(Tearlier, tau.earlierTc)T(1800,)
}
# Midyear total: other countries, 3-level
for (i in 1:n.notrich){
  T.c[cnotrich.index[i]] ~ dnorm(T.subreg[subreg.c[cnotrich.index[i]]], tau.Tc)T(1800,)
}

for (subreg in 1:n.subreg){
  unmet.subreg[subreg] ~ dnorm(0, tau.unmetworld)
	w.subreg[subreg] ~ dnorm(w.reg[reg.subreg[subreg]] , tau.wsubreg)#T(-4.5, 0)
	T.subreg[subreg] ~ dnorm(T.reg[reg.subreg[subreg]] , tau.Tsubreg)#T(1800,)
	Rw.subreg[subreg] ~ dnorm(Rw.reg[reg.subreg[subreg]] , tau.Rwsubreg)#T(-4.5, 0)
	RT.subreg[subreg] ~ dnorm(RT.reg[reg.subreg[subreg]] , tau.RTsubreg)#T(1800,)	
}

for (reg in 1:n.reg){
	w.reg[reg] ~ dnorm(w.world, tau.wreg)#T(-4.5, 0)
	T.reg[reg] ~ dnorm(T.world, tau.Treg)#T(1800,)
	Rw.reg[reg] ~ dnorm(Rw.world, tau.Rwreg)#T(-4.5, 0)
	RT.reg[reg] ~ dnorm(RT.world, tau.RTreg)#T(1800,)
}

lp.world ~ dnorm(0,0.01)
lr.world ~ dnorm(0,0.01)

# omegas
#Rw.world ~ dnorm(-2, 0.01)
#w.world ~ dnorm(-2, 0.01)
Rw.world ~ dnorm(-1, 0.01)
w.world ~ dnorm(-1, 0.01)
tau.wreg <- pow(sigma.wreg, -2)
tau.wsubreg <- pow(sigma.wsubreg, -2)
tau.Rwreg <- pow(sigma.Rwreg, -2)
tau.Rwsubreg <- pow(sigma.Rwsubreg, -2)
sigma.wsubreg ~ dunif(0,sigmawregsubreg.upper)
sigma.wreg ~ dunif(0,sigmawregsubreg.upper)
sigma.Rwreg ~ dunif(0,sigmawregsubreg.upper)
sigma.Rwsubreg ~ dunif(0,sigmawregsubreg.upper)

# Midpoints
RT.world ~ dnorm(mean.RTworld, tau0.T) 
T.world ~ dnorm(mean.Tworld, tau0.T)
Tearlier ~ dnorm(mean.Tearlier, tau0.T) 

tau.RTreg <- pow(sigma.RTreg, -2)
tau.RTsubreg <- pow(sigma.RTsubreg, -2)
tau.Treg <- pow(sigma.Treg, -2)
tau.Tsubreg <- pow(sigma.Tsubreg, -2)

sigma.RTreg ~ dunif(0,sigmaTregsubreg.upper)
sigma.RTsubreg ~ dunif(0,sigmaTregsubreg.upper)
sigma.Treg ~ dunif(0,sigmaTregsubreg.upper)
sigma.Tsubreg ~ dunif(0,sigmaTregsubreg.upper)

#--------------------------------------------------------------
# Likelihood
# gett.j gives (year - year.start+1) for obs j
# getc.j gives gives c for obs j
# ind.j refers to counter (1 if not applicable), ind1 refers to yes/no (1/0).

for (j in 1:J){ 
  trad.j[j] <- p.ci[getc.j[j], geti.j[j]] * (1-R.ci[getc.j[j], geti.j[j]])
  modern.j[j] <- p.ci[getc.j[j], geti.j[j]] * (R.ci[getc.j[j], geti.j[j]])
  unmet.j[j] <- (1-p.ci[getc.j[j], geti.j[j]])/(1+exp(-logitZ.j[j]))
  logitZstar.j[j] <- (
        unmet.intercept.c[getc.j[j]]
        + a.unmet 
        + b.unmet * (p.ci[getc.j[j], geti.j[j]] - pmid.for.unmet) 
        + c.unmet*pow(p.ci[getc.j[j], geti.j[j]] - pmid.for.unmet,2)) 
  # logitZ.j defined in AR-loop, just adding theta

  sump.j[j] <- (trad.j[j]*Vtrad.j[j] 
      + modern.j[j]* Vmodern.j[j] 
      + (1- p.ci[getc.j[j], geti.j[j]]))
  p.perturb.ij[1,j] <-  trad.j[j]*Vtrad.j[j]/sump.j[j]
  p.perturb.ij[2,j] <-  modern.j[j]* Vmodern.j[j]/sump.j[j]
  p.perturb.ij[3,j] <- unmet.j[j]/sump.j[j]
  p.perturb.ij[4,j] <- (1- trad.j[j] - modern.j[j] - unmet.j[j])/sump.j[j]
  
  folkbias.j[j] <- step(folk.ind1.j[j]-0.5)*v.folk* p.perturb.ij[3,j]
  micsbias.j[j] <- step(source.MICS.ind1.j[j]-0.5)* v.mics * p.perturb.ij[1,j]
  modposbias.j[j] <- step(mpos.ind1.j[j]-0.5)*v.mpos* p.perturb.ij[4,j]
  modnegbias.j[j] <- step(mneg.ind1.j[j]-0.5)*v.mneg * p.perturb.ij[2,j]

  q.ij[1,j] <- p.perturb.ij[1,j] - micsbias.j[j] + folkbias.j[j]
  q.ij[2,j] <- p.perturb.ij[2,j] + modposbias.j[j] - modnegbias.j[j]
  q.ij[3,j] <- p.perturb.ij[3,j] + micsbias.j[j] - folkbias.j[j]
  q.ij[4,j] <- p.perturb.ij[4,j] - modposbias.j[j] + modnegbias.j[j]
  none.adj.j[j] <- max(0.01,q.ij[3,j]) + q.ij[4,j]
  
  mu.jn[j,1] <- log(max(0.01, q.ij[1,j])/none.adj.j[j])
  mu.jn[j,2] <- log(max(0.01, q.ij[2,j])/none.adj.j[j])

  Vtrad.j[j] <- (
        V.geo.12i[1,geo.ind.j[j]] 
        * V.age.12i[1,age.ind.j[j]] 
        * V.hw.12i[1,hw.ind.j[j]] 
        * V.emal.12i[1,emal.ind.j[j]] 
        * V.sa.12i[1,sa.ind.j[j]]
        * V.posbias.12i[1,posbias.ind.j[j]]
        * V.posage.12i[1, posage.ind.j[j]]
        * V.negage.12i[1, negage.ind.j[j]] )
  
  Vmodern.j[j] <- (
        V.geo.12i[2,geo.ind.j[j]]
        * V.age.12i[2,age.ind.j[j]] 
        * V.hw.12i[2,hw.ind.j[j]] 
        * V.emal.12i[2,emal.ind.j[j]] 
        * V.sa.12i[2,sa.ind.j[j]]
        * V.posbias.12i[2,posbias.ind.j[j]]
        * V.posage.12i[2, posage.ind.j[j]]
        * V.negage.12i[2, negage.ind.j[j]] )

  T.j[j,1:2, 1:2] <- T.s[source.ind.j[j],,]
} # end loop J observations

# Unmet observations, get mean, skip NAs
for(k in 1:N.unmet) { 
 logitratio.yunmet.hat.j[getj.unmet.k[k]] <-  
        logit(max(0.01,q.ij[3,getj.unmet.k[k]])/none.adj.j[getj.unmet.k[k]])
}

#--------------------------------------------------------------
# Training sets
for (k in 1:n.training.breakdown){
  ratios.trad.modern.jn[getj.training.k[k],1:2] ~ 
    dmnorm(mu.jn[getj.training.k[k], ],T.j[getj.training.k[k],,])
}

# unmet training
for (k in 1:n.training.unmet){
 logitratio.yunmet.j[getj.training.unmet.k[k]] ~ dnorm(
      logitratio.yunmet.hat.j[getj.training.unmet.k[k]],    
      tau.unmet.source.s[source.ind.unmet.j[getj.training.unmet.k[k]]])
} 

# total training set
for (k in 1:n.training.tot){
##  logit.ytothat.j[j] <- logit(1-none.adj.j[j])
### temp
  logit.ytothat.j[getj.training.tot.k[k]] <- logit(max(0.01, 1-none.adj.j[getj.training.tot.k[k]]))
  logit.ytot.j[getj.training.tot.k[k]] ~ dnorm(
    logit.ytothat.j[getj.training.tot.k[k]], tau.sourcetot)
} 

#--------------------------------------------------------------
# data model parameters

# Biases: Folk, MICS, sterilization included/excluded (M+/M-)
v.mics ~ dunif(0,1)
v.mneg ~ dunif(0,1)
v.folk ~ dunif(0,1)
v.mpos ~ dunif(0,1)

# Multipliers V in [0,inf): geo, emal, hw, age other, sa (for trad only)
V.sa.12i[1,1] <- 1
V.geo.12i[1,1] <- 1 
V.geo.12i[2,1] <- 1 
V.age.12i[1,1] <- 1 
V.age.12i[2,1] <- 1 
V.hw.12i[1,1] <- 1 
V.hw.12i[2,1] <- 1 
V.emal.12i[1,1] <- 1 
V.emal.12i[2,1] <- 1 
for (m in 1:2){
  for (i in 2:ncat.emal){ 
    V.emal.12i[m,i] ~ dlnorm(0, tau.geo.m[m]) 
  }
  for (i in 2:ncat.hw){ 
    V.hw.12i[m,i] ~ dlnorm(0, tau.geo.m[m]) 
  }
  for (i in 2:ncat.geo){ 
    V.geo.12i[m,i] ~ dlnorm(0, tau.geo.m[m]) 
  }
  for (i in 2:ncat.age){ 
    V.age.12i[m,i] ~ dlnorm(0, tau.geo.m[m]) 
  }
  tau.geo.m[m] <- pow(sigma.geo.m[m], -2)
  sigma.geo.m[m] ~ dunif(0.01,2)
} # end m-loop
    
for (i in 2:ncat.sa){ 
  V.sa.12i[1,i] ~ dlnorm(0, tau.geo.m[1])
}

# multipliers V > 1: sa (for modern only), posbias, age+
# multipliers V < 1: age-

V.sa.12i[2,1] <- 1 
V.posbias.12i[1,1] <- 1 
V.posbias.12i[2,1] <- 1 
V.posage.12i[1,1] <- 1 
V.posage.12i[2,1] <- 1
V.negage.12i[1,1] <- 1 
V.negage.12i[2,1] <- 1 

# m = 2:
for (i in 2:ncat.sa){ 
  W.sa.12i[2,i] ~ dlnorm(mu.pos.m[2], tau.pos)
  V.sa.12i[2,i] <- 1+W.sa.12i[2,i]
}
for (i in 2:ncat.posbias){ 
  W.posbias.12i[2,i] ~ dlnorm(mu.pos.m[2], tau.pos) 
  V.posbias.12i[2,i] <- 1+W.posbias.12i[2,i]
}
for (i in 2:ncat.posage){ 
  W.posage.12i[2,i] ~ dlnorm(mu.pos.m[2], tau.pos) 
  V.posage.12i[2,i] <-1+W.posage.12i[2,i]
}
for (i in 2:ncat.negage){ 
  W.negage.12i[2,i] ~ dlnorm(mu.pos.m[2], tau.pos) 
  V.negage.12i[2,i] <- 1/(1+W.negage.12i[2,i])
}
tau.pos <- pow(sigma.pos, -2)
sigma.pos ~ dunif(0.01,2)
mu.pos.m[2] ~ dnorm(-2, 0.64)# 1/1.25^2 ) # 0.01)

# m=1
# note: could simplify code and throw out these V's
for (i in 2:ncat.posbias){ 
  V.posbias.12i[1,i] <- 1+exp(mu.pos.m[1])
}
for (i in 2:ncat.posage){ 
  V.posage.12i[1,i] <- 1+exp(mu.pos.m[1])
}
for (i in 2:ncat.negage){ 
  V.negage.12i[1,i] <- 1/(1+exp(mu.pos.m[1]))
}
#tau.pos.m[1] <- pow(sigma.pos.m[1], -2)
#sigma.pos.m[1] ~ dunif(0,2)
mu.pos.m[1] ~ dnorm(-2, 0.64)#~ dnorm(-2, 0.01)#T(-10,) # in Bugs: use I()


#--------------------------------------------------------------
# Source variances
tau.sourcetot ~ dgamma(0.5, halfsigma2.sourcetot0)
sigma.sourcetot <- 1/sqrt(tau.sourcetot)
for (s in 1:4){
  T1.source.s[s] <- T.s[s,1,1] 
  T2.source.s[s] <- T.s[s,2,2]
  T12.source.s[s] <- T.s[s,1,2]
  T.s[s,1:2,1:2] ~ dwish(R[1:2,1:2], 3) #change JR, 20131113: from 4 to 3 # change JR, 20131031: from 2 to 4
}

#--------------------------------------------------------------
a.unmet ~ dnorm(a0.unmet, tau.a0)
b.unmet ~ dnorm(b0.unmet, tau.b0) 
c.unmet ~ dunif(-10,0)

tau.unmet.source.s[1] <- pow(sigma.unmet.dhs,-2)
tau.unmet.source.s[2] <- pow(sigma.unmet.other,-2)
sigma.unmet.other ~ dunif(0.01,2)
sigma.unmet.dhs ~ dunif(0.01,2)

sigma.unmetworld ~ dunif(0,5)
tau.unmetworld <- pow(sigma.unmetworld,-2) 

#end main model code, see validation and change priors below
", file = file.path(mcmc.meta$general$output.dir, "/model.txt"), fill = TRUE)

#--------------------------------------------------------------
# Validation
if (!is.null(mcmc.meta$validation.list)){
# validation: same as above, but now with test data (and save the output!)
cat("
for (k in 1:n.test.unmet){
 pred.logitratio.yunmet.j[getj.test.unmet.k[k]] ~ dnorm(
  logitratio.yunmet.hat.j[getj.test.unmet.k[k]], 
  tau.unmet.source.s[source.ind.unmet.j[getj.test.unmet.k[k]]])
  q.unmet.j[getj.test.unmet.k[k]] <- q.ij[3,getj.test.unmet.k[k]]
} 
", file = file.path(mcmc.meta$general$output.dir, "/model.txt"), fill = TRUE, append = T)

if (mcmc.meta$validation.list$exclude.unmet.only){

cat("
for (k in 1:n.test.unmet){
  # get total too
  pred.ratios.trad.modern.jn[getj.test.unmet.k[k],1:2] ~ dmnorm(
    mu.jn[getj.test.unmet.k[k], ], T.j[getj.test.unmet.k[k],,])
  pred.logratio.ytrad.j[getj.test.unmet.k[k]] <- pred.ratios.trad.modern.jn[getj.test.unmet.k[k],1]
  pred.logratio.ymodern.j[getj.test.unmet.k[k]] <- pred.ratios.trad.modern.jn[getj.test.unmet.k[k],2]
} 
", file = file.path(mcmc.meta$general$output.dir, "/model.txt"), fill = TRUE, append = T)

} else {
#if (!validation.list$exclude.unmet.only){
cat("
for (k in 1:n.test.breakdown){
	pred.ratios.trad.modern.jn[getj.test.k[k],1:2] ~ dmnorm(
    mu.jn[getj.test.k[k], ], T.j[getj.test.k[k],,])
  pred.logratio.ytrad.j[getj.test.k[k]] <- pred.ratios.trad.modern.jn[getj.test.k[k],1]
  pred.logratio.ymodern.j[getj.test.k[k]] <- pred.ratios.trad.modern.jn[getj.test.k[k],2]
  q.trad.j[getj.test.k[k]] <- q.ij[1,getj.test.k[k]]
  q.modern.j[getj.test.k[k]] <- q.ij[2,getj.test.k[k]]
} 
", file = file.path(mcmc.meta$general$output.dir, "/model.txt"), fill = TRUE, append = T)

} # end if-!exclude.unmet.only loop
} # end if-validation loop

#------------------------------------------------------------------------------------------
# Priors on country variances kappa
if (mcmc.meta$general$change.priors.to.zerolower){
cat("
sigma.lpc ~ dunif(0,5)
sigma.lrc ~ dunif(0,5)
sigma.wc ~ dunif(0,2)
sigma.Rwc ~ dunif(0,2)
sigma.Tc ~ dunif(0,30)
sigma.RTc ~ dunif(0,30)
sigma.earlierTc ~ dunif(0,70)
sigma.unmetc ~ dunif(0,5) 
tau.lrc <- pow(sigma.lrc,-2)
tau.lpc <- pow(sigma.lpc,-2)
tau.wc <- pow(sigma.wc, -2)
tau.Rwc <- pow(sigma.Rwc, -2)
tau.Tc <- pow(sigma.Tc, -2)
tau.RTc <- pow(sigma.RTc, -2)
tau.earlierTc <- pow(sigma.earlierTc, -2)
tau.unmetc <- pow(sigma.unmetc,-2) 

", file = file.path(mcmc.meta$general$output.dir, "/model.txt"), fill = TRUE, append = T)

} else {
cat("
tau.unmetc ~ dgamma(0.5,halfsigma2.unmetc0)
tau.earlierTc ~ dgamma(halfnu0_rich,halfnu0_rich_sigma2.earlierTc0)
tau.Tc ~ dgamma(halfnu0_poor,halfnu0_poor_sigma2.Tc0)
#tau.earlierTc ~ dgamma(halfnu0,halfnu0sigma2.earlierTc0)
#tau.Tc ~ dgamma(halfnu0,halfnu0sigma2.Tc0)
tau.lpc  ~ dgamma(halfnu0,halfnu0sigma2.lpc0)
tau.lrc  ~ dgamma(halfnu0,halfnu0sigma2.lrc0)
tau.wc ~ dgamma(halfnu0,halfnu0sigma2.wc0)
tau.Rwc  ~ dgamma(halfnu0,halfnu0sigma2.Rwc0)
tau.RTc  ~ dgamma(halfnu0,halfnu0sigma2.RTc0)

sigma.unmetc <- 1/sqrt(tau.unmetc) 
sigma.earlierTc <- 1/sqrt(tau.earlierTc)
sigma.lpc <- 1/sqrt(tau.lpc)
sigma.lrc <- 1/sqrt(tau.lrc)
sigma.wc <- 1/sqrt(tau.wc)
sigma.Rwc <- 1/sqrt(tau.Rwc)
sigma.Tc <- 1/sqrt(tau.Tc)
sigma.RTc <- 1/sqrt(tau.RTc)
 ", file = file.path(mcmc.meta$general$output.dir, "/model.txt"), fill = TRUE, append = T)

} # end change priors
#--------------------------------------------------------------

#--------------------------------------------------------------
if (mcmc.meta$include.AR){
cat("
for (c in 1:C){
  for (i in 1:N.unique.c[c]){
    p.ci[c, i] <- 1/(1+exp(-( logit(pstar.ci[c,i]) + eps.ci[c,i])))
    R.ci[c, i] <- 1/(1+exp(-( logit(Rstar.ci[c,i]) + eta.ci[c,i])))
  }
}
for (j in 1:J){ 
  logitZ.j[j] <- logitZstar.j[j] + theta.ci[getc.j[j], geti.j[j]]
}
  ", file = file.path(mcmc.meta$general$output.dir, "/model.txt"), fill = TRUE, append = T)
} else {
cat("
for (c in 1:C){
  for (i in 1:N.unique.c[c]){
    p.ci[c, i] <- pstar.ci[c,i]
    R.ci[c, i] <- Rstar.ci[c,i]
  }
}   
for (j in 1:J){ 
  logitZ.j[j] <- logitZstar.j[j] 
}
  ", file = file.path(mcmc.meta$general$output.dir, "/model.txt"), fill = TRUE, append = T)

} # end extra AR part

#--------------------------------------------------------------
cat("} # end model", file = file.path(mcmc.meta$general$output.dir, "/model.txt"), fill = TRUE, append = T)
# THE END!
#--------------------------------------------------------------

#if (){
# cat("
#", file = file.path(mcmc.meta$general$output.dir, "/model.txt"), fill = TRUE, append = T)
#} else {
# cat("
# ", file = file.path(mcmc.meta$general$output.dir, "/model.txt"), fill = TRUE, append = T)
#} # end extra part
 ##value<< NULL
} # end function
#----------------------------------------------------------------------   
# The End!
