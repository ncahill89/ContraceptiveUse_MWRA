WriteCountryModelForTargets <- function(# Writes BUGS model
  ### Writes BUGS model to output.dir (stored in mcmc.meta)
  mcmc.meta ##<< Object
) {
  filename.model <- paste0(ifelse(mcmc.meta$general$do.SS.run.first.pass, "model_pre.txt", "model.txt")) # change JR, 20140414
#------ MODEL ----------
cat("
#--------------------------------------------------------------
# Model for contraceptive use (country-specific for targets)
# Leontine Alkema, 2011 & Jin Rou New, 2015
#--------------------------------------------------------------

model{

tau.tot.st <- tau.tot*(1-pow(rho.tot,2))
rho.tot <- rho.tot0
tau.tot <- pow(sigma.tot, -2)
sigma.tot <- sigma.tot0

tau.rat.st <- tau.rat*(1-pow(rho.rat,2))
rho.rat <- rho.rat0
tau.rat <- pow(sigma.rat, -2)
sigma.rat <- sigma.rat0

tau.unmet.st <- tau.unmet*(1-pow(rho.unmet,2))
rho.unmet <- rho.unmet0
tau.unmet <- pow(sigma.ar.unmet, -2)
sigma.ar.unmet <- sigma.ar.unmet0

for (c in 1:C) {
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

# dummy
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
  w.reg[reg] <- w.reg0[reg]
	T.reg[reg] <- T.reg0[reg]
	Rw.reg[reg] <- Rw.reg0[reg]
	RT.reg[reg] <- RT.reg0[reg]
}

tau.unmetworld <- pow(sigma.unmetworld, -2)
sigma.unmetworld <- sigma.unmetworld0
tau.wsubreg <- pow(sigma.wsubreg, -2)
tau.Rwsubreg <- pow(sigma.Rwsubreg, -2)
sigma.wsubreg <- sigma.wsubreg0
sigma.Rwsubreg <- sigma.Rwsubreg0

tau.RTsubreg <- pow(sigma.RTsubreg, -2)
tau.Tsubreg <- pow(sigma.Tsubreg, -2)
sigma.RTsubreg <- sigma.RTsubreg0
sigma.Tsubreg <- sigma.Tsubreg0

Tearlier <- Tearlier0

lp.world <- lp.world0
lr.world <- lr.world0
#--------------------------------------------------------------
# Likelihood
# gett.j gives (year - year.start+1) for obs j
# getc.j gives gives c for obs j

for (j in 1:J) { 
  trad.j[j] <- p.ci[getc.j[j], geti.j[j]] * (1-R.ci[getc.j[j], geti.j[j]])
  modern.j[j] <- p.ci[getc.j[j], geti.j[j]] * (R.ci[getc.j[j], geti.j[j]])
  unmet.j[j] <- (1-p.ci[getc.j[j], geti.j[j]])/(1+exp(-logitZ.j[j]))
  none.j[j] <- 1 - trad.j[j] - modern.j[j] - unmet.j[j]
  logitZstar.j[j] <- (
        unmet.intercept.c[getc.j[j]]
        + a.unmet 
        + b.unmet * (p.ci[getc.j[j], geti.j[j]] - pmid.for.unmet) 
        + c.unmet * pow(p.ci[getc.j[j], geti.j[j]] - pmid.for.unmet, 2)) 
  # logitZ.j defined in AR-loop, just adding theta  
  logratio.y.jn[j, 1:3] ~ dmnorm(mu.logratios.jn[j, 1:3], Tau.logratios[1:3, 1:3])
  mu.logratios.jn[j, 1] <- log(max(0.01, trad.j[j])/none.j[j])
  mu.logratios.jn[j, 2] <- log(max(0.01, modern.j[j])/none.j[j])
  mu.logratios.jn[j, 3] <- log(max(0.01, unmet.j[j])/none.j[j])
} # end loop J observations

# dummy
# mu.logratios.jn[J+1, 1] <- 0

a.unmet <- a.unmet0
b.unmet <- b.unmet0
c.unmet <- c.unmet0

#end main model code, change priors below
", file = file.path(mcmc.meta$general$output.dir, filename.model), fill = TRUE)
#--------------------------------------------------------------
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
", file = file.path(mcmc.meta$general$output.dir, filename.model), fill = TRUE, append = T)
} else {
cat("
tau.unmetc ~ dgamma(0.5,halfsigma2.unmetc0)
tau.earlierTc ~ dgamma(halfnu0_rich,halfnu0_rich_sigma2.earlierTc0)
tau.Tc ~ dgamma(halfnu0_poor,halfnu0_poor_sigma2.Tc0)
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
 ", file = file.path(mcmc.meta$general$output.dir, filename.model), fill = TRUE, append = T)
} # end change priors
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
  ", file = file.path(mcmc.meta$general$output.dir, filename.model), fill = TRUE, append = T)
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
  ", file = file.path(mcmc.meta$general$output.dir, filename.model), fill = TRUE, append = T)

} # end extra AR part

#--------------------------------------------------------------
cat("} # end model", file = file.path(mcmc.meta$general$output.dir, filename.model), fill = TRUE, append = T)
# THE END!
#--------------------------------------------------------------

#if (){
# cat("
#", file = file.path(mcmc.meta$general$output.dir, filename.model), fill = TRUE, append = T)
#} else {
# cat("
# ", file = file.path(mcmc.meta$general$output.dir, filename.model), fill = TRUE, append = T)
#} # end extra part
 ##value<< NULL
} # end function
#----------------------------------------------------------------------   
# The End!
