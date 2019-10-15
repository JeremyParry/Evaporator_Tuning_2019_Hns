#Written By Jeremy Parry

#Note
#this package was developed in order to perform various
#process control functions on an evaporator with one specific data structure
#this package has not been tested on other data structures

#################################################
#This function calculates tunings established in Chen et. al. (2006) when passes a process gain, time constant, a, b and deadtime.
#PID calculations are supported.
#note that the robustness factor (lambda must be set)

#Chen, P., Zhang, W., & Zhu, L. (2006). Design and tuning method of PID controller for a class of inverse response processes. Paper presented at the 2006 American Control Conference.
Chen <- function(Kp,
                 tau1,
                 tau2,
                 tauz,
                 taud,
                 lambda=4,
                 cname="tune"){
  #Summary
  h=taud-(2*tauz)
  
  Kc = (1/Kp)*((tau1+tau2)/(lambda+(2*tauz)+h))
  tauI=tau1+tau2
  tauD=(tau1*tau2)/(tau1+tau2)
  N=((lambda-h)*tauz)/(lambda+(2*tauz)+h)
  tauF=tauD/N
  
  output = data.frame(c(Kc,tauI,tauD,tauF))
  row.names(output) = c("Kc","taui","taud","tauf")
  colnames(output) = cname
  return(output)
}
#################################################
#this function returns a seris of response criteria including rise time, overshoot, controller overshoot
# and various ITAE criteria
Criteria = function(tm,
                    Pv,
                    Mv,
                    setpoint,
                    PVnm = "PV",
                    MVnm = "MV",
                    pltrisetime = TRUE,
                    pltovershoot = TRUE,
                    pltcovershoot = TRUE,
                    signf = 3,
                    normalise = FALSE,
                    legoff = (max(Pv)-min(Pv))/2.5,
                    rtlabjust=c(-0.03,0),
                    oslabjust=c(0.6,0),
                    ITAEmxtm=max(tm)){
  
  par(mar = c(3, 3, 3, 3))
  par(xpd = TRUE)
  #Plot PV
  plot(tm,Pv,type = 'l',xlab = "Time (min)",ylab = " ",ylim = c(min(Pv),max(Pv)+(setpoint-Pv[1])*0.07))
  mtext("Time (min)", side=1, line=2)
  mtext(PVnm, side=2, line=2)
  legend(tm[1],max(Pv)+legoff,c(PVnm,MVnm),lty = c(1,2),ncol=2)
  par(xpd = FALSE)
  abline(h=setpoint, lty=3)
  
  #Rise Time
  RiseTime= signif(min(tm[Pv>=setpoint-(setpoint-Pv[1])*0.05]),signf)
  if (pltrisetime==TRUE){
    abline(v=RiseTime,lty=3)
    text(RiseTime,min(Pv),paste("Rise Time = ",RiseTime,"min"),rtlabjust)
  }
  
  #Overshoot
  a = max(Pv)[1]-setpoint
  c = setpoint-Pv[1]
  Overshoot= signif(100*(a/c),signf)
  if (pltovershoot==TRUE){
    arrows(tm[Pv==max(Pv)][1],setpoint,tm[Pv==max(Pv)][1],setpoint+a,0.05,code = 3)
    text(tm[Pv==max(Pv)][1],max(Pv)[1]+(setpoint-Pv[1])*0.01,paste("Overshoot = ",Overshoot,"%"),oslabjust)
  }
  
  #Plot MV
  par(new = TRUE)
  plot(tm, Mv, type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "",lty=2,ylim = c(min(Mv),max(Mv)+((max(Mv)-min(Mv))*0.05)))
  axis(side=4, at = pretty(range(Mv)))
  mtext(MVnm, side=4, line=2)
  
  #Controller Overshoot
  b=(max(Mv)-Mv[tm==max(tm)])
  d=(Mv[tm==max(tm)]-Mv[1])
  COvershoot = signif(100*(b/d),signf)
  if (pltcovershoot==TRUE){
    text(min(tm),max(Mv),paste("Controller overshoot =",COvershoot,"%"),c(0,-0.2))
    arrows(tm[Mv==max(Mv)][1],max(Mv),tm[Mv==max(Mv)][1],Mv[length(Mv)],0.05,code = 3)
  }
  #peaks
  peak=max(Pv)[1]
  Cpeak=max(Mv)[1]
  
  #IAE
  ts = tm[2]-tm[1]
  e = abs(setpoint-Pv)
  e=e[tm<=ITAEmxtm]
  tm=tm[tm<=ITAEmxtm]
  ne = e/e[1]
  
  IAE = signif(ts*sum(e),signf)
  ISE = signif(ts*sum(e^2),signf)
  ITAE = signif(ts*sum(e*tm),signf)
  ITSE = signif(ts*sum((e^2)*tm),signf)
  
  NIAE = signif(ts*sum(ne),signf)
  NISE = signif(ts*sum(ne^2),signf)
  NITAE = signif(ts*sum(ne*tm),signf)
  NITSE = signif(ts*sum((ne^2)*tm),signf)
  
  if(normalise == TRUE){
    result = data.frame(RiseTime,Overshoot,peak,COvershoot,Cpeak,NIAE,NISE,NITAE,NITSE)
  }else{
    result = data.frame(RiseTime,Overshoot,peak,COvershoot,Cpeak,IAE,ISE,ITAE,ITSE)
  }
  par(mar = c(5.1,4.1,4.1,2.1))
  return(result)
}
##################################################
#This function is required for peakmax and returns the mode of a vector of numbers V
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
#################################################
#This function calculates ITAE tunings when passes a process gain, time constant, and deadtime.
#P, PI, or PID calculations are supported
ITAE_Tune <- function(Kp,
                      tau,
                      taud,
                      Mode = "PI") {

  if (Mode == "P") {
    numtune = 1
    Value = matrix(c(0.49,-1.084), ncol =  numtune)
  }else if (Mode == "PI") {
    numtune = 2
    Value = matrix(c(0.859,-0.977,(1/0.674),0.680), ncol =  numtune)
  }else if (Mode == "PID") {
    numtune = 3
    Value = matrix(c(1.357,-0.947,(1/0.842),0.738,0.381,0.995), ncol =  numtune)
  }else{
    stop("Mode must be P, PI, or PID")
  }
  
  x = c((1/Kp), tau, tau)
  
  ITAE = rep(NA, numtune)
  for (i in 1:numtune) {
    ITAE[i] = (x[i]*Value[1,i])*((taud/tau)^Value[2,i])
  }
  output = data.frame(ITAE)
  row.names(output) = c("Kc", "taui", "taud")
  return(output)
}
#################################################
#This function calculates Tau tunings when passes a process gain, time constant, and deadtime.
#PI and PID calculations are supported

#Åström, K. J., & Hägglund, T. (2004). Revisiting the Ziegler-Nichols step response method for PID control. Journal of Process Control, 14(6), 635-650.
Kappa_Tau <- function(Kp,
                      tau,
                      taud,
                      Mode = "PI") {
  #Explain this function here
  
  a = Kp * (taud / tau)
  Tau = taud / (taud + tau)
  
  if (Mode == "PI") {
    Value = matrix(
      c(
        0.29,-2.7,3.7,0.78,-4.1,5.7,8.9,-6.6,
        3.0,8.9,-6.6,3.0,0.79,-1.4,2.4,0.79,
        -1.4, 2.4,0.81,0.73,1.9,0.44,0.78,-0.45
      ),
      nrow = 6
    )
  } else if (Mode == "PID") {
    Value = matrix(
      c(
        3.8,-8.4,7.3,8.6,-9.6,9.8,5.2,-2.5,
        -1.4,3.2,-1.5,-0.93,0.077,5,-4.8,0.076,
        3.4,-1.1,0.4,0.18,2.8,0.22,6.5,0.051
      ),
      nrow = 6
    )
  } else{stop('Mode must be PI or PID')}
  
  x = c((1 / a), taud, tau, 1)
  
  M1.4 = rep(NA, 4)
  M2.0 = rep(NA, 4)
  for (i in 1:4) {
    M1.4[i] = x[i] * (Value[1, i] * exp((Value[2, i] * Tau)+(Value[3, i] * (Tau^2))))
    M2.0[i] = x[i] * (Value[4, i] * exp((Value[5, i] * Tau)+(Value[6, i] * (Tau^2))))
  }
  output = data.frame(M1.4, M2.0)
  row.names(output) = c("Kc", "taui", "taud", "beta")
  return(output)
}
#################################################
#This function returns a numerically fitted process reaction curve alongside a plot showing resulting fit
# currently supported process reaction curves are First order, First orderplus dead time, Second Order, 
# Second Order plus dead time, second order plus right hand plane zero, second order plus right hand plane 
# zero plus dead time.
Lagfit <- function(t,
                   y,
                   mdty = "FOPDT",#model type
                   M = 1,#step change magnitude
                   y0 = y[t == min(t)], #steady state level of Process variable before step change
                   Kp = 1,#inital guess for process gain
                   tau = 1,#inital guess for first order time consatints
                   tau1 = 1,#inital guess for second order time consatint1
                   tau2 = 0,#inital guess for second order time consatint2
                   taud = 1,#inital guess for deadtime
                   A = 0.5,#inital guess for a
                   B = 0.5,#inital guess for b
                   mv = NA,#the index of a point to force the curve through
                   mvw = length(t)*5,#the number of points worth of weight to give this point
                   plt = TRUE,#plot a graph?
                   marg = c(3, 3, 1, 1),#plot margins
                   ttl = paste("Fitted", mdty, "curve"),#Plot heading
                   xax = "Time (s)",#x axis label
                   yax = "%",#y axis label
                   r2 = TRUE,#display rsquared in plot?
                   pos=ifelse(M>0,"topleft","bottomleft"),#where to put the legend
                   Datcol="grey", 
                   Lncol="black",
                   ...) {
  #this function fits a response delay model to experimental data.
  library(minpack.lm)
  
  we=rep(1,length(t))
  if(is.numeric(mv)==TRUE){
    we[mv]=mvw
  }
  
  #models
  if (mdty == "FO") { #first order model
    mdl = nls(y ~ y0 + (kp * M) * (1 - exp(-t / Tau)),
              start = list(kp = Kp, Tau = tau),weights=we)
  } else if (mdty == "FOPDT") { #first order plus dead time model
    mdl = nlsLM(
      y ~ ifelse(t >= Taud, y0 + (kp * M) * (1 - exp(-(t - Taud) / Tau)), y0),
      start = list(kp = Kp, Tau = tau, Taud = taud),
      lower = c(-Inf, 0.1, 0.1),weights=we
    )
  } else if (mdty == "SO") { #second order model
    mdl = nlsLM(
      y ~ y0 + (kp * M) * (1 - (((Tau1 * exp(-t / Tau1)) - (Tau2 * exp(-t / Tau2))) / (Tau1 - Tau2))),
      start = list(
        kp = Kp,
        Tau1 = tau1,
        Tau2 = tau2
      ),
      lower = c(-Inf, 0.1, 0.1),weights=we
    )
  } else if (mdty == "SOPDT") { #second order plus dead time model
    mdl = nlsLM(
      y ~ ifelse(t >= Taud, y0 + (kp * M) * (1 - (((Tau1 * exp(-(t - Taud) / Tau1)) - (Tau2 * exp(-(t - Taud) / Tau2))) / (Tau1 - Tau2)
      )), y0),
      start = list(
        kp = Kp,
        Tau1 = tau1,
        Tau2 = tau2,
        Taud = taud
      ),
      lower = c(-Inf, 0.1, 0.1, 0.1),weights=we
    )
  } else if (mdty == "SOPRHPZ") { #second order plus right hand plane zero model
    mdl = nlsLM(y ~ y0+(kp*M*(1+((((1+b)*exp(-t/Tau))-((a+b)*exp(-t/(Tau*a))))/(a-1)))),
                start = list(
                  kp = Kp, 
                  Tau = tau,
                  b=B,
                  a=A
                ),
                lower = c(-Inf, 0.1, 0.001, 0.001),
                weights=we)
  }else if (mdty == "SOPRHPZPDT") { #second order plus right hand plane zero plus dead time model
    mdl = nlsLM(y ~ ifelse(t >= Taud,y0+(kp*M*(1+((((1+b)*exp(-(t-Taud)/Tau))-((a+b)*exp(-(t-Taud)/(Tau*a))))/(a-1)))),y0),
                start = list(
                  kp = Kp, 
                  Tau = tau,
                  b=B,
                  a=A,
                  Taud=taud
                ),
                lower = c(-Inf, 0.1, 0.13, 0.01, 0.1),
                weights=we)
  } else {stop("Invalid model code")}
  
  if (r2 == TRUE) { #should r squared be displayed on plot
    rsq = signif(1 - (sum(summary(mdl)$residuals ^ 2) / sum((
      predict(mdl, t) - mean(y)
    ) ^ 2)), 3)
  }
  
  if (plt == TRUE) { #should a plot be generated
    par(mar = marg)
    plot(t, y, col = Datcol, main = ttl,...)
    lines(t, predict(mdl),col = Lncol)
    
    #legend
    legend(
      pos,
      legend = c("data", "fitted"),
      col = c(Datcol, Lncol),
      pch = c("o", "-")
    )
    mtext(text = xax,
          side = 1,
          line = 2)
    mtext(text = yax,
          side = 2,
          line = 2)
    if (r2 == TRUE) {
      text(
        x = max(t / 2),
        y = (max(y) + min(y)) / 2,
        labels = paste("Rsquared=", rsq)
      )
    }
    par(mar = c(5.1,4.1,4.1,2.1)) #return plot borders to defult
  }
  return(mdl)
}
#################################################
#This function calculates tunings established in Alfaro & Vilanova (2013) when passes a process gain, time constant, a and b.
#PID calculations are supported.

#Alfaro, V. M., & Vilanova, R. (2013). Robust tuning of 2DoF five-parameter PID controllers for inverse response controlled processes. Journal of Process Control, 23(4), 453-462. https://doi.org/10.1016/j.jprocont.2013.01.005
MoReRT_RR <- function(a,
                      b,
                      kp,
                      tau) {
  x = c((1/kp), tau, tau, tau, 1)
  
  c2.0 = matrix(
    c(
      3.808,0.8608, 0.2728,0.316,0.4172,0.4111,
      0.988,0.4044,0.1914,-0.08662,-6.217,1.599,
      0.4831,1.082,0.2333,-0.3489,0.03134,0.09845,
      -0.01138,-0.05651,4.25,-0.9303,-0.3138,
      -0.4785,-0.07367,0.09239,-0.06749,-0.0327,
      -0.03542,0.01223,-0.9784,0.2252,0.06335,
      0.1128,0.01114
    ),
    nrow = 5
  )
  colnames(c2.0) = c("c0", "c1", "c2", "c3", "c4", "c5", "c6")
  row.names(c2.0) = c("kp*", "ti*", "td*", "tf", "b*")
  M2.0 = rep(NA, 5)
  for (i in 1:5) {
    M2.0[i] = (
      c2.0[i, 1] + (c2.0[i, 2] * a) + (c2.0[i, 3] * b) + (c2.0[i, 4] * a * b) +
        (c2.0[i, 5] * b ^ 2) + (c2.0[i, 6] * a * b ^ 2) + (c2.0[i, 7] * b ^ 3)
    ) * x[i]
  }
  
  c1.6 = matrix(
    c(
      3.042, 1.079,0.3174,0.4545,0.4418,-0.07878,
      1.121,0.4253, 0.2903, -0.08582, -4.409,
      1.355, 0.4505, 1.042,0.2538, 0.1679,
      -0.26, -0.03513,-0.189,-0.07043, 2.153,
      -0.4896,-0.3367,-0.1696,-0.0276
    ),
    nrow = 5
  )
  colnames(c1.6) = c("c0", "c1", "c2", "c3", "c4")
  row.names(c1.6) = c("kp*", "ti*", "td*", "tf", "b*")
  M1.6 = rep(NA, 5)
  for (i in 1:5) {
    M1.6[i] = (c1.6[i, 1] + (c1.6[i, 2] * a) + (c1.6[i, 3] * b) + (c1.6[i, 4] * a * b) + (c1.6[i, 5] * b ^ 2)) * x[i]
  }
  M2.0[4]=M2.0[3]/M2.0[4]
  M1.6[4]=M1.6[3]/M1.6[4]
  
  output = data.frame(M1.6,M2.0)
  row.names(output) = c("Kc","taui","taud","tauf","beta")
  return(output)
}
#################################################
#this function finds peaks for the relay method
maxpeak <- function(x, y, w=1, ...) {
  require(zoo)
  n <- length(y)
  y.smooth <- loess(y ~ x, ...)$fitted
  y.max <- rollapply(zoo(y.smooth), 2*w+1, max, align="center")
  delta <- y.max - y.smooth[-c(1:w, n+1-1:w)]
  i.max <- which(delta <= 0) + w
  list(x=x[i.max], i=i.max, y.hat=y.smooth)
}
#################################################
#this function applies the relay method when passed the process data and establishes kappa tunings
# currently only supports PI tunings
Relay <-
  function(t,#time vector in seconds mins or hours
           MV, #vector of Manipulated variabl square wave
           PV, #vector of Process variable occilations
           Kp, #Process gain
           marg = c(3, 3, 1, 1), #plot margins
           xax = "Time (s)", #x axis label
           yax = "%",#y axis label
           off = 1, #peak label offsets
           ABoff= max(t)/50, #a and b bar offset
           PIDmode = 'PI', #controller mode
           Ms = '1.4', #robustness
           timeunit = 's', #the unit time is for Pu adjustment
           dispsmth = TRUE, #display the smoothed line?
           displab = TRUE, #display the peak labels?
           dispAB = TRUE, #display AB bars?
           PVw = 100, #PV width parameter for peak finding
           MVw = 90, #MV width parameter for peak finding
           PVspan = 0.05, #PV degree of smoothing for peak finding
           MVspan = 0.3, #MV degree of smoothing for peak finding
           plt = TRUEn,
           ttl = "Relay Method"
  ) {
    source('peakmax.R')
    PVpeak = maxpeak(t, PV, w = PVw, span = PVspan)
    PVtr = maxpeak(t,-PV, w = PVw, span = PVspan)
    MVpeak = maxpeak(t, MV, w = MVw, span = MVspan)
    MVtr = maxpeak(t,-MV, w = MVw, span = MVspan)
    A = mean(PV[PVpeak$i]) - mean(PV[PVtr$i])
    B = mean(MV[MVpeak$i]) - mean(MV[MVtr$i])
    Ku = (4 * B) / (pi * A)
    Pu = mean(sort(dist(t[PVpeak$i]))[1:length(t[PVpeak$i]) - 1])
    k = 1 / (Ku * Kp)
    
    if (timeunit == 's') {
      Pu = Pu / 60
    } else if (timeunit == 'min') {
      Pu = Pu
    } else if (timeunit == 'hr') {
      Pu = Pu * 60
    } else{
      stop('The time unit must be s, min or hr')
    }
    #calculating Kc taui and b
    if (PIDmode == 'PI') {
      if (Ms == '1.4') {
        a0 = c(0.053, 0.9, 1.1)
        a1 = c(2.9, -4.4, -0.0061)
        a2 = c(-2.6, 2.7, 1.8)
      } else if (Ms == '2.0') {
        a0 = c(0.13, 0.9, 0.48)
        a1 = c(1.9, -4.4, 0.4)
        a2 = c(-1.3, 2.7, -0.17)
      } else{
        stop('MS must equal 1.4 or 2.0')
      }
      Kc = (a0[1] * exp((a1[1] * k) + (a2[1] * k ^ 2))) * Ku
      ti = (a0[2] * exp((a1[2] * k) + (a2[2] * k ^ 2))) * (Pu)
      beta = (a0[3] * exp((a1[3] * k) + (a2[3] * k ^ 2)))
    } else{
      stop('only PI is currently supported')
    }
    if(plt==TRUE){
      #plotting data
      par(mar = marg)
      max(c(PV, MV))
      plot(t, PV, ylim = c(min(c(PV, MV)), max(c(PV, MV))), main = ttl,pch=20)
      points(t, MV, pch = 4)
      leg.nam = c("PV data", "MV data")
      leg.col = c("Black", "Black")
      leg.sym = c(20, 4)
      if (dispsmth == TRUE) {
        #Smoothed lines
        points(t, PVpeak$y.hat, col = 'grey')
        points(t, MVpeak$y.hat, col = 'grey', pch = 4)
        leg.nam = c(leg.nam, "PV smoothed", "MV smoothed")
        leg.col = c(leg.col, "grey", "grey")
        leg.sym = c(leg.sym, 1, 4)
      }
      if (displab == TRUE) {
        #peak labels
        text(x = t[PVpeak$i], y = PV[PVpeak$i] + off)
        text(x = t[PVtr$i], y = PV[PVtr$i] + off)
        text(x = t[MVpeak$i], y = MV[MVpeak$i] + off)
        text(x = t[MVtr$i], y = MV[MVtr$i] + off)
      }
      if (dispAB == TRUE) {
        #display magnitude of A and B on plot
        arrows(max(t), mean(MV[MVpeak$i]), max(t), mean(MV[MVtr$i]), length=0.05, angle=90, code=3)
        text(x=max(t)+ABoff,y=mean((MV[MVpeak$i])+mean(MV[MVtr$i]))/2,labels="B")
        arrows(max(t), mean(PV[PVpeak$i]), max(t), mean(PV[PVtr$i]), length=0.05, angle=90, code=3)
        text(x=max(t)+ABoff,y=mean((PV[PVpeak$i])+mean(PV[PVtr$i]))/2,labels="A")
      }
      #legend
      legend("topleft",
             legend = leg.nam,
             col = leg.col,
             pch = leg.sym)
      
      mtext(text = xax,
            side = 1,
            line = 2)
      mtext(text = yax,
            side = 2,
            line = 2)
    }
    if (PIDmode == 'PI') {
      tunings = as.data.frame(c(Kc, ti, beta,Ku,Pu,k))
      row.names(tunings)=c("Kc","taui (min)","beta","Ku","Pu (min)","kappa")
      colnames(tunings)=c("value")
    }
    return(tunings)
  }