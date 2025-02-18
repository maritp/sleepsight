### SETTING UP ENVIRONMENT
Sys.setenv(KMP_DUPLICATE_LIB_OK=TRUE)
library('reticulate')
library('beeswarm')
library('nloptr')
library('car')
library('zoo')
library('plotrix')
library('lme4')
library('RColorBrewer')
library('lsr')
library('prettyGraphs')
library('scales')
library('matrixStats')
library('BFpack')
np = import("numpy")
basedir = '/Users/loewe/Documents/Science/Projects/SleepSight/SleepSight_Upload'

# adds se shadows to data x
se_shadows = function(x, y, se, ccol='#66666666', border = NA) {
  polygon(c(x, rev(x)), c(y + se, rev(y - se)), col = ccol, border = border)
}

# define colour schemes
sleep_cols = c(brewer.pal(n = 9, name = "YlGn")[6],brewer.pal(n = 11, name = "Spectral")[9], brewer.pal(n = 11, name = "BrBG")[10])
insight_cols = c(brewer.pal(n = 9, name = "YlGnBu")[7], brewer.pal(n = 9, name = "Blues")[4])

###########################
### LOADING DATA
###########################

data = data.frame(ACC = NULL, COH = NULL, TRIALS = NULL, ID = NULL, BLOCK = NULL, CUP = NULL, SLEEP = NULL, OG_ID = NULL)
rt_data = data.frame(RT = NULL, COH = NULL, TRIALS = NULL, ID = NULL, BLOCK = NULL, CUP = NULL, SLEEP = NULL, OG_ID = NULL)
pvt_data = data.frame(PVT = NULL, TRIALS = NULL, ID = NULL)
cacc = np$load(paste(basedir,'/Data/', 'eeg_sleep', '_accuracy.npy', sep = ''))
crt = np$load(paste(basedir,'/Data/', 'eeg_sleep', '_reaction.npy', sep = ''))
ccoh = np$load(paste(basedir,'/Data/', 'eeg_sleep', '_coherence.npy', sep = ''))
cpvt = np$load(paste(basedir,'/Data/', 'eeg_pvt.npy', sep = ''))
nids = dim(cacc)[1]
ntrials = dim(cacc)[2]
ptrials = dim(cpvt)[2]
blocksize = 100
nblocks = ntrials/50
#nblocks = ntrials/5
ccup = csleep = csleeprep = og_id = array(NA, dim=c(nids,ntrials))
for (cid in 1:nids){
  ccup[cid,] = rep(read.csv(paste(basedir,'/Data/Sleep/','cup_drop_table_incl_EEGscoring_corr_2.csv', sep = ''), header=FALSE)[cid,3],ntrials)
  csleep[cid,] = rep(read.csv(paste(basedir,'/Data/Sleep/','cup_drop_table_incl_EEGscoring_corr_2.csv', sep = ''), header=FALSE)[cid,6],ntrials)
  csleeprep[cid,] = rep(read.csv(paste(basedir,'/Data/Sleep/','cup_drop_table_incl_EEGscoring_corr_2.csv', sep = ''), header=FALSE)[cid,5],ntrials)
  og_id[cid,] = rep(read.csv(paste(basedir,'/Data/Sleep/','cup_drop_table_incl_EEGscoring_corr_2.csv', sep = ''), header=FALSE)[cid,2],ntrials)
}
cdf = data.frame(ACC = c(cacc), COH = as.factor(c(ccoh)), TRIALS = rep(1:ntrials, each = nids), ID = gl(nids, k = 1, length = nids*ntrials), CUP = as.factor(c(ccup)), SLEEP = as.factor(c(csleep)), SLEEPREP = as.factor(c(csleeprep)), OG_ID = c(og_id))
cdf$BLOCK = ceiling(cdf$TRIALS/50)
#cdf$BLOCK = ceiling(cdf$TRIALS/5)
rtcdf = data.frame(RT = c(crt), COH = as.factor(c(ccoh)), TRIALS = rep(1:ntrials, each = nids), ID = gl(nids, k = 1, length = nids*ntrials), CUP = as.factor(c(ccup)), SLEEP = as.factor(c(csleep)), SLEEPREP = as.factor(c(csleeprep)), OG_ID = c(og_id))
rtcdf$BLOCK = ceiling(cdf$TRIALS/50)
#rtcdf$BLOCK = ceiling(cdf$TRIALS/5)
pvtcdf = data.frame(PVT = c(cpvt), TRIALS = rep(1:ptrials, each = nids), ID = gl(nids, k = 1, length = nids*ptrials))
data = rbind(data, cdf)
rt_data = rbind(rt_data, rtcdf)
pvt_data = rbind(pvt_data, pvtcdf)
summary(data)
summary(rt_data)



###########################
### INSIGHT CLASSIFICATION
###########################

# define generalised logistic function
switchfun = function(fmin, fmax, steepness, inflection, time) {
  x = ((fmax-fmin)*(1/(1 + exp(-steepness*(time-inflection)))) + fmin)
  return(x)
}

# define model evaluation function
switchfun_eval = function(x, cdata = cdata) {
  time = 1:length(cdata)
  fmin = mean(cdata[1:3])
  y = switchfun(fmin, x[1], x[3], x[2], time)
  loss = sum((cdata - y)^2)
  return(loss)
}

# define bounds + starting values
clbs = c(0.35, 0, 0)
cubs = c(1, 12, 9)
x0 = c(0.5, 3, 7)

# fit generalised logistic function to data
opts = list('algorithm'='NLOPT_LN_COBYLA',"xtol_rel" = 1.0e-8, 'maxeval' = 1e4)
csteep = cfit = clevel = cswitch = cscore = matrix(NA, nids)
cparams = array(NA, dim = c(4, nids))
pred = array(NA, dim = c(12, nids))
lo_acc = array(NA, dim = c(12, nids))

accuracyblocks = tapply(data$ACC, list(data$COH, data$BLOCK, data$ID), mean)[1,5:16,]
baseblocks = tapply(data$ACC, list(data$COH, data$BLOCK, data$ID), mean)[3:5,6:7,]
baseblocks = apply(baseblocks, 3, mean)


for (cid in 1:nids) {
  cdata = accuracyblocks[,cid]
  k = nloptr(x0 = x0, eval_f = switchfun_eval, cdata = cdata, ub = cubs, lb = clbs, opts = opts)
  csteep[cid] = k$solution[3] * ((k$solution[1] - mean(cdata[1:3]))/4)
  cswitch[cid] = k$solution[2]
  cfit[cid] = k$objective*1
  cscore[cid] = csteep[cid] - cfit[cid]
  clevel[cid] = k$solution[1]
  cparams[,cid] = c(mean(cdata[1:3]), k$solution)
  lo_acc[,cid] = cdata
  pred[1:length(cdata),cid] = switchfun(mean(accuracyblocks[1:3,cid]), k$solution[1], k$solution[3], k$solution[2], 1:length(cdata))
}


# define insight subjects
swi_sleep = which(colMeans(accuracyblocks[11:12,]) > 0.85)

# subjects passing performance criterion
low_ids = which(baseblocks < 0.8)
in_ids = setdiff(1:nids, low_ids)

crit = array(NA, dim=length(nids))
for (cid in 1:nids) {
  if (cid %in% in_ids) {
    crit[cid] = 1
  } else {
    crit[cid] = 0
  }
}

# insight class
switch_class = array(NA, dim=length(nids))
for (cid in 1:nids) {
  if (cid %in% as.numeric(swi_sleep)) {
    switch_class[cid] = 1
  } else {
    switch_class[cid] = 0
  }
}


###########################
### FIG 1C
###########################

# plot accuracy of all subjects
accuracyblocks = tapply(data$ACC, list(data$COH, data$BLOCK, data$ID), mean)[,,in_ids]
cmeans = apply(accuracyblocks, c(2, 1), mean, na.rm = TRUE)[5:18,]
csds = apply(accuracyblocks, c(2, 1), std.error, na.rm = TRUE)[5:18,]
par(mar=c(5,5,4,2))
matplot(cmeans, lwd =4, lty = 1, type = 'o', pch = 16, col = c(brewer.pal(n = 9, name = "YlGnBu")[7:2],'black'), bty = 'n', cex.axis = 1.75, cex.lab = 2, ylim = c(0.4, 1), yaxt='n', xlim = c(1,15),ylab = '% correct', xlab = 'Block')
axis(side = 2,at = c(0.4,0.5,0.6,0.7,0.8,0.9,1), labels = c(40,50,60,70,80,90,100), cex.axis = 1.75)
abline(h = 0.5, col = '#666666', lwd = 1, lty = 2)
for (ccoh in 1:5) {
  se_shadows(1:dim(cmeans)[1], cmeans[,ccoh],csds[,ccoh], alpha(c(brewer.pal(n = 9, name = "YlGnBu")[7:2][ccoh],'black'), 0.4))
  legend("bottomright", legend = c("5", "23", "41", "59", "76", "100"), col = c(brewer.pal(n = 9, name = "YlGnBu")[7:2],'black'), pch = (19), bty = "n", cex = 1, title = expression(atop("Motion", "coh. (%)")))
  abline(v = 3, col = '#666666', lwd = 1, lty = 2)
  abline(v = 12, col = '#666666', lwd = 1, lty = 2)
}

# plot RT of all subjects 
rtblocks = tapply(rt_data$RT, list(rt_data$COH, rt_data$BLOCK, rt_data$ID), mean)[,,]
cmeans = apply(rtblocks, c(2, 1), mean, na.rm = TRUE)[5:18,]
csds = apply(rtblocks, c(2, 1), std.error, na.rm = TRUE)[5:18,]
par(mar=c(5,5,4,2))
matplot(cmeans, lwd =4, lty = 1, type = 'o', pch = 16, col = c(brewer.pal(n = 9, name = "YlGnBu")[7:2],'black'), bty = 'n', cex.axis = 1.75, cex.lab = 2, xlim = c(1,15),ylab = 'RT (ms)', xlab = 'Block')
for (ccoh in 1:5) {
  se_shadows(1:dim(cmeans)[1], cmeans[,ccoh],csds[,ccoh], alpha(c(brewer.pal(n = 9, name = "YlGnBu")[7:2][ccoh],'black'), 0.4))
  legend("topright", legend = c("5", "23", "41", "59", "76", "100"), col = c(brewer.pal(n = 9, name = "YlGnBu")[7:2],'black'), pch = (19), bty = "n", cex = 1.2, title = 'coh. (%)')
  abline(v = 3, col = '#666666', lwd = 1, lty = 2)
  abline(v = 12, col = '#666666', lwd = 1, lty = 2)
}


###########################
### FIG 1D
###########################

# accuracy by motion coherence in last block before colour correlation
accuracyblocks = tapply(data$ACC, list(data$COH, data$BLOCK, data$ID), mean)[,5:7,in_ids]
cmeans = apply(accuracyblocks, c(2, 1), mean, na.rm = TRUE)[,]
cmeans = colMeans(cmeans)
csds = apply(accuracyblocks, c(2, 1), std.error, na.rm = TRUE)[,]
csds = colMeans(csds)
par(mar=c(5,5,4,2))
plot(cmeans, pch = 19, cex = 3, col = brewer.pal(n = 9, name = "YlGnBu")[7:2], bty = 'n', cex.axis = 1.75, cex.lab = 2, ylim = c(0.4,1), xlim = c(0.7,5), xaxt = 'n', yaxt = 'n', ylab = '% correct', xlab = 'Motion coherence (%)')
points(cmeans, pch = 21, cex = 0.5, bg = 'black')
axis(side = 1,at = c(1,2,3,4,5), labels = c("5", "23", "41", "59", "76"), cex.axis = 1.75)
axis(side = 2,at = c(0.4,0.5,0.6,0.7,0.8,0.9,1), labels = c(40,50,60,70,80,90,100), cex.axis = 1.75)
for (i in 1:length(cmeans)){
  arrows(i, (cmeans[i]-csds[i]),i,(cmeans[i]+csds[i]), lwd = 2, length=0, angle=90, code=3, col = 'black')
}

accblocks = tapply(data$ACC, list(data$TRIALS, data$ID), mean)[201:300, in_ids]
cohblocks = data[data$TRIALS %in% c(201:300) & data$ID %in% in_ids, "COH"]
pred_df <- data.frame(
  acc = as.vector(accblocks),
  coh = as.vector(as.numeric(cohblocks)),
  id = as.vector(in_ids)
)
summary(pred_df)

glcontrol <- glmerControl(optCtrl = list(algorithm = "NLOPT_LN_NELDERMEAD", xtol_abs = 1e-12, ftol_abs = 1e-12), calc.derivs = FALSE)
acc_mod <- glmer(formula = acc ~ coh + (1 + coh | id), data = pred_df, family = binomial, control = glcontrol, na.action = na.exclude)
summary(acc_mod)
Anova(acc_mod)


###########################
### LOADING EEG SLEEP DATA
###########################

# sleep scored EEG data
sleep_all = data$SLEEP[1:nids]

# sleep self reports
rep_all = data$SLEEPREP[1:nids]

# cup drops
cup_all = data$CUP[1:nids]

# data exclusions
firstacc = colMeans(tapply(data$ACC, list(data$COH, data$BLOCK, data$ID), mean)[1,9:10,])
exc = array(1, dim=length(1:nids))
# performance criterion not met
perf = which(crit == 0)
# insight before nap
presw = which(cswitch < 4 & firstacc > 0.85)
exc[perf] = 2
exc[presw] = 3

# data frame for modelling
sleepdata = data.frame(ID = NULL, CUP = NULL, SLEEP = NULL, SLEEPREP = NULL, ACC = NULL, FACC = NULL, INSIGHT = NULL, SWITCH = NULL, CRIT = NULL, EX = NULL, OG_ID = NULL)
lastacc = colMeans(tapply(data$ACC, list(data$COH, data$BLOCK, data$ID), mean)[1,15:16,])
cdf = data.frame(ID = gl(nids, k = 1, length = nids), CUP = as.factor(cup_all), SLEEP = as.factor(sleep_all), SLEEPREP = as.factor(rep_all), ACC = lastacc, FACC = firstacc, INSIGHT = as.factor(switch_class), SWITCH = cswitch, CRIT = as.factor(crit), EX = as.factor(exc), OG_ID = og_id[1:nids])
sleepdata = rbind(sleepdata, cdf)
summary(sleepdata)

# data frame with exclusions
DATA = sleepdata[sleepdata$CRIT == 1 & sleepdata$SLEEP != 'NaN' & sleepdata$EX == 1 & sleepdata$OG_ID != 46 & sleepdata$OG_ID != 93, ]
levels(DATA$SLEEP) = c('W', 'N1', 'N2', 'NaN')

# define insight subjects
switchids = DATA$ID[DATA$INSIGHT == 1]
nswitchids = DATA$ID[DATA$INSIGHT == 0]
preswitchids = sleepdata$ID[sleepdata$INSIGHT == 1 & sleepdata$SWITCH < 4 & sleepdata$CRIT == 1 & sleepdata$FACC > 0.85]

# sleep and insight subjects
w_all = sleepdata$ID[sleepdata$SLEEP == 0]
w_in = DATA$ID[DATA$SLEEP == 'W' & DATA$INSIGHT == 1]
w_no = DATA$ID[DATA$SLEEP == 'W' & DATA$INSIGHT == 0]
w = c(w_in, w_no)
n1_all = sleepdata$ID[sleepdata$SLEEP == 1]
n1_in = DATA$ID[DATA$SLEEP == 'N1' & DATA$INSIGHT == 1]
n1_no = DATA$ID[DATA$SLEEP == 'N1' & DATA$INSIGHT == 0]
n1 = c(n1_in, n1_no)
n2_all = sleepdata$ID[sleepdata$SLEEP == 2]
n2_in = DATA$ID[DATA$SLEEP == 'N2' & DATA$INSIGHT == 1]
n2_no = DATA$ID[DATA$SLEEP == 'N2' & DATA$INSIGHT == 0]
n2 = c(n2_in, n2_no)

# sleep report and insight subjects
rep_w_all = sleepdata$ID[sleepdata$SLEEPREP == 0]
rep_w_in = DATA$ID[DATA$SLEEPREP == 0 & DATA$INSIGHT == 1]
rep_w_no = DATA$ID[DATA$SLEEPREP == 0 & DATA$INSIGHT == 0]
rep_w = c(rep_w_in, rep_w_no)
rep_n1_all = sleepdata$ID[sleepdata$SLEEPREP == 1]
rep_n1_in = DATA$ID[DATA$SLEEPREP == 1 & DATA$INSIGHT == 1]
rep_n1_no = DATA$ID[DATA$SLEEPREP == 1 & DATA$INSIGHT == 0]
rep_n1 = c(rep_n1_in, rep_n1_no)
rep_n2_all = sleepdata$ID[sleepdata$SLEEPREP == 2]
rep_n2_in = DATA$ID[DATA$SLEEPREP == 2 & DATA$INSIGHT == 1]
rep_n2_no = DATA$ID[DATA$SLEEPREP == 2 & DATA$INSIGHT == 0]
rep_n2 = c(rep_n2_in, rep_n2_no)



###########################
### FIG 2B
###########################
par(mar=c(5,5,2,2), lwd = 2)
in_perc = c('Base' = 49.5,'Wake'= length(w_in)*100/length(w), 'N1'= length(n1_in)*100/length(n1), 'N2'= length(n2_in)*100/length(n2))
barplot(in_perc, col = c('grey',sleep_cols), border = 'black', ylab = '% insight', width = 0.2, space = 0.7, xlab = '', main = '', ylim = c(0,100), cex.lab = 2, cex.axis = 1.5, cex.names = 1.5, lwd = 1.5)
abline(h = 49.5, lty = 2, lwd = 1.5)

###########################
### FIG 5A
###########################
par(mar=c(5,5,4,2), lwd = 2)
rep_in_perc = c('Base' = 49.5,'Wake'= length(rep_w_in)*100/length(rep_w), 'N1'= length(rep_n1_in)*100/length(rep_n1), 'N2'= length(rep_n2_in)*100/length(rep_n2))
barplot(rep_in_perc, col = c('grey',sleep_cols), border = 'black', ylab = '% insight', width= 0.2, space = 0.7, xlab = '', main = '', ylim = c(0,100), cex.lab = 2, cex.axis = 1.5, cex.names = 1.5, lwd = 1.5)
abline(h = 49.5, lty = 2, lwd = 1.5)

# Fisher test
# all sleep groups
sleep_effect = matrix(c(length(w_in), length(n1_in), length(n2_in), length(w_no), length(n1_no), length(n2_no)), nrow = 3, dimnames = list(sleep = c("W", "N1", "N2"), insight = c("yes", "no")))
fisher.test(sleep_effect)
# wake + n2
sleep_effect = matrix(c(length(w_in), length(w_no), length(n2_in), length(n2_no)), nrow = 2, dimnames = list(insight = c("yes", "no"), sleep = c("W", "N2")))
fisher.test(sleep_effect)
# wake + n1
sleep_effect = matrix(c(length(w_in), length(w_no), length(n1_in), length(n1_no)), nrow = 2, dimnames = list(insight = c("yes", "no"), sleep = c("W", "N1")))
fisher.test(sleep_effect)
# n1 + n2
sleep_effect = matrix(c(length(n1_in), length(n1_no), length(n2_in), length(n2_no)), nrow = 2, dimnames = list(insight = c("yes", "no"), sleep = c("N1", "N2")))
fisher.test(sleep_effect)
# wake + sleep
sleep_effect = matrix(c(length(w_in), length(w_no), length(c(n2_in,n1_in)), length(c(n2_no,n1_no))), nrow = 2, dimnames = list(insight = c("yes", "no"), sleep = c("W", "N2")))
fisher.test(sleep_effect)

# insight sleep + baseline
sleep_effect = matrix(c(length(switchids), length(nswitchids), 49, 50), nrow = 2, dimnames = list(insight = c("yes", "no"), sleep = c("sleep", "base")))
fisher.test(sleep_effect)

# pre-nap insight + baseline
pre_effect = matrix(c(length(preswitchids), length(switchids)+length(preswitchids), 15, 49), nrow = 2, dimnames = list(insight = c("yes", "no"), sleep = c("sleep", "base")))
fisher.test(pre_effect)


# Insight classification glm
glm0 = glm(INSIGHT ~ 1, data = DATA, family = 'binomial')
glm1 = glm(INSIGHT ~ 1 + SLEEP, data = DATA, family = 'binomial')

# Switch point prediction glm
glm0 = glm(SWITCH ~ 1, data = DATA, family = 'gaussian')
glm1 = glm(SWITCH ~ 1 + SLEEP, data = DATA, family = 'gaussian')

# Bayes Factor
# H1
BF(glm1, hypothesis = 'SLEEPN1 > Intercept; SLEEPN1 = Intercept')
# H2
BF(glm1, hypothesis = 'SLEEPN1 > SLEEPN2; SLEEPN1 = SLEEPN2')
# H3
BF(glm1, hypothesis = 'SLEEPN2 > Intercept; SLEEPN2 = Intercept')


### SINGLE TRIAL
# adjust nblocks in line 44 to ntrials/5
# number of blocks before switch
cstart = 6
# number of blocks after switch
cend = 1
nblocks <- 12*10
cblocks <- array(data = NA, dim = c(nids, nblocks))
switch_align  <- array(data = NA, dim = c(nids, nblocks*2))

accuracyblocks = tapply(data$ACC, list(data$COH, data$BLOCK, data$ID), mean)[,,]
accuracyse = tapply(data$ACC, list(data$COH, data$BLOCK, data$ID), std.error)[,,]
accuracyblocks = accuracyblocks[1,1:(12*10),1:nids]
accuracyse = accuracyse[1,1:(12*10),1:nids]
for (cid in 1:nids){
  cblocks[cid,] = (1:nblocks) + (nblocks - cswitch[cid]*10)
  switch_align[cid, as.integer(min(cblocks[cid,])):as.integer(max(cblocks[cid,]))] = 
    accuracyblocks[,cid]
}

###########################
### FIG 1H
###########################
ex_ins = 58
par(mar=c(5,5,4,2))
plot(accuracyblocks[81:86,ex_ins], type = 'o', lty = 1, lwd = 6, pch = 16, cex = 2, col = insight_cols[1], bty = 'n', cex.axis = 1.75, cex.lab = 2, ylim = c(0,1), ylab = 'Accuracy', xlab = 'Trial', xaxt = 'n', yaxt = 'n')
axis(1, at = c(1,2,3,4,5,6), labels = c(-2,-1,0,1,2,3), cex.lab = 1.75, cex.axis = 1.75)
axis(2, at = c(0,1), labels = c(0,1), cex.lab = 1.75, cex.axis = 1.75)
abline(v=3.5, lty = 2, lwd = 2)


###########################
### INSIGHT-ALIGNED ACCURACY
###########################

# number of blocks before switch
cstart = 3
# number of blocks after switch
cend = 2
nblocks <- 12
cblocks <- array(data = NA, dim = c(nids, nblocks))
switch_align  <- rt_align <- array(data = NA, dim = c(nids, nblocks*2))

## align data to individual switch points
accuracyblocks = tapply(data$ACC, list(data$COH, data$BLOCK, data$ID), mean)[,,]
accuracyse = tapply(data$ACC, list(data$COH, data$BLOCK, data$ID), std.error)[,,]
rtblocks = tapply(rt_data$RT, list(rt_data$COH, rt_data$BLOCK, rt_data$ID), mean)[,,]
rtse = tapply(rt_data$RT, list(rt_data$COH, rt_data$BLOCK, rt_data$ID), std.error)[,,]

accuracyblocks = accuracyblocks[1,5:16,]
accuracyse = accuracyse[1,5:16,]
rtblocks = rtblocks[1,5:16,]
rtse = rtse[1,5:16,]

for (cid in 1:nids){
  cblocks[cid,] = (1:nblocks) + (nblocks - cswitch[cid])
  switch_align[cid, as.integer(min(cblocks[cid,])):as.integer(max(cblocks[cid,]))] <- 
    accuracyblocks[ ,cid]
  rt_align[cid, as.integer(min(cblocks[cid,])):as.integer(max(cblocks[cid,]))] <- 
    rtblocks[ ,cid]
}

# mean insight/no insight switch aligned accuracy
switchmean = apply(switch_align[switchids, ], 2, mean, na.rm = TRUE)[(nblocks-cstart):(nblocks+cend)]
switchse = apply(switch_align[switchids, ], 2, std.error, na.rm = TRUE)[(nblocks-cstart):(nblocks+cend)]
nswitchmean = apply(switch_align[nswitchids, ], 2, mean, na.rm = TRUE)[(nblocks-cstart):(nblocks+cend)]
nswitchse = apply(switch_align[nswitchids, ], 2, std.error, na.rm = TRUE)[(nblocks-cstart):(nblocks+cend)]

###########################
### FIG 1G
###########################
par(mar=c(5,5,4,2))
matplot(-cstart:cend, switchmean, type = 'o', lty = 1, pch = c(16), cex = 2, lwd = 4, bty = 'n', xlab = 'Block', ylab = '% correct', ylim = c(0.4, 1), col = insight_cols[1], cex.lab = 2, cex.axis = 1.75, yaxt = 'n', xaxt = 'n')
se_shadows(-cstart:cend, switchmean, switchse, alpha(insight_cols[1],0.2))
lines(-cstart:cend, nswitchmean, type = 'o', lty = 1, pch = c(16), cex = 2, lwd = 4, bty = 'n', xlab = 'Block', ylab = '% correct', ylim = c(0.4, 1), col = insight_cols[2], cex.lab = 2, cex.axis = 1.75, yaxt = 'n', xaxt = 'n')
se_shadows(-cstart:cend, nswitchmean, nswitchse, alpha(insight_cols[2],0.2))
axis(2, at = c(0.4,0.5,0.6,0.7,0.8,0.9,1), labels = c(40,50,60,70,80,90,100), cex.lab = 1.75, cex.axis = 1.75)
axis(1, at = c(-3,-2,-1,0,1,2), labels = c(-2,-1,0,1,2,3), cex.lab = 1.75, cex.axis = 1.75)
abline(v = -1, col = 1, lty = 3, lwd = 2)
legend("bottomright", legend = c("Insight", "No insight"), col = insight_cols, lty = 1,lwd = 4, bty = "n", cex = 1.3)


# mean insight/no insight switch aligned RT
rtswitchmean = apply(rt_align[switchids, ], 2, mean, na.rm = TRUE)[(nblocks-cstart):(nblocks+cend)]
rtswitchse = apply(rt_align[switchids, ], 2, std.error, na.rm = TRUE)[(nblocks-cstart):(nblocks+cend)]
rtnswitchmean = apply(rt_align[nswitchids, ], 2, mean, na.rm = TRUE)[(nblocks-cstart):(nblocks+cend)]
rtnswitchse = apply(rt_align[nswitchids, ], 2, std.error, na.rm = TRUE)[(nblocks-cstart):(nblocks+cend)]

# plot mean switch aligned reaction times
par(mar=c(5,5,4,2))
matplot(-cstart:cend, rtswitchmean, type = 'o', lty = 1, pch = c(16), cex = 2, lwd = 4, bty = 'n', xlab = 'Block', ylab = 'RT(ms)', ylim = c(500,1300), col = insight_cols[1], cex.lab = 2, cex.axis = 1.75, xaxt = 'n')
se_shadows(-cstart:cend, rtswitchmean, rtswitchse, alpha(insight_cols[1],0.2))
lines(-cstart:cend, rtnswitchmean, type = 'o', lty = 1, pch = c(16), cex = 2, lwd = 4, bty = 'n', xlab = 'Block', ylab = 'RT(ms)', col = insight_cols[2], cex.lab = 2, cex.axis = 1.75, xaxt = 'n')
se_shadows(-cstart:cend, rtnswitchmean, rtnswitchse, alpha(insight_cols[2],0.2))
axis(1, at = c(-3,-2,-1,0,1,2), labels = c(-2,-1,0,1,2,3), cex.lab = 1.75, cex.axis = 1.75)
abline(v = -1, col = 1, lty = 3, lwd = 2)
legend("bottomleft", legend = c("insight", "no insight"), col = insight_cols, lty = 1,lwd = 4, bty = "n", cex = 1.4)

###########################
### FIG 1F
###########################
shist = hist(cswitch[switchids],breaks = seq(2, 7.5, by=0.5), plot=FALSE)
plot(shist, border = insight_cols[1], col = insight_cols[1], bty = 'n', ylab = 'Frequency',xlab = 'Block', main = '', cex.lab = 2, cex.axis = 1.75, xlim = c(3,12))
abline(v = 3.5, lty=2, lwd = 1.25)
abline(v = 4, lty=2, lwd = 1.25)


# mean sleep/no sleep switch aligned accuracy (EEG)
sleepmean = apply(switch_align[n2, ], 2, mean, na.rm = TRUE)[(nblocks-cstart):(nblocks+cend)]
sleepse = apply(switch_align[n2, ], 2, std.error, na.rm = TRUE)[(nblocks-cstart):(nblocks+cend)]
nsleepmean = apply(switch_align[w, ], 2, mean, na.rm = TRUE)[(nblocks-cstart):(nblocks+cend)]
nsleepse = apply(switch_align[w, ], 2, std.error, na.rm = TRUE)[(nblocks-cstart):(nblocks+cend)]
bsleepmean = apply(switch_align[n1, ], 2, mean, na.rm = TRUE)[(nblocks-cstart):(nblocks+cend)]
bsleepse = apply(switch_align[n1, ], 2, std.error, na.rm = TRUE)[(nblocks-cstart):(nblocks+cend)]

# mean sleep/no sleep switch aligned accuracy (reports)
sleepmean = apply(switch_align[rep_n2, ], 2, mean, na.rm = TRUE)[(nblocks-cstart):(nblocks+cend)]
sleepse = apply(switch_align[rep_n2, ], 2, std.error, na.rm = TRUE)[(nblocks-cstart):(nblocks+cend)]
nsleepmean = apply(switch_align[rep_w, ], 2, mean, na.rm = TRUE)[(nblocks-cstart):(nblocks+cend)]
nsleepse = apply(switch_align[rep_w, ], 2, std.error, na.rm = TRUE)[(nblocks-cstart):(nblocks+cend)]
bsleepmean = apply(switch_align[rep_n1, ], 2, mean, na.rm = TRUE)[(nblocks-cstart):(nblocks+cend)]
bsleepse = apply(switch_align[rep_n1, ], 2, std.error, na.rm = TRUE)[(nblocks-cstart):(nblocks+cend)]

# plot mean sleep aligned accuracy
par(mar=c(5,5,4,2))
matplot(-cstart:cend, nsleepmean, type = 'o', lty = 1, pch = c(16), cex = 2, lwd = 4, bty = 'n', xlab = 'Block', ylab = '% correct', ylim = c(0.4, 1), col = sleep_cols[1], cex.lab = 2, cex.axis = 1.75, yaxt = 'n', xaxt = 'n')
se_shadows(-cstart:cend, nsleepmean, nsleepse, alpha(sleep_cols[1],0.2))
lines(-cstart:cend, bsleepmean, type = 'o', lty = 1, pch = c(16), cex = 2, lwd = 4, bty = 'n', xlab = 'Block', ylab = '% correct', ylim = c(0.4, 1), col = sleep_cols[2], cex.lab = 2, cex.axis = 1.75, yaxt = 'n', xaxt = 'n')
se_shadows(-cstart:cend, bsleepmean, bsleepse, alpha(sleep_cols[2],0.2))
lines(-cstart:cend, sleepmean, type = 'o', lty = 1, pch = c(16), cex = 2, lwd = 4, bty = 'n', xlab = 'Block', ylab = '% correct', ylim = c(0.4, 1), col = sleep_cols[3], cex.lab = 2, cex.axis = 1.75, yaxt = 'n', xaxt = 'n')
se_shadows(-cstart:cend, sleepmean, sleepse, alpha(sleep_cols[3],0.2))
axis(2, at = c(0.4,0.5,0.6,0.7,0.8,0.9,1), labels = c(40,50,60,70,80,90,100), cex.lab = 1.75, cex.axis = 1.75)
axis(1, at = c(-3,-2,-1,0,1,2), labels = c(-2,-1,0,1,2,3), cex.lab = 1.75, cex.axis = 1.75)
abline(v = -1, col = 1, lty = 3, lwd = 2)
legend("bottomright", legend = c('Wake',"N1", 'N2'), col = sleep_cols, lty = 1, lwd = 4, bty = "n", cex = 1.1)


###########################
### FIG 4A
###########################
accuracyblocks = tapply(data$ACC, list(data$COH, data$BLOCK, data$ID), mean)[1,5:18,]
sleepmean = apply(accuracyblocks[, n2], 1, mean, na.rm = TRUE)
sleepse = apply(accuracyblocks[, n2], 1, std.error, na.rm = TRUE)
bsleepmean = apply(accuracyblocks[, n1], 1, mean, na.rm = TRUE)
bsleepse = apply(accuracyblocks[, n1], 1, std.error, na.rm = TRUE)
nsleepmean = apply(accuracyblocks[, w], 1, mean, na.rm = TRUE)
nsleepse = apply(accuracyblocks[, w], 1, std.error, na.rm = TRUE)

par(mar=c(5,5,4,2))
matplot(nsleepmean, type = 'o', lty = 1, pch = c(16), cex = 2, lwd = 4, bty = 'n', xlab = 'Block', ylab = '% correct', ylim = c(0.4, 1), col = sleep_cols[1], cex.lab = 2, cex.axis = 1.75, yaxt = 'n')
se_shadows(1:14,nsleepmean, nsleepse, alpha(sleep_cols[1],0.2))
lines(bsleepmean, type = 'o', lty = 1, pch = c(16), cex = 2, lwd = 4, bty = 'n', xlab = 'Block', ylab = '% correct', ylim = c(0.4, 1), col = sleep_cols[2], cex.lab = 2, cex.axis = 1.75)
se_shadows(1:14,bsleepmean, bsleepse, alpha(sleep_cols[2],0.2))
lines(sleepmean, type = 'o', lty = 1, pch = c(16), cex = 2, lwd = 4, bty = 'n', xlab = 'Block', ylab = '% correct', ylim = c(0.4, 1), col = sleep_cols[3], cex.lab = 2, cex.axis = 1.75)
se_shadows(1:14,sleepmean, sleepse, alpha(sleep_cols[3],0.2))
axis(2, at = c(0.4,0.5,0.6,0.7,0.8,0.9,1), labels = c(40,50,60,70,80,90,100), cex.lab = 1.75, cex.axis = 1.75)
legend("bottomright", legend = c('Wake',"N1", 'N2'), col = sleep_cols, lty = 1, lwd = 4, bty = "n", cex = 1.1)
abline(v = 3.5, col = '#666666', lwd = 1.75, lty = 2)
abline(v = 4, col = '#666666', lwd = 1.75, lty = 2)
abline(v = 12, col = '#666666', lwd = 1.75, lty = 2)


###########################
### FIG 4B
###########################
accuracyblocks = tapply(rt_data$RT, list(rt_data$COH, rt_data$BLOCK, rt_data$ID), mean)[1,5:18,]
sleepmean = apply(accuracyblocks[, n2], 1, mean, na.rm = TRUE)
sleepse = apply(accuracyblocks[, n2], 1, std.error, na.rm = TRUE)
bsleepmean = apply(accuracyblocks[, n1], 1, mean, na.rm = TRUE)
bsleepse = apply(accuracyblocks[, n1], 1, std.error, na.rm = TRUE)
nsleepmean = apply(accuracyblocks[, w], 1, mean, na.rm = TRUE)
nsleepse = apply(accuracyblocks[, w], 1, std.error, na.rm = TRUE)

par(mar=c(5,5,4,2))
matplot(nsleepmean, type = 'o', lty = 1, pch = c(16), cex = 2, lwd = 4, bty = 'n', xlab = 'Block', ylab = 'RT (ms)', ylim = c(400, 1200), col = sleep_cols[1], cex.lab = 2, cex.axis = 1.75)
se_shadows(1:14,nsleepmean, nsleepse, alpha(sleep_cols[1],0.2))
lines(bsleepmean, type = 'o', lty = 1, pch = c(16), cex = 2, lwd = 4, bty = 'n', xlab = 'Block', ylab = 'RT (ms)', ylim = c(400, 1200), col = sleep_cols[2], cex.lab = 2, cex.axis = 1.75)
se_shadows(1:14,bsleepmean, bsleepse, alpha(sleep_cols[2],0.2))
lines(sleepmean, type = 'o', lty = 1, pch = c(16), cex = 2, lwd = 4, bty = 'n', xlab = 'Block', ylab = 'RT (ms)', ylim = c(400, 1200), col = sleep_cols[3], cex.lab = 2, cex.axis = 1.75)
se_shadows(1:14,sleepmean, sleepse, alpha(sleep_cols[3],0.2))
#axis(2, at = c(0.4,0.5,0.6,0.7,0.8,0.9,1), labels = c(40,50,60,70,80,90,100), cex.lab = 1.75, cex.axis = 1.75)
legend("topright", legend = c('Wake',"N1", 'N2'), col = sleep_cols, lty = 1, lwd = 4, bty = "n", cex = 1.1)
abline(v = 3.5, col = '#666666', lwd = 1.75, lty = 2)
abline(v = 4, col = '#666666', lwd = 1.75, lty = 2)
abline(v = 12, col = '#666666', lwd = 1.75, lty = 2)


###########################
### FIG 2C
###########################
DATA$SLEEP = droplevels(DATA$SLEEP)
DATA$SLEEPREP = droplevels(DATA$SLEEPREP)
DATA$X = ((DATA$SWITCH + 5)/2)-3.5
switch_dists = tapply(DATA$X[DATA$INSIGHT == 1], DATA$SLEEP[DATA$INSIGHT == 1], mean)[c('W', 'N1', 'N2')]
dists_se = tapply(DATA$X[DATA$INSIGHT == 1], DATA$SLEEP[DATA$INSIGHT == 1], std.error)[c('W', 'N1', 'N2')]
par(mar=c(5,5,4,2), lwd = 3) 
k = barplot(switch_dists, col = 'white', border = sleep_cols, xlab = 'Block', xlim = c(0,4), xaxt = 'n', ylab = '', main = '', width = 0.8, space = 0.5, cex.lab = 1.75, cex.names = 1.75, lwd = 1.5, horiz = TRUE, names = c('Wake', 'N1', 'N2'), las = 2)
beeswarm(DATA$X[DATA$INSIGHT == 1] ~ DATA$SLEEP[DATA$INSIGHT == 1], at = k,  pch = 21, col = sleep_cols, method = 'swarm', horiz = TRUE, xlim = c(0,8), axes = FALSE, add = TRUE)
#beeswarm(DATA$X[DATA$INSIGHT == 1] ~ DATA$SLEEPREP[DATA$INSIGHT == 1], at = k,  pch = 21, col = sleep_cols, method = 'swarm', horiz = TRUE, xlim = c(0,8), axes = FALSE, add = TRUE)
axis(1, at = c(0,0.5,1.5,2.5,3.5), labels = c('',4,5,6,7), cex.axis = 1.75, lwd = 1.5) 
for (i in 1:length(switch_dists)){
  arrows(switch_dists[i] - dists_se[i], k[i], switch_dists[i] + dists_se[i], k[i], lwd = 2.5, length=0.1, angle=90, code=3, col = c(c(brewer.pal(n = 9, name = "YlGn")[6],brewer.pal(n = 11, name = "Spectral")[9], brewer.pal(n = 11, name = "BrBG")[10]))[i])
}
abline(v = 0.05, col = '#666666', lwd = 1.75, lty = 2)
abline(v = 0.5, col = '#666666', lwd = 1.75, lty = 2)


# KS test switch point distribution
ks.test(cswitch[w], cswitch[n2])
ks.test(cswitch[w], cswitch[n1])
ks.test(cswitch[n1], cswitch[n2])

ks.test(cswitch[w_in], cswitch[n2_in])
ks.test(cswitch[w_in], cswitch[n1_in])
ks.test(cswitch[n1_in], cswitch[n2_in])


# mean sleep/no sleep insight accuracy
n2_sleepmean = apply(switch_align[n2_in, ], 2, mean, na.rm = TRUE)[(nblocks-cstart):(nblocks+cend)]
n2_sleepse = apply(switch_align[n2_in, ], 2, std.error, na.rm = TRUE)[(nblocks-cstart):(nblocks+cend)]
n1_sleepmean = apply(switch_align[n1_in, ], 2, mean, na.rm = TRUE)[(nblocks-cstart):(nblocks+cend)]
n1_sleepse = apply(switch_align[n1_in, ], 2, std.error, na.rm = TRUE)[(nblocks-cstart):(nblocks+cend)]
w_sleepmean = apply(switch_align[w_in, ], 2, mean, na.rm = TRUE)[(nblocks-cstart):(nblocks+cend)]
w_sleepse = apply(switch_align[w_in, ], 2, std.error, na.rm = TRUE)[(nblocks-cstart):(nblocks+cend)]

# mean sleep/no sleep report insight accuracy
n2_sleepmean = apply(switch_align[rep_n2_in, ], 2, mean, na.rm = TRUE)[(nblocks-cstart):(nblocks+cend)]
n2_sleepse = apply(switch_align[rep_n2_in, ], 2, std.error, na.rm = TRUE)[(nblocks-cstart):(nblocks+cend)]
n1_sleepmean = apply(switch_align[rep_n1_in, ], 2, mean, na.rm = TRUE)[(nblocks-cstart):(nblocks+cend)]
n1_sleepse = apply(switch_align[rep_n1_in, ], 2, std.error, na.rm = TRUE)[(nblocks-cstart):(nblocks+cend)]
w_sleepmean = apply(switch_align[rep_w_in, ], 2, mean, na.rm = TRUE)[(nblocks-cstart):(nblocks+cend)]
w_sleepse = apply(switch_align[rep_w_in, ], 2, std.error, na.rm = TRUE)[(nblocks-cstart):(nblocks+cend)]

# mean sleep/no sleep no insight accuracy
n2_sleepmean = apply(switch_align[n2_no, ], 2, mean, na.rm = TRUE)[(nblocks-cstart):(nblocks+cend)]
n2_sleepse = apply(switch_align[n2_no, ], 2, std.error, na.rm = TRUE)[(nblocks-cstart):(nblocks+cend)]
n1_sleepmean = apply(switch_align[n1_no, ], 2, mean, na.rm = TRUE)[(nblocks-cstart):(nblocks+cend)]
n1_sleepse = apply(switch_align[n1_no, ], 2, std.error, na.rm = TRUE)[(nblocks-cstart):(nblocks+cend)]
w_sleepmean = apply(switch_align[w_no, ], 2, mean, na.rm = TRUE)[(nblocks-cstart):(nblocks+cend)]
w_sleepse = apply(switch_align[w_no, ], 2, std.error, na.rm = TRUE)[(nblocks-cstart):(nblocks+cend)]


###########################
### FIG 2D
###########################
par(mar=c(5,5,4,2))
matplot(-cstart:cend, n2_sleepmean, type = 'o', lty = 1, pch = c(16), cex = 2, lwd = 4, bty = 'n', xlab = 'Block', ylab = '% correct', ylim = c(0.4, 1), col = sleep_cols[3], cex.lab = 2, cex.axis = 1.75, yaxt = 'n', xaxt = 'n')
se_shadows(-cstart:cend, n2_sleepmean, n2_sleepse, alpha(sleep_cols[3],0.2))
lines(-cstart:cend, n1_sleepmean, type = 'o', lty = 1, pch = c(16), cex = 2, lwd = 4, bty = 'n', xlab = 'Block', ylab = '% correct', ylim = c(0.4, 1), col = sleep_cols[2], cex.lab = 2, cex.axis = 1.75, yaxt = 'n', xaxt = 'n')
se_shadows(-cstart:cend, n1_sleepmean, n1_sleepse, alpha(sleep_cols[2],0.2))
lines(-cstart:cend, w_sleepmean, type = 'o', lty = 1, pch = c(16), cex = 2, lwd = 4, bty = 'n', xlab = 'Block', ylab = '% correct', ylim = c(0.4, 1), col = sleep_cols[1], cex.lab = 2, cex.axis = 1.75, yaxt = 'n', xaxt = 'n')
se_shadows(-cstart:cend, w_sleepmean, w_sleepse, alpha(sleep_cols[1],0.2))
axis(2, at = c(0.4,0.5,0.6,0.7,0.8,0.9,1), labels = c(40,50,60,70,80,90,100), cex.lab = 1.75, cex.axis = 1.75)
axis(1, at = c(-3,-2,-1,0,1,2), labels = c(-2,-1,0,1,2,3), cex.lab = 1.75, cex.axis = 1.75)
abline(v = -1, col = 1, lty = 3, lwd = 2)
legend("bottomright", legend = c("Wake", 'N1', 'N2'), col = sleep_cols, lty = 1, lwd = 4, bty = "n", cex = 1.1)


# mean sleep insight switch aligned RT
rt_n2_switchmean = apply(rt_align[n2_in, ], 2, mean, na.rm = TRUE)[(nblocks-cstart):(nblocks+cend)]
rt_n2_switchse = apply(rt_align[n2_in, ], 2, std.error, na.rm = TRUE)[(nblocks-cstart):(nblocks+cend)]
rt_n1_switchmean = apply(rt_align[n1_in, ], 2, mean, na.rm = TRUE)[(nblocks-cstart):(nblocks+cend)]
rt_n1_switchse = apply(rt_align[n1_in, ], 2, std.error, na.rm = TRUE)[(nblocks-cstart):(nblocks+cend)]
rt_w_switchmean = apply(rt_align[w_in, ], 2, mean, na.rm = TRUE)[(nblocks-cstart):(nblocks+cend)]
rt_w_switchse = apply(rt_align[w_in, ], 2, std.error, na.rm = TRUE)[(nblocks-cstart):(nblocks+cend)]

# mean sleep report insight switch aligned RT
rt_n2_switchmean = apply(rt_align[rep_n2, ], 2, mean, na.rm = TRUE)[(nblocks-cstart):(nblocks+cend)]
rt_n2_switchse = apply(rt_align[rep_n2, ], 2, std.error, na.rm = TRUE)[(nblocks-cstart):(nblocks+cend)]
rt_n1_switchmean = apply(rt_align[rep_n1, ], 2, mean, na.rm = TRUE)[(nblocks-cstart):(nblocks+cend)]
rt_n1_switchse = apply(rt_align[rep_n1, ], 2, std.error, na.rm = TRUE)[(nblocks-cstart):(nblocks+cend)]
rt_w_switchmean = apply(rt_align[rep_w, ], 2, mean, na.rm = TRUE)[(nblocks-cstart):(nblocks+cend)]
rt_w_switchse = apply(rt_align[rep_w, ], 2, std.error, na.rm = TRUE)[(nblocks-cstart):(nblocks+cend)]

# mean sleep no insight switch aligned RT
rt_n2_switchmean = apply(rt_align[n2_no, ], 2, mean, na.rm = TRUE)[(nblocks-cstart):(nblocks+cend)]
rt_n2_switchse = apply(rt_align[n2_no, ], 2, std.error, na.rm = TRUE)[(nblocks-cstart):(nblocks+cend)]
rt_n1_switchmean = apply(rt_align[n1_no, ], 2, mean, na.rm = TRUE)[(nblocks-cstart):(nblocks+cend)]
rt_n1_switchse = apply(rt_align[n1_no, ], 2, std.error, na.rm = TRUE)[(nblocks-cstart):(nblocks+cend)]
rt_w_switchmean = apply(rt_align[w_no, ], 2, mean, na.rm = TRUE)[(nblocks-cstart):(nblocks+cend)]
rt_w_switchse = apply(rt_align[w_no, ], 2, std.error, na.rm = TRUE)[(nblocks-cstart):(nblocks+cend)]


###########################
### FIG 2E
###########################
par(mar=c(5,5,4,2))
matplot(-cstart:cend, rt_n2_switchmean, type = 'o', lty = 1, pch = c(16), cex = 2, lwd = 4, bty = 'n', xlab = 'Block', ylab = 'RT(ms)', ylim = c(500,1300), col = sleep_cols[3], cex.lab = 2, cex.axis = 1.75, xaxt = 'n')
se_shadows(-cstart:cend, rt_n2_switchmean, rt_n2_switchse, alpha(sleep_cols[3],0.2))
lines(-cstart:cend, rt_n1_switchmean, type = 'o', lty = 1, pch = c(16), cex = 2, lwd = 4, bty = 'n', xlab = 'Block', ylab = 'RT(ms)', col = sleep_cols[2], cex.lab = 2, cex.axis = 1.75, xaxt = 'n')
se_shadows(-cstart:cend, rt_n1_switchmean, rt_n1_switchse, alpha(sleep_cols[2],0.2))
lines(-cstart:cend, rt_w_switchmean, type = 'o', lty = 1, pch = c(16), cex = 2, lwd = 4, bty = 'n', xlab = 'Block', ylab = 'RT(ms)', col = sleep_cols[1], cex.lab = 2, cex.axis = 1.75, xaxt = 'n')
se_shadows(-cstart:cend, rt_w_switchmean, rt_w_switchse, alpha(sleep_cols[1],0.2))
axis(1, at = c(-3,-2,-1,0,1,2), labels = c(-2,-1,0,1,2,3), cex.lab = 1.75, cex.axis = 1.75)
abline(v = -1, col = 1, lty = 3, lwd = 2)
legend("topright", legend = c('Wake', 'N1','N2'), col = sleep_cols, lty = 1, lwd = 4, bty = "n", cex = 1.1)


# consciousness reports (about colour rule)
consc = np$load(paste(basedir, '/Data/eeg_conscious.npy', sep = ''))
consc_time = np$load(paste(basedir, '/Data/eeg_conscious_time.npy', sep = ''))
consc_yes = which(consc == 0)
consc_no = which(consc == 1)


###########################
### FIG 1E
###########################
ratio = c(length(switchids),length(nswitchids))
pie(ratio, labels = '', col = insight_cols, border = 'black')

