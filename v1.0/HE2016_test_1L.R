#' @title Detecting and quantifying natural selection at two linked loci from time series data of allele frequencies
#' @author Zhangyi He, Xiaoyang Dai, Mark Beaumont and Feng Yu

#' version 1.0

# set the directory
setwd("~/Dropbox/Jeffery He/iResearch/Publications/2016/HE2019-WFD-2L-MMSE-PMMH-Genetics")

source("./Code/Code v2.1/HE2016_rfun_1L.R")

#install.packages("RColorBrewer")
library("RColorBrewer")

#install.packages("viridis")
library("viridis")

#install.packages("ggplot2")  
library("ggplot2")

#install.packages("plot3D")
library("plot3D")

################################################################################

#' Simulate the mutant allele frequency trajectory according to the one-locus Wright-Fisher model with selection
#' Parameter setting
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param pop_siz the number of the diploid individuals in the population
#' @param int_frq the initial mutant allele frequency of the population
#' @param int_gen the first generation of the simulated mutant allele frequency trajectory
#' @param lst_gen the last generation of the simulated mutant allele frequency trajectory

sel_cof <- 5e-03
dom_par <- 5e-01
pop_siz <- 5e+03
int_frq <- 2e-01
int_gen <- 0
lst_gen <- 500

frq_pth <- cmpsimulateOLWFMS(sel_cof, dom_par, pop_siz, int_frq, int_gen, lst_gen)

k <- int_gen:lst_gen
plot(k, frq_pth, type = 'l', lwd = 1.5, 
     xlab = "Generation", ylab = "Allele frequency", 
     main = "A frequency trajectory of the mutant allele generated with the Wright-Fisher model")

########################################

#' Simulate the mutant allele frequency trajectory according to the one-locus Wright-Fisher diffusion with selection using the Euler-Maruyama method
#' Parameter settings
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param pop_siz the number of the diploid individuals in the population
#' @param int_frq the initial mutant allele frequency of the population
#' @param int_gen the first generation of the simulated mutant allele frequency trajectory
#' @param lst_gen the last generation of the simulated mutant allele frequency trajectory
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param data_augmentation = TRUE/FALSE (return the simulated sample trajectory with data augmentation or not)

sel_cof <- 5e-03
dom_par <- 5e-01
pop_siz <- 5e+03
int_frq <- 2e-01
int_gen <- 0
lst_gen <- 500
ptn_num <- 5e+00

frq_pth <- cmpsimulateOLWFDS(sel_cof, dom_par, pop_siz, int_frq, int_gen, lst_gen, ptn_num, data_augmentation = TRUE)

t <- (int_gen:(int_gen + (lst_gen - int_gen) * ptn_num)) / 2 / pop_siz
plot(t, frq_pth, type = 'l', lwd = 1.5, 
     xlab = "Generation", ylab = "Allele frequency", 
     main = "A frequency trajectory of the mutant allele generated with the Wright-Fisher diffusion")

########################################

#' Compare the simulation generated with the Wright-Fisher model and the Wright-Fisher diffusion
sel_cof <- 5e-03
dom_par <- 5e-01
pop_siz <- 5e+03
int_frq <- 2e-01
int_gen <- 0
lst_gen <- 500
ptn_num <- 5e+00
sim_num <- 1e+06

sim_frq_WFM <- numeric(sim_num)
sim_frq_WFD <- numeric(sim_num)
for (i in 1:sim_num) {
  print(i)
  sim_frq_WFM[i] <- cmpsimulateOLWFMS(sel_cof, dom_par, pop_siz, int_frq, int_gen, lst_gen)[(lst_gen - int_gen) + 1]
  sim_frq_WFD[i] <- cmpsimulateOLWFDS(sel_cof, dom_par, pop_siz, int_frq, int_gen, lst_gen, ptn_num, data_augmentation = FALSE)[(lst_gen - int_gen) + 1]
}

save(sel_cof, sel_cof, pop_siz, int_frq, int_gen, lst_gen, ptn_num, sim_num, sim_frq_WFM, sim_frq_WFD, 
     file = "./Output/Output v2.1/TEST_1L_WFM_vs_WFD.rda")

load("./Output/Output v2.1/TEST_1L_WFM_vs_WFD.rda")

pdf(file = "./Output/Output v2.1/TEST_1L_WFM_vs_WFD.rda", width = 20, height = 10)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
hist(sim_frq_WFM, breaks = seq(min(sim_frq_WFM, sim_frq_WFD), max(sim_frq_WFM, sim_frq_WFD), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5), 
     xlim = c(min(sim_frq_WFM, sim_frq_WFD), max(sim_frq_WFM, sim_frq_WFD)), 
     xlab = "Allele frequency", 
     main = paste("Histograms of the mutant allele frequency in generation", lst_gen, "under the Wright-Fisher model and the Wright-Fisher diffusion"))
hist(sim_frq_WFD, breaks = seq(min(sim_frq_WFM, sim_frq_WFD), max(sim_frq_WFM, sim_frq_WFD), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
dev.off()

################################################################################

#' Simulate the hidden Markov model
#' Parameter setting
#' @param model = "WFM"/"WFD" (return the observations from the underlying population evolving according to the WFM or the WFD)
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param pop_siz the number of the diploid individuals in the population
#' @param int_frq the initial mutant allele frequency of the population
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the chromosomes drawn from the population at all sampling time points
#' @param ptn_num the number of the subintervals divided per generation in the Euler-Maruyama method for the WFD

#' Simulate the dataset under the Wright-Fisher model 
model <- "WFM"
sel_cof <- 5e-03
dom_par <- 5e-01
pop_siz <- 5e+03
int_frq <- 2e-01
smp_gen <- (0:10) * 50
smp_siz <- rep(50, 11)

sim_HMM_WFM <- cmpsimulateHMM(model, sel_cof, dom_par, pop_siz, int_frq, smp_gen, smp_siz)
smp_gen <- sim_HMM_WFM$smp_gen
smp_siz <- sim_HMM_WFM$smp_siz
smp_cnt <- sim_HMM_WFM$smp_cnt
pop_frq <- sim_HMM_WFM$pop_frq

k <- min(smp_gen):max(smp_gen)
smp_frq <- smp_cnt / smp_siz
plot(k, pop_frq, type = 'l', lwd = 1.5, 
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq, pop_frq), max(smp_frq, pop_frq)), 
     xlab = "Generation", ylab = "Allele frequency", 
     main = "A simulated dataset of the mutant allele generated with the Wright-Fisher model")
points(smp_gen, smp_frq, col = 'red', pch = 17, cex = 1)

####################

#' Simulate the dataset under the Wright-Fisher diffusion
model <- "WFD"
sel_cof <- 5e-03
dom_par <- 5e-01
pop_siz <- 5e+03
int_frq <- 2e-01
smp_gen <- (0:10) * 50
smp_siz <- rep(50, 11)
ptn_num <- 5e+00

sim_HMM_WFD <- cmpsimulateHMM(model, sel_cof, dom_par, pop_siz, int_frq, smp_gen, smp_siz, ptn_num)
smp_gen <- sim_HMM_WFD$smp_gen
smp_siz <- sim_HMM_WFD$smp_siz
smp_cnt <- sim_HMM_WFD$smp_cnt
pop_frq <- sim_HMM_WFD$pop_frq

k <- min(smp_gen):max(smp_gen)
smp_frq <- smp_cnt / smp_siz
plot(k, pop_frq, type = 'l', lwd = 1.5, 
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq, pop_frq), max(smp_frq, pop_frq)), 
     xlab = "Generation", ylab = "Allele frequency", 
     main = "A simulated dataset of the mutant allele generated with the Wright-Fisher diffusion")
points(smp_gen, smp_frq, col = 'red', pch = 17, cex = 1)

################################################################################

#' Generate a simulated dataset under the Wright-Fisher model
test_seed <- 7
set.seed(test_seed)

model <- "WFM"
sel_cof <- 5e-03
dom_par <- 5e-01
pop_siz <- 5e+03
int_frq <- 2e-01
smp_gen <- (0:10) * 50
smp_siz <- rep(50, 11)

SimData <- cmpsimulateHMM(model, sel_cof, dom_par, pop_siz, int_frq, smp_gen, smp_siz)
smp_gen <- SimData$smp_gen
smp_siz <- SimData$smp_siz
smp_cnt <- SimData$smp_cnt
pop_frq <- SimData$pop_frq

save(sel_cof, dom_par, pop_siz, smp_gen, smp_siz, smp_cnt, pop_frq, 
     file = "./Output/Output v2.1/TEST_1L_SimData.rda")

load("./Output/Output v2.1/TEST_1L_SimData.rda")

k <- min(smp_gen):max(smp_gen)
smp_frq <- smp_cnt / smp_siz

pdf(file = "./Output/Output v2.1/TEST_1L_SimData.pdf", width = 10, height = 10)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
plot(k, pop_frq, type = 'l', lwd = 1.5, 
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_frq, pop_frq), max(smp_frq, pop_frq)), 
     xlab = "Generation", ylab = "Allele frequency", 
     main = "A simulated dataset without missing values generated with the Wright-Fisher model")
points(smp_gen, smp_frq, col = 'red', pch = 17, cex = 1)
dev.off()

########################################

#' Run the bootstrap particle filter (BPF) with the one-locus Wright-Fisher diffusion with selection
#' Parameter settings
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param pop_siz the number of the diploid individuals in the population
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the chromosomes drawn from the population at all sampling time points
#' @param smp_cnt the count of the mutant alleles observed in the sample at all sampling time points
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter

load("./Output/Output v2.1/TEST_1L_SimData.rda")

set.seed(test_seed)

sel_cof
dom_par
pop_siz
smp_gen
smp_siz
smp_cnt
ptn_num <- 5e+00
pcl_num <- 5e+04

system.time(BPF <- cmprunBPF(sel_cof, dom_par, pop_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num))

save(sel_cof, dom_par, pop_siz, smp_gen, smp_siz, smp_cnt, pop_frq, ptn_num, pcl_num, BPF, 
     file = "./Output/Output v2.1/TEST_1L_BPF.rda")

load("./Output/Output v2.1/TEST_1L_BPF.rda")

lik <- rep(1, pcl_num)
wght <- BPF$wght
for (k in 1:length(smp_gen)) {
  lik <- lik * (cumsum(wght[, k]) / (1:pcl_num))
}
pdf(file = "./Output/Output v2.1/TEST_1L_BPF_Likelihood.pdf", width = 20, height = 10)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
plot(1:pcl_num, log(lik), type = 'l', 
     xlab = "Number of particles", ylab = "Log likelihood", 
     main = "Log likelihood estimated with the bootstrap particle filter")
dev.off()

smp_frq <- smp_cnt / smp_siz
pop_ale_frq_pre_resmp <- BPF$pop_frq_pre_resmp
pop_ale_frq_pst_resmp <- BPF$pop_frq_pst_resmp

pdf(file = "./Output/Output v2.1/TEST_1L_BPF_Particle.pdf", width = 20, height = 10)
par(mfrow = c(3, 4), oma = c(0, 0, 3, 0), mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
for (k in 1:length(smp_gen)) {
  hist(pop_ale_frq_pst_resmp[, k], breaks = seq(min(pop_ale_frq_pst_resmp[, k], pop_ale_frq_pre_resmp[, k]), max(pop_ale_frq_pst_resmp[, k], pop_ale_frq_pre_resmp[, k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5), 
       xlim = c(min(pop_ale_frq_pst_resmp[, k], pop_ale_frq_pre_resmp[, k]), max(pop_ale_frq_pst_resmp[, k], pop_ale_frq_pre_resmp[, k])), 
       xlab = "Allele frequency", 
       main = paste("Allele A1 in generation", smp_gen[k]))
  hist(pop_ale_frq_pre_resmp[, k], breaks = seq(min(pop_ale_frq_pst_resmp[, k], pop_ale_frq_pre_resmp[, k]), max(pop_ale_frq_pst_resmp[, k], pop_ale_frq_pre_resmp[, k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_frq[k], col = 'red', lty = 2, lwd = 2)
}
title("Histograms of the pre- and post-resampling particles", outer = TRUE, cex = 2)
dev.off()

########################################

#' Calculate the optimal particle number in the particle marginal Metropolis-Hastings (PMMH)
#' Parameter settings
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param pop_siz the number of the diploid individuals in the population
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the chromosomes drawn from the population at all sampling time points
#' @param smp_cnt the count of the mutant alleles observed in the sample at all sampling time points
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param gap_num the number of particles increased or decreased in the optimal particle number search

load("./Output/Output v2.1/TEST_1L_SimData.rda")

set.seed(test_seed)

sel_cof
dom_par
pop_siz
smp_gen
smp_siz
smp_cnt
ptn_num <- 5e+00
pcl_num <- 1e+03
gap_num <- 1e+02

system.time(OptNum <- calculateOptimalParticleNum(sel_cof, dom_par, pop_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, gap_num))

save(sel_cof, dom_par, pop_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, gap_num, OptNum,
     file = "./Output/Output v2.1/TEST_1L_OptNum.rda")

load("./Output/Output v2.1/TEST_1L_OptNum.rda")

opt_pcl_num <- OptNum$opt_pcl_num
log_lik_sdv <- OptNum$log_lik_sdv

pdf(file = "./Output/Output v2.1/TEST_1L_OptNum.pdf", width = 12, height = 9)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(opt_pcl_num, log_lik_sdv, type = 'b', lwd = 2, 
     xlab = "Particle number", ylab = "Log-likelihood standard deviation", 
     main = "Optimal particle number in the PMMH")
abline(h = 1.7, col = 'red', lty = 2, lwd = 2)
abline(h = 1.0, col = 'red', lty = 2, lwd = 2)
dev.off()

########################################


#' Run the particle marginal Metropolis-Hastings (PMMH)
#' Parameter settings
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param pop_siz the number of the diploid individuals in the population
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the chromosomes drawn from the population at all sampling time points
#' @param smp_cnt the count of the mutant alleles observed in the sample at all sampling time points
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param itn_num the number of the iterations carried out in the particle marginal Metropolis-Hastings

load("./Output/Output v2.1/TEST_1L_SimData.rda")

set.seed(test_seed)

sel_cof
dom_par
pop_siz
smp_gen
smp_siz
smp_cnt
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 5e+04

system.time(PMMH <- cmprunPMMH(sel_cof, dom_par, pop_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num))

save(sel_cof, dom_par, pop_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, PMMH, 
     file = "./Output/Output v2.1/TEST_1L_PMMH.rda")

load("./Output/Output v2.1/TEST_1L_PMMH.rda")

load("./Output/Output v2.1/TEST_1L_SimData.rda")

sel_cof_chn <- PMMH$sel_cof_chn
pdf(file = "./Output/Output v2.1/TEST_1L_PMMH_Traceplot.pdf", width = 20, height = 5)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
plot(1:itn_num, sel_cof_chn[1:itn_num], type = 'l', 
     xlab = "Iteration", ylab = "Selection coefficient", 
     main = "Trace plot of the selection coefficient")
abline(h = sel_cof, col = 'red', lty = 2, lwd = 2)
dev.off()

brn_num <- 1e+04
sel_cof_chn <- sel_cof_chn[brn_num:length(sel_cof_chn)]

thn_num <- 8e+00
sel_cof_chn <- sel_cof_chn[(1:round(length(sel_cof_chn) / thn_num)) * thn_num]

grd_num <- 1e+03
sel_cof_pdf <- density(sel_cof_chn, n = grd_num)
pdf(file = "./Output/Output v2.1/TEST_1L_PMMH_posterior.pdf", width = 10, height = 10)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
hist(sel_cof_chn, breaks = seq(min(sel_cof_chn), max(sel_cof_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient", 
     main = "Posterior for the selection coefficient")
lines(sel_cof_pdf, lwd = 2, col = 'black')
abline(v = sel_cof, col = 'red', lty = 2, lwd = 2)
dev.off()

########################################

#' Run the Bayesian Procedure for the inference of the selection coefficients
#' Parameter settings
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param pop_siz the number of the diploid individuals in the population
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the chromosomes drawn from the population at all sampling time points
#' @param smp_cnt the count of the mutant alleles observed in the sample at all sampling time points
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param itn_num the number of the iterations carried out in the particle marginal Metropolis-Hastings
#' @param brn_num the number of the iterations for burn-in
#' @param thn_num the number of the iterations for thinning
#' @param grd_num the number of the grids in the kernel density estimation

load("./Output/Output v2.1/TEST_1L_SimData.rda")

set.seed(test_seed)

sel_cof
dom_par
pop_siz
smp_gen
smp_siz
smp_cnt
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 5e+04
brn_num <- 1e+04
thn_num <- 8e+00
grd_num <- 1e+03

system.time(BayesianProcedure <- cmprunBayesianProcedure(sel_cof, dom_par, pop_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, brn_num, thn_num, grd_num))

save(sel_cof, dom_par, pop_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, brn_num, thn_num, grd_num, BayesianProcedure, 
     file = "./Output/Output v2.1/TEST_1L_BayesianProcedure.rda")

load("./Output/Output v2.1/TEST_1L_BayesianProcedure.rda")

load("./Output/Output v2.1/TEST_1L_SimData.rda")

sel_cof_chn <- BayesianProcedure$sel_cof_chn

sel_cof_pdf <- BayesianProcedure$sel_cof_pdf

sel_cof_map <- BayesianProcedure$sel_cof_map

sel_cof_mmse <- BayesianProcedure$sel_cof_mmse

sel_cof_hpd <- BayesianProcedure$sel_cof_hpd

pdf(file = "./Output/Output v2.1/TEST_1L_BayesianProcedure_Posterior.pdf", width = 10, height = 10)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
hist(sel_cof_chn, breaks = seq(min(sel_cof_chn), max(sel_cof_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient", 
     main = "Posterior for the selection coefficient")
lines(sel_cof_pdf, lwd = 2, col = 'black')
abline(v = sel_cof, col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_map, col = 'black', lty = 4, lwd = 2)
abline(v = sel_cof_mmse, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

################################################################################
