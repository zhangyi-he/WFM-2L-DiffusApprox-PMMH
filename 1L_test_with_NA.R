#' @title An MCMC-based method for Bayesian inference of natural selection from time series DNA data across linked loci
#' @author Zhangyi He, Xiaoyang Dai, Mark Beaumont and Feng Yu

#' this version is able to handle missing values in DNA data

source("1L_rfun_with_NA.R")

#install.packages("RColorBrewer")
library("RColorBrewer")

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
#' @param int_ale_frq the initial mutant allele frequency of the population
#' @param int_gen the first generation of the simulated mutant allele frequency trajectory
#' @param lst_gen the last generation of the simulated mutant allele frequency trajectory

sel_cof <- 5e-03
dom_par <- 5e-01
pop_siz <- 5e+03
int_ale_frq <- 2e-01
int_gen <- 0
lst_gen <- 500

ale_frq_pth <- cmpsimulateOLWFMS(sel_cof, dom_par, pop_siz, int_ale_frq, int_gen, lst_gen)

k <- int_gen:lst_gen
plot(k, ale_frq_pth, type = 'l', lwd = 1.5, 
     xlab = "Generation", ylab = "Allele frequency", 
     main = "A mutant allele frequency generated with the Wright-Fisher model")

########################################

#' Simulate the mutant allele frequency trajectory according to the one-locus Wright-Fisher diffusion with selection using the Euler-Maruyama method
#' Parameter settings
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param pop_siz the number of the diploid individuals in the population
#' @param int_ale_frq the initial mutant allele frequency of the population
#' @param int_gen the first generation of the simulated mutant allele frequency trajectory
#' @param lst_gen the last generation of the simulated mutant allele frequency trajectory
#' @param ptn_num the number of the subintervals divided per generation in the Euler-Maruyama method
#' @param data_augmentation = TRUE/FALSE (return the simulated sample trajectory with data augmentation or not)

sel_cof <- 5e-03
dom_par <- 5e-01
pop_siz <- 5e+03
int_ale_frq <- 2e-01
int_gen <- 0
lst_gen <- 500
ptn_num <- 1e+01

ale_frq_pth <- cmpsimulateOLWFDS(sel_cof, dom_par, pop_siz, int_ale_frq, int_gen, lst_gen, ptn_num, data_augmentation = TRUE)

t <- (int_gen:(int_gen + (lst_gen - int_gen) * ptn_num)) / 2 / pop_siz
plot(t, ale_frq_pth, type = 'l', lwd = 1.5, 
     xlab = "Time", ylab = "Allele frequency", 
     main = "A mutant allele frequency generated with the Wright-Fisher diffusion")

########################################

#' Compare the simulation generated with the Wright-Fisher model and the Wright-Fisher diffusion
sel_cof <- 5e-03
dom_par <- 5e-01
pop_siz <- 5e+03
int_ale_frq <- 2e-01
int_gen <- 0
lst_gen <- 500
ptn_num <- 1e+01
sim_num <- 1e+06

sim_ale_frq_WFM <- numeric(sim_num)
sim_ale_frq_WFD <- numeric(sim_num)
for (i in 1:sim_num) {
  print(i)
  sim_ale_frq_WFM[i] <- cmpsimulateOLWFMS(sel_cof, dom_par, pop_siz, int_ale_frq, int_gen, lst_gen)[(lst_gen - int_gen) + 1]
  sim_ale_frq_WFD[i] <- cmpsimulateOLWFDS(sel_cof, dom_par, pop_siz, int_ale_frq, int_gen, lst_gen, ptn_num, data_augmentation = FALSE)[(lst_gen - int_gen) + 1]
}

save(sel_cof, sel_cof, pop_siz, int_ale_frq, int_gen, lst_gen, ptn_num, sim_num, sim_ale_frq_WFM, sim_ale_frq_WFD, 
     file = "TEST_1L_comparison_WFM_and_WFD.rda")

load("TEST_1L_comparison_WFM_and_WFD.rda")

pdf(file = "TEST_1L_comparison_WFM_and_WFD.pdf", width = 16, height = 9)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
hist(sim_ale_frq_WFM, breaks = seq(min(sim_ale_frq_WFM, sim_ale_frq_WFD), max(sim_ale_frq_WFM, sim_ale_frq_WFD), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5), 
     xlim = c(min(sim_ale_frq_WFM, sim_ale_frq_WFD), max(sim_ale_frq_WFM, sim_ale_frq_WFD)), 
     xlab = "Allele frequency", main = paste("Histograms of the allele frequency in generation", lst_gen, "under the W-F model and the W-F diffusion"))
hist(sim_ale_frq_WFD, breaks = seq(min(sim_ale_frq_WFM, sim_ale_frq_WFD), max(sim_ale_frq_WFM, sim_ale_frq_WFD), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
dev.off()

################################################################################

#' Simulate the hidden Markov model
#' Parameter setting
#' @param model = "WFM"/"WFD" (return the observations from the underlying population evolving according to the WFM or the WFD)
#' @param missing = TRUE/FALSE (return the observations with missing values or not)
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param pop_siz the number of the diploid individuals in the population
#' @param int_ale_frq the initial mutant allele frequency of the population
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_chr_cnt the count of the chromosomes drawn from the population at all sampling time points
#' @param ptn_num the number of the subintervals divided per generation in the Euler-Maruyama method for the WFD

#' Simulate the dataset under the Wright-Fisher model 
model <- "WFM"
missing <- TRUE
sel_cof <- 5e-03
dom_par <- 5e-01
pop_siz <- 5e+03
int_ale_frq <- 2e-01
smp_gen <- (0:5) * 100
smp_chr_cnt <- rep(100, 6)

sim_HMM_WFM <- cmpsimulateHMM(model, missing, sel_cof, dom_par, pop_siz, int_ale_frq, smp_gen, smp_chr_cnt)
smp_gen <- sim_HMM_WFM$smp_gen
smp_chr_cnt <- sim_HMM_WFM$smp_chr_cnt
smp_ale_cnt <- sim_HMM_WFM$smp_ale_cnt
mis_ale_cnt <- sim_HMM_WFM$mis_ale_cnt
smp_ale_frq <- sim_HMM_WFM$smp_ale_cnt / sim_HMM_WFM$smp_chr_cnt
pop_ale_frq <- sim_HMM_WFM$pop_ale_frq

k <- min(smp_gen):max(smp_gen)
plot(k, pop_ale_frq, type = 'l', lwd = 1.5, 
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_ale_frq, pop_ale_frq), max(smp_ale_frq, pop_ale_frq)), 
     xlab = "Generation", ylab = "Allele frequency", 
     main = "A simulated dataset of the mutant allele generated with the Wright-Fisher model")
points(smp_gen, smp_ale_frq, col = 'red', pch = 17, cex = 1)

####################

#' Simulate the dataset under the Wright-Fisher diffusion
model <- "WFD"
missing <- TRUE
sel_cof <- 5e-03
dom_par <- 5e-01
pop_siz <- 5e+03
int_ale_frq <- 2e-01
smp_gen <- (0:5) * 100
smp_chr_cnt <- rep(100, 6)
ptn_num <- 1e+01

sim_HMM_WFD <- cmpsimulateHMM(model, missing, sel_cof, dom_par, pop_siz, int_ale_frq, smp_gen, smp_chr_cnt, ptn_num)
smp_gen <- sim_HMM_WFD$smp_gen
smp_chr_cnt <- sim_HMM_WFD$smp_chr_cnt
smp_ale_cnt <- sim_HMM_WFD$smp_ale_cnt
mis_ale_cnt <- sim_HMM_WFD$mis_ale_cnt
smp_ale_frq <- sim_HMM_WFD$smp_ale_cnt / sim_HMM_WFD$smp_chr_cnt
pop_ale_frq <- sim_HMM_WFD$pop_ale_frq

k <- min(smp_gen):max(smp_gen)
plot(k, pop_ale_frq, type = 'l', lwd = 1.5, 
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_ale_frq, pop_ale_frq), max(smp_ale_frq, pop_ale_frq)), 
     xlab = "Generation", ylab = "Allele frequency", 
     main = "A simulated dataset of the mutant allele generated with the Wright-Fisher diffusion")
points(smp_gen, smp_ale_frq, col = 'red', pch = 17, cex = 1)

################################################################################

#' Generate a simulated dataset under the Wright-Fisher model
test_seed <- 21
set.seed(test_seed)

model <- "WFM"
missing <- TRUE
sel_cof <- 5e-03
dom_par <- 5e-01
pop_siz <- 5e+03
int_ale_frq <- 2e-01
smp_gen <- (0:5) * 100
smp_chr_cnt <- rep(100, 6)

sim_HMM_WFM <- cmpsimulateHMM(model, missing, sel_cof, dom_par, pop_siz, int_ale_frq, smp_gen, smp_chr_cnt)
smp_gen <- sim_HMM_WFM$smp_gen
smp_chr_cnt <- sim_HMM_WFM$smp_chr_cnt
smp_ale_cnt <- sim_HMM_WFM$smp_ale_cnt
mis_ale_cnt <- sim_HMM_WFM$mis_ale_cnt
smp_ale_frq <- sim_HMM_WFM$smp_ale_cnt / sim_HMM_WFM$smp_chr_cnt
pop_ale_frq <- sim_HMM_WFM$pop_ale_frq

save(sel_cof, dom_par, pop_siz, smp_gen, smp_chr_cnt, smp_ale_cnt, mis_ale_cnt, smp_ale_frq, pop_ale_frq, 
     file = "TEST_1L_simulated_dataset.rda")

load("TEST_1L_simulated_dataset.rda")

pdf(file = "TEST_1L_simulated_dataset.pdf", width = 16, height = 9)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
k <- min(smp_gen):max(smp_gen)
plot(k, pop_ale_frq, type = 'l', lwd = 1.5, 
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_ale_frq, pop_ale_frq), max(smp_ale_frq, pop_ale_frq)), 
     xlab = "Generation", ylab = "Allele frequency", main = "A simulated dataset generated with the W-F model")
points(smp_gen, smp_ale_frq, col = 'red', pch = 17, cex = 1)
dev.off()

########################################

#' Run the bootstrap particle filter (BPF) with the one-locus Wright-Fisher diffusion with selection
#' Parameter settings
#' @param gen_par the population genetic quantities of interest i.e., gen_par = (s, h, N)
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_chr_cnt the count of the chromosomes drawn from the population at all sampling time points
#' @param smp_ale_cnt the count of the mutant alleles observed in the sample at all sampling time points
#' @param mis_ale_cnt the count of the unknown alleles observed in the sample at all sampling time points
#' @param ptn_num the number of the subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of the particles generated in the bootstrap particle filter

load("TEST_1L_simulated_dataset.rda")

set.seed(test_seed)

gen_par <- c(sel_cof, dom_par, pop_siz)
smp_gen
smp_chr_cnt
smp_ale_cnt
mis_ale_cnt
ptn_num <- 1e+01
pcl_num <- 1e+06

system.time(BPF <- cmprunBPF(gen_par, smp_gen, smp_chr_cnt, smp_ale_cnt, mis_ale_cnt, ptn_num, pcl_num))

save(gen_par, smp_gen, smp_chr_cnt, smp_ale_cnt, mis_ale_cnt, smp_ale_frq, pop_ale_frq, ptn_num, pcl_num, BPF, 
     file = "TEST_1L_BPF.rda")

load("TEST_1L_BPF.rda")

lik <- rep(1, pcl_num)
wght <- BPF$wght
for (k in 1:length(smp_gen)) {
  lik <- lik * (cumsum(wght[, k]) / (1:pcl_num))
}
pdf(file = "TEST_1L_BPF_likelihood.pdf", width = 16, height = 9)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
plot(1:pcl_num, lik, type = 'l', 
     xlab = "Number of particles", ylab = "Likelihood", main = "Likelihood estimated with the bootstrap particle filter")
dev.off()

pop_ale_frq_pre_resmp <- BPF$pop_frq_pre_resmp
pop_ale_frq_pst_resmp <- BPF$pop_frq_pst_resmp
pdf(file = "TEST_1L_BPF_particle.pdf", width = 16, height = 9)
par(mfrow = c(2, 3), oma = c(0, 0, 3, 0), mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
for (k in 1:length(smp_gen)) {
  hist(pop_ale_frq_pst_resmp[, k], breaks = seq(min(pop_ale_frq_pst_resmp[, k], pop_ale_frq_pre_resmp[, k]), max(pop_ale_frq_pst_resmp[, k], pop_ale_frq_pre_resmp[, k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5), 
       xlim = c(min(pop_ale_frq_pst_resmp[, k], pop_ale_frq_pre_resmp[, k]), max(pop_ale_frq_pst_resmp[, k], pop_ale_frq_pre_resmp[, k])), 
       xlab = "Allele frequency", main = paste("A1 allele in generation", smp_gen[k]))
  hist(pop_ale_frq_pre_resmp[, k], breaks = seq(min(pop_ale_frq_pst_resmp[, k], pop_ale_frq_pre_resmp[, k]), max(pop_ale_frq_pst_resmp[, k], pop_ale_frq_pre_resmp[, k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_ale_frq[k], col = 'red', lty = 2, lwd = 2)
}
title("Histograms of the pre- and post-resampling particles", outer = TRUE, cex = 2)
dev.off()

########################################

#' Run the particle marginal Metropolis-Hastings (PMMH)
#' Parameter settings
#' @param int_gen_par the initial values of the population genetic quantities of interest i.e., int_gen_par = (s, h, N)
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_chr_cnt the count of the chromosomes drawn from the population at all sampling time points
#' @param smp_ale_cnt the count of the mutant alleles observed in the sample at all sampling time points
#' @param mis_ale_cnt the count of the unknown alleles observed in the sample at all sampling time points
#' @param ptn_num the number of the subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of the particles generated in the bootstrap particle filter
#' @param itn_num the number of the iterations carried out in the particle marginal Metropolis-Hastings

load("TEST_1L_simulated_dataset.rda")

set.seed(test_seed)

gen_par <- c(sel_cof, dom_par, pop_siz)

int_gen_par <- c(0e+00, dom_par, pop_siz)
smp_gen
smp_chr_cnt
smp_ale_cnt
mis_ale_cnt
ptn_num <- 1e+01
pcl_num <- 1e+03
itn_num <- 1e+05

system.time(PMMH <- cmprunPMMH(int_gen_par, smp_gen, smp_chr_cnt, smp_ale_cnt, mis_ale_cnt, ptn_num, pcl_num, itn_num))

save(gen_par, smp_gen, smp_chr_cnt, smp_ale_cnt, mis_ale_cnt, smp_ale_frq, pop_ale_frq, ptn_num, pcl_num, itn_num, PMMH, 
     file = "TEST_1L_PMMH.rda")

load("TEST_1L_PMMH.rda")

sel_cof_chn <- PMMH$sel_cof_chn
pdf(file = "TEST_1L_PMMH_trace_plot.pdf", width = 16, height = 9)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
plot(1:itn_num, sel_cof_chn[1:itn_num], type = 'l', 
     xlab = "Iteration", ylab = "Selection coefficient", main = "Trace plot of the selection coefficient")
abline(h = gen_par[1], col = 'red', lty = 2, lwd = 2)
dev.off()

brn_num <- 1e+04
sel_cof_chn <- sel_cof_chn[brn_num:length(sel_cof_chn)]

thn_num <- 9e+00
sel_cof_chn <- sel_cof_chn[(1:round(length(sel_cof_chn) / thn_num)) * thn_num]

grd_num <- 1e+03
sel_cof_pdf <- density(sel_cof_chn, n = grd_num)
pdf(file = "TEST_1L_PMMH_posterior.pdf", width = 16, height = 9)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
hist(sel_cof_chn, breaks = seq(min(sel_cof_chn), max(sel_cof_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient", main = "Posterior for the selection coefficient")
lines(sel_cof_pdf, lwd = 2, col = 'black')
abline(v = gen_par[1], col = 'red', lty = 2, lwd = 2)
dev.off()

########################################

#' Run the Bayesian Procedure for the inference of the selection coefficients
#' Parameter settings
#' @param int_gen_par the initial values of the population genetic quantities of interest i.e., int_gen_par = (s, h, N)
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_chr_cnt the count of the chromosomes drawn from the population at all sampling time points
#' @param smp_ale_cnt the count of the mutant alleles observed in the sample at all sampling time points
#' @param mis_ale_cnt the count of the unknown alleles observed in the sample at all sampling time points
#' @param ptn_num the number of the subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of the particles generated in the bootstrap particle filter
#' @param itn_num the number of the iterations carried out in the particle marginal Metropolis-Hastings
#' @param brn_num the number of the iterations for burn-in
#' @param thn_num the number of the iterations for thinning
#' @param grd_num the number of the grids in the kernel density estimation

load("TEST_1L_simulated_dataset.rda")

set.seed(test_seed)

gen_par <- c(sel_cof, dom_par, pop_siz)

int_gen_par <- c(0e+00, dom_par, pop_siz)
smp_gen
smp_chr_cnt
smp_ale_cnt
mis_ale_cnt
ptn_num <- 1e+01
pcl_num <- 1e+03
itn_num <- 1e+05
brn_num <- 1e+04
thn_num <- 9e+00
grd_num <- 1e+03

system.time(BayesianProcedure <- cmprunBayesianProcedure(int_gen_par, smp_gen, smp_chr_cnt, smp_ale_cnt, mis_ale_cnt, ptn_num, pcl_num, itn_num, brn_num, thn_num, grd_num))

save(gen_par, smp_gen, smp_chr_cnt, smp_ale_cnt, mis_ale_cnt, smp_ale_frq, pop_ale_frq, ptn_num, pcl_num, itn_num, brn_num, thn_num, grd_num, BayesianProcedure, 
     file = "TEST_1L_BayesianProcedure.rda")

load("TEST_1L_BayesianProcedure.rda")

sel_cof_chn <- BayesianProcedure$sel_cof_chn
sel_cof_pdf <- BayesianProcedure$sel_cof_pdf
sel_cof_map <- BayesianProcedure$sel_cof_map
sel_cof_mmse <- BayesianProcedure$sel_cof_mmse
sel_cof_hpd <- BayesianProcedure$sel_cof_hpd
pdf(file = "TEST_1L_BayesianProcedure_posterior.pdf", width = 16, height = 9)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
hist(sel_cof_chn, breaks = seq(min(sel_cof_chn), max(sel_cof_chn), length.out = 50), freq = FALSE,
     xlab = "Selection coefficient", main = "Posterior for the selection coefficient")
lines(sel_cof_pdf, lwd = 2, col = 'black')
abline(v = gen_par[1], col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_map, col = 'black', lty = 4, lwd = 2)
abline(v = sel_cof_mmse, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

################################################################################
