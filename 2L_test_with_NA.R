#' @title An MCMC-based method for Bayesian inference of natural selection from time series DNA data across linked loci
#' @author Zhangyi He, Xiaoyang Dai, Mark Beaumont and Feng Yu

#' this version is able to handle missing values in DNA data

source("2L_rfun_with_NA.R")

#install.packages("RColorBrewer")
library("RColorBrewer")

#install.packages("ggplot2")  
library("ggplot2")

#install.packages("plot3D")
library("plot3D")

################################################################################

#' Simulate the haplotype frequency trajectories according to the two-locus Wright-Fisher model with selection
#' Parameter setting
#' @param sel_cof_A the selection coefficient at locus A
#' @param dom_par_A the dominance parameter at locus A
#' @param sel_cof_B the selection coefficient at locus B
#' @param dom_par_B the dominance parameter at locus B
#' @param rec_rat the recombination rate between loci A and B
#' @param pop_siz the number of the diploid individuals in the population
#' @param int_hap_frq the initial haplotype frequencies of the population
#' @param int_gen the first generation of the simulated haplotype frequency trajectories
#' @param lst_gen the last generation of the simulated haplotype frequency trajectories

sel_cof_A <- 1e-02
dom_par_A <- 5e-01
sel_cof_B <- 5e-03
dom_par_B <- 5e-01
rec_rat <- 1e-03
pop_siz <- 5e+03
int_hap_frq <- c(1e-01, 2e-01, 3e-01, 4e-01)
int_gen <- 0
lst_gen <- 500

hap_frq_pth <- cmpsimulateTLWFMS(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_hap_frq, int_gen, lst_gen)

k <- int_gen:lst_gen
plot(k, hap_frq_pth[1, ], type = "l", lwd = 1.5, 
     xlab = "Generation", ylab = "Haplotype frequency", 
     main = "A haplotype frequency trajectory of A1B1 generated with the Wright-Fisher model")
plot(k, hap_frq_pth[2, ], type = "l", lwd = 1.5, 
     xlab = "Generation", ylab = "Haplotype frequency", 
     main = "A haplotype frequency trajectory of A1B2 generated with the Wright-Fisher model")
plot(k, hap_frq_pth[3, ], type = "l", lwd = 1.5, 
     xlab = "Generation", ylab = "Haplotype frequency", 
     main = "A haplotype frequency trajectory of A2B1 generated with the Wright-Fisher model")
plot(k, hap_frq_pth[4, ], type = "l", lwd = 1.5, 
     xlab = "Generation", ylab = "Haplotype frequency", 
     main = "A haplotype frequency trajectory of A2B2 generated with the Wright-Fisher model")

########################################

#' Simulate the haplotype frequency trajectories according to the two-locus Wright-Fisher diffusion with selection using the Euler-Maruyama method
#' Parameter setting
#' @param sel_cof_A the selection coefficient at locus A
#' @param dom_par_A the dominance parameter at locus A
#' @param sel_cof_B the selection coefficient at locus B
#' @param dom_par_B the dominance parameter at locus B
#' @param rec_rat the recombination rate between loci A and B
#' @param pop_siz the number of the diploid individuals in the population
#' @param int_hap_frq the initial haplotype frequencies of the population
#' @param int_gen the first generation of the simulated haplotype frequency trajectories
#' @param lst_gen the last generation of the simulated haplotype frequency trajectories
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param data_augmentation = TRUE/FALSE (return the simulated sample trajectory with data augmentation or not)

sel_cof_A <- 1e-02
dom_par_A <- 5e-01
sel_cof_B <- 5e-03
dom_par_B <- 5e-01
rec_rat <- 1e-03
pop_siz <- 5e+03
int_hap_frq <- c(1e-01, 2e-01, 3e-01, 4e-01)
int_gen <- 0
lst_gen <- 500
ptn_num <- 1e+01

hap_frq_pth <- cmpsimulateTLWFDS(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_hap_frq, int_gen, lst_gen, ptn_num, data_augmentation = TRUE)

t <- (int_gen:(int_gen + (lst_gen - int_gen) * ptn_num)) / 2 / pop_siz
plot(t, hap_frq_pth[1, ], type = "l", lwd = 1.5, 
     xlab = "Time", ylab = "Haplotype frequency", 
     main = "A haplotype frequency trajectory of A1B1 generated with the Wright-Fisher diffusion")
plot(t, hap_frq_pth[2, ], type = "l", lwd = 1.5, 
     xlab = "Time", ylab = "Haplotype frequency", 
     main = "A haplotype frequency trajectory of A1B2 generated with the Wright-Fisher diffusion")
plot(t, hap_frq_pth[3, ], type = "l", lwd = 1.5, 
     xlab = "Time", ylab = "Haplotype frequency", 
     main = "A haplotype frequency trajectory of A2B1 generated with the Wright-Fisher diffusion")
plot(t, hap_frq_pth[4, ], type = "l", lwd = 1.5, 
     xlab = "Time", ylab = "Haplotype frequency", 
     main = "A haplotype frequency trajectory of A2B2 generated with the Wright-Fisher diffusion")

########################################

#' Compare the simulation generated with the Wright-Fisher model and the Wright-Fisher diffusion
sel_cof_A <- 1e-02
dom_par_A <- 5e-01
sel_cof_B <- 5e-03
dom_par_B <- 5e-01
rec_rat <- 1e-03
pop_siz <- 5e+03
int_hap_frq <- c(1e-01, 2e-01, 3e-01, 4e-01)
int_gen <- 0
lst_gen <- 500
ptn_num <- 1e+01
sim_num <- 1e+06

sim_hap_frq_WFM <- matrix(NA, nrow = 4, ncol = sim_num)
sim_hap_frq_WFD <- matrix(NA, nrow = 4, ncol = sim_num)
for (i in 1:sim_num) {
  print(i)
  sim_hap_frq_WFM[, i] <- cmpsimulateTLWFMS(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_hap_frq, int_gen, lst_gen)[, (lst_gen - int_gen) + 1]
  sim_hap_frq_WFD[, i] <- cmpsimulateTLWFDS(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_hap_frq, int_gen, lst_gen, ptn_num, data_augmentation = FALSE)[, (lst_gen - int_gen) + 1]
}

save(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_hap_frq, int_gen, lst_gen, ptn_num, sim_num, sim_hap_frq_WFM, sim_hap_frq_WFD, 
     file = "TEST_2L_comparison_WFM_and_WFD.rda")

load("TEST_2L_comparison_WFM_and_WFD.rda")

pdf(file = "TEST_2L_comparison_WFM_and_WFD.pdf", width = 16, height = 9)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), oma = c(0, 0, 3, 0), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
hist(sim_hap_frq_WFM[1, ], breaks = seq(min(sim_hap_frq_WFM[1, ], sim_hap_frq_WFD[1, ]), max(sim_hap_frq_WFM[1, ], sim_hap_frq_WFD[1, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5), 
     xlim = c(min(sim_hap_frq_WFM[1, ], sim_hap_frq_WFD[1, ]), max(sim_hap_frq_WFM[1, ], sim_hap_frq_WFD[1, ])), 
     xlab = "Haplotype frequency", main = "A1B1 haplotype")
hist(sim_hap_frq_WFD[1, ], breaks = seq(min(sim_hap_frq_WFM[1, ], sim_hap_frq_WFD[1, ]), max(sim_hap_frq_WFM[1, ], sim_hap_frq_WFD[1, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
hist(sim_hap_frq_WFM[2, ], breaks = seq(min(sim_hap_frq_WFM[2, ], sim_hap_frq_WFD[2, ]), max(sim_hap_frq_WFM[2, ], sim_hap_frq_WFD[2, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5), 
     xlim = c(min(sim_hap_frq_WFM[2, ], sim_hap_frq_WFD[2, ]), max(sim_hap_frq_WFM[2, ], sim_hap_frq_WFD[2, ])), 
     xlab = "Haplotype frequency", main = "A1B2 haplotype")
hist(sim_hap_frq_WFD[2, ], breaks = seq(min(sim_hap_frq_WFM[2, ], sim_hap_frq_WFD[2, ]), max(sim_hap_frq_WFM[2, ], sim_hap_frq_WFD[2, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
hist(sim_hap_frq_WFM[3, ], breaks = seq(min(sim_hap_frq_WFM[3, ], sim_hap_frq_WFD[3, ]), max(sim_hap_frq_WFM[3, ], sim_hap_frq_WFD[3, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5), 
     xlim = c(min(sim_hap_frq_WFM[3, ], sim_hap_frq_WFD[3, ]), max(sim_hap_frq_WFM[3, ], sim_hap_frq_WFD[3, ])), 
     xlab = "Haplotype frequency", main = "A2B1 haplotype")
hist(sim_hap_frq_WFD[3, ], breaks = seq(min(sim_hap_frq_WFM[3, ], sim_hap_frq_WFD[3, ]), max(sim_hap_frq_WFM[3, ], sim_hap_frq_WFD[3, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
hist(sim_hap_frq_WFM[4, ], breaks = seq(min(sim_hap_frq_WFM[4, ], sim_hap_frq_WFD[4, ]), max(sim_hap_frq_WFM[4, ], sim_hap_frq_WFD[4, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5), 
     xlim = c(min(sim_hap_frq_WFM[4, ], sim_hap_frq_WFD[4, ]), max(sim_hap_frq_WFM[4, ], sim_hap_frq_WFD[4, ])), 
     xlab = "Haplotype frequency", main = "A2B2 haplotype")
hist(sim_hap_frq_WFD[4, ], breaks = seq(min(sim_hap_frq_WFM[4, ], sim_hap_frq_WFD[4, ]), max(sim_hap_frq_WFM[4, ], sim_hap_frq_WFD[4, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
title(paste("Histograms of the haplotype frequencies in generation", lst_gen, "under the W-F model and the W-F diffusion"), outer = TRUE)
dev.off()

################################################################################

#' Simulate the hidden Markov model
#' Parameter setting
#' @param model = "WFM"/"WFD" (return the observations from the underlying population evolving according to the WFM or the WFD)
#' @param phased = TRUE/FALSE (return the observations with phased or unphased chromosomes)
#' @param missing = TRUE/FALSE (return the observations with missing values or not)
#' @param sel_cof_A the selection coefficient at locus A
#' @param dom_par_A the dominance parameter at locus A
#' @param sel_cof_B the selection coefficient at locus B
#' @param dom_par_B the dominance parameter at locus B
#' @param rec_rat the recombination rate between loci A and B
#' @param pop_siz the number of individuals in the population
#' @param int_hap_frq the initial frequencies of the four haplotypes in the population
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_chr_cnt the count of the chromosomes drawn from the population at all sampling time points
#' @param ptn_num the number of the subintervals divided per generation in the Euler-Maruyama method for the WFD

#' Simulate the dataset under the Wright-Fisher model 
model <- "WFM"
missing <- TRUE
sel_cof_A <- 1e-02
dom_par_A <- 5e-01
sel_cof_B <- 5e-03
dom_par_B <- 5e-01
rec_rat <- 1e-03
pop_siz <- 5e+03
int_hap_frq <- c(1e-01, 2e-01, 3e-01, 4e-01)
smp_gen <- (0:5) * 100
smp_chr_cnt <- rep(100, 6)

phased <- TRUE
sim_HMM_WFM <- cmpsimulateHMM(model, phased, missing, sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_hap_frq, smp_gen, smp_chr_cnt)
smp_gen <- sim_HMM_WFM$smp_gen
smp_chr_cnt <- sim_HMM_WFM$smp_chr_cnt
smp_hap_cnt <- sim_HMM_WFM$smp_hap_cnt
mis_hap_cnt <- sim_HMM_WFM$mis_hap_cnt
smp_hap_frq <- sim_HMM_WFM$smp_hap_cnt / sim_HMM_WFM$smp_chr_cnt
pop_hap_frq <- sim_HMM_WFM$pop_hap_frq

k <- min(smp_gen):max(smp_gen)
plot(k, pop_hap_frq[1, ], type = 'l', lwd = 1.5, 
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_hap_frq[1, ], pop_hap_frq[1, ]), max(smp_hap_frq[1, ], pop_hap_frq[1, ])), 
     xlab = "Generation", ylab = "Haplotype frequency", 
     main = "A simulated dataset of the A1B1 haplotype generated with the Wright-Fisher model")
points(smp_gen, smp_hap_frq[1, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_hap_frq[2, ], type = 'l', lwd = 1.5, 
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_hap_frq[2, ], pop_hap_frq[2, ]), max(smp_hap_frq[2, ], pop_hap_frq[2, ])), 
     xlab = "Generation", ylab = "Haplotype frequency", 
     main = "A simulated dataset of the A1B2 haplotype generated with the Wright-Fisher model")
points(smp_gen, smp_hap_frq[2, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_hap_frq[3, ], type = 'l', lwd = 1.5, 
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_hap_frq[3, ], pop_hap_frq[3, ]), max(smp_hap_frq[3, ], pop_hap_frq[3, ])), 
     xlab = "Generation", ylab = "Haplotype frequency", 
     main = "A simulated dataset of the A2B1 haplotype generated with the Wright-Fisher model")
points(smp_gen, smp_hap_frq[3, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_hap_frq[4, ], type = 'l', lwd = 1.5, 
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_hap_frq[4, ], pop_hap_frq[4, ]), max(smp_hap_frq[4, ], pop_hap_frq[4, ])), 
     xlab = "Generation", ylab = "Haplotype frequency", 
     main = "A simulated dataset of the A2B2 haplotype generated with the Wright-Fisher model")
points(smp_gen, smp_hap_frq[4, ], col = 'red', pch = 17, cex = 1)

phased <- FALSE
sim_HMM_WFM <- cmpsimulateHMM(model, phased, missing, sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_hap_frq, smp_gen, smp_chr_cnt)
smp_gen <- sim_HMM_WFM$smp_gen
smp_chr_cnt <- sim_HMM_WFM$smp_chr_cnt
smp_ale_cnt <- sim_HMM_WFM$smp_ale_cnt
mis_ale_cnt <- sim_HMM_WFM$mis_ale_cnt
smp_ale_frq <- sim_HMM_WFM$smp_ale_cnt / sim_HMM_WFM$smp_chr_cnt
pop_ale_frq <- sim_HMM_WFM$pop_ale_frq

k <- min(smp_gen):max(smp_gen)
plot(k, pop_ale_frq[1, ], type = 'l', lwd = 1.5, 
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_ale_frq[1, ], pop_ale_frq[1, ]), max(smp_ale_frq[1, ], pop_ale_frq[1, ])), 
     xlab = "Generation", ylab = "Allele frequency", 
     main = "A simulated dataset of the A1 allele generated with the Wright-Fisher model")
points(smp_gen, smp_ale_frq[1, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_ale_frq[2, ], type = 'l', lwd = 1.5, 
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_ale_frq[2, ], pop_ale_frq[2, ]), max(smp_ale_frq[2, ], pop_ale_frq[2, ])), 
     xlab = "Generation", ylab = "Allele frequency", 
     main = "A simulated dataset of the B1 allele generated with the Wright-Fisher model")
points(smp_gen, smp_ale_frq[2, ], col = 'red', pch = 17, cex = 1)

####################

#' Simulate the dataset under the Wright-Fisher diffusion
model <- "WFD"
missing <- TRUE
sel_cof_A <- 1e-02
dom_par_A <- 5e-01
sel_cof_B <- 5e-03
dom_par_B <- 5e-01
rec_rat <- 1e-03
pop_siz <- 5e+03
int_hap_frq <- c(1e-01, 2e-01, 3e-01, 4e-01)
smp_gen <- (0:5) * 100
smp_chr_cnt <- rep(100, 6)
ptn_num <- 1e+01

phased <- TRUE
sim_HMM_WFD <- cmpsimulateHMM(model, phased, missing, sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_hap_frq, smp_gen, smp_chr_cnt, ptn_num)
smp_gen <- sim_HMM_WFD$smp_gen
smp_chr_cnt <- sim_HMM_WFD$smp_chr_cnt
smp_hap_cnt <- sim_HMM_WFD$smp_hap_cnt
mis_hap_cnt <- sim_HMM_WFD$mis_hap_cnt
smp_hap_frq <- sim_HMM_WFD$smp_hap_cnt / sim_HMM_WFD$smp_chr_cnt
pop_hap_frq <- sim_HMM_WFD$pop_hap_frq

k <- min(smp_gen):max(smp_gen)
plot(k, pop_hap_frq[1, ], type = 'l', lwd = 1.5, 
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_hap_frq[1, ], pop_hap_frq[1, ]), max(smp_hap_frq[1, ], pop_hap_frq[1, ])), 
     xlab = "Generation", ylab = "Haplotype frequency", 
     main = "A simulated dataset of the A1B1 haplotype generated with the Wright-Fisher diffusion")
points(smp_gen, smp_hap_frq[1, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_hap_frq[2, ], type = 'l', lwd = 1.5, 
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_hap_frq[2, ], pop_hap_frq[2, ]), max(smp_hap_frq[2, ], pop_hap_frq[2, ])), 
     xlab = "Generation", ylab = "Haplotype frequency", 
     main = "A simulated dataset of the A1B2 haplotype generated with the Wright-Fisher diffusion")
points(smp_gen, smp_hap_frq[2, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_hap_frq[3, ], type = 'l', lwd = 1.5, 
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_hap_frq[3, ], pop_hap_frq[3, ]), max(smp_hap_frq[3, ], pop_hap_frq[3, ])), 
     xlab = "Generation", ylab = "Haplotype frequency", 
     main = "A simulated dataset of the A2B1 haplotype generated with the Wright-Fisher diffusion")
points(smp_gen, smp_hap_frq[3, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_hap_frq[4, ], type = 'l', lwd = 1.5, 
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_hap_frq[4, ], pop_hap_frq[4, ]), max(smp_hap_frq[4, ], pop_hap_frq[4, ])), 
     xlab = "Generation", ylab = "Haplotype frequency", 
     main = "A simulated dataset of the A2B2 haplotype generated with the Wright-Fisher diffusion")
points(smp_gen, smp_hap_frq[4, ], col = 'red', pch = 17, cex = 1)

phased <- FALSE
sim_HMM_WFD <- cmpsimulateHMM(model, phased, missing, sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_hap_frq, smp_gen, smp_chr_cnt, ptn_num)
smp_gen <- sim_HMM_WFD$smp_gen
smp_chr_cnt <- sim_HMM_WFD$smp_chr_cnt
smp_ale_cnt <- sim_HMM_WFD$smp_ale_cnt
mis_ale_cnt <- sim_HMM_WFD$mis_ale_cnt
smp_ale_frq <- sim_HMM_WFD$smp_ale_cnt / sim_HMM_WFD$smp_chr_cnt
pop_ale_frq <- sim_HMM_WFD$pop_ale_frq

k <- min(smp_gen):max(smp_gen)
plot(k, pop_ale_frq[1, ], type = 'l', lwd = 1.5, 
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_ale_frq[1, ], pop_ale_frq[1, ]), max(smp_ale_frq[1, ], pop_ale_frq[1, ])), 
     xlab = "Generation", ylab = "Allele frequency", 
     main = "A simulated dataset of the A1 allele generated with the Wright-Fisher diffusion")
points(smp_gen, smp_ale_frq[1, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_ale_frq[2, ], type = 'l', lwd = 1.5, 
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_ale_frq[2, ], pop_ale_frq[2, ]), max(smp_ale_frq[2, ], pop_ale_frq[2, ])), 
     xlab = "Generation", ylab = "Allele frequency", 
     main = "A simulated dataset of the B1 allele generated with the Wright-Fisher diffusion")
points(smp_gen, smp_ale_frq[2, ], col = 'red', pch = 17, cex = 1)

################################################################################

#' Generate a simulated dataset with phased chromosomes under the Wright-Fisher model 

test_seed <- 21
set.seed(test_seed)

model <- "WFM"
phased <- TRUE
missing <- TRUE
sel_cof_A <- 1e-02
dom_par_A <- 5e-01
sel_cof_B <- 5e-03
dom_par_B <- 5e-01
rec_rat <- 1e-03
pop_siz <- 5e+03
int_hap_frq <- c(1e-01, 2e-01, 3e-01, 4e-01)
smp_gen <- (0:5) * 100
smp_chr_cnt <- rep(100, 6)

sim_HMM_WFM <- cmpsimulateHMM(model, phased, missing, sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_hap_frq, smp_gen, smp_chr_cnt)
smp_gen <- sim_HMM_WFM$smp_gen
smp_chr_cnt <- sim_HMM_WFM$smp_chr_cnt
smp_hap_cnt <- sim_HMM_WFM$smp_hap_cnt
mis_hap_cnt <- sim_HMM_WFM$mis_hap_cnt
smp_hap_frq <- sim_HMM_WFM$smp_hap_cnt / sim_HMM_WFM$smp_chr_cnt
pop_hap_frq <- sim_HMM_WFM$pop_hap_frq

save(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, smp_gen, smp_chr_cnt, smp_hap_cnt, mis_hap_cnt, smp_hap_frq, pop_hap_frq, 
     file = "TEST_2L_PhasedChr_simulated_dataset.rda")

load("TEST_2L_PhasedChr_simulated_dataset.rda")

pdf(file = "TEST_2L_PhasedChr_simulated_dataset.pdf", width = 16, height = 9)
par(mfrow = c(2, 2), oma = c(0, 0, 3, 0), mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
k <- min(smp_gen):max(smp_gen)
plot(k, pop_hap_frq[1, ], type = 'l', lwd = 1.5, 
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(pop_hap_frq[1, ], smp_hap_frq[1, ]), max(pop_hap_frq[1, ], smp_hap_frq[1, ])), 
     xlab = "Generation", ylab = "Haplotype frequency", 
     main = "A1B1 haplotype")
points(smp_gen, smp_hap_frq[1, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_hap_frq[2, ], type = 'l', lwd = 1.5, 
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(pop_hap_frq[2, ], smp_hap_frq[2, ]), max(pop_hap_frq[2, ], smp_hap_frq[2, ])), 
     xlab = "Generation", ylab = "Haplotype frequency", 
     main = "A1B2 haplotype")
points(smp_gen, smp_hap_frq[2, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_hap_frq[3, ], type = 'l', lwd = 1.5, 
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(pop_hap_frq[3, ], smp_hap_frq[3, ]), max(pop_hap_frq[3, ], smp_hap_frq[3, ])), 
     xlab = "Generation", ylab = "Haplotype frequency", 
     main = "A2B1 haplotype")
points(smp_gen, smp_hap_frq[3, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_hap_frq[4, ], type = 'l', lwd = 1.5, 
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(pop_hap_frq[4, ], smp_hap_frq[4, ]), max(pop_hap_frq[4, ], smp_hap_frq[4, ])), 
     xlab = "GenerationG", ylab = "Haplotype frequency", 
     main = "A2B2 haplotype")
points(smp_gen, smp_hap_frq[4, ], col = 'red', pch = 17, cex = 1)
title("A simulated dataset with phased chromosomes generated with the W-F model", outer = TRUE)
dev.off()

####################

#' Generate a simulated dataset with unphased chromosomes under the Wright-Fisher model

test_seed <- 21
set.seed(test_seed)

model <- "WFM"
phased <- FALSE
missing <- TRUE
sel_cof_A <- 1e-02
dom_par_A <- 5e-01
sel_cof_B <- 5e-03
dom_par_B <- 5e-01
rec_rat <- 1e-03
pop_siz <- 5e+03
int_hap_frq <- c(1e-01, 2e-01, 3e-01, 4e-01)
smp_gen <- (0:5) * 100
smp_chr_cnt <- rep(100, 6)

sim_HMM_WFM <- cmpsimulateHMM(model, phased, missing, sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_hap_frq, smp_gen, smp_chr_cnt)
smp_gen <- sim_HMM_WFM$smp_gen
smp_chr_cnt <- sim_HMM_WFM$smp_chr_cnt
smp_ale_cnt <- sim_HMM_WFM$smp_ale_cnt
mis_ale_cnt <- sim_HMM_WFM$mis_ale_cnt
smp_ale_frq <- sim_HMM_WFM$smp_ale_cnt / smp_chr_cnt
pop_ale_frq <- sim_HMM_WFM$pop_ale_frq
smp_hap_cnt <- sim_HMM_WFM$smp_hap_cnt
smp_hap_frq <- sim_HMM_WFM$smp_hap_cnt / smp_chr_cnt
pop_hap_frq <- sim_HMM_WFM$pop_hap_frq

save(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, smp_gen, smp_chr_cnt, smp_ale_cnt, mis_ale_cnt, smp_ale_frq, pop_ale_frq, smp_hap_cnt, smp_hap_frq, pop_hap_frq, 
     file = "TEST_2L_UnphasedChr_simulated_dataset.rda")

load("TEST_2L_UnphasedChr_simulated_dataset.rda")

pdf(file = "TEST_2L_UnphasedChr_simulated_dataset.pdf", width = 16, height = 9)
par(mfrow = c(1, 2), oma = c(0, 0, 3, 0), mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
k <- min(smp_gen):max(smp_gen)
plot(k, pop_ale_frq[1, ], type = 'l', lwd = 1.5, 
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(pop_ale_frq[1, ], smp_ale_frq[1, ]), max(pop_ale_frq[1, ], smp_ale_frq[1, ])), 
     xlab = "Generation", ylab = "Allele frequency", 
     main = "A1 allele")
points(smp_gen, smp_ale_frq[1, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_ale_frq[2, ], type = 'l', lwd = 1.5, 
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(pop_ale_frq[2, ], smp_ale_frq[2, ]), max(pop_ale_frq[2, ], smp_ale_frq[2, ])), 
     xlab = "Generation", ylab = "Allele frequency", 
     main = "B1 allele")
points(smp_gen, smp_ale_frq[2, ], col = 'red', pch = 17, cex = 1)
title("A simulated dataset with unphased chromosomes generated with the W-F model", outer = TRUE)
dev.off()

########################################

#' Run the bootstrap particle filter (BPF) with the two-locus Wright-Fisher diffusion with selection
#' Parameter setting
#' @param gen_par the population genetic quantities of interest i.e., gen_par = (s_A, h_A, s_B, h_B, r, N)
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_chr_cnt the count of the chromosomes drawn from the population at all sampling time points
#' @param smp_cnt the count of the two mutant alleles or the four haplotypes observed in the sample at all sampling time points
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter

load("TEST_2L_PhasedChr_simulated_dataset.rda")

set.seed(test_seed)

gen_par <- c(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz)
smp_gen
smp_chr_cnt
smp_hap_cnt
mis_hap_cnt
ptn_num <- 1e+01
pcl_num <- 1e+06

system.time(BPF <- cmprunBPF(gen_par, smp_gen, smp_chr_cnt, smp_hap_cnt, mis_hap_cnt, ptn_num, pcl_num))

save(gen_par, smp_gen, smp_chr_cnt, smp_hap_cnt, mis_hap_cnt, smp_hap_frq, pop_hap_frq, ptn_num, pcl_num, BPF, 
     file = "TEST_2L_PhasedChr_BPF.rda")

load("TEST_2L_PhasedChr_BPF.rda")

lik <- rep(1, pcl_num)
wght <- BPF$wght
for (k in 1:length(smp_gen)) {
  lik <- lik * (cumsum(wght[, k]) / (1:pcl_num))
}
pdf(file = "TEST_2L_PhasedChr_BPF_likelihood.pdf", width = 16, height = 9)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
plot(1:pcl_num, lik, type = 'l', 
     xlab = "Number of particles", ylab = "Likelihood", main = "Likelihood estimated with the bootstrap particle filter")
dev.off()

pop_hap_frq_pre_resmp <- BPF$pop_frq_pre_resmp
pop_hap_frq_pst_resmp <- BPF$pop_frq_pst_resmp
pdf(file = "TEST_2L_PhasedChr_BPF_particle.pdf", width = 16, height = 27)
par(mfrow = c(6, 4), oma = c(0, 0, 3, 0), mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
for (k in 1:length(smp_gen)) {
  hist(pop_hap_frq_pst_resmp[1, , k], breaks = seq(min(pop_hap_frq_pst_resmp[1, , k], pop_hap_frq_pre_resmp[1, , k]), max(pop_hap_frq_pst_resmp[1, , k], pop_hap_frq_pre_resmp[1, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5), 
       xlim = c(min(pop_hap_frq_pst_resmp[1, , k], pop_hap_frq_pre_resmp[1, , k], smp_hap_frq[1, k]), max(pop_hap_frq_pst_resmp[1, , k], pop_hap_frq_pre_resmp[1, , k], smp_hap_frq[1, k])), 
       xlab = "Haplotype frequency", main = paste("A1B1 haplotype in generation", smp_gen[k]))
  hist(pop_hap_frq_pre_resmp[1, , k], breaks = seq(min(pop_hap_frq_pst_resmp[1, , k], pop_hap_frq_pre_resmp[1, , k]), max(pop_hap_frq_pst_resmp[1, , k], pop_hap_frq_pre_resmp[1, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_hap_frq[1, k], col = 'red', lty = 2, lwd = 2)
  hist(pop_hap_frq_pst_resmp[2, , k], breaks = seq(min(pop_hap_frq_pst_resmp[2, , k], pop_hap_frq_pre_resmp[2, , k]), max(pop_hap_frq_pst_resmp[2, , k], pop_hap_frq_pre_resmp[2, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5), 
       xlim = c(min(pop_hap_frq_pst_resmp[2, , k], pop_hap_frq_pre_resmp[2, , k], smp_hap_frq[2, k]), max(pop_hap_frq_pst_resmp[2, , k], pop_hap_frq_pre_resmp[2, , k], smp_hap_frq[2, k])), 
       xlab = "Haplotype frequency", main = paste("A1B2 haplotype in generation", smp_gen[k]))
  hist(pop_hap_frq_pre_resmp[2, , k], breaks = seq(min(pop_hap_frq_pst_resmp[2, , k], pop_hap_frq_pre_resmp[2, , k]), max(pop_hap_frq_pst_resmp[2, , k], pop_hap_frq_pre_resmp[2, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_hap_frq[2, k], col = 'red', lty = 2, lwd = 2)
  hist(pop_hap_frq_pst_resmp[3, , k], breaks = seq(min(pop_hap_frq_pst_resmp[3, , k], pop_hap_frq_pre_resmp[3, , k]), max(pop_hap_frq_pst_resmp[3, , k], pop_hap_frq_pre_resmp[3, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5), 
       xlim = c(min(pop_hap_frq_pst_resmp[3, , k], pop_hap_frq_pre_resmp[3, , k], smp_hap_frq[3, k]), max(pop_hap_frq_pst_resmp[3, , k], pop_hap_frq_pre_resmp[3, , k], smp_hap_frq[3, k])), 
       xlab = "Haplotype frequency", main = paste("A2B1 haplotype in generation", smp_gen[k]))
  hist(pop_hap_frq_pre_resmp[3, , k], breaks = seq(min(pop_hap_frq_pst_resmp[3, , k], pop_hap_frq_pre_resmp[3, , k]), max(pop_hap_frq_pst_resmp[3, , k], pop_hap_frq_pre_resmp[3, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_hap_frq[3, k], col = 'red', lty = 2, lwd = 2)
  hist(pop_hap_frq_pst_resmp[4, , k], breaks = seq(min(pop_hap_frq_pst_resmp[4, , k], pop_hap_frq_pre_resmp[4, , k]), max(pop_hap_frq_pst_resmp[4, , k], pop_hap_frq_pre_resmp[4, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5), 
       xlim = c(min(pop_hap_frq_pst_resmp[4, , k], pop_hap_frq_pre_resmp[4, , k], smp_hap_frq[4, k]), max(pop_hap_frq_pst_resmp[4, , k], pop_hap_frq_pre_resmp[4, , k], smp_hap_frq[4, k])), 
       xlab = "Haplotype frequency", main = paste("A2B2 haplotype in generation", smp_gen[k]))
  hist(pop_hap_frq_pre_resmp[4, , k], breaks = seq(min(pop_hap_frq_pst_resmp[4, , k], pop_hap_frq_pre_resmp[4, , k]), max(pop_hap_frq_pst_resmp[4, , k], pop_hap_frq_pre_resmp[4, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_hap_frq[4, k], col = 'red', lty = 2, lwd = 2)
}
title("Histograms of the pre- and post-resampling particles", outer = TRUE)
dev.off()

####################

load("TEST_2L_UnphasedChr_simulated_dataset.rda")

set.seed(test_seed)

gen_par <- c(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz)
smp_gen
smp_chr_cnt
smp_ale_cnt
mis_ale_cnt
ptn_num <- 1e+01
pcl_num <- 1e+06

system.time(BPF <- cmprunBPF(gen_par, smp_gen, smp_chr_cnt, smp_ale_cnt, mis_ale_cnt, ptn_num, pcl_num))

save(gen_par, smp_gen, smp_chr_cnt, smp_hap_cnt, smp_hap_frq, pop_hap_frq, smp_ale_cnt, mis_ale_cnt, smp_ale_frq, pop_ale_frq, ptn_num, pcl_num, BPF, 
     file = "TEST_2L_UnphasedChr_BPF.rda")

load("TEST_2L_UnphasedChr_BPF.rda")

lik <- rep(1, pcl_num)
wght <- BPF$wght
for (k in 1:length(smp_gen)) {
  lik <- lik * (cumsum(wght[, k]) / (1:pcl_num))
}
pdf(file = "TEST_2L_UnphasedChr_BPF_likelihood.pdf", width = 16, height = 9)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
plot(1:pcl_num, lik, type = 'l', 
     xlab = "Number of particles", ylab = "Likelihood", main = "Likelihood estimated with the bootstrap particle filter")
dev.off()

pop_ale_frq_pre_resmp <- BPF$pop_frq_pre_resmp
pop_ale_frq_pst_resmp <- BPF$pop_frq_pst_resmp
pdf(file = "TEST_2L_UnphasedChr_BPF_particle.pdf", width = 16, height = 27)
par(mfrow = c(6, 2), oma = c(0, 0, 3, 0), mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
for (k in 1:length(smp_gen)) {
  hist(pop_ale_frq_pst_resmp[1, , k], breaks = seq(min(pop_ale_frq_pst_resmp[1, , k], pop_ale_frq_pre_resmp[1, , k]), max(pop_ale_frq_pst_resmp[1, , k], pop_ale_frq_pre_resmp[1, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5), 
       xlim = c(min(pop_ale_frq_pst_resmp[1, , k], pop_ale_frq_pre_resmp[1, , k], smp_ale_frq[1, k]), max(pop_ale_frq_pst_resmp[1, , k], pop_ale_frq_pre_resmp[1, , k], smp_ale_frq[1, k])), 
       xlab = "Allele frequency", main = paste("A1 allele in generation", smp_gen[k]))
  hist(pop_ale_frq_pre_resmp[1, , k], breaks = seq(min(pop_ale_frq_pst_resmp[1, , k], pop_ale_frq_pre_resmp[1, , k]), max(pop_ale_frq_pst_resmp[1, , k], pop_ale_frq_pre_resmp[1, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_ale_frq[1, k], col = 'red', lty = 2, lwd = 2)
  hist(pop_ale_frq_pst_resmp[2, , k], breaks = seq(min(pop_ale_frq_pst_resmp[2, , k], pop_ale_frq_pre_resmp[2, , k]), max(pop_ale_frq_pst_resmp[2, , k], pop_ale_frq_pre_resmp[2, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5), 
       xlim = c(min(pop_ale_frq_pst_resmp[2, , k], pop_ale_frq_pre_resmp[2, , k], smp_ale_frq[2, k]), max(pop_ale_frq_pst_resmp[2, , k], pop_ale_frq_pre_resmp[2, , k], smp_ale_frq[2, k])), 
       xlab = "Allele frequency", main = paste("B1 allele in generation", smp_gen[k]))
  hist(pop_ale_frq_pre_resmp[2, , k], breaks = seq(min(pop_ale_frq_pst_resmp[2, , k], pop_ale_frq_pre_resmp[2, , k]), max(pop_ale_frq_pst_resmp[2, , k], pop_ale_frq_pre_resmp[2, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_ale_frq[2, k], col = 'red', lty = 2, lwd = 2)
}
title("Histograms of the pre- and post-resampling particles", outer = TRUE)
dev.off()

########################################

#' Run the particle marginal Metropolis-Hastings (PMMH)
#' Parameter settings
#' @param int_gen_par the initial values of the population genetic quantities of interest i.e., int_gen_par = (s_A, h_A, s_B, h_B, r, N)
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_chr_cnt the count of the chromosomes drawn from the population at all sampling time points
#' @param smp_cnt the count of the two mutant alleles or the four haplotypes observed in the sample at all sampling time points
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param itn_num the number of the iterations carried out in the particle marginal Metropolis-Hastings

load("TEST_2L_PhasedChr_simulated_dataset.rda")

set.seed(test_seed)

gen_par <- c(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz)

int_gen_par <- c(0e+00, dom_par_A, 0e+00, dom_par_B, rec_rat, pop_siz)
smp_gen
smp_chr_cnt
smp_hap_cnt
mis_hap_cnt
ptn_num <- 1e+01
pcl_num <- 1e+03
itn_num <- 2e+05

system.time(PMMH <- cmprunPMMH(int_gen_par, smp_gen, smp_chr_cnt, smp_hap_cnt, mis_hap_cnt, ptn_num, pcl_num, itn_num))

save(gen_par, smp_gen, smp_chr_cnt, smp_hap_cnt, mis_hap_cnt, smp_hap_frq, pop_hap_frq, ptn_num, pcl_num, itn_num, PMMH, 
     file = "TEST_2L_PhasedChr_PMMH.rda")

load("TEST_2L_PhasedChr_PMMH.rda")

sel_cof_A_chn <- PMMH$sel_cof_A_chn
sel_cof_B_chn <- PMMH$sel_cof_B_chn
pdf(file = "TEST_2L_PhasedChr_PMMH_trace_plot.pdf", width = 16, height = 9)
par(mfrow = c(2, 1), mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
plot(1:itn_num, sel_cof_A_chn[1:itn_num], type = 'l', 
     xlab = "Iteration", ylab = "Selection coefficient", main = "Trace plot of the selection coefficient at locus A")
abline(h = gen_par[1], col = 'red', lty = 2, lwd = 2)
plot(1:itn_num, sel_cof_B_chn[1:itn_num], type = 'l', 
     xlab = "Iteration", ylab = "Selection coefficient", main = "Trace plot of the selection coefficient at locus B")
abline(h = gen_par[3], col = 'red', lty = 2, lwd = 2)
dev.off()

brn_num <- 2e+04
sel_cof_A_chn <- sel_cof_A_chn[brn_num:length(sel_cof_A_chn)]
sel_cof_B_chn <- sel_cof_B_chn[brn_num:length(sel_cof_B_chn)]

thn_num <- 9e+00
sel_cof_A_chn <- sel_cof_A_chn[(1:round(length(sel_cof_A_chn) / thn_num)) * thn_num]
sel_cof_B_chn <- sel_cof_B_chn[(1:round(length(sel_cof_B_chn) / thn_num)) * thn_num]

grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_A_chn, sel_cof_B_chn, n = grd_num)
pdf(file = "TEST_2L_PhasedChr_PMMH_posterior.pdf", width = 16, height = 9)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32), 
      xlab = "Selection coefficient at locus A", ylab = "Selection coefficient at locus B", 
      main = "Posterior for the selection coefficients at loci A and B")
abline(v = gen_par[1], col = 'red', lty = 2, lwd = 2)
abline(h = gen_par[3], col = 'red', lty = 2, lwd = 2)
dev.off()

####################

load("TEST_2L_UnphasedChr_simulated_dataset.rda")

set.seed(test_seed)

int_gen_par <- c(0e+00, dom_par_A, 0e+00, dom_par_B, rec_rat, pop_siz)
smp_gen
smp_chr_cnt
smp_ale_cnt
mis_ale_cnt
ptn_num <- 1e+01
pcl_num <- 1e+03
itn_num <- 2e+05

system.time(PMMH <- cmprunPMMH(int_gen_par, smp_gen, smp_chr_cnt, smp_ale_cnt, mis_ale_cnt, ptn_num, pcl_num, itn_num))

save(gen_par, smp_gen, smp_chr_cnt, smp_hap_cnt, smp_hap_frq, pop_hap_frq, smp_ale_cnt, mis_ale_cnt, smp_ale_frq, pop_ale_frq, ptn_num, pcl_num, itn_num, PMMH, 
     file = "TEST_2L_UnphasedChr_PMMH.rda")

load("TEST_2L_UnphasedChr_PMMH.rda")

sel_cof_A_chn <- PMMH$sel_cof_A_chn
sel_cof_B_chn <- PMMH$sel_cof_B_chn
pdf(file = "TEST_2L_UnphasedChr_PMMH_trace_plot.pdf", width = 16, height = 9)
par(mfrow = c(2, 1), mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
plot(1:itn_num, sel_cof_A_chn[1:itn_num], type = 'l', 
     xlab = "Iteration", ylab = "Selection coefficient", main = "Trace plot of the selection coefficient at locus A")
abline(h = gen_par[1], col = 'red', lty = 2, lwd = 2)
plot(1:itn_num, sel_cof_B_chn[1:itn_num], type = 'l', 
     xlab = "Iteration", ylab = "Selection coefficient", main = "Trace plot of the selection coefficient at locus B")
abline(h = gen_par[3], col = 'red', lty = 2, lwd = 2)
dev.off()

brn_num <- 2e+04
sel_cof_A_chn <- sel_cof_A_chn[brn_num:length(sel_cof_A_chn)]
sel_cof_B_chn <- sel_cof_B_chn[brn_num:length(sel_cof_B_chn)]

thn_num <- 9e+00
sel_cof_A_chn <- sel_cof_A_chn[(1:round(length(sel_cof_A_chn) / thn_num)) * thn_num]
sel_cof_B_chn <- sel_cof_B_chn[(1:round(length(sel_cof_B_chn) / thn_num)) * thn_num]

grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_A_chn, sel_cof_B_chn, n = grd_num)
pdf(file = "TEST_2L_UnphasedChr_PMMH_posterior.pdf", width = 16, height = 9)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32), 
      xlab = "Selection coefficient at locus A", ylab = "Selection coefficient at locus B", 
      main = "Posterior for the selection coefficients at loci A and B")
abline(v = gen_par[1], col = 'red', lty = 2, lwd = 2)
abline(h = gen_par[3], col = 'red', lty = 2, lwd = 2)
dev.off()

########################################

#' Run the Bayesian procedure for the inference of natural selection
#' Parameter settings
#' @param int_gen_par the initial values of the population genetic quantities of interest i.e., int_gen_par = (s_A, h_A, s_B, h_B, r, N)
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_chr_cnt the count of the chromosomes drawn from the population at all sampling time points
#' @param smp_cnt the count of the two mutant alleles or the four haplotypes observed in the sample at all sampling time points
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of the particles generated in the bootstrap particle filter
#' @param itn_num the number of the iterations carried out in the particle marginal Metropolis-Hastings
#' @param brn_num the number of the iterations for burn-in
#' @param thn_num the number of the iterations for thinning
#' @param grd_num the number of the grids in the kernel density estimation

load("TEST_2L_PhasedChr_simulated_dataset.rda")

set.seed(test_seed)

gen_par <- c(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz)

int_gen_par <- c(0e+00, dom_par_A, 0e+00, dom_par_B, rec_rat, pop_siz)
smp_gen
smp_chr_cnt
smp_hap_cnt
mis_hap_cnt
ptn_num <- 1e+01
pcl_num <- 1e+03
itn_num <- 2e+05
brn_num <- 2e+04
thn_num <- 9e+00
grd_num <- 1e+03

system.time(BayesianProcedure <- cmprunBayesianProcedure(int_gen_par, smp_gen, smp_chr_cnt, smp_hap_cnt, mis_hap_cnt, ptn_num, pcl_num, itn_num, brn_num, thn_num, grd_num))

save(gen_par, smp_gen, smp_chr_cnt, smp_hap_cnt, mis_hap_cnt, smp_hap_frq, pop_hap_frq, ptn_num, pcl_num, itn_num, brn_num, thn_num, grd_num, BayesianProcedure, 
     file = "TEST_2L_PhasedChr_BayesianProcedure.rda")

load("TEST_2L_PhasedChr_BayesianProcedure.rda")

sel_cof_A_chn <- BayesianProcedure$sel_cof_A_chn
sel_cof_B_chn <- BayesianProcedure$sel_cof_B_chn
sel_cof_pdf <- BayesianProcedure$sel_cof_pdf
sel_cof_A_map <- BayesianProcedure$sel_cof_A_map
sel_cof_B_map <- BayesianProcedure$sel_cof_B_map
sel_cof_A_mmse <- BayesianProcedure$sel_cof_A_mmse
sel_cof_B_mmse <- BayesianProcedure$sel_cof_B_mmse
sel_cof_B_hpd <- BayesianProcedure$sel_cof_B_hpd
sel_cof_A_hpd <- BayesianProcedure$sel_cof_A_hpd
pdf(file = "TEST_2L_PhasedChr_BayesianProcedure_posterior.pdf", width = 16, height = 9)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
layout(matrix(c(1, 1, 2, 3), nrow = 2, ncol = 2))
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32), 
      xlab = "Selection coefficient at locus A", ylab = "Selection coefficient at locus B", 
      main = "Joint posterior for the selection coefficients at loci A and B")
abline(v = gen_par[1], col = 'red', lty = 2, lwd = 2)
abline(h = gen_par[3], col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_A_map, col = 'black', lty = 4, lwd = 2)
abline(h = sel_cof_B_map, col = 'black', lty = 4, lwd = 2)
abline(v = sel_cof_A_mmse, col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_B_mmse, col = 'black', lty = 2, lwd = 2)

hist(sel_cof_A_chn, breaks = seq(min(sel_cof_A_chn), max(sel_cof_A_chn), length.out = 50), freq = FALSE, 
     xlab = "Selection coefficient", main = "Marginal posterior for the selection coefficients at locus A")
lines(density(sel_cof_A_chn), lwd = 2, col = 'black')
abline(v = sel_cof_A_map, col = 'black', lty = 4, lwd = 2)
abline(v = sel_cof_A_mmse, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_A_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_A_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_B_chn, breaks = seq(min(sel_cof_B_chn), max(sel_cof_B_chn), length.out = 50), freq = FALSE, 
     xlab = "Selection coefficient", main = "Marginal posterior for the selection coefficients at locus B")
lines(density(sel_cof_B_chn), lwd = 2, col = 'black')
abline(v = sel_cof_B_map, col = 'black', lty = 4, lwd = 2)
abline(v = sel_cof_B_mmse, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_B_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_B_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

####################

load("TEST_2L_UnphasedChr_simulated_dataset.rda")

set.seed(test_seed)

gen_par <- c(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz)

int_gen_par <- c(0e+00, dom_par_A, 0e+00, dom_par_B, rec_rat, pop_siz)
smp_gen
smp_chr_cnt
smp_ale_cnt
mis_ale_cnt
ptn_num <- 1e+01
pcl_num <- 1e+03
itn_num <- 2e+05
brn_num <- 2e+04
thn_num <- 9e+00
grd_num <- 1e+03

system.time(BayesianProcedure <- cmprunBayesianProcedure(int_gen_par, smp_gen, smp_chr_cnt, smp_ale_cnt, mis_ale_cnt, ptn_num, pcl_num, itn_num, brn_num, thn_num, grd_num))

save(gen_par, smp_gen, smp_chr_cnt, smp_hap_cnt, smp_hap_frq, pop_hap_frq, smp_ale_cnt, mis_ale_cnt, smp_ale_frq, pop_ale_frq, ptn_num, pcl_num, itn_num, brn_num, thn_num, grd_num, BayesianProcedure, 
     file = "TEST_2L_UnphasedChr_BayesianProcedure.rda")

load("TEST_2L_UnphasedChr_BayesianProcedure.rda")

sel_cof_A_chn <- BayesianProcedure$sel_cof_A_chn
sel_cof_B_chn <- BayesianProcedure$sel_cof_B_chn
sel_cof_pdf <- BayesianProcedure$sel_cof_pdf
sel_cof_A_map <- BayesianProcedure$sel_cof_A_map
sel_cof_B_map <- BayesianProcedure$sel_cof_B_map
sel_cof_A_mmse <- BayesianProcedure$sel_cof_A_mmse
sel_cof_B_mmse <- BayesianProcedure$sel_cof_B_mmse
sel_cof_B_hpd <- BayesianProcedure$sel_cof_B_hpd
sel_cof_A_hpd <- BayesianProcedure$sel_cof_A_hpd
pdf(file = "TEST_2L_UnphasedChr_BayesianProcedure_posterior.pdf", width = 16, height = 9)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
layout(matrix(c(1, 1, 2, 3), nrow = 2, ncol = 2))
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32), 
      xlab = "Selection coefficient at locus A", ylab = "Selection coefficient at locus B", 
      main = "Joint posterior for the selection coefficients at loci A and B")
abline(v = gen_par[1], col = 'red', lty = 2, lwd = 2)
abline(h = gen_par[3], col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_A_map, col = 'black', lty = 4, lwd = 2)
abline(h = sel_cof_B_map, col = 'black', lty = 4, lwd = 2)
abline(v = sel_cof_A_mmse, col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_B_mmse, col = 'black', lty = 2, lwd = 2)

hist(sel_cof_A_chn, breaks = seq(min(sel_cof_A_chn), max(sel_cof_A_chn), length.out = 50), freq = FALSE, 
     xlab = "Selection coefficient", main = "Marginal posterior for the selection coefficients at locus A")
lines(density(sel_cof_A_chn), lwd = 2, col = 'black')
abline(v = sel_cof_A_map, col = 'black', lty = 4, lwd = 2)
abline(v = sel_cof_A_mmse, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_A_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_A_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_B_chn, breaks = seq(min(sel_cof_B_chn), max(sel_cof_B_chn), length.out = 50), freq = FALSE, 
     xlab = "Selection coefficient", main = "Marginal posterior for the selection coefficients at locus B")
lines(density(sel_cof_B_chn), lwd = 2, col = 'black')
abline(v = sel_cof_B_map, col = 'black', lty = 4, lwd = 2)
abline(v = sel_cof_B_mmse, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_B_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_B_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

################################################################################
