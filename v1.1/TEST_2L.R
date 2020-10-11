#' @title Detecting and quantifying natural selection at two linked loci from time series data of allele frequencies with forward-in-time simulations
#' @author Zhangyi He, Xiaoyang Dai, Mark Beaumont and Feng Yu

#' version 1.1
#' Two-loucs case (N/A is allowed)

#install.packages("RColorBrewer")
library("RColorBrewer")

#install.packages("viridis")
library("viridis")

#install.packages("ggplot2")  
library("ggplot2")

#install.packages("plot3D")
library("plot3D")

# call R functions
source("./RFUN_2L.R")

################################################################################

#' Simulate the haplotype frequency trajectories according to the two-locus Wright-Fisher model with selection
#' Parameter setting
#' @param sel_cof the selection coefficients at loci A and B
#' @param dom_par the dominance parameters at loci A and B
#' @param rec_rat the recombination rate between loci A and B
#' @param pop_siz the number of the diploid individuals in the population
#' @param int_frq the initial haplotype frequencies of the population
#' @param int_gen the first generation of the simulated haplotype frequency trajectories
#' @param lst_gen the last generation of the simulated haplotype frequency trajectories

sel_cof <- c(1e-02, 5e-03)
dom_par <- c(5e-01, 5e-01)
rec_rat <- 1e-03
pop_siz <- 5e+03
int_frq <- c(1e-01, 2e-01, 3e-01, 4e-01)
int_gen <- 0
lst_gen <- 500

frq_pth <- cmpsimulateTLWFMS(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen)

k <- int_gen:lst_gen
plot(k, frq_pth[1, ], type = "l", lwd = 1.5, 
     xlab = "Generation", ylab = "Haplotype frequency", 
     main = "A frequency trajectory of the A1B1 haplotype generated with the Wright-Fisher model")
plot(k, frq_pth[2, ], type = "l", lwd = 1.5, 
     xlab = "Generation", ylab = "Haplotype frequency", 
     main = "A frequency trajectory of the A1B2 haplotype generated with the Wright-Fisher model")
plot(k, frq_pth[3, ], type = "l", lwd = 1.5, 
     xlab = "Generation", ylab = "Haplotype frequency", 
     main = "A frequency trajectory of the A2B1 haplotype generated with the Wright-Fisher model")
plot(k, frq_pth[4, ], type = "l", lwd = 1.5, 
     xlab = "Generation", ylab = "Haplotype frequency", 
     main = "A frequency trajectory of the A2B2 haplotype generated with the Wright-Fisher model")

########################################

#' Simulate the haplotype frequency trajectories according to the two-locus Wright-Fisher diffusion with selection using the Euler-Maruyama method
#' Parameter setting
#' @param sel_cof the selection coefficients at loci A and B
#' @param dom_par the dominance parameters at loci A and B
#' @param rec_rat the recombination rate between loci A and B
#' @param pop_siz the number of the diploid individuals in the population
#' @param int_frq the initial haplotype frequencies of the population
#' @param int_gen the first generation of the simulated haplotype frequency trajectories
#' @param lst_gen the last generation of the simulated haplotype frequency trajectories
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param data_augmentation = TRUE/FALSE (return the simulated sample trajectory with data augmentation or not)

sel_cof <- c(1e-02, 5e-03)
dom_par <- c(5e-01, 5e-01)
rec_rat <- 1e-03
pop_siz <- 5e+03
int_frq <- c(1e-01, 2e-01, 3e-01, 4e-01)
int_gen <- 0
lst_gen <- 500
ptn_num <- 5e+00

frq_pth <- cmpsimulateTLWFDS(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, ptn_num, data_augmentation = TRUE)

t <- (int_gen:(int_gen + (lst_gen - int_gen) * ptn_num)) / 2 / pop_siz
plot(t, frq_pth[1, ], type = "l", lwd = 1.5, 
     xlab = "Time", ylab = "Haplotype frequency", 
     main = "A frequency trajectory of the A1B1 haplotype generated with the Wright-Fisher diffusion")
plot(t, frq_pth[2, ], type = "l", lwd = 1.5, 
     xlab = "Time", ylab = "Haplotype frequency", 
     main = "A frequency trajectory of the A1B2 haplotype generated with the Wright-Fisher diffusion")
plot(t, frq_pth[3, ], type = "l", lwd = 1.5, 
     xlab = "Time", ylab = "Haplotype frequency", 
     main = "A frequency trajectory of the A2B1 haplotype generated with the Wright-Fisher diffusion")
plot(t, frq_pth[4, ], type = "l", lwd = 1.5, 
     xlab = "Time", ylab = "Haplotype frequency", 
     main = "A frequency trajectory of the A2B2 haplotype generated with the Wright-Fisher diffusion")

########################################

#' Compare the simulation generated with the Wright-Fisher model and the Wright-Fisher diffusion
sel_cof <- c(1e-02, 5e-03)
dom_par <- c(5e-01, 5e-01)
rec_rat <- 1e-03
pop_siz <- 5e+03
int_frq <- c(1e-01, 2e-01, 3e-01, 4e-01)
int_gen <- 0
lst_gen <- 500
ptn_num <- 5e+00
sim_num <- 1e+06

sim_frq_WFM <- matrix(NA, nrow = 4, ncol = sim_num)
sim_frq_WFD <- matrix(NA, nrow = 4, ncol = sim_num)
for (i in 1:sim_num) {
  print(i)
  sim_frq_WFM[, i] <- cmpsimulateTLWFMS(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen)[, (lst_gen - int_gen) + 1]
  sim_frq_WFD[, i] <- cmpsimulateTLWFDS(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, ptn_num, data_augmentation = FALSE)[, (lst_gen - int_gen) + 1]
}

save(sel_cof, dom_par, rec_rat, pop_siz, int_frq, int_gen, lst_gen, ptn_num, sim_num, sim_frq_WFM, sim_frq_WFD, 
     file = "./Output/Output v2.1/Test v1.1/TEST_2L_WFM_vs_WFD.rda")

load("./Output/Output v2.1/Test v1.1/TEST_2L_WFM_vs_WFD.rda")

pdf(file = "./Output/Output v2.1/Test v1.1/TEST_2L_WFM_vs_WFD.pdf", width = 20, height = 10)
par(mfrow = c(2, 2), mar = c(5.5, 5, 5.5, 2.5), oma = c(0, 0, 3, 0), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
hist(sim_frq_WFM[1, ], breaks = seq(min(sim_frq_WFM[1, ], sim_frq_WFD[1, ]), max(sim_frq_WFM[1, ], sim_frq_WFD[1, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5), 
     xlim = c(min(sim_frq_WFM[1, ], sim_frq_WFD[1, ]), max(sim_frq_WFM[1, ], sim_frq_WFD[1, ])), 
     xlab = "Haplotype frequency", main = "Haplotype A1B1")
hist(sim_frq_WFD[1, ], breaks = seq(min(sim_frq_WFM[1, ], sim_frq_WFD[1, ]), max(sim_frq_WFM[1, ], sim_frq_WFD[1, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

hist(sim_frq_WFM[2, ], breaks = seq(min(sim_frq_WFM[2, ], sim_frq_WFD[2, ]), max(sim_frq_WFM[2, ], sim_frq_WFD[2, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5), 
     xlim = c(min(sim_frq_WFM[2, ], sim_frq_WFD[2, ]), max(sim_frq_WFM[2, ], sim_frq_WFD[2, ])), 
     xlab = "Haplotype frequency", main = "Haplotype A1B2")
hist(sim_frq_WFD[2, ], breaks = seq(min(sim_frq_WFM[2, ], sim_frq_WFD[2, ]), max(sim_frq_WFM[2, ], sim_frq_WFD[2, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

hist(sim_frq_WFM[3, ], breaks = seq(min(sim_frq_WFM[3, ], sim_frq_WFD[3, ]), max(sim_frq_WFM[3, ], sim_frq_WFD[3, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5), 
     xlim = c(min(sim_frq_WFM[3, ], sim_frq_WFD[3, ]), max(sim_frq_WFM[3, ], sim_frq_WFD[3, ])), 
     xlab = "Haplotype frequency", main = "Haplotype A2B1")
hist(sim_frq_WFD[3, ], breaks = seq(min(sim_frq_WFM[3, ], sim_frq_WFD[3, ]), max(sim_frq_WFM[3, ], sim_frq_WFD[3, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)

hist(sim_frq_WFM[4, ], breaks = seq(min(sim_frq_WFM[4, ], sim_frq_WFD[4, ]), max(sim_frq_WFM[4, ], sim_frq_WFD[4, ]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5), 
     xlim = c(min(sim_frq_WFM[4, ], sim_frq_WFD[4, ]), max(sim_frq_WFM[4, ], sim_frq_WFD[4, ])), 
     xlab = "Haplotype frequency", main = "Haplotype A2B2")
hist(sim_frq_WFD[4, ], breaks = seq(min(sim_frq_WFM[4, ], sim_frq_WFD[4, ]), max(sim_frq_WFM[4, ], sim_frq_WFD[4, ]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
title(paste("Histograms of the haplotype frequencies in generation", lst_gen, "under the Wright-Fisher model and the Wright-Fisher diffusion"), outer = TRUE)
dev.off()

################################################################################

#' Simulate the hidden Markov model
#' Parameter setting
#' @param model = "WFM"/"WFD" (return the observations from the underlying population evolving according to the WFM or the WFD)
#' @param missing = TRUE/FALSE (return the observations with missing values or not)
#' @param mis_rat the rate of the sampled chromosomes that contain variants with unknown alleles at each locus
#' @param sel_cof the selection coefficients at loci A and B
#' @param dom_par the dominance parameters at loci A and B
#' @param rec_rat the recombination rate between loci A and B
#' @param pop_siz the number of the diploid individuals in the population
#' @param int_frq the initial haplotype frequencies of the population
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the chromosomes drawn from the population at all sampling time points
#' @param ptn_num the number of the subintervals divided per generation in the Euler-Maruyama method for the WFD

#' Simulate the dataset under the Wright-Fisher model 
model <- "WFM"
missing <- TRUE
mis_rat <- 2e-02
sel_cof <- c(1e-02, 5e-03)
dom_par <- c(5e-01, 5e-01)
rec_rat <- 1e-03
pop_siz <- 5e+03
int_frq <- c(1e-01, 2e-01, 3e-01, 4e-01)
smp_gen <- (0:10) * 50
smp_siz <- rep(50, 11)

sim_HMM_WFM <- cmpsimulateHMM(model, missing, mis_rat, sel_cof, dom_par, rec_rat, pop_siz, int_frq, smp_gen, smp_siz)
smp_gen <- sim_HMM_WFM$smp_gen
smp_siz <- sim_HMM_WFM$smp_siz
smp_hap_cnt <- sim_HMM_WFM$smp_hap_cnt
pop_hap_frq <- sim_HMM_WFM$pop_hap_frq
smp_ale_cnt <- sim_HMM_WFM$smp_ale_cnt
pop_ale_frq <- sim_HMM_WFM$pop_ale_frq

k <- min(smp_gen):max(smp_gen)
smp_ale_frq <- smp_ale_cnt %*% diag(1 / smp_siz)
plot(k, pop_ale_frq[1, ], type = 'l', lwd = 1.5, 
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_ale_frq[1, ], pop_ale_frq[1, ]), max(smp_ale_frq[1, ], pop_ale_frq[1, ])), 
     xlab = "Generation", ylab = "Allele frequency", 
     main = "A simulated dataset of the mutant allele at locus A generated with the Wright-Fisher model")
points(smp_gen, smp_ale_frq[1, ], col = 'red', pch = 17, cex = 1)
plot(k, 1 - pop_ale_frq[1, ], type = 'l', lwd = 1.5, 
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_ale_frq[2, ], 1 - pop_ale_frq[1, ]), max(smp_ale_frq[2, ], 1 - pop_ale_frq[1, ])), 
     xlab = "Generation", ylab = "Allele frequency", 
     main = "A simulated dataset of the ancestral allele at locus A generated with the Wright-Fisher model")
points(smp_gen, smp_ale_frq[2, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_ale_frq[2, ], type = 'l', lwd = 1.5, 
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_ale_frq[3, ], pop_ale_frq[2, ]), max(smp_ale_frq[3, ], pop_ale_frq[2, ])), 
     xlab = "Generation", ylab = "Allele frequency", 
     main = "A simulated dataset of the mutant allele at locus B generated with the Wright-Fisher model")
points(smp_gen, smp_ale_frq[3, ], col = 'red', pch = 17, cex = 1)
plot(k, 1 - pop_ale_frq[2, ], type = 'l', lwd = 1.5, 
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_ale_frq[4, ], 1 - pop_ale_frq[2, ]), max(smp_ale_frq[4, ], 1 - pop_ale_frq[2, ])), 
     xlab = "Generation", ylab = "Allele frequency", 
     main = "A simulated dataset of the ancestral allele at locus B generated with the Wright-Fisher model")
points(smp_gen, smp_ale_frq[4, ], col = 'red', pch = 17, cex = 1)

####################

#' Simulate the dataset under the Wright-Fisher diffusion 
model <- "WFD"
missing <- TRUE
mis_rat <- 2e-02
sel_cof <- c(1e-02, 5e-03)
dom_par <- c(5e-01, 5e-01)
rec_rat <- 1e-03
pop_siz <- 5e+03
int_frq <- c(1e-01, 2e-01, 3e-01, 4e-01)
smp_gen <- (0:10) * 50
smp_siz <- rep(50, 11)
ptn_num <- 5e+00

sim_HMM_WFD <- cmpsimulateHMM(model, missing, mis_rat, sel_cof, dom_par, rec_rat, pop_siz, int_frq, smp_gen, smp_siz)
smp_gen <- sim_HMM_WFD$smp_gen
smp_siz <- sim_HMM_WFD$smp_siz
smp_hap_cnt <- sim_HMM_WFD$smp_hap_cnt
pop_hap_frq <- sim_HMM_WFD$pop_hap_frq
smp_ale_cnt <- sim_HMM_WFD$smp_ale_cnt
pop_ale_frq <- sim_HMM_WFD$pop_ale_frq

k <- min(smp_gen):max(smp_gen)
smp_ale_frq <- smp_ale_cnt %*% diag(1 / smp_siz)
plot(k, pop_ale_frq[1, ], type = 'l', lwd = 1.5, 
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_ale_frq[1, ], pop_ale_frq[1, ]), max(smp_ale_frq[1, ], pop_ale_frq[1, ])), 
     xlab = "Generation", ylab = "Allele frequency", 
     main = "A simulated dataset of the mutant allele at locus A generated with the Wright-Fisher diffusion")
points(smp_gen, smp_ale_frq[1, ], col = 'red', pch = 17, cex = 1)
plot(k, 1 - pop_ale_frq[1, ], type = 'l', lwd = 1.5, 
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_ale_frq[2, ], 1 - pop_ale_frq[1, ]), max(smp_ale_frq[2, ], 1 - pop_ale_frq[1, ])), 
     xlab = "Generation", ylab = "Allele frequency", 
     main = "A simulated dataset of the ancestral allele at locus A generated with the Wright-Fisher diffusion")
points(smp_gen, smp_ale_frq[2, ], col = 'red', pch = 17, cex = 1)
plot(k, pop_ale_frq[2, ], type = 'l', lwd = 1.5, 
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_ale_frq[3, ], pop_ale_frq[2, ]), max(smp_ale_frq[3, ], pop_ale_frq[2, ])), 
     xlab = "Generation", ylab = "Allele frequency", 
     main = "A simulated dataset of the mutant allele at locus B generated with the Wright-Fisher diffusion")
points(smp_gen, smp_ale_frq[3, ], col = 'red', pch = 17, cex = 1)
plot(k, 1 - pop_ale_frq[2, ], type = 'l', lwd = 1.5, 
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_ale_frq[4, ], 1 - pop_ale_frq[2, ]), max(smp_ale_frq[4, ], 1 - pop_ale_frq[2, ])), 
     xlab = "Generation", ylab = "Allele frequency", 
     main = "A simulated dataset of the ancestral allele at locus B generated with the Wright-Fisher diffusion")
points(smp_gen, smp_ale_frq[4, ], col = 'red', pch = 17, cex = 1)

################################################################################

#' Generate a simulated dataset under the Wright-Fisher model 
test_seed <- 7
set.seed(test_seed)

model <- "WFM"
missing <- TRUE
mis_rat <- 2e-02
sel_cof <- c(1e-02, 5e-03)
dom_par <- c(5e-01, 5e-01)
rec_rat <- 1e-03
pop_siz <- 5e+03
int_frq <- c(1e-01, 2e-01, 3e-01, 4e-01)
smp_gen <- (0:10) * 50
smp_siz <- rep(50, 11)

SimData <- cmpsimulateHMM(model, missing, mis_rat, sel_cof, dom_par, rec_rat, pop_siz, int_frq, smp_gen, smp_siz)
smp_gen <- SimData$smp_gen
smp_siz <- SimData$smp_siz
smp_hap_cnt <- SimData$smp_hap_cnt
pop_hap_frq <- SimData$pop_hap_frq
smp_ale_cnt <- SimData$smp_ale_cnt
pop_ale_frq <- SimData$pop_ale_frq

save(sel_cof, dom_par, rec_rat, pop_siz, smp_gen, smp_siz, smp_hap_cnt, pop_hap_frq, smp_ale_cnt, pop_ale_frq, 
     file = "./Output/Output v2.1/Test v1.1/TEST_2L_SimData.rda")

load("./Output/Output v2.1/Test v1.1/TEST_2L_SimData.rda")

k <- min(smp_gen):max(smp_gen)
smp_ale_frq <- smp_ale_cnt %*% diag(1 / smp_siz)

pdf(file = "./Output/Output v2.1/Test v1.1/TEST_2L_SimData.pdf", width = 20, height = 20)
par(mfrow = c(2, 2), oma = c(0, 0, 3, 0), mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
plot(k, pop_ale_frq[1, ], type = 'l', lwd = 1.5, 
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_ale_frq[1, ], pop_ale_frq[1, ]), max(smp_ale_frq[1, ], pop_ale_frq[1, ])), 
     xlab = "Generation", ylab = "Allele frequency", 
     main = "Mutant allele at locus A")
points(smp_gen, smp_ale_frq[1, ], col = 'red', pch = 17, cex = 1)

plot(k, 1 - pop_ale_frq[1, ], type = 'l', lwd = 1.5, 
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_ale_frq[2, ], 1 - pop_ale_frq[1, ]), max(smp_ale_frq[2, ], 1 - pop_ale_frq[1, ])), 
     xlab = "Generation", ylab = "Allele frequency", 
     main = "Ancestral allele at locus A")
points(smp_gen, smp_ale_frq[2, ], col = 'red', pch = 17, cex = 1)

plot(k, pop_ale_frq[2, ], type = 'l', lwd = 1.5, 
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_ale_frq[3, ], pop_ale_frq[2, ]), max(smp_ale_frq[3, ], pop_ale_frq[2, ])), 
     xlab = "Generation", ylab = "Allele frequency", 
     main = "Mutant allele at locus B")
points(smp_gen, smp_ale_frq[3, ], col = 'red', pch = 17, cex = 1)

plot(k, 1 - pop_ale_frq[2, ], type = 'l', lwd = 1.5, 
     xlim = c(min(smp_gen), max(smp_gen)), ylim = c(min(smp_ale_frq[4, ], 1 - pop_ale_frq[2, ]), max(smp_ale_frq[4, ], 1 - pop_ale_frq[2, ])), 
     xlab = "Generation", ylab = "Allele frequency", 
     main = "Ancestral allele at locus B")
points(smp_gen, smp_ale_frq[4, ], col = 'red', pch = 17, cex = 1)
title("A simulated dataset with missing values generated with the Wright-Fisher model", outer = TRUE)
dev.off()

########################################

#' Run the bootstrap particle filter (BPF) with the two-locus Wright-Fisher diffusion with selection
#' Parameter setting
#' @param sel_cof the selection coefficients at loci A and B
#' @param dom_par the dominance parameters at loci A and B
#' @param rec_rat the recombination rate between loci A and B
#' @param pop_siz the number of the diploid individuals in the population
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the chromosomes drawn from the population at all sampling time points
#' @param smp_cnt the count of the mutant alleles and ancestral alleles observed in the sample at all sampling time points
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter

load("./Output/Output v2.1/Test v1.1/TEST_2L_SimData.rda")

set.seed(test_seed)

sel_cof
dom_par
rec_rat
pop_siz
smp_gen
smp_siz
smp_ale_cnt
ptn_num <- 5e+00
pcl_num <- 5e+04

system.time(BPF <- cmprunBPF(sel_cof, dom_par, rec_rat, pop_siz, smp_gen, smp_siz, smp_ale_cnt, ptn_num, pcl_num))

save(sel_cof, dom_par, rec_rat, pop_siz, smp_gen, smp_siz, smp_ale_cnt, ptn_num, pcl_num, BPF, 
     file = "./Output/Output v2.1/Test v1.1/TEST_2L_BPF.rda")

load("./Output/Output v2.1/Test v1.1/TEST_2L_BPF.rda")

lik <- rep(1, pcl_num)
wght <- BPF$wght
for (k in 1:length(smp_gen)) {
  lik <- lik * (cumsum(wght[, k]) / (1:pcl_num))
}

pdf(file = "./Output/Output v2.1/Test v1.1/TEST_2L_BPF_Likelihood.pdf", width = 10, height = 10)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(1:pcl_num, log(lik), type = 'l', 
     xlab = "Number of particles", ylab = "Log likelihood", main = "Log likelihood estimated with the bootstrap particle filter")
dev.off()

smp_hap_frq <- smp_hap_cnt %*% diag(1 / smp_siz)
pop_hap_frq_pre_resmp <- BPF$pop_frq_pre_resmp
pop_hap_frq_pst_resmp <- BPF$pop_frq_pst_resmp

pdf(file = "./Output/Output v2.1/Test v1.1/TEST_2L_BPF_Particle.pdf", width = 20, height = 55)
par(mfrow = c(11, 4), oma = c(0, 0, 3, 0), mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
for (k in 1:length(smp_gen)) {
  hist_pst_resmp <- hist(pop_hap_frq_pst_resmp[1, , k], breaks = seq(min(pop_hap_frq_pst_resmp[1, , k], pop_hap_frq_pre_resmp[1, , k]), max(pop_hap_frq_pst_resmp[1, , k], pop_hap_frq_pre_resmp[1, , k]), length.out = 50), plot = FALSE)
  hist_pre_resmp <- hist(pop_hap_frq_pre_resmp[1, , k], breaks = seq(min(pop_hap_frq_pst_resmp[1, , k], pop_hap_frq_pre_resmp[1, , k]), max(pop_hap_frq_pst_resmp[1, , k], pop_hap_frq_pre_resmp[1, , k]), length.out = 50), plot = FALSE)
  hist(pop_hap_frq_pst_resmp[1, , k], breaks = seq(min(pop_hap_frq_pst_resmp[1, , k], pop_hap_frq_pre_resmp[1, , k]), max(pop_hap_frq_pst_resmp[1, , k], pop_hap_frq_pre_resmp[1, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5), 
       xlim = c(min(pop_hap_frq_pst_resmp[1, , k], pop_hap_frq_pre_resmp[1, , k], smp_hap_frq[1, k]), max(pop_hap_frq_pst_resmp[1, , k], pop_hap_frq_pre_resmp[1, , k], smp_hap_frq[1, k])), ylim = c(0, max(hist_pst_resmp$density, hist_pre_resmp$density)), 
       xlab = "Haplotype frequency", main = paste("Haplotype A1B1 in generation", smp_gen[k]))
  hist(pop_hap_frq_pre_resmp[1, , k], breaks = seq(min(pop_hap_frq_pst_resmp[1, , k], pop_hap_frq_pre_resmp[1, , k]), max(pop_hap_frq_pst_resmp[1, , k], pop_hap_frq_pre_resmp[1, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_hap_frq[1, k], col = 'red', lty = 2, lwd = 2)
  
  hist_pst_resmp <- hist(pop_hap_frq_pst_resmp[2, , k], breaks = seq(min(pop_hap_frq_pst_resmp[2, , k], pop_hap_frq_pre_resmp[2, , k]), max(pop_hap_frq_pst_resmp[2, , k], pop_hap_frq_pre_resmp[2, , k]), length.out = 50), plot = FALSE)
  hist_pre_resmp <- hist(pop_hap_frq_pre_resmp[2, , k], breaks = seq(min(pop_hap_frq_pst_resmp[2, , k], pop_hap_frq_pre_resmp[2, , k]), max(pop_hap_frq_pst_resmp[2, , k], pop_hap_frq_pre_resmp[2, , k]), length.out = 50), plot = FALSE)  
  hist(pop_hap_frq_pst_resmp[2, , k], breaks = seq(min(pop_hap_frq_pst_resmp[2, , k], pop_hap_frq_pre_resmp[2, , k]), max(pop_hap_frq_pst_resmp[2, , k], pop_hap_frq_pre_resmp[2, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5), 
       xlim = c(min(pop_hap_frq_pst_resmp[2, , k], pop_hap_frq_pre_resmp[2, , k], smp_hap_frq[2, k]), max(pop_hap_frq_pst_resmp[2, , k], pop_hap_frq_pre_resmp[2, , k], smp_hap_frq[2, k])), ylim = c(0, max(hist_pst_resmp$density, hist_pre_resmp$density)), 
       xlab = "Haplotype frequency", main = paste("Haplotype A1B2 in generation", smp_gen[k]))
  hist(pop_hap_frq_pre_resmp[2, , k], breaks = seq(min(pop_hap_frq_pst_resmp[2, , k], pop_hap_frq_pre_resmp[2, , k]), max(pop_hap_frq_pst_resmp[2, , k], pop_hap_frq_pre_resmp[2, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_hap_frq[2, k], col = 'red', lty = 2, lwd = 2)
  
  hist_pst_resmp <- hist(pop_hap_frq_pst_resmp[3, , k], breaks = seq(min(pop_hap_frq_pst_resmp[3, , k], pop_hap_frq_pre_resmp[3, , k]), max(pop_hap_frq_pst_resmp[3, , k], pop_hap_frq_pre_resmp[3, , k]), length.out = 50), plot = FALSE)
  hist_pre_resmp <- hist(pop_hap_frq_pre_resmp[3, , k], breaks = seq(min(pop_hap_frq_pst_resmp[3, , k], pop_hap_frq_pre_resmp[3, , k]), max(pop_hap_frq_pst_resmp[3, , k], pop_hap_frq_pre_resmp[3, , k]), length.out = 50), plot = FALSE)
  hist(pop_hap_frq_pst_resmp[3, , k], breaks = seq(min(pop_hap_frq_pst_resmp[3, , k], pop_hap_frq_pre_resmp[3, , k]), max(pop_hap_frq_pst_resmp[3, , k], pop_hap_frq_pre_resmp[3, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5), 
       xlim = c(min(pop_hap_frq_pst_resmp[3, , k], pop_hap_frq_pre_resmp[3, , k], smp_hap_frq[3, k]), max(pop_hap_frq_pst_resmp[3, , k], pop_hap_frq_pre_resmp[3, , k], smp_hap_frq[3, k])), ylim = c(0, max(hist_pst_resmp$density, hist_pre_resmp$density)), 
       xlab = "Haplotype frequency", main = paste("Haplotype A2B1 in generation", smp_gen[k]))
  hist(pop_hap_frq_pre_resmp[3, , k], breaks = seq(min(pop_hap_frq_pst_resmp[3, , k], pop_hap_frq_pre_resmp[3, , k]), max(pop_hap_frq_pst_resmp[3, , k], pop_hap_frq_pre_resmp[3, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_hap_frq[3, k], col = 'red', lty = 2, lwd = 2)
  
  hist_pst_resmp <- hist(pop_hap_frq_pst_resmp[4, , k], breaks = seq(min(pop_hap_frq_pst_resmp[4, , k], pop_hap_frq_pre_resmp[4, , k]), max(pop_hap_frq_pst_resmp[4, , k], pop_hap_frq_pre_resmp[4, , k]), length.out = 50), plot = FALSE)
  hist_pre_resmp <- hist(pop_hap_frq_pre_resmp[4, , k], breaks = seq(min(pop_hap_frq_pst_resmp[4, , k], pop_hap_frq_pre_resmp[4, , k]), max(pop_hap_frq_pst_resmp[4, , k], pop_hap_frq_pre_resmp[4, , k]), length.out = 50), plot = FALSE)
  hist(pop_hap_frq_pst_resmp[4, , k], breaks = seq(min(pop_hap_frq_pst_resmp[4, , k], pop_hap_frq_pre_resmp[4, , k]), max(pop_hap_frq_pst_resmp[4, , k], pop_hap_frq_pre_resmp[4, , k]), length.out = 50), freq = FALSE, col = rgb(0.1, 0.1, 0.1, 0.5), 
       xlim = c(min(pop_hap_frq_pst_resmp[4, , k], pop_hap_frq_pre_resmp[4, , k], smp_hap_frq[4, k]), max(pop_hap_frq_pst_resmp[4, , k], pop_hap_frq_pre_resmp[4, , k], smp_hap_frq[4, k])), ylim = c(0, max(hist_pst_resmp$density, hist_pre_resmp$density)), 
       xlab = "Haplotype frequency", main = paste("Haplotype A2B2 in generation", smp_gen[k]))
  hist(pop_hap_frq_pre_resmp[4, , k], breaks = seq(min(pop_hap_frq_pst_resmp[4, , k], pop_hap_frq_pre_resmp[4, , k]), max(pop_hap_frq_pst_resmp[4, , k], pop_hap_frq_pre_resmp[4, , k]), length.out = 50), freq = FALSE, col = rgb(0.8, 0.8, 0.8, 0.5), add = TRUE)
  abline(v = smp_hap_frq[4, k], col = 'red', lty = 2, lwd = 2)
}
title("Histograms of the pre- and post-resampling particles", outer = TRUE)
dev.off()

########################################

#' Calculate the optimal particle number in the particle marginal Metropolis-Hastings (PMMH)
#' Parameter settings
#' @param sel_cof the selection coefficients at loci A and B
#' @param dom_par the dominance parameters at loci A and B
#' @param rec_rat the recombination rate between loci A and B
#' @param pop_siz the number of the diploid individuals in the population
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the chromosomes drawn from the population at all sampling time points
#' @param smp_cnt the count of the mutant alleles and ancestral alleles observed in the sample at all sampling time points
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param gap_num the number of particles increased or decreased in the optimal particle number search

load("./Output/Output v2.1/Test v1.1/TEST_2L_SimData.rda")

set.seed(test_seed)

sel_cof
dom_par
rec_rat
pop_siz
smp_gen
smp_siz
smp_ale_cnt
ptn_num <- 5e+00
pcl_num <- 1e+03
gap_num <- 1e+02

system.time(OptNum <- calculateOptimalParticleNum(sel_cof, dom_par, rec_rat, pop_siz, smp_gen, smp_siz, smp_ale_cnt, ptn_num, pcl_num, gap_num))

save(sel_cof, dom_par, rec_rat, pop_siz, smp_gen, smp_siz, smp_ale_cnt, ptn_num, pcl_num, gap_num, OptNum,
     file = "./Output/Output v2.1/Test v1.1/TEST_2L_OptNum.rda")

load("./Output/Output v2.1/Test v1.1/TEST_2L_OptNum.rda")

opt_pcl_num <- OptNum$opt_pcl_num
log_lik_sdv <- OptNum$log_lik_sdv

pdf(file = "./Output/Output v2.1/Test v1.1/TEST_2L_OptNum.pdf", width = 12, height = 9)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 1.75, cex.sub = 1.5, cex.axis = 1.5, cex.lab = 1.5)
plot(opt_pcl_num, log_lik_sdv, type = 'b', lwd = 2, 
     xlab = "Particle number", ylab = "Log-likelihood standard deviation", main = "Optimal particle number in the PMMH")
abline(h = 1.7, col = 'red', lty = 2, lwd = 2)
abline(h = 1.0, col = 'red', lty = 2, lwd = 2)
dev.off()

########################################

#' Run the particle marginal Metropolis-Hastings (PMMH)
#' Parameter settings
#' @param sel_cof the selection coefficients at loci A and B
#' @param dom_par the dominance parameters at loci A and B
#' @param rec_rat the recombination rate between loci A and B
#' @param pop_siz the number of the diploid individuals in the population
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the chromosomes drawn from the population at all sampling time points
#' @param smp_cnt the count of the mutant alleles and ancestral alleles observed in the sample at all sampling time points
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param itn_num the number of the iterations carried out in the particle marginal Metropolis-Hastings

load("./Output/Output v2.1/Test v1.1/TEST_2L_SimData.rda")

set.seed(test_seed)

sel_cof <- c(0e+00, 0e+00)
dom_par <- c(5e-01, 5e-01)
rec_rat <- 1e-03
pop_siz <- 5e+03
smp_gen
smp_siz
smp_ale_cnt
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 5e+04

system.time(PMMH <- cmprunPMMH(sel_cof, dom_par, rec_rat, pop_siz, smp_gen, smp_siz, smp_ale_cnt, ptn_num, pcl_num, itn_num))

save(sel_cof, dom_par, rec_rat, pop_siz, smp_gen, smp_siz, smp_ale_cnt, ptn_num, pcl_num, itn_num, PMMH, 
     file = "./Output/Output v2.1/Test v1.1/TEST_2L_PMMH.rda")

load("./Output/Output v2.1/Test v1.1/TEST_2L_PMMH.rda")

load("./Output/Output v2.1/Test v1.1/TEST_2L_SimData.rda")

sel_cof_A_chn <- PMMH$sel_cof_A_chn
sel_cof_B_chn <- PMMH$sel_cof_B_chn
pdf(file = "./Output/Output v2.1/Test v1.1/TEST_2L_PMMH_Traceplot.pdf", width = 20, height = 10)
par(mfrow = c(2, 1), mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
plot(1:itn_num, sel_cof_A_chn[1:itn_num], type = 'l', 
     xlab = "Iteration", ylab = "Selection coefficient", main = "Trace plot of the selection coefficient at locus A")
abline(h = sel_cof[1], col = 'red', lty = 2, lwd = 2)
plot(1:itn_num, sel_cof_B_chn[1:itn_num], type = 'l', 
     xlab = "Iteration", ylab = "Selection coefficient", main = "Trace plot of the selection coefficient at locus B")
abline(h = sel_cof[2], col = 'red', lty = 2, lwd = 2)
dev.off()

brn_num <- 1e+04
sel_cof_A_chn <- sel_cof_A_chn[brn_num:length(sel_cof_A_chn)]
sel_cof_B_chn <- sel_cof_B_chn[brn_num:length(sel_cof_B_chn)]

thn_num <- 8e+00
sel_cof_A_chn <- sel_cof_A_chn[(1:round(length(sel_cof_A_chn) / thn_num)) * thn_num]
sel_cof_B_chn <- sel_cof_B_chn[(1:round(length(sel_cof_B_chn) / thn_num)) * thn_num]

grd_num <- 1e+03
sel_cof_pdf <- kde2d(sel_cof_A_chn, sel_cof_B_chn, n = grd_num)
pdf(file = "./Output/Output v2.1/Test v1.1/TEST_2L_PMMH_Posterior.pdf", width = 10, height = 10)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32), 
      xlab = "Selection coefficient at locus A", ylab = "Selection coefficient at locus B", 
      main = "Posterior for the selection coefficients at loci A and B")
abline(v = sel_cof[1], col = 'red', lty = 2, lwd = 2)
abline(h = sel_cof[2], col = 'red', lty = 2, lwd = 2)
dev.off()

########################################

#' Run the Bayesian procedure for the inference of natural selection
#' Parameter settings
#' @param sel_cof the selection coefficients at loci A and B
#' @param dom_par the dominance parameters at loci A and B
#' @param rec_rat the recombination rate between loci A and B
#' @param pop_siz the number of the diploid individuals in the population
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_siz the count of the chromosomes drawn from the population at all sampling time points
#' @param smp_cnt the count of the mutant alleles and ancestral alleles observed in the sample at all sampling time points
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter
#' @param itn_num the number of the iterations carried out in the particle marginal Metropolis-Hastings
#' @param brn_num the number of the iterations for burn-in
#' @param thn_num the number of the iterations for thinning
#' @param grd_num the number of the grids in the kernel density estimation

load("./Output/Output v2.1/Test v1.1/TEST_2L_SimData.rda")

set.seed(test_seed)

sel_cof <- c(0e+00, 0e+00)
dom_par <- c(5e-01, 5e-01)
rec_rat <- 1e-03
pop_siz <- 5e+03
smp_gen
smp_siz
smp_ale_cnt
ptn_num <- 5e+00
pcl_num <- 1e+03
itn_num <- 5e+04
brn_num <- 1e+04
thn_num <- 8e+00
grd_num <- 1e+03

system.time(BayesianProcedure <- cmprunBayesianProcedure(sel_cof, dom_par, rec_rat, pop_siz, smp_gen, smp_siz, smp_ale_cnt, ptn_num, pcl_num, itn_num, brn_num, thn_num, grd_num))

save(sel_cof, dom_par, rec_rat, pop_siz, smp_gen, smp_siz, smp_ale_cnt, ptn_num, pcl_num, itn_num, brn_num, thn_num, grd_num, BayesianProcedure, 
     file = "./Output/Output v2.1/Test v1.1/TEST_2L_BayesProc.rda")

load("./Output/Output v2.1/Test v1.1/TEST_2L_BayesProc.rda")

load("./Output/Output v2.1/Test v1.1/TEST_2L_SimData.rda")

sel_cof_A_chn <- BayesianProcedure$sel_cof_A_chn
sel_cof_B_chn <- BayesianProcedure$sel_cof_B_chn

sel_cof_pdf <- BayesianProcedure$sel_cof_pdf

sel_cof_A_map <- BayesianProcedure$sel_cof_A_map
sel_cof_B_map <- BayesianProcedure$sel_cof_B_map

sel_cof_A_mmse <- BayesianProcedure$sel_cof_A_mmse
sel_cof_B_mmse <- BayesianProcedure$sel_cof_B_mmse

sel_cof_B_hpd <- BayesianProcedure$sel_cof_B_hpd
sel_cof_A_hpd <- BayesianProcedure$sel_cof_A_hpd

pdf(file = "./Output/Output v2.1/Test v1.1/TEST_2L_BayesProc_Posterior.pdf", width = 20, height = 10)
par(mar = c(5.5, 5, 5.5, 2.5), cex.main = 2, cex.sub = 1.75, cex.axis = 1.75, cex.lab = 1.75)
layout(matrix(c(1, 1, 2, 3), nrow = 2, ncol = 2))
image(sel_cof_pdf, col = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32), 
      xlab = "Selection coefficient at locus A", ylab = "Selection coefficient at locus B", 
      main = "Joint posterior for the selection coefficients at loci A and B")
abline(v = sel_cof[1], col = 'red', lty = 2, lwd = 2)
abline(h = sel_cof[2], col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_A_map, col = 'black', lty = 4, lwd = 2)
abline(h = sel_cof_B_map, col = 'black', lty = 4, lwd = 2)
abline(v = sel_cof_A_mmse, col = 'black', lty = 2, lwd = 2)
abline(h = sel_cof_B_mmse, col = 'black', lty = 2, lwd = 2)

hist(sel_cof_A_chn, breaks = seq(min(sel_cof_A_chn), max(sel_cof_A_chn), length.out = 50), freq = FALSE, 
     xlab = "Selection coefficient", main = "Marginal posterior for the selection coefficients at locus A")
lines(density(sel_cof_A_chn), lwd = 2, col = 'black')
abline(v = sel_cof[1], col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_A_map, col = 'black', lty = 4, lwd = 2)
abline(v = sel_cof_A_mmse, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_A_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_A_hpd[2], col = 'blue', lty = 2, lwd = 2)

hist(sel_cof_B_chn, breaks = seq(min(sel_cof_B_chn), max(sel_cof_B_chn), length.out = 50), freq = FALSE, 
     xlab = "Selection coefficient", main = "Marginal posterior for the selection coefficients at locus B")
lines(density(sel_cof_B_chn), lwd = 2, col = 'black')
abline(v = sel_cof[2], col = 'red', lty = 2, lwd = 2)
abline(v = sel_cof_B_map, col = 'black', lty = 4, lwd = 2)
abline(v = sel_cof_B_mmse, col = 'black', lty = 2, lwd = 2)
abline(v = sel_cof_B_hpd[1], col = 'blue', lty = 2, lwd = 2)
abline(v = sel_cof_B_hpd[2], col = 'blue', lty = 2, lwd = 2)
dev.off()

################################################################################
