#' @title Detecting and quantifying natural selection at two linked loci from time series data of allele frequencies
#' @author Zhangyi He, Xiaoyang Dai, Mark Beaumont and Feng Yu

#' version 1.0

#' R functions

#install.packages("MASS")
library("MASS")

#install.packages("coda")
library("coda")

#install.packages("inline")
library("inline")
#install.packages("Rcpp")
library("Rcpp")
#install.packages("RcppArmadillo")
library("RcppArmadillo")

#install.packages("compiler")
library("compiler")
#enableJIT(1)

# call C++ functions
sourceCpp("1L_CFUN.cpp")

################################################################################

#' Simulate the mutant allele frequency trajectory according to the one-locus Wright-Fisher model with selection
#' Parameter setting
#' @param sel_cof the selection coefficient
#' @param dom_par the dominance parameter
#' @param pop_siz the number of the diploid individuals in the population
#' @param int_frq the initial mutant allele frequency of the population
#' @param int_gen the first generation of the simulated mutant allele frequency trajectory
#' @param lst_gen the last generation of the simulated mutant allele frequency trajectory

#' Standard version
simulateOLWFMS <- function(sel_cof, dom_par, pop_siz, int_frq, int_gen, lst_gen) {
  frq_pth <- simulateOLWFMS_arma(sel_cof, dom_par, pop_siz, int_frq, int_gen, lst_gen)
  frq_pth <- as.vector(frq_pth)
  
  return(frq_pth)
}
#' Compiled version
cmpsimulateOLWFMS <- cmpfun(simulateOLWFMS)

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

#' Standard version
simulateOLWFDS <- function(sel_cof, dom_par, pop_siz, int_frq, int_gen, lst_gen, ptn_num, data_augmentation = TRUE) {
  frq_pth <- simulateOLWFDS_arma(sel_cof, dom_par, pop_siz, int_frq, int_gen, lst_gen, ptn_num)
  frq_pth <- as.vector(frq_pth)
  
  if (data_augmentation == FALSE) {
    return(frq_pth[(0:(lst_gen - int_gen)) * ptn_num + 1])
  } else {
    return(frq_pth)
  }  
}
#' Compiled version
cmpsimulateOLWFDS <- cmpfun(simulateOLWFDS)

########################################

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

#' Standard version
simulateHMM <- function(model, sel_cof, dom_par, pop_siz, int_frq, smp_gen, smp_siz, ...) {
  int_gen <- min(smp_gen)
  lst_gen <- max(smp_gen)
  
  # generate the population mutant allele frequency trajectory
  if (model == "WFM") {
    pop_frq <- cmpsimulateOLWFMS(sel_cof, dom_par, pop_siz, int_frq, int_gen, lst_gen)
  }
  if (model == "WFD") {
    pop_frq <- cmpsimulateOLWFDS(sel_cof, dom_par, pop_siz, int_frq, int_gen, lst_gen, ptn_num, data_augmentation = FALSE)
  }
  
  # generate the sample allele counts at all sampling time points
  smp_cnt <- numeric(length(smp_gen))
  for (k in 1:length(smp_gen)) {
    smp_cnt[k] <- rbinom(1, size = smp_siz[k], prob = pop_frq[smp_gen[k] - int_gen + 1])
  }
  
  return(list(smp_gen = smp_gen, 
              smp_siz = smp_siz, 
              smp_cnt = smp_cnt, 
              pop_frq = pop_frq))
}
#' Compiled version
cmpsimulateHMM <- cmpfun(simulateHMM)

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

#' Standard version
runBPF <- function(sel_cof, dom_par, pop_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num) {
  # run the BPF
  BPF <- runBPF_arma(sel_cof, dom_par, pop_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num)
  
  return(list(lik = BPF$lik, 
              wght = BPF$wght, 
              pop_frq_pre_resmp = BPF$part_pre_resmp, 
              pop_frq_pst_resmp = BPF$part_pst_resmp))
}
#' Compiled version
cmprunBPF <- cmpfun(runBPF)

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

#' Standard version
calculateOptimalParticleNum <- function(sel_cof, dom_par, pop_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, gap_num) {
  OptNum <- calculateOptimalParticleNum_arma(sel_cof, dom_par, pop_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, gap_num)
  
  return(list(opt_pcl_num = as.vector(OptNum$opt_pcl_num), 
              log_lik_sdv = as.vector(OptNum$log_lik_sdv)))
}
#' Compiled version
cmpcalculateOptimalParticleNum <- cmpfun(calculateOptimalParticleNum)

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

#' Standard version
runPMMH <- function(sel_cof, dom_par, pop_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num) {
  # run the PMMH
  PMMH <- runPMMH_arma(sel_cof, dom_par, pop_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num)
  
  return(list(sel_cof_chn = as.vector(PMMH$sel_cof_chn)))
}
#' Compiled version
cmprunPMMH <- cmpfun(runPMMH)

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

#' Standard version
runBayesianProcedure <- function(sel_cof, dom_par, pop_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num, brn_num, thn_num, grd_num) {
  # run the PMMH
  PMMH <- runPMMH_arma(sel_cof, dom_par, pop_siz, smp_gen, smp_siz, smp_cnt, ptn_num, pcl_num, itn_num)
  
  # burn-in and thinning
  sel_cof_chn <- as.vector(PMMH$sel_cof_chn)
  sel_cof_chn <- sel_cof_chn[brn_num:length(sel_cof_chn)]
  sel_cof_chn <- sel_cof_chn[(1:round(length(sel_cof_chn) / thn_num)) * thn_num]
  
  # MAP estimate for the selection coefficient
  if (length(sel_cof_chn) < 1e+05) {
    sel_cof_pdf <- density(sel_cof_chn, n = grd_num)
    sel_cof_grd <- sel_cof_pdf$x
  } else {
    sel_cof_pdf <- density(tail(sel_cof_chn, 1e+05), n = grd_num)
    sel_cof_grd <- sel_cof_pdf$x
  }
  sel_cof_map <- sel_cof_grd[which(sel_cof_pdf$y == max(sel_cof_pdf$y))]
  
  # MMSE estimate for the selection coefficient
  sel_cof_mmse <- mean(sel_cof_chn)
  
  # 95% HPD interval for the selection coefficient
  sel_cof_hpd <- HPDinterval(as.mcmc(sel_cof_chn), prob = 0.95)
  
  return(list(sel_cof_pdf = sel_cof_pdf, 
              sel_cof_map = sel_cof_map, 
              sel_cof_mmse = sel_cof_mmse, 
              sel_cof_hpd = sel_cof_hpd, 
              sel_cof_chn = sel_cof_chn))
}
#' Compiled version
cmprunBayesianProcedure <- cmpfun(runBayesianProcedure)

################################################################################
