#' @title An MCMC-based method for Bayesian inference of natural selection from time series DNA data across linked loci
#' @author Zhangyi He, Xiaoyang Dai, Mark Beaumont and Feng Yu

#' this version is unable to handle missing values in DNA data

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
sourceCpp("2L_cfun_without_NA.cpp")

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

#' Standard version
simulateTLWFMS <- function(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_hap_frq, int_gen, lst_gen) {
  hap_frq_pth <- simulateTLWFMS_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_hap_frq, int_gen, lst_gen)
  
  return(hap_frq_pth)
}
#' Compiled version
cmpsimulateTLWFMS <- cmpfun(simulateTLWFMS)

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

#' Standard version
simulateTLWFDS <- function(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_hap_frq, int_gen, lst_gen, ptn_num, data_augmentation = TRUE) {
  hap_frq_pth <- simulateTLWFDS_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_hap_frq, int_gen, lst_gen, ptn_num)
  
  # return the simulated sample trajectories without data augmentation
  if (data_augmentation == FALSE) {
    hap_frq_pth <- hap_frq_pth[, (0:(lst_gen - int_gen)) * ptn_num + 1]
  }
  
  return(hap_frq_pth)
}
#' Compiled version
cmpsimulateTLWFDS <- cmpfun(simulateTLWFDS)

########################################

#' Simulate the hidden Markov model
#' Parameter setting
#' @param model = "WFM"/"WFD" (return the observations from the underlying population evolving according to the WFM or the WFD)
#' @param phased = TRUE/FALSE (return the time serial sample with phased or unphased chromosomes)
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

#' Standard version
simulateHMM <- function(model, phased, sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_hap_frq, smp_gen, smp_chr_cnt, ...) {
  int_gen <- min(smp_gen)
  lst_gen <- max(smp_gen)
  
  if (phased == FALSE) {
    # generate the population haplotype frequency trajectories and the population allele frequency trajectories
    if (model == "WFM") {
      pop_hap_frq <- cmpsimulateTLWFMS(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_hap_frq, int_gen, lst_gen)
    }
    if (model == "WFD") {
      pop_hap_frq <- cmpsimulateTLWFDS(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_hap_frq, int_gen, lst_gen, ptn_num, data_augmentation = FALSE)
    }
    pop_ale_frq <- matrix(NA, nrow = 2, ncol = (lst_gen - int_gen) + 1)
    pop_ale_frq[1, ] <- pop_hap_frq[1, ] + pop_hap_frq[2, ]
    pop_ale_frq[2, ] <- pop_hap_frq[1, ] + pop_hap_frq[3, ]
    
    # generate the sample haplotype counts and the sample allele counts at all sampling time points
    smp_hap_cnt <- matrix(NA, nrow = 4, ncol = length(smp_gen))
    smp_ale_cnt <- matrix(NA, nrow = 2, ncol = length(smp_gen))
    for (k in 1:length(smp_gen)) {
      smp_hap_cnt[, k] <- rmultinom(1, size = smp_chr_cnt[k], prob = pop_hap_frq[, smp_gen[k] - int_gen + 1])
      smp_ale_cnt[1, k] <- smp_hap_cnt[1, k] + smp_hap_cnt[2, k]
      smp_ale_cnt[2, k] <- smp_hap_cnt[1, k] + smp_hap_cnt[3, k]
    }
    
    return(list(smp_gen = smp_gen, 
                smp_chr_cnt = smp_chr_cnt, 
                smp_hap_cnt = smp_hap_cnt,
                smp_ale_cnt = smp_ale_cnt, 
                pop_hap_frq = pop_hap_frq, 
                pop_ale_frq = pop_ale_frq))
  } else {
    # generate the population haplotype frequency trajectories
    if (model == "WFM") {
      pop_hap_frq <- cmpsimulateTLWFMS(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_hap_frq, int_gen, lst_gen)
    }
    if (model == "WFD") {
      pop_hap_frq <- cmpsimulateTLWFDS(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, int_hap_frq, int_gen, lst_gen, ptn_num, data_augmentation = FALSE)
    }
    
    # generate the sample haplotype counts at all sampling time points
    smp_hap_cnt <- matrix(NA, nrow = 4, ncol = length(smp_gen))
    for (k in 1:length(smp_gen)) {
      smp_hap_cnt[, k] <- rmultinom(1, size = smp_chr_cnt[k], prob = pop_hap_frq[, smp_gen[k] - int_gen + 1])
    }
    
    return(list(smp_gen = smp_gen, 
                smp_chr_cnt = smp_chr_cnt, 
                smp_hap_cnt = smp_hap_cnt,
                pop_hap_frq = pop_hap_frq))
  }
}
#' Compiled version
cmpsimulateHMM <- cmpfun(simulateHMM)

########################################

#' Run the bootstrap particle filter (BPF) with the two-locus Wright-Fisher diffusion with selection
#' Parameter setting
#' @param gen_par the population genetic quantities of interest i.e., gen_par = (s_A, h_A, s_B, h_B, r, N)
#' @param smp_gen the sampling time points measured in one generation
#' @param smp_chr_cnt the count of the chromosomes drawn from the population at all sampling time points
#' @param smp_cnt the count of the two mutant alleles or the four haplotypes observed in the sample at all sampling time points
#' @param ptn_num the number of subintervals divided per generation in the Euler-Maruyama method
#' @param pcl_num the number of particles generated in the bootstrap particle filter

#' Standard version
runBPF <- function(gen_par, smp_gen, smp_chr_cnt, smp_cnt, ptn_num, pcl_num) {
  sel_cof_A <- gen_par[1]
  dom_par_A <- gen_par[2]
  sel_cof_B <- gen_par[3]
  dom_par_B <- gen_par[4]
  rec_rat <- gen_par[5]
  pop_siz <- gen_par[6]
  
  # run the BPF
  phased <- ifelse(nrow(smp_cnt) == 4, TRUE, FALSE)
  if (phased) {
    BPF <- runBPF_PhasedChr_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, smp_gen, smp_chr_cnt, smp_cnt, ptn_num, pcl_num)
    
    lik = BPF$lik
    wght = BPF$wght
    pop_hap_frq_pre_resmp = BPF$part_pre_resmp
    pop_hap_frq_pst_resmp = BPF$part_pst_resmp
    
    return(list(lik = lik, 
                wght = wght, 
                pop_frq_pre_resmp = pop_hap_frq_pre_resmp, 
                pop_frq_pst_resmp = pop_hap_frq_pst_resmp))
  } else {
    BPF <- runBPF_UnphasedChr_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, smp_gen, smp_chr_cnt, smp_cnt, ptn_num, pcl_num)
    
    lik = BPF$lik
    wght = BPF$wght
    pop_hap_frq_pre_resmp = BPF$part_pre_resmp
    pop_hap_frq_pst_resmp = BPF$part_pst_resmp
    
    pop_ale_frq_pre_resmp <- array(NA, dim = c(2, dim(pop_hap_frq_pre_resmp)[2], dim(pop_hap_frq_pre_resmp)[3]))
    pop_ale_frq_pre_resmp[1, , ] <- pop_hap_frq_pre_resmp[1, , ] + pop_hap_frq_pre_resmp[2, , ]
    pop_ale_frq_pre_resmp[2, , ] <- pop_hap_frq_pre_resmp[1, , ] + pop_hap_frq_pre_resmp[3, , ]
    pop_ale_frq_pst_resmp <- array(NA, dim = c(2, dim(pop_hap_frq_pst_resmp)[2], dim(pop_hap_frq_pst_resmp)[3]))
    pop_ale_frq_pst_resmp[1, , ] <- pop_hap_frq_pst_resmp[1, , ] + pop_hap_frq_pst_resmp[2, , ]
    pop_ale_frq_pst_resmp[2, , ] <- pop_hap_frq_pst_resmp[1, , ] + pop_hap_frq_pst_resmp[3, , ]
    
    return(list(lik = lik, 
                wght = wght, 
                pop_frq_pre_resmp = pop_ale_frq_pre_resmp, 
                pop_frq_pst_resmp = pop_ale_frq_pst_resmp))
  }
}
#' Compiled version
cmprunBPF <- cmpfun(runBPF)

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

#' Standard version
runPMMH <- function(int_gen_par, smp_gen, smp_chr_cnt, smp_cnt, ptn_num, pcl_num, itn_num) {
  int_sel_cof_A <- int_gen_par[1]
  int_dom_par_A <- int_gen_par[2]
  int_sel_cof_B <- int_gen_par[3]
  int_dom_par_B <- int_gen_par[4]
  int_rec_rat <- int_gen_par[5]
  int_pop_siz <- int_gen_par[6]
  
  # run the PMMH
  phased <- ifelse(nrow(smp_cnt) == 4, TRUE, FALSE)
  if (phased) {
    PMMH <- runPMMH_PhasedChr_arma(int_sel_cof_A, int_dom_par_A, int_sel_cof_B, int_dom_par_B, int_rec_rat, int_pop_siz, smp_gen, smp_chr_cnt, smp_cnt, ptn_num, pcl_num, itn_num)
  } else {
    PMMH <- runPMMH_UnphasedChr_arma(int_sel_cof_A, int_dom_par_A, int_sel_cof_B, int_dom_par_B, int_rec_rat, int_pop_siz, smp_gen, smp_chr_cnt, smp_cnt, ptn_num, pcl_num, itn_num)
  }
  
  return(list(sel_cof_A_chn = as.vector(PMMH$sel_cof_A_chn), 
              dom_par_A_chn = as.vector(PMMH$dom_par_A_chn), 
              sel_cof_B_chn = as.vector(PMMH$sel_cof_B_chn), 
              dom_par_B_chn = as.vector(PMMH$dom_par_B_chn), 
              rec_rat_chn = as.vector(PMMH$rec_rat_chn), 
              pop_siz_chn = as.vector(PMMH$pop_siz_chn)))
}
#' Compiled version
cmprunPMMH <- cmpfun(runPMMH)

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

#' Standard version
runBayesianProcedure <- function(int_gen_par, smp_gen, smp_chr_cnt, smp_cnt, ptn_num, pcl_num, itn_num, brn_num, thn_num, grd_num) {
  int_sel_cof_A <- int_gen_par[1]
  int_dom_par_A <- int_gen_par[2]
  int_sel_cof_B <- int_gen_par[3]
  int_dom_par_B <- int_gen_par[4]
  int_rec_rat <- int_gen_par[5]
  int_pop_siz <- int_gen_par[6]
  
  # run the PMMH
  phased <- ifelse(nrow(smp_cnt) == 4, TRUE, FALSE)
  if (phased) {
    PMMH <- runPMMH_PhasedChr_arma(int_sel_cof_A, int_dom_par_A, int_sel_cof_B, int_dom_par_B, int_rec_rat, int_pop_siz, smp_gen, smp_chr_cnt, smp_cnt, ptn_num, pcl_num, itn_num)
  } else {
    PMMH <- runPMMH_UnphasedChr_arma(int_sel_cof_A, int_dom_par_A, int_sel_cof_B, int_dom_par_B, int_rec_rat, int_pop_siz, smp_gen, smp_chr_cnt, smp_cnt, ptn_num, pcl_num, itn_num)
  }
  
  # burn-in and thinning
  sel_cof_A_chn <- as.vector(PMMH$sel_cof_A_chn)
  sel_cof_A_chn <- sel_cof_A_chn[brn_num:length(sel_cof_A_chn)]
  sel_cof_A_chn <- sel_cof_A_chn[(1:round(length(sel_cof_A_chn) / thn_num)) * thn_num]
  sel_cof_B_chn <- as.vector(PMMH$sel_cof_B_chn)
  sel_cof_B_chn <- sel_cof_B_chn[brn_num:length(sel_cof_B_chn)]
  sel_cof_B_chn <- sel_cof_B_chn[(1:round(length(sel_cof_B_chn) / thn_num)) * thn_num]
  
  dom_par_A_chn <- as.vector(PMMH$dom_par_A_chn)
  dom_par_A_chn <- dom_par_A_chn[brn_num:length(dom_par_A_chn)]
  dom_par_A_chn <- dom_par_A_chn[(1:round(length(dom_par_A_chn) / thn_num)) * thn_num]
  dom_par_B_chn <- as.vector(PMMH$dom_par_B_chn)
  dom_par_B_chn <- dom_par_B_chn[brn_num:length(dom_par_B_chn)]
  dom_par_B_chn <- dom_par_B_chn[(1:round(length(dom_par_B_chn) / thn_num)) * thn_num]
  
  rec_rat_chn <- as.vector(PMMH$rec_rat_chn)
  rec_rat_chn <- rec_rat_chn[brn_num:length(rec_rat_chn)]
  rec_rat_chn <- rec_rat_chn[(1:round(length(rec_rat_chn) / thn_num)) * thn_num]
  
  pop_siz_chn <- as.vector(PMMH$pop_siz_chn)
  pop_siz_chn <- pop_siz_chn[brn_num:length(pop_siz_chn)]
  pop_siz_chn <- pop_siz_chn[(1:round(length(pop_siz_chn) / thn_num)) * thn_num]
  
  # MAP estimates for the selection coefficients
  if (length(sel_cof_A_chn) < 1e+05) {
    sel_cof_pdf <- kde2d(sel_cof_A_chn, sel_cof_B_chn, n = grd_num)
    sel_cof_A_grd <- sel_cof_pdf$x
    sel_cof_B_grd <- sel_cof_pdf$y
  } else {
    sel_cof_pdf <- kde2d(tail(sel_cof_A_grd, 1e+05), tail(sel_cof_B_grd, 1e+05), n = grd_num)
    sel_cof_A_grd <- sel_cof_pdf$x
    sel_cof_B_grd <- sel_cof_pdf$y
  }
  sel_cof_A_map <- sel_cof_A_grd[which(sel_cof_pdf$z == max(sel_cof_pdf$z), arr.ind = TRUE)[1]]
  sel_cof_B_map <- sel_cof_B_grd[which(sel_cof_pdf$z == max(sel_cof_pdf$z), arr.ind = TRUE)[2]]
  
  # MMSE estimates for the selection coefficients
  sel_cof_A_mmse <- mean(sel_cof_A_chn)
  sel_cof_B_mmse <- mean(sel_cof_B_chn)
  
  # 95% HPD intervals for the selection coefficients
  sel_cof_A_hpd <- HPDinterval(as.mcmc(sel_cof_A_chn), prob = 0.95)
  sel_cof_B_hpd <- HPDinterval(as.mcmc(sel_cof_B_chn), prob = 0.95)
  
  return(list(sel_cof_pdf = sel_cof_pdf, 
              sel_cof_A_map = sel_cof_A_map, 
              sel_cof_B_map = sel_cof_B_map, 
              sel_cof_A_mmse = sel_cof_A_mmse, 
              sel_cof_B_mmse = sel_cof_B_mmse, 
              sel_cof_A_hpd = sel_cof_A_hpd, 
              sel_cof_B_hpd = sel_cof_B_hpd, 
              sel_cof_A_chn = sel_cof_A_chn, 
              dom_par_A_chn = dom_par_A_chn, 
              sel_cof_B_chn = sel_cof_B_chn, 
              dom_par_B_chn = dom_par_B_chn, 
              rec_rat_chn = rec_rat_chn, 
              pop_siz_chn = pop_siz_chn))
}
#' Compiled version
cmprunBayesianProcedure <- cmpfun(runBayesianProcedure)

################################################################################
