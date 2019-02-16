// An MCMC-based method for Bayesian inference of natural selection from time series DNA data across linked loci
// Zhangyi He, Xiaoyang Dai, Mark Beaumont and Feng Yu

// C functions

#define ARMA_64BIT_WORD 1

#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::plugins(cpp11)]]
#include <math.h>

using namespace Rcpp;
using namespace std;
// using namespace arma;

/********** WFM **********/
// Simulate the haplotype frequency trajectories according to the two-locus Wright-Fisher model with selection
// [[Rcpp::export]]
arma::dmat simulateTLWFMS_arma(const double& sel_cof_A, const double& dom_par_A, const double& sel_cof_B, const double& dom_par_B, const double& rec_rat, const int& pop_siz, const arma::dcolvec& int_frq, const int& int_gen, const int& lst_gen) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  // declare the haplotype frequency trajectories
  arma::dmat hap_frq_pth(4, arma::uword(lst_gen - int_gen) + 1);
  
  // initialise the haplotype frequencies in generation 0
  hap_frq_pth.col(0) = int_frq;
  
  // declare the fitness
  arma::dmat fts(4, 4);
  fts(0, 0) = 1;
  fts(1, 0) = (1 - dom_par_B * sel_cof_B);
  fts(2, 0) = (1 - dom_par_A * sel_cof_A);
  fts(3, 0) = (1 - dom_par_A * sel_cof_A) * (1 - dom_par_B * sel_cof_B);
  fts(0, 1) = (1 - dom_par_B * sel_cof_B);
  fts(1, 1) = (1 - sel_cof_B);
  fts(2, 1) = (1 - dom_par_A * sel_cof_A) * (1 - dom_par_B * sel_cof_B);
  fts(3, 1) = (1 - dom_par_A * sel_cof_A) * (1 - sel_cof_B);
  fts(0, 2) = (1 - dom_par_A * sel_cof_A);
  fts(1, 2) = (1 - dom_par_A * sel_cof_A) * (1 - dom_par_B * sel_cof_B);
  fts(2, 2) = (1 - sel_cof_A);
  fts(3, 2) = (1 - sel_cof_A) * (1 - dom_par_B * sel_cof_B);
  fts(0, 3) = (1 - dom_par_A * sel_cof_A) * (1 - dom_par_B * sel_cof_B);
  fts(1, 3) = (1 - dom_par_A * sel_cof_A) * (1 - sel_cof_B);
  fts(2, 3) = (1 - sel_cof_A) * (1 - dom_par_B * sel_cof_B);
  fts(3, 3) = (1 - sel_cof_A) * (1 - sel_cof_B);
  
  // declare and initialise the haplotype frequencies during a single generation of the life cycle
  arma::dcolvec hap_frq = int_frq;
  
  // declare and initialise the genotype frequencies during a single generation of the life cycle
  arma::dmat gen_frq = hap_frq * hap_frq.t();
  
  // declare eta
  arma::dcolvec eta(4);
  eta(0) = -1;
  eta(1) = 1;
  eta(2) = 1;
  eta(3) = -1;
  
  for(arma::uword t = 1; t < arma::uword(lst_gen - int_gen) + 1; t++) {
    // random union of gametes
    gen_frq = hap_frq * hap_frq.t();
    
    // viability selection
    gen_frq = (fts % gen_frq) / arma::as_scalar(sum(sum(fts % gen_frq, 0), 1));
    
    // meiosis (calculate the sampling probability)
    arma::dcolvec prob = arma::zeros<arma::dcolvec>(4);
    for(arma::uword i = 0; i < 4; i++) {
      prob(i) = (sum(gen_frq.row(i)) + sum(gen_frq.col(i))) / 2;
    }
    prob = prob + eta * rec_rat * (prob(0) * prob(3) - prob(1) * prob(2));
    
    // reproduction (the Wright-Fisher sampling)
    IntegerVector hap_cnt(4);
    R::rmultinom(2 * pop_siz, prob.begin(), 4, hap_cnt.begin());
    hap_frq = as<arma::dcolvec>(hap_cnt) / 2 / pop_siz;
    
    hap_frq_pth.col(t) = hap_frq;
  }
  
  return hap_frq_pth;
}
/*************************/


/********** WFD **********/
// Simulate the haplotype frequency trajectories according to the two-locus Wright-Fisher diffusion with selection using the Euler-Maruyama method
// [[Rcpp::export]]
arma::dmat simulateTLWFDS_arma(const double& sel_cof_A, const double& dom_par_A, const double& sel_cof_B, const double& dom_par_B, const double& rec_rat, const int& pop_siz, const arma::dcolvec& int_frq, const int& int_gen, const int& lst_gen, const arma::uword& ptn_num) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  // rescale the selection coefficients and the recombination rate
  double scl_sel_cof_A = 2 * pop_siz * sel_cof_A;
  double scl_sel_cof_B = 2 * pop_siz * sel_cof_B;
  double scl_rec_rat = 4 * pop_siz * rec_rat;
  
  // declare delta t
  double dt = 1.0 / (2 * pop_siz) / ptn_num;  
  // declare delta W
  arma::dmat dW = pow(dt, 0.5) * arma::randn<arma::dmat>(6, arma::uword(lst_gen - int_gen) * ptn_num);
  
  // declare the haplotype frequency trajectories
  arma::dmat hap_frq_pth(4, arma::uword(lst_gen - int_gen) * ptn_num + 1);
  
  // initialise the haplotype frequencies at generation 0
  hap_frq_pth.col(0) = int_frq;
  
  for (arma::uword t = 1; t < arma::uword(lst_gen - int_gen) * ptn_num + 1; t++) {
    // calculate the drift coefficient vector
    arma::dcolvec mu(4);
    mu(0) =  scl_sel_cof_A * hap_frq_pth(0, t - 1) * (hap_frq_pth(2, t - 1) + hap_frq_pth(3, t - 1)) * ((hap_frq_pth(0, t - 1) + hap_frq_pth(1, t - 1)) * dom_par_A + (hap_frq_pth(2, t - 1) + hap_frq_pth(3, t - 1)) * (1 - dom_par_A)) + scl_sel_cof_B * hap_frq_pth(0, t - 1) * (hap_frq_pth(1, t - 1) + hap_frq_pth(3, t - 1)) * ((hap_frq_pth(0, t - 1) + hap_frq_pth(2, t - 1)) * dom_par_B + (hap_frq_pth(1, t - 1) + hap_frq_pth(3, t - 1)) * (1 - dom_par_B)) - 0.5 * scl_rec_rat * (hap_frq_pth(0, t - 1) * hap_frq_pth(3, t - 1) - hap_frq_pth(1, t - 1) * hap_frq_pth(2, t - 1));
    mu(1) =  scl_sel_cof_A * hap_frq_pth(1, t - 1) * (hap_frq_pth(2, t - 1) + hap_frq_pth(3, t - 1)) * ((hap_frq_pth(0, t - 1) + hap_frq_pth(1, t - 1)) * dom_par_A + (hap_frq_pth(2, t - 1) + hap_frq_pth(3, t - 1)) * (1 - dom_par_A)) - scl_sel_cof_B * hap_frq_pth(1, t - 1) * (hap_frq_pth(0, t - 1) + hap_frq_pth(2, t - 1)) * ((hap_frq_pth(0, t - 1) + hap_frq_pth(2, t - 1)) * dom_par_B + (hap_frq_pth(1, t - 1) + hap_frq_pth(3, t - 1)) * (1 - dom_par_B)) + 0.5 * scl_rec_rat * (hap_frq_pth(0, t - 1) * hap_frq_pth(3, t - 1) - hap_frq_pth(1, t - 1) * hap_frq_pth(2, t - 1));
    mu(2) = -scl_sel_cof_A * hap_frq_pth(2, t - 1) * (hap_frq_pth(0, t - 1) + hap_frq_pth(1, t - 1)) * ((hap_frq_pth(0, t - 1) + hap_frq_pth(1, t - 1)) * dom_par_A + (hap_frq_pth(2, t - 1) + hap_frq_pth(3, t - 1)) * (1 - dom_par_A)) + scl_sel_cof_B * hap_frq_pth(2, t - 1) * (hap_frq_pth(1, t - 1) + hap_frq_pth(3, t - 1)) * ((hap_frq_pth(0, t - 1) + hap_frq_pth(2, t - 1)) * dom_par_B + (hap_frq_pth(1, t - 1) + hap_frq_pth(3, t - 1)) * (1 - dom_par_B)) + 0.5 * scl_rec_rat * (hap_frq_pth(0, t - 1) * hap_frq_pth(3, t - 1) - hap_frq_pth(1, t - 1) * hap_frq_pth(2, t - 1));
    mu(3) = -scl_sel_cof_A * hap_frq_pth(3, t - 1) * (hap_frq_pth(0, t - 1) + hap_frq_pth(1, t - 1)) * ((hap_frq_pth(0, t - 1) + hap_frq_pth(1, t - 1)) * dom_par_A + (hap_frq_pth(2, t - 1) + hap_frq_pth(3, t - 1)) * (1 - dom_par_A)) - scl_sel_cof_B * hap_frq_pth(3, t - 1) * (hap_frq_pth(0, t - 1) + hap_frq_pth(2, t - 1)) * ((hap_frq_pth(0, t - 1) + hap_frq_pth(2, t - 1)) * dom_par_B + (hap_frq_pth(1, t - 1) + hap_frq_pth(3, t - 1)) * (1 - dom_par_B)) - 0.5 * scl_rec_rat * (hap_frq_pth(0, t - 1) * hap_frq_pth(3, t - 1) - hap_frq_pth(1, t - 1) * hap_frq_pth(2, t - 1));
    
    // calculate the diffusion coefficient matrix
    arma::dmat sigma(4, 6);
    sigma(0, 0) = pow(hap_frq_pth(0, t - 1) * hap_frq_pth(1, t - 1), 0.5);
    sigma(0, 1) = pow(hap_frq_pth(0, t - 1) * hap_frq_pth(2, t - 1), 0.5);
    sigma(0, 2) = pow(hap_frq_pth(0, t - 1) * hap_frq_pth(3, t - 1), 0.5);
    sigma(0, 3) = 0;
    sigma(0, 4) = 0;
    sigma(0, 5) = 0;
    sigma(1, 0) = -pow(hap_frq_pth(1, t - 1) * hap_frq_pth(0, t - 1), 0.5);
    sigma(1, 1) = 0;
    sigma(1, 2) = 0;
    sigma(1, 3) = pow(hap_frq_pth(1, t - 1) * hap_frq_pth(2, t - 1), 0.5);
    sigma(1, 4) = pow(hap_frq_pth(1, t - 1) * hap_frq_pth(3, t - 1), 0.5);
    sigma(1, 5) = 0;
    sigma(2, 0) = 0;
    sigma(2, 1) = -pow(hap_frq_pth(2, t - 1) * hap_frq_pth(0, t - 1), 0.5);
    sigma(2, 2) = 0;
    sigma(2, 3) = -pow(hap_frq_pth(2, t - 1) * hap_frq_pth(1, t - 1), 0.5);
    sigma(2, 4) = 0;
    sigma(2, 5) = pow(hap_frq_pth(2, t - 1) * hap_frq_pth(3, t - 1), 0.5);
    sigma(3, 0) = 0;
    sigma(3, 1) = 0;
    sigma(3, 2) = -pow(hap_frq_pth(3, t - 1) * hap_frq_pth(0, t - 1), 0.5);
    sigma(3, 3) = 0;
    sigma(3, 4) = -pow(hap_frq_pth(3, t - 1) * hap_frq_pth(1, t - 1), 0.5);
    sigma(3, 5) = -pow(hap_frq_pth(3, t - 1) * hap_frq_pth(2, t - 1), 0.5);
    
    // proceed the Euler-Maruyama scheme
    hap_frq_pth.col(t) = hap_frq_pth.col(t - 1) + mu * dt + sigma * dW.col(t - 1);
    
    // remove the noise from the numerical techniques
    for(arma::uword i = 0; i < 4; i++) {
      if(hap_frq_pth(i, t) < 0) {
        hap_frq_pth(i, t) = 0;
      } 
      if(hap_frq_pth(i, t) > 1) {
        hap_frq_pth(i, t) = 1;
      }
    }
    hap_frq_pth.col(t) = hap_frq_pth.col(t) / sum(hap_frq_pth.col(t));
  }
  
  return hap_frq_pth;
}
/*************************/


/********** BPF **********/
// Calculate the multinomial probabilities
// [[Rcpp::export]]
double calculateMultinomProb_arma(const arma::icolvec& smp_hap_cnt, const int& smp_chr_cnt, const arma::dcolvec& pop_hap_frq) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  if (arma::any(pop_hap_frq == 0)) {
    if (arma::any(smp_hap_cnt.elem(arma::find(pop_hap_frq == 0)) != 0)) {
      return 0;
    }
    
    arma::icolvec smp_hap_cnt_nonzero = smp_hap_cnt.elem(arma::find(pop_hap_frq != 0));
    arma::dcolvec pop_hap_frq_nonzero = pop_hap_frq.elem(arma::find(pop_hap_frq != 0));
    
    return exp(lgamma(smp_chr_cnt + 1) + sum(arma::conv_to<arma::dcolvec>::from(smp_hap_cnt_nonzero) % log(pop_hap_frq_nonzero) - lgamma(arma::conv_to<arma::dcolvec>::from(smp_hap_cnt_nonzero) + 1)));
  } else {
    return exp(lgamma(smp_chr_cnt + 1) + sum(arma::conv_to<arma::dcolvec>::from(smp_hap_cnt) % log(pop_hap_frq) - lgamma(arma::conv_to<arma::dcolvec>::from(smp_hap_cnt) + 1)));
  }
}

/* Observations with unphased chromosomes */
// Calculate the possible unobserved haplotype counts from the observed allele counts in the sample
// [[Rcpp::export]]
arma::imat calculateHaploCount_arma(const int& smp_chr_cnt, const arma::icolvec& smp_ale_cnt) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  if (smp_ale_cnt(0) + smp_ale_cnt(1) - smp_chr_cnt > 0) {
    int smp_hap_cnt_min = smp_ale_cnt(0) + smp_ale_cnt(1) - smp_chr_cnt;
    int smp_hap_cnt_max = min(smp_ale_cnt(0), smp_ale_cnt(1));
    arma::irowvec smp_hap_cnt_rng = arma::linspace<arma::irowvec>(smp_hap_cnt_min, smp_hap_cnt_max, smp_hap_cnt_max - smp_hap_cnt_min + 1);
    
    arma::imat smp_hap_cnt = arma::zeros<arma::imat>(4, smp_hap_cnt_rng.n_elem);
    smp_hap_cnt.row(0) = smp_hap_cnt_rng;
    smp_hap_cnt.row(1) = smp_ale_cnt(0) - smp_hap_cnt.row(0);
    smp_hap_cnt.row(2) = smp_ale_cnt(1) - smp_hap_cnt.row(0);
    smp_hap_cnt.row(3) = smp_chr_cnt - smp_hap_cnt.row(0) - smp_hap_cnt.row(1) - smp_hap_cnt.row(2);
    
    return smp_hap_cnt;
  } else {
    int smp_hap_cnt_min = 0;
    int smp_hap_cnt_max = min(smp_ale_cnt(0), smp_ale_cnt(1));
    arma::irowvec smp_hap_cnt_rng = arma::linspace<arma::irowvec>(smp_hap_cnt_min, smp_hap_cnt_max, smp_hap_cnt_max - smp_hap_cnt_min + 1);
    
    arma::imat smp_hap_cnt = arma::zeros<arma::imat>(4, smp_hap_cnt_rng.n_elem);
    smp_hap_cnt.row(0) = smp_hap_cnt_rng;
    smp_hap_cnt.row(1) = smp_ale_cnt(0) - smp_hap_cnt.row(0);
    smp_hap_cnt.row(2) = smp_ale_cnt(1) - smp_hap_cnt.row(0);
    smp_hap_cnt.row(3) = smp_chr_cnt - smp_hap_cnt.row(0) - smp_hap_cnt.row(1) - smp_hap_cnt.row(2);
    
    return smp_hap_cnt;
  }
}

// Run the particle filter for the observations with unphased chromosomes
// [[Rcpp::export]]
List runBPF_UnphasedChr_arma(const double& sel_cof_A, const double& dom_par_A, const double& sel_cof_B, const double& dom_par_B, const double& rec_rat, const int& pop_siz, const arma::irowvec& smp_gen, const arma::irowvec& smp_chr_cnt, const arma::imat& smp_ale_cnt, const arma::uword& ptn_num, const arma::uword& pcl_num) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  double lik = 1;
  
  arma::dmat wght = arma::zeros<arma::dmat>(pcl_num, smp_gen.n_elem);
  arma::dcube part_pre = arma::zeros<arma::dcube>(4, pcl_num, smp_gen.n_elem);
  arma::dcube part_pst = arma::zeros<arma::dcube>(4, pcl_num, smp_gen.n_elem);
  
  arma::dcolvec wght_tmp = arma::zeros<arma::dcolvec>(pcl_num);
  arma::dmat part_tmp = arma::zeros<arma::dmat>(4, pcl_num);
  
  // initialise the particles
  cout << "generation: " << smp_gen(0) << endl;
  part_tmp = arma::normalise(arma::randu<arma::dmat>(4, pcl_num), 1, 0);
  arma::imat smp_hap_cnt = calculateHaploCount_arma(smp_chr_cnt(0), smp_ale_cnt.col(0));
  for (arma::uword i = 0; i < pcl_num; i++) {
    for (arma::uword j = 0; j < smp_hap_cnt.n_cols; j++) {
      wght_tmp(i) = wght_tmp(i) + calculateMultinomProb_arma(smp_hap_cnt.col(j), smp_chr_cnt(0), part_tmp.col(i));
    }
  }
  
  if (arma::sum(wght_tmp) > 0) {
    arma::dcolvec prob = arma::normalise(wght_tmp, 1);
    arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
    arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, TRUE, prob);
    
    lik = lik * (arma::sum(wght_tmp) / pcl_num);
    wght.col(0) = wght_tmp;
    part_pre.slice(0) = part_tmp;
    part_pst.slice(0) = part_tmp.cols(indx);
  } else {
    lik = 0;
    wght.shed_cols(0, smp_gen.n_elem - 1);
    part_pre.shed_slices(0, smp_gen.n_elem - 1);
    part_pst.shed_slices(0, smp_gen.n_elem - 1);
    
    return List::create(Named("lik", lik), 
                        Named("wght", wght), 
                        Named("part_pre_resmp", part_pre), 
                        Named("part_pst_resmp", part_pst));
  }
  
  // run the bootstrap particle filter
  for (arma::uword k = 1; k < smp_gen.n_elem; k++) {
    cout << "generation: " << smp_gen(k) << endl;
    wght_tmp = arma::zeros<arma::dcolvec>(pcl_num);
    part_tmp = part_pst.slice(k - 1);
    arma::imat smp_hap_cnt = calculateHaploCount_arma(smp_chr_cnt(k), smp_ale_cnt.col(k));
    for (arma::uword i = 0; i < pcl_num; i++) {
      arma::dmat path = simulateTLWFDS_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, part_tmp.col(i), smp_gen(k - 1), smp_gen(k), ptn_num);
      part_tmp.col(i) = arma::vectorise(path.tail_cols(1), 0);
      for (arma::uword j = 0; j < smp_hap_cnt.n_cols; j++) {
        wght_tmp(i) = wght_tmp(i) + calculateMultinomProb_arma(smp_hap_cnt.col(j), smp_chr_cnt(k), part_tmp.col(i));
      }
    }
    
    if (arma::sum(wght_tmp) > 0) {
      arma::dcolvec prob = arma::normalise(wght_tmp, 1);
      arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
      arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, TRUE, prob);
      
      lik = lik * (arma::sum(wght_tmp) / pcl_num);
      wght.col(k) = wght_tmp;
      part_pre.slice(k) = part_tmp;
      part_pst.slice(k) = part_tmp.cols(indx);
    } else {
      lik = 0;
      wght.shed_cols(k, smp_gen.n_elem - 1);
      part_pre.shed_slices(k, smp_gen.n_elem - 1);
      part_pst.shed_slices(k, smp_gen.n_elem - 1);
      break;
    }
  }
  
  return List::create(Named("lik", lik),
                      Named("wght", wght), 
                      Named("part_pre_resmp", part_pre), 
                      Named("part_pst_resmp", part_pst));
}

/* Observations with phased chromosomes */
// Run the particle filter for the observations with phased chromosomes
// [[Rcpp::export]]
List runBPF_PhasedChr_arma(const double& sel_cof_A, const double& dom_par_A, const double& sel_cof_B, const double& dom_par_B, const double& rec_rat, const int& pop_siz, const arma::irowvec& smp_gen, const arma::irowvec& smp_chr_cnt, const arma::imat& smp_hap_cnt, const arma::uword& ptn_num, const arma::uword& pcl_num) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  double lik = 1;
  
  arma::dmat wght = arma::zeros<arma::dmat>(pcl_num, smp_gen.n_elem);
  arma::dcube part_pre = arma::zeros<arma::dcube>(4, pcl_num, smp_gen.n_elem);
  arma::dcube part_pst = arma::zeros<arma::dcube>(4, pcl_num, smp_gen.n_elem);
  
  arma::dcolvec wght_tmp = arma::zeros<arma::dcolvec>(pcl_num);
  arma::dmat part_tmp = arma::zeros<arma::dmat>(4, pcl_num);
  
  // initialise the particles
  cout << "generation: " << smp_gen(0) << endl;
  part_tmp = arma::normalise(arma::randu<arma::dmat>(4, pcl_num), 1, 0);
  for (arma::uword i = 0; i < pcl_num; i++) {
    wght_tmp(i) = calculateMultinomProb_arma(smp_hap_cnt.col(0), smp_chr_cnt(0), part_tmp.col(i));
  }
  
  if (arma::sum(wght_tmp) > 0) {
    arma::dcolvec prob = arma::normalise(wght_tmp, 1);
    arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
    arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, TRUE, prob);
    
    lik = lik * (arma::sum(wght_tmp) / pcl_num);
    wght.col(0) = wght_tmp;
    part_pre.slice(0) = part_tmp;
    part_pst.slice(0) = part_tmp.cols(indx);
  } else {
    lik = 0;
    wght.shed_cols(0, smp_gen.n_elem - 1);
    part_pre.shed_slices(0, smp_gen.n_elem - 1);
    part_pst.shed_slices(0, smp_gen.n_elem - 1);
    
    return List::create(Named("lik", lik), 
                        Named("wght", wght), 
                        Named("part_pre_resmp", part_pre), 
                        Named("part_pst_resmp", part_pst));
  }
  
  // run the bootstrap particle filter  
  for (arma::uword k = 1; k < smp_gen.n_elem; k++) {
    cout << "generation: " << smp_gen(k) << endl;
    part_tmp = part_pst.slice(k - 1);
    for (arma::uword i = 0; i < pcl_num; i++) {
      arma::dmat path = simulateTLWFDS_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, part_tmp.col(i), smp_gen(k - 1), smp_gen(k), ptn_num);
      part_tmp.col(i) = arma::vectorise(path.tail_cols(1), 0);
      wght_tmp(i) = calculateMultinomProb_arma(smp_hap_cnt.col(k), smp_chr_cnt(k), part_tmp.col(i));
    }
    
    if (arma::sum(wght_tmp) > 0) {
      arma::dcolvec prob = arma::normalise(wght_tmp, 1);
      arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
      arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, TRUE, prob);
      
      lik = lik * (arma::sum(wght_tmp) / pcl_num);
      wght.col(k) = wght_tmp;
      part_pre.slice(k) = part_tmp;
      part_pst.slice(k) = part_tmp.cols(indx);
    } else {
      lik = 0;
      wght.shed_cols(k, smp_gen.n_elem - 1);
      part_pre.shed_slices(k, smp_gen.n_elem - 1);
      part_pst.shed_slices(k, smp_gen.n_elem - 1);
      break;
    }
  }
  
  return List::create(Named("lik", lik),
                      Named("wght", wght), 
                      Named("part_pre_resmp", part_pre), 
                      Named("part_pst_resmp", part_pst));
}
/*************************/


/********** PMMH **********/
/* Observations with unphased chromosomes */
// Calculate the log-likelihood using the bootstrap particle filter for the observations with unphased chromosomes
// [[Rcpp::export]]
void calculateLogLikelihood_UnphasedChr_arma(double& log_lik, const double& sel_cof_A, const double& dom_par_A, const double& sel_cof_B, const double& dom_par_B, const double& rec_rat, const int& pop_siz, const arma::irowvec& smp_gen, const arma::irowvec& smp_chr_cnt, const arma::imat& smp_ale_cnt, const arma::uword& ptn_num, const arma::uword& pcl_num) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  arma::dcolvec wght = arma::zeros<arma::dcolvec>(pcl_num);
  arma::dmat part_pre = arma::zeros<arma::dmat>(4, pcl_num);
  arma::dmat part_pst = arma::zeros<arma::dmat>(4, pcl_num);
  
  // initialise the particles
  part_pre = arma::normalise(arma::randu<arma::dmat>(4, pcl_num), 1, 0);
  arma::imat smp_hap_cnt = calculateHaploCount_arma(smp_chr_cnt(0), smp_ale_cnt.col(0));
  for (arma::uword i = 0; i < pcl_num; i++) {
    for (arma::uword j = 0; j < smp_hap_cnt.n_cols; j++) {
      wght(i) = wght(i) + calculateMultinomProb_arma(smp_hap_cnt.col(j), smp_chr_cnt(0), part_pre.col(i));
    }
  }
  
  if (arma::sum(wght) > 0) {
    log_lik = log_lik + log(arma::sum(wght)) - log(pcl_num);
    arma::dcolvec prob = arma::normalise(wght, 1);
    arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
    arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, TRUE, prob);
    part_pst = part_pre.cols(indx);
  } else {
    log_lik = -(arma::datum::inf);
    return;
  }
  
  // run the bootstrap particle filter  
  for (arma::uword k = 1; k < smp_gen.n_elem; k++) {
    wght = arma::zeros<arma::dcolvec>(pcl_num);
    arma::imat smp_hap_cnt = calculateHaploCount_arma(smp_chr_cnt(k), smp_ale_cnt.col(k));
    for (arma::uword i = 0; i < pcl_num; i++) {
      arma::dmat path = simulateTLWFDS_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, part_pst.col(i), smp_gen(k - 1), smp_gen(k), ptn_num);
      part_pre.col(i) = arma::vectorise(path.tail_cols(1), 0);
      for (arma::uword j = 0; j < smp_hap_cnt.n_cols; j++) {
        wght(i) = wght(i) + calculateMultinomProb_arma(smp_hap_cnt.col(j), smp_chr_cnt(k), part_pre.col(i));
      }
    }
    
    if (arma::sum(wght) > 0) {
      log_lik = log_lik + (log(arma::sum(wght)) - log(pcl_num));
      arma::dcolvec prob = arma::normalise(wght, 1);
      arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
      arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, TRUE, prob);
      part_pst = part_pre.cols(indx);
    } else {
      log_lik = -(arma::datum::inf);
      break;
    }
  }
}

// Run the particle marginal Metropolis-Hastings for the observations with unphased chromosomes
//[[Rcpp::export]]
List runPMMH_UnphasedChr_arma(const double& int_sel_cof_A, const double& int_dom_par_A, const double& int_sel_cof_B, const double& int_dom_par_B, const double& int_rec_rat, const int& int_pop_siz, const arma::irowvec& smp_gen, const arma::irowvec& smp_chr_cnt, const arma::imat& smp_ale_cnt, const arma::uword& ptn_num, const arma::uword& pcl_num, const arma::uword& itn_num) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  arma::drowvec sel_cof_A_chn = arma::zeros<arma::drowvec>(itn_num);
  arma::drowvec dom_par_A_chn = arma::zeros<arma::drowvec>(itn_num);
  arma::drowvec sel_cof_B_chn = arma::zeros<arma::drowvec>(itn_num);
  arma::drowvec dom_par_B_chn = arma::zeros<arma::drowvec>(itn_num);
  arma::drowvec rec_rat_chn = arma::zeros<arma::drowvec>(itn_num);
  arma::irowvec pop_siz_chn = arma::zeros<arma::irowvec>(itn_num);
  
  // suppose that the dominance parameters, the recombination rate and the population size are known
  dom_par_A_chn.fill(int_dom_par_A);
  dom_par_B_chn.fill(int_dom_par_B);
  rec_rat_chn.fill(int_rec_rat);
  pop_siz_chn.fill(int_pop_siz);
  
  arma::drowvec log_pri_chn = arma::zeros<arma::drowvec>(itn_num);
  arma::drowvec log_lik_chn = arma::zeros<arma::drowvec>(itn_num);
  
  double sel_cof_sd = 1e-02;
  //double rec_rat_sd = 1e-02;
  
  // initialise the population genetic parameters
  cout << "iteration: " << 1 << endl;
  // take the uniform prior and fix the initial value of the selection coefficient to zero
  // or take the beta prior with alpha = 1 and beta = 3
  sel_cof_A_chn(0) = int_sel_cof_A;
  sel_cof_B_chn(0) = int_sel_cof_B;
  // take the uniform prior and fix the initial value of the recombination rate to zero
  // or take the beta prior with alpha = 1 and beta = 3  
  //rec_rat_chn(0) = int_rec_rat;
  
  calculateLogLikelihood_UnphasedChr_arma(log_lik_chn(0), sel_cof_A_chn(0), dom_par_A_chn(0), sel_cof_B_chn(0), dom_par_B_chn(0), rec_rat_chn(0), pop_siz_chn(0), smp_gen, smp_chr_cnt, smp_ale_cnt, ptn_num, pcl_num);
  
  double apt_rto = 0;
  for (arma::uword i = 1; i < itn_num; i++) {
    cout << "iteration: " << i + 1 << endl;
    apt_rto = 1;
    
    // draw the candidates of the population genetic parameters from the proposals
    sel_cof_A_chn(i) = sel_cof_A_chn(i - 1) + sel_cof_sd * arma::randn();
    if (sel_cof_A_chn(i) > 1) {
      apt_rto = 0;
    }
    sel_cof_B_chn(i) = sel_cof_B_chn(i - 1) + sel_cof_sd * arma::randn();
    if (sel_cof_B_chn(i) > 1) {
      apt_rto = 0;
    }
    /*
    rec_rat_chn(i) = rec_rat_chn(i - 1) + rec_rat_sd * arma::randn();
    if (rec_rat_chn(i) < 0 || rec_rat_chn(i) > 0.5) {
      apt_rto = 0;
    }
    */
    
    if (apt_rto == 1) {
      //double log_psl_old_new = log(arma::normpdf(sel_cof_A_chn(i - 1), sel_cof_A_chn(i), sel_cof_sd)) + log(arma::normpdf(sel_cof_B_chn(i - 1), sel_cof_B_chn(i), sel_cof_sd)) + log(arma::normpdf(rec_rat_chn(i - 1), rec_rat_chn(i), rec_rat_sd));
      //double log_psl_new_old = log(arma::normpdf(sel_cof_A_chn(i), sel_cof_A_chn(i - 1), sel_cof_sd)) + log(arma::normpdf(sel_cof_B_chn(i), sel_cof_B_chn(i - 1), sel_cof_sd)) + log(arma::normpdf(rec_rat_chn(i), rec_rat_chn(i - 1), rec_rat_sd));
      
      calculateLogLikelihood_UnphasedChr_arma(log_lik_chn(i), sel_cof_A_chn(i), dom_par_A_chn(i), sel_cof_B_chn(i), dom_par_B_chn(i), rec_rat_chn(i), pop_siz_chn(i), smp_gen, smp_chr_cnt, smp_ale_cnt, ptn_num, pcl_num);
      
      // calculate the acceptance ratio
      //apt_rto = exp((log_pri_chn(i) + log_lik_chn(i) + log_psl_old_new) - (log_pri_chn(i - 1) + log_lik_chn(i - 1) + log_psl_new_old));
      apt_rto = exp((log_pri_chn(i) + log_lik_chn(i)) - (log_pri_chn(i - 1) + log_lik_chn(i - 1)));
      if (arma::randu() > apt_rto) {
        sel_cof_A_chn(i) = sel_cof_A_chn(i - 1);
        sel_cof_B_chn(i) = sel_cof_B_chn(i - 1);
        //rec_rat_chn(i) = rec_rat_chn(i - 1);
        log_lik_chn(i) = log_lik_chn(i - 1);
      }
    } else {
      sel_cof_A_chn(i) = sel_cof_A_chn(i - 1);
      sel_cof_B_chn(i) = sel_cof_B_chn(i - 1);
      //rec_rat_chn(i) = rec_rat_chn(i - 1);
      log_lik_chn(i) = log_lik_chn(i - 1);
    }
  }
  
  return List::create(Named("sel_cof_A_chn", sel_cof_A_chn), 
                      Named("dom_par_A_chn", dom_par_A_chn), 
                      Named("sel_cof_B_chn", sel_cof_B_chn), 
                      Named("dom_par_B_chn", dom_par_B_chn), 
                      Named("rec_rat_chn", rec_rat_chn), 
                      Named("pop_siz_chn", pop_siz_chn));
}

/* Observations with phased chromosomes */
// Calculate the log-likelihood using the bootstrap particle filter for the observations with phased chromosomes
// [[Rcpp::export]]
void calculateLogLikelihood_PhasedChr_arma(double& log_lik, const double& sel_cof_A, const double& dom_par_A, const double& sel_cof_B, const double& dom_par_B, const double& rec_rat, const int& pop_siz, const arma::irowvec& smp_gen, const arma::irowvec& smp_chr_cnt, const arma::imat& smp_hap_cnt, const arma::uword& ptn_num, const arma::uword& pcl_num) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  arma::dcolvec wght = arma::zeros<arma::dcolvec>(pcl_num);
  arma::dmat part_pre = arma::zeros<arma::dmat>(4, pcl_num);
  arma::dmat part_pst = arma::zeros<arma::dmat>(4, pcl_num);
  
  // initialise the particles
  part_pre = arma::normalise(arma::randu<arma::dmat>(4, pcl_num), 1, 0);
  for (arma::uword i = 0; i < pcl_num; i++) {
    wght(i) = calculateMultinomProb_arma(smp_hap_cnt.col(0), smp_chr_cnt(0), part_pre.col(i));
  }

  if (arma::sum(wght) > 0) {
    log_lik = log_lik + log(arma::sum(wght)) - log(pcl_num);
    arma::dcolvec prob = arma::normalise(wght, 1);
    arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
    arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, TRUE, prob);
    part_pst = part_pre.cols(indx);
  } else {
    log_lik = -(arma::datum::inf);
    return;
  }
  
  // run the bootstrap particle filter  
  for (arma::uword k = 1; k < smp_gen.n_elem; k++) {
    for (arma::uword i = 0; i < pcl_num; i++) {
      arma::dmat path = simulateTLWFDS_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, part_pst.col(i), smp_gen(k - 1), smp_gen(k), ptn_num);
      part_pre.col(i) = arma::vectorise(path.tail_cols(1), 0);
      wght(i) = calculateMultinomProb_arma(smp_hap_cnt.col(k), smp_chr_cnt(k), part_pre.col(i));
    }

    if (arma::sum(wght) > 0) {
      log_lik = log_lik + (log(arma::sum(wght)) - log(pcl_num));
      arma::dcolvec prob = arma::normalise(wght, 1);
      arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
      arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, TRUE, prob);
      part_pst = part_pre.cols(indx);
    } else {
      log_lik = -(arma::datum::inf);
      break;
    }
  }
}

// Run the particle marginal Metropolis-Hastings for the observations with phased chromosomes
//[[Rcpp::export]]
List runPMMH_PhasedChr_arma(const double& int_sel_cof_A, const double& int_dom_par_A, const double& int_sel_cof_B, const double& int_dom_par_B, const double& int_rec_rat, const int& int_pop_siz, const arma::irowvec& smp_gen, const arma::irowvec& smp_chr_cnt, const arma::imat& smp_hap_cnt, const arma::uword& ptn_num, const arma::uword& pcl_num, const arma::uword& itn_num) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  arma::drowvec sel_cof_A_chn = arma::zeros<arma::drowvec>(itn_num);
  arma::drowvec dom_par_A_chn = arma::zeros<arma::drowvec>(itn_num);
  arma::drowvec sel_cof_B_chn = arma::zeros<arma::drowvec>(itn_num);
  arma::drowvec dom_par_B_chn = arma::zeros<arma::drowvec>(itn_num);
  arma::drowvec rec_rat_chn = arma::zeros<arma::drowvec>(itn_num);
  arma::irowvec pop_siz_chn = arma::zeros<arma::irowvec>(itn_num);
  
  // suppose that the dominance parameters, the recombination rate and the population size are known
  dom_par_A_chn.fill(int_dom_par_A);
  dom_par_B_chn.fill(int_dom_par_B);
  rec_rat_chn.fill(int_rec_rat);
  pop_siz_chn.fill(int_pop_siz);
  
  arma::drowvec log_pri_chn = arma::zeros<arma::drowvec>(itn_num);
  arma::drowvec log_lik_chn = arma::zeros<arma::drowvec>(itn_num);
  
  double sel_cof_sd = 1e-02;
  //double rec_rat_sd = 1e-02;
  
  // initialise the population genetic parameters
  cout << "iteration: " << 1 << endl;
  // take the uniform prior and fix the initial value of the selection coefficient to zero
  // or take the beta prior with alpha = 1 and beta = 3
  sel_cof_A_chn(0) = int_sel_cof_A;
  sel_cof_B_chn(0) = int_sel_cof_B;
  // take the uniform prior and fix the initial value of the recombination rate to zero
  // or take the beta prior with alpha = 1 and beta = 3  
  //rec_rat_chn(0) = int_rec_rat;
  
  calculateLogLikelihood_PhasedChr_arma(log_lik_chn(0), sel_cof_A_chn(0), dom_par_A_chn(0), sel_cof_B_chn(0), dom_par_B_chn(0), rec_rat_chn(0), pop_siz_chn(0), smp_gen, smp_chr_cnt, smp_hap_cnt, ptn_num, pcl_num);
  
  double apt_rto = 0;
  for (arma::uword i = 1; i < itn_num; i++) {
    cout << "iteration: " << i + 1 << endl;
    apt_rto = 1;
    
    // draw the candidates of the population genetic parameters from the proposals
    sel_cof_A_chn(i) = sel_cof_A_chn(i - 1) + sel_cof_sd * arma::randn();
    if (sel_cof_A_chn(i) > 1) {
      apt_rto = 0;
    }
    sel_cof_B_chn(i) = sel_cof_B_chn(i - 1) + sel_cof_sd * arma::randn();
    if (sel_cof_B_chn(i) > 1) {
      apt_rto = 0;
    }
    /*
    rec_rat_chn(i) = rec_rat_chn(i - 1) + rec_rat_sd * arma::randn();
    if (rec_rat_chn(i) < 0 || rec_rat_chn(i) > 0.5) {
      apt_rto = 0;
    }
    */
    
    if (apt_rto == 1) {
      //double log_psl_old_new = log(arma::normpdf(sel_cof_A_chn(i - 1), sel_cof_A_chn(i), sel_cof_sd)) + log(arma::normpdf(sel_cof_B_chn(i - 1), sel_cof_B_chn(i), sel_cof_sd)) + log(arma::normpdf(rec_rat_chn(i - 1), rec_rat_chn(i), rec_rat_sd));
      //double log_psl_new_old = log(arma::normpdf(sel_cof_A_chn(i), sel_cof_A_chn(i - 1), sel_cof_sd)) + log(arma::normpdf(sel_cof_B_chn(i), sel_cof_B_chn(i - 1), sel_cof_sd)) + log(arma::normpdf(rec_rat_chn(i), rec_rat_chn(i - 1), rec_rat_sd));
      
      calculateLogLikelihood_PhasedChr_arma(log_lik_chn(i), sel_cof_A_chn(i), dom_par_A_chn(i), sel_cof_B_chn(i), dom_par_B_chn(i), rec_rat_chn(i), pop_siz_chn(i), smp_gen, smp_chr_cnt, smp_hap_cnt, ptn_num, pcl_num);
      
      // calculate the acceptance ratio
      //apt_rto = exp((log_pri_chn(i) + log_lik_chn(i) + log_psl_old_new) - (log_pri_chn(i - 1) + log_lik_chn(i - 1) + log_psl_new_old));
      apt_rto = exp((log_pri_chn(i) + log_lik_chn(i)) - (log_pri_chn(i - 1) + log_lik_chn(i - 1)));
      if (arma::randu() > apt_rto) {
        sel_cof_A_chn(i) = sel_cof_A_chn(i - 1);
        sel_cof_B_chn(i) = sel_cof_B_chn(i - 1);
        //rec_rat_chn(i) = rec_rat_chn(i - 1);
        log_lik_chn(i) = log_lik_chn(i - 1);
      }
    } else {
      sel_cof_A_chn(i) = sel_cof_A_chn(i - 1);
      sel_cof_B_chn(i) = sel_cof_B_chn(i - 1);
      //rec_rat_chn(i) = rec_rat_chn(i - 1);
      log_lik_chn(i) = log_lik_chn(i - 1);
    }
  }
  
  return List::create(Named("sel_cof_A_chn", sel_cof_A_chn), 
                      Named("dom_par_A_chn", dom_par_A_chn), 
                      Named("sel_cof_B_chn", sel_cof_B_chn), 
                      Named("dom_par_B_chn", dom_par_B_chn), 
                      Named("rec_rat_chn", rec_rat_chn), 
                      Named("pop_siz_chn", pop_siz_chn));
}
/*************************/
