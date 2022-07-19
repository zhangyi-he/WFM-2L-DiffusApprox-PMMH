// Detecting and quantifying natural selection at two linked loci from time series data of allele frequencies with forward-in-time simulations
// Zhangyi He, Xiaoyang Dai, Mark Beaumont and Feng Yu

// version 1.1
// Two linked loci with missing genotypes

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
  
  for (arma::uword t = 1; t < arma::uword(lst_gen - int_gen) + 1; t++) {
    // random union of gametes
    gen_frq = hap_frq * hap_frq.t();
    
    // viability selection
    gen_frq = (fts % gen_frq) / arma::as_scalar(sum(sum(fts % gen_frq, 0), 1));
    
    // meiosis (calculate the sampling probability)
    arma::dcolvec prob = arma::zeros<arma::dcolvec>(4);
    for (arma::uword i = 0; i < 4; i++) {
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
    for (arma::uword i = 0; i < 4; i++) {
      if (hap_frq_pth(i, t) < 0) {
        hap_frq_pth(i, t) = 0;
      } 
      if (hap_frq_pth(i, t) > 1) {
        hap_frq_pth(i, t) = 1;
      }
    }
    hap_frq_pth.col(t) = hap_frq_pth.col(t) / sum(hap_frq_pth.col(t));
  }
  
  return hap_frq_pth;
}
/*************************/


/********** BPF **********/
// Calculate the possible haplotype counts in the sample
// [[Rcpp::export]]
arma::imat calculateHaploCnt_arma(const int& smp_siz, const arma::icolvec& smp_cnt) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  arma::imat smp_hap_cnt = arma::zeros<arma::imat>(4, 1);
  
  for (int i = 0; i <= min(smp_siz - smp_cnt(1), smp_siz - smp_cnt(3)); i++) {
    for (int j = 0; j <= smp_siz - smp_cnt(1) - i; j++) {
      if (i + j >= smp_cnt(0)) {
        for (int k = 0; k <= smp_siz - smp_cnt(3) - i; k++) {
          if (i + k >= smp_cnt(2)) {
            if (i + j + k <= smp_siz) {
              smp_hap_cnt(0, 0) = i;
              smp_hap_cnt(1, 0) = j;
              smp_hap_cnt(2, 0) = k;
              smp_hap_cnt(3, 0) = smp_siz - i - j - k;
              smp_hap_cnt.insert_cols(0, 1);
            }
          }
        }
      }
    }
  }
  smp_hap_cnt.shed_cols(0, 0);
  
  return smp_hap_cnt;
}

// Calculate the multinomial probabilities
// [[Rcpp::export]]
double calculateMultinomProb_arma(const arma::icolvec& smp_cnt, const int& smp_siz, const arma::dcolvec& pop_frq) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  if (arma::any(pop_frq == 0)) {
    if (arma::any(smp_cnt.elem(arma::find(pop_frq == 0)) != 0)) {
      return 0;
    }
    
    arma::icolvec smp_cnt_nonzero = smp_cnt.elem(arma::find(pop_frq != 0));
    arma::dcolvec pop_frq_nonzero = pop_frq.elem(arma::find(pop_frq != 0));
    
    return exp(lgamma(smp_siz + 1) + sum(arma::conv_to<arma::dcolvec>::from(smp_cnt_nonzero) % log(pop_frq_nonzero) - lgamma(arma::conv_to<arma::dcolvec>::from(smp_cnt_nonzero) + 1)));
  } else {
    return exp(lgamma(smp_siz + 1) + sum(arma::conv_to<arma::dcolvec>::from(smp_cnt) % log(pop_frq) - lgamma(arma::conv_to<arma::dcolvec>::from(smp_cnt) + 1)));
  }
}

// Calculate the emission probabilities
// [[Rcpp::export]]
double calculateEmissionProb_arma(const arma::icolvec& smp_cnt, const arma::icolvec& smp_ale_cnt, const arma::icolvec& smp_hap_cnt, const int& smp_siz, const arma::dcolvec& pop_hap_frq, const double& phi) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  double emn_prob = calculateMultinomProb_arma(smp_hap_cnt, smp_siz, pop_hap_frq);
  emn_prob = emn_prob * R::dbinom(smp_cnt(0), smp_ale_cnt(0), 1 - phi, false);
  emn_prob = emn_prob * R::dbinom(smp_cnt(1), smp_ale_cnt(1), 1 - phi, false);
  emn_prob = emn_prob * R::dbinom(smp_cnt(2), smp_ale_cnt(2), 1 - phi, false);
  emn_prob = emn_prob * R::dbinom(smp_cnt(3), smp_ale_cnt(3), 1 - phi, false);
  
  return emn_prob;
}

// Initialise the particles in the particle filter (uniform generation from the flat Dirichlet distribution)
// [[Rcpp::export]]
arma::dmat initialiseParticle_arma(const arma::uword& pcl_num) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  NumericMatrix part(pcl_num, 4);
  for (int j = 0; j < 4; j++) {
    part(_, j) = rgamma(pcl_num, 1.0, 1.0);
  }
  for (int i = 0; i < pcl_num; i++) {
    part(i, _) = part(i, _) / sum(part(i, _));
  }
  
  return as<arma::dmat>(transpose(part));
}

// Run the bootstrap particle filter
// [[Rcpp::export]]
List runBPF_arma(const double& sel_cof_A, const double& dom_par_A, const double& sel_cof_B, const double& dom_par_B, const double& rec_rat, const int& pop_siz, const arma::irowvec& smp_gen, const arma::irowvec& smp_siz, const arma::imat& smp_cnt, const arma::uword& ptn_num, const arma::uword& pcl_num) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  double lik = 1;
  
  arma::dmat wght = arma::zeros<arma::dmat>(pcl_num, smp_gen.n_elem);
  arma::dcube part_pre = arma::zeros<arma::dcube>(4, pcl_num, smp_gen.n_elem);
  arma::dcube part_pst = arma::zeros<arma::dcube>(4, pcl_num, smp_gen.n_elem);
  
  arma::dcolvec wght_tmp = arma::zeros<arma::dcolvec>(pcl_num);
  arma::dmat part_tmp = arma::zeros<arma::dmat>(4, pcl_num);
  
  double phi = 1.0 - (double) arma::accu(smp_cnt) / (double) arma::accu(smp_siz) / 2;
  
  if (phi == 0) {
    // initialise the particles
    cout << "generation: " << smp_gen(0) << endl;
    wght_tmp = arma::zeros<arma::dcolvec>(pcl_num);
    part_tmp = initialiseParticle_arma(pcl_num);
    arma::imat smp_hap_cnt = calculateHaploCnt_arma(smp_siz(0), smp_cnt.col(0));
    for (arma::uword i = 0; i < pcl_num; i++) {
      for (arma::uword j = 0; j < smp_hap_cnt.n_cols; j++) {
        wght_tmp(i) = wght_tmp(i) + calculateMultinomProb_arma(smp_hap_cnt.col(j), smp_siz(0), part_tmp.col(i));
      }
    }
    
    if (arma::sum(wght_tmp) > 0) {
      arma::dcolvec prob = arma::normalise(wght_tmp, 1);
      arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
      arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, TRUE, prob);
      
      lik = lik * arma::mean(wght_tmp);
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
      arma::imat smp_hap_cnt = calculateHaploCnt_arma(smp_siz(k), smp_cnt.col(k));
      for (arma::uword i = 0; i < pcl_num; i++) {
        arma::dmat path = simulateTLWFDS_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, part_tmp.col(i), smp_gen(k - 1), smp_gen(k), ptn_num);
        part_tmp.col(i) = arma::vectorise(path.tail_cols(1), 0);
        for (arma::uword j = 0; j < smp_hap_cnt.n_cols; j++) {
          wght_tmp(i) = wght_tmp(i) + calculateMultinomProb_arma(smp_hap_cnt.col(j), smp_siz(k), part_tmp.col(i));
        }
      }
      
      if (arma::sum(wght_tmp) > 0) {
        arma::dcolvec prob = arma::normalise(wght_tmp, 1);
        arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
        arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, TRUE, prob);
        
        lik = lik * arma::mean(wght_tmp);
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
  } else {
    // initialise the particles
    cout << "generation: " << smp_gen(0) << endl;
    wght_tmp = arma::zeros<arma::dcolvec>(pcl_num);
    part_tmp = initialiseParticle_arma(pcl_num);
    arma::imat smp_hap_cnt = calculateHaploCnt_arma(smp_siz(0), smp_cnt.col(0));
    arma::imat smp_ale_cnt = smp_hap_cnt;
    smp_ale_cnt.row(0) = smp_hap_cnt.row(0) + smp_hap_cnt.row(1);
    smp_ale_cnt.row(1) = smp_hap_cnt.row(2) + smp_hap_cnt.row(3);
    smp_ale_cnt.row(2) = smp_hap_cnt.row(0) + smp_hap_cnt.row(2);
    smp_ale_cnt.row(3) = smp_hap_cnt.row(1) + smp_hap_cnt.row(3);
    for (arma::uword i = 0; i < pcl_num; i++) {
      for (arma::uword j = 0; j < smp_hap_cnt.n_cols; j++) {
        wght_tmp(i) = wght_tmp(i) + calculateEmissionProb_arma(smp_cnt.col(0), smp_ale_cnt.col(j), smp_hap_cnt.col(j), smp_siz(0), part_tmp.col(i), phi);
      }
    }
    
    if (arma::sum(wght_tmp) > 0) {
      arma::dcolvec prob = arma::normalise(wght_tmp, 1);
      arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
      arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, TRUE, prob);
      
      lik = lik * arma::mean(wght_tmp);
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
      arma::imat smp_hap_cnt = calculateHaploCnt_arma(smp_siz(k), smp_cnt.col(k));
      arma::imat smp_ale_cnt = smp_hap_cnt;
      smp_ale_cnt.row(0) = smp_hap_cnt.row(0) + smp_hap_cnt.row(1);
      smp_ale_cnt.row(1) = smp_hap_cnt.row(2) + smp_hap_cnt.row(3);
      smp_ale_cnt.row(2) = smp_hap_cnt.row(0) + smp_hap_cnt.row(2);
      smp_ale_cnt.row(3) = smp_hap_cnt.row(1) + smp_hap_cnt.row(3);
      for (arma::uword i = 0; i < pcl_num; i++) {
        arma::dmat path = simulateTLWFDS_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, part_tmp.col(i), smp_gen(k - 1), smp_gen(k), ptn_num);
        part_tmp.col(i) = arma::vectorise(path.tail_cols(1), 0);
        for (arma::uword j = 0; j < smp_hap_cnt.n_cols; j++) {
          wght_tmp(i) = wght_tmp(i) + calculateEmissionProb_arma(smp_cnt.col(k), smp_ale_cnt.col(j), smp_hap_cnt.col(j), smp_siz(k), part_tmp.col(i), phi);
        }
      }
      
      if (arma::sum(wght_tmp) > 0) {
        arma::dcolvec prob = arma::normalise(wght_tmp, 1);
        arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
        arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, TRUE, prob);
        
        lik = lik * arma::mean(wght_tmp);
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
  }
  
  return List::create(Named("lik", lik),
                      Named("wght", wght), 
                      Named("part_pre_resmp", part_pre), 
                      Named("part_pst_resmp", part_pst));
}
/*************************/


/********** PMMH **********/
// Calculate the log-likelihood using the bootstrap particle filter
// [[Rcpp::export]]
double calculateLogLikelihood_arma(const double& sel_cof_A, const double& dom_par_A, const double& sel_cof_B, const double& dom_par_B, const double& rec_rat, const int& pop_siz, const arma::irowvec& smp_gen, const arma::irowvec& smp_siz, const arma::imat& smp_cnt, const arma::field<arma::imat>& ptl_ale_cnt, const arma::field<arma::imat>& ptl_hap_cnt, const double& phi, const arma::uword& ptn_num, const arma::uword& pcl_num) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  double log_lik = 0;
  
  arma::dcolvec wght = arma::zeros<arma::dcolvec>(pcl_num);
  arma::dmat part_pre = arma::zeros<arma::dmat>(4, pcl_num);
  arma::dmat part_pst = arma::zeros<arma::dmat>(4, pcl_num);
  
  if (phi == 0) {
    // initialise the particles
    part_pre = initialiseParticle_arma(pcl_num);
    arma::imat smp_hap_cnt = ptl_hap_cnt(0);
    for (arma::uword i = 0; i < pcl_num; i++) {
      for (arma::uword j = 0; j < smp_hap_cnt.n_cols; j++) {
        wght(i) = wght(i) + calculateMultinomProb_arma(smp_hap_cnt.col(j), smp_siz(0), part_pre.col(i));
      }
    }
    
    if (arma::mean(wght) > 0) {
      log_lik = log_lik + log(arma::mean(wght));
      arma::dcolvec prob = arma::normalise(wght, 1);
      arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
      arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, TRUE, prob);
      part_pst = part_pre.cols(indx);
    } else {
      log_lik = -(arma::datum::inf);
      return log_lik;
    }
    
    // run the bootstrap particle filter
    for (arma::uword k = 1; k < smp_gen.n_elem; k++) {
      wght = arma::zeros<arma::dcolvec>(pcl_num);
      arma::imat smp_hap_cnt = ptl_hap_cnt(k);
      for (arma::uword i = 0; i < pcl_num; i++) {
        arma::dmat path = simulateTLWFDS_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, part_pst.col(i), smp_gen(k - 1), smp_gen(k), ptn_num);
        part_pre.col(i) = arma::vectorise(path.tail_cols(1), 0);
        for (arma::uword j = 0; j < smp_hap_cnt.n_cols; j++) {
          wght(i) = wght(i) + calculateMultinomProb_arma(smp_hap_cnt.col(j), smp_siz(k), part_pre.col(i));
        }
      }
      
      if (arma::mean(wght) > 0) {
        log_lik = log_lik + log(arma::mean(wght));
        arma::dcolvec prob = arma::normalise(wght, 1);
        arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
        arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, TRUE, prob);
        part_pst = part_pre.cols(indx);
      } else {
        log_lik = -(arma::datum::inf);
        return log_lik;
      }
    }
  } else {
    // initialise the particles
    part_pre = initialiseParticle_arma(pcl_num);
    arma::imat smp_hap_cnt = ptl_hap_cnt(0);
    arma::imat smp_ale_cnt = ptl_ale_cnt(0);
    for (arma::uword i = 0; i < pcl_num; i++) {
      for (arma::uword j = 0; j < smp_hap_cnt.n_cols; j++) {
        wght(i) = wght(i) + calculateEmissionProb_arma(smp_cnt.col(0), smp_ale_cnt.col(j), smp_hap_cnt.col(j), smp_siz(0), part_pre.col(i), phi);
      }
    }
    
    if (arma::mean(wght) > 0) {
      log_lik = log_lik + log(arma::mean(wght));
      arma::dcolvec prob = arma::normalise(wght, 1);
      arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
      arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, TRUE, prob);
      part_pst = part_pre.cols(indx);
    } else {
      log_lik = -(arma::datum::inf);
      return log_lik;
    }
    
    // run the bootstrap particle filter
    for (arma::uword k = 1; k < smp_gen.n_elem; k++) {
      wght = arma::zeros<arma::dcolvec>(pcl_num);
      arma::imat smp_hap_cnt = ptl_hap_cnt(k);
      arma::imat smp_ale_cnt = ptl_ale_cnt(k);
      for (arma::uword i = 0; i < pcl_num; i++) {
        arma::dmat path = simulateTLWFDS_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, part_pst.col(i), smp_gen(k - 1), smp_gen(k), ptn_num);
        part_pre.col(i) = arma::vectorise(path.tail_cols(1), 0);
        for (arma::uword j = 0; j < smp_hap_cnt.n_cols; j++) {
          wght(i) = wght(i) + calculateEmissionProb_arma(smp_cnt.col(k), smp_ale_cnt.col(j), smp_hap_cnt.col(j), smp_siz(k), part_pre.col(i), phi);
        }
      }
      
      if (arma::mean(wght) > 0) {
        log_lik = log_lik + log(arma::mean(wght));
        arma::dcolvec prob = arma::normalise(wght, 1);
        arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
        arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, TRUE, prob);
        part_pst = part_pre.cols(indx);
      } else {
        log_lik = -(arma::datum::inf);
        return log_lik;
      }
    }
  }
  
  return log_lik;
}

// Calculate the optimal particle number in the particle marginal Metropolis-Hastings
// [[Rcpp::export]]
List calculateOptimalParticleNum_arma(const double& sel_cof_A, const double& dom_par_A, const double& sel_cof_B, const double& dom_par_B, const double& rec_rat, const int& pop_siz, const arma::irowvec& smp_gen, const arma::irowvec& smp_siz, const arma::imat& smp_cnt, const arma::uword& ptn_num, const arma::uword& pcl_num, const arma::uword& gap_num) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  double phi = 1.0 - (double) arma::accu(smp_cnt) / (double) arma::accu(smp_siz) / 2;
  
  arma::field<arma::imat> ptl_hap_cnt(smp_gen.n_elem);
  arma::field<arma::imat> ptl_ale_cnt(smp_gen.n_elem);
  for (arma::uword k = 0; k < smp_gen.n_elem; k++) {
    arma::imat smp_hap_cnt = calculateHaploCnt_arma(smp_siz(k), smp_cnt.col(k));
    ptl_hap_cnt(k) = smp_hap_cnt;
    
    arma::imat smp_ale_cnt = smp_hap_cnt;
    smp_ale_cnt.row(0) = smp_hap_cnt.row(0) + smp_hap_cnt.row(1);
    smp_ale_cnt.row(1) = smp_hap_cnt.row(2) + smp_hap_cnt.row(3);
    smp_ale_cnt.row(2) = smp_hap_cnt.row(0) + smp_hap_cnt.row(2);
    smp_ale_cnt.row(3) = smp_hap_cnt.row(1) + smp_hap_cnt.row(3);
    ptl_ale_cnt(k) = smp_ale_cnt;
  }
  
  arma::drowvec log_lik(300);
  for (arma::uword i = 0; i < 300; i++) {
    log_lik(i) = calculateLogLikelihood_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, smp_gen, smp_siz, smp_cnt, ptl_ale_cnt, ptl_hap_cnt, phi, ptn_num, pcl_num);
  }
  
  arma::drowvec log_lik_sdv(1);
  log_lik_sdv(0) = arma::stddev(log_lik);
  log_lik_sdv.print();
  arma::urowvec opt_pcl_num(1);
  opt_pcl_num(0) = pcl_num;
  
  if (log_lik_sdv(0) > 1.7) {
    while (log_lik_sdv(0) > 1.0) {
      opt_pcl_num.insert_cols(0, 1);
      opt_pcl_num(0) = opt_pcl_num(1) + gap_num;
      for (arma::uword i = 0; i < 300; i++) {
        log_lik(i) = calculateLogLikelihood_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, smp_gen, smp_siz, smp_cnt, ptl_ale_cnt, ptl_hap_cnt, phi, ptn_num, opt_pcl_num(0));
      }
      log_lik_sdv.insert_cols(0, 1);
      log_lik_sdv(0) = arma::stddev(log_lik);
      log_lik_sdv.print();
    }
    opt_pcl_num = arma::reverse(opt_pcl_num);
    log_lik_sdv = arma::reverse(log_lik_sdv);
  } else if (log_lik_sdv(0) < 1.0) {
    while (log_lik_sdv(0) < 1.7 && opt_pcl_num(0) > gap_num) {
      opt_pcl_num.insert_cols(0, 1);
      opt_pcl_num(0) = opt_pcl_num(1) - gap_num;
      for (arma::uword i = 0; i < 300; i++) {
        log_lik(i) = calculateLogLikelihood_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, smp_gen, smp_siz, smp_cnt, ptl_ale_cnt, ptl_hap_cnt, phi, ptn_num, opt_pcl_num(0));
      }
      log_lik_sdv.insert_cols(0, 1);
      log_lik_sdv(0) = arma::stddev(log_lik);
      log_lik_sdv.print();
    }
  } else {
    while (log_lik_sdv(0) > 1.0) {
      opt_pcl_num.insert_cols(0, 1);
      opt_pcl_num(0) = opt_pcl_num(1) + gap_num;
      for (arma::uword i = 0; i < 300; i++) {
        log_lik(i) = calculateLogLikelihood_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, smp_gen, smp_siz, smp_cnt, ptl_ale_cnt, ptl_hap_cnt, phi, ptn_num, opt_pcl_num(0));
      }
      log_lik_sdv.insert_cols(0, 1);
      log_lik_sdv(0) = arma::stddev(log_lik);
      log_lik_sdv.print();
    }
    opt_pcl_num = arma::reverse(opt_pcl_num);
    log_lik_sdv = arma::reverse(log_lik_sdv);
    
    while (log_lik_sdv(0) < 1.7 && opt_pcl_num(0) > gap_num) {
      opt_pcl_num.insert_cols(0, 1);
      opt_pcl_num(0) = opt_pcl_num(1) - gap_num;
      for (arma::uword i = 0; i < 300; i++) {
        log_lik(i) = calculateLogLikelihood_arma(sel_cof_A, dom_par_A, sel_cof_B, dom_par_B, rec_rat, pop_siz, smp_gen, smp_siz, smp_cnt, ptl_ale_cnt, ptl_hap_cnt, phi, ptn_num, opt_pcl_num(0));
      }
      log_lik_sdv.insert_cols(0, 1);
      log_lik_sdv(0) = arma::stddev(log_lik);
      log_lik_sdv.print();
    }
  }
  
  return List::create(Named("opt_pcl_num", opt_pcl_num), 
                      Named("log_lik_sdv", log_lik_sdv));
}

// Run the particle marginal Metropolis-Hastings
//[[Rcpp::export]]
List runPMMH_arma(const double& sel_cof_A, const double& dom_par_A, const double& sel_cof_B, const double& dom_par_B, const double& rec_rat, const int& pop_siz, const arma::irowvec& smp_gen, const arma::irowvec& smp_siz, const arma::imat& smp_cnt, const arma::uword& ptn_num, const arma::uword& pcl_num, const arma::uword& itn_num) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  double phi = 1.0 - (double) arma::accu(smp_cnt) / (double) arma::accu(smp_siz) / 2;
  
  arma::field<arma::imat> ptl_hap_cnt(smp_gen.n_elem);
  arma::field<arma::imat> ptl_ale_cnt(smp_gen.n_elem);
  for (arma::uword k = 0; k < smp_gen.n_elem; k++) {
    arma::imat smp_hap_cnt = calculateHaploCnt_arma(smp_siz(k), smp_cnt.col(k));
    ptl_hap_cnt(k) = smp_hap_cnt;
    
    arma::imat smp_ale_cnt = smp_hap_cnt;
    smp_ale_cnt.row(0) = smp_hap_cnt.row(0) + smp_hap_cnt.row(1);
    smp_ale_cnt.row(1) = smp_hap_cnt.row(2) + smp_hap_cnt.row(3);
    smp_ale_cnt.row(2) = smp_hap_cnt.row(0) + smp_hap_cnt.row(2);
    smp_ale_cnt.row(3) = smp_hap_cnt.row(1) + smp_hap_cnt.row(3);
    ptl_ale_cnt(k) = smp_ale_cnt;
  }
  
  arma::drowvec sel_cof_A_chn = arma::zeros<arma::drowvec>(itn_num);
  arma::drowvec sel_cof_B_chn = arma::zeros<arma::drowvec>(itn_num);
  
  //arma::drowvec log_pri_chn = arma::zeros<arma::drowvec>(itn_num);
  arma::drowvec log_lik_chn = arma::zeros<arma::drowvec>(itn_num);
  
  double sel_cof_sd = 1e-02;
  
  // initialise the population genetic parameters
  cout << "iteration: " << 1 << endl;
  // take the uniform prior and fix the initial value of the selection coefficient to zero
  // or take the beta prior with alpha = 1 and beta = 3
  sel_cof_A_chn(0) = sel_cof_A;
  sel_cof_B_chn(0) = sel_cof_B;
  
  log_lik_chn(0) = calculateLogLikelihood_arma(sel_cof_A_chn(0), dom_par_A, sel_cof_B_chn(0), dom_par_B, rec_rat, pop_siz, smp_gen, smp_siz, smp_cnt, ptl_ale_cnt, ptl_hap_cnt, phi, ptn_num, pcl_num);
  
  double apt_rto = 0;
  for (arma::uword i = 1; i < itn_num; i++) {
    cout << "iteration: " << i + 1 << endl;
    
    // draw the candidates of the selection coefficients from the random walk proposal
    sel_cof_A_chn(i) = sel_cof_A_chn(i - 1) + sel_cof_sd * arma::randn();
    sel_cof_B_chn(i) = sel_cof_B_chn(i - 1) + sel_cof_sd * arma::randn();
    
    if (sel_cof_A_chn(i) > 1 || sel_cof_B_chn(i) > 1) {
      sel_cof_A_chn(i) = sel_cof_A_chn(i - 1);
      sel_cof_B_chn(i) = sel_cof_B_chn(i - 1);
      log_lik_chn(i) = log_lik_chn(i - 1);
    } else {
      // calculate the proposal
      //double log_psl_old_new = log(arma::normpdf(sel_cof_A_chn(i - 1), sel_cof_A_chn(i), sel_cof_sd)) + log(arma::normpdf(sel_cof_B_chn(i - 1), sel_cof_B_chn(i), sel_cof_sd));
      //double log_psl_new_old = log(arma::normpdf(sel_cof_A_chn(i), sel_cof_A_chn(i - 1), sel_cof_sd)) + log(arma::normpdf(sel_cof_B_chn(i), sel_cof_B_chn(i - 1), sel_cof_sd));
      
      // calculate the likelihood
      log_lik_chn(i) = calculateLogLikelihood_arma(sel_cof_A_chn(i), dom_par_A, sel_cof_B_chn(i), dom_par_B, rec_rat, pop_siz, smp_gen, smp_siz, smp_cnt, ptl_ale_cnt, ptl_hap_cnt, phi, ptn_num, pcl_num);
      
      // calculate the acceptance ratio
      apt_rto = exp(log_lik_chn(i) - log_lik_chn(i - 1));
      //apt_rto = exp((log_pri_chn(i) + log_lik_chn(i) + log_psl_old_new) - (log_pri_chn(i - 1) + log_lik_chn(i - 1) + log_psl_new_old));
      
      if (arma::randu() > apt_rto) {
        sel_cof_A_chn(i) = sel_cof_A_chn(i - 1);
        sel_cof_B_chn(i) = sel_cof_B_chn(i - 1);
        log_lik_chn(i) = log_lik_chn(i - 1);
      }
    }
  }
  
  return List::create(Named("sel_cof_A_chn", sel_cof_A_chn), 
                      Named("sel_cof_B_chn", sel_cof_B_chn));
}
/*************************/
