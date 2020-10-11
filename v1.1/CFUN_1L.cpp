// Detecting and quantifying natural selection at two linked loci from time series data of allele frequencies with forward-in-time simulations
// Zhangyi He, Xiaoyang Dai, Mark Beaumont and Feng Yu

// version 1.1
// Single-loucs case (N/A is allowed)

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
// Simulate the mutant allele frequency trajectory according to the one-locus Wright-Fisher model with selection
// [[Rcpp::export]]
arma::drowvec simulateOLWFMS_arma(const double& sel_cof, const double& dom_par, const int& pop_siz, const double& int_frq, const int& int_gen, const int& lst_gen) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  // declare the mutant allele frequency trajectory
  arma::drowvec ale_frq_pth(arma::uword(lst_gen - int_gen) + 1);
  
  // initialise the mutant allele frequency in generation 0
  ale_frq_pth(0) = int_frq;
  
  // declare the fitness
  arma::dcolvec fts = arma::zeros<arma::dcolvec>(3);
  fts(0) = 1.0;
  fts(1) = 1.0 - sel_cof * dom_par;
  fts(2) = 1.0 - sel_cof;
  
  // declare and initialise the genotype frequencies during a single generation of the life cycle
  arma::dcolvec gen_frq = arma::zeros<arma::dcolvec>(3);
  
  for (arma::uword t = 1; t < arma::uword(lst_gen - int_gen) + 1; t++) {
    // random union of gametes
    gen_frq(0) = ale_frq_pth(t - 1) * ale_frq_pth(t - 1);
    gen_frq(1) = 2 * ale_frq_pth(t - 1) * (1 - ale_frq_pth(t - 1));
    gen_frq(2) = (1 - ale_frq_pth(t - 1)) * (1 - ale_frq_pth(t - 1));
    
    // viability selection
    gen_frq = (fts % gen_frq) / arma::accu(fts % gen_frq);
    
    // meiosis (calculate the sampling probability)
    double prob = gen_frq(0) + gen_frq(1) / 2;
    
    // reproduction (the Wright-Fisher sampling)
    ale_frq_pth(t) = R::rbinom(2 * pop_siz, prob) / 2 / pop_siz;
  }
  
  return ale_frq_pth;
}
/*************************/


/********** WFD **********/
// Simulate the mutant allele frequency trajectory according to the one-locus Wright-Fisher diffusion with selection using the Euler-Maruyama method
// [[Rcpp::export]]
arma::drowvec simulateOLWFDS_arma(const double& sel_cof, const double& dom_par, const int& pop_siz, const double& int_frq, const int& int_gen, const int& lst_gen, const arma::uword& ptn_num) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  // rescale the selection coefficient
  double scl_sel_cof = 2 * pop_siz * sel_cof;
  
  // declare delta t
  double dt = 1.0 / (2 * pop_siz) / ptn_num;
  // declare delta W
  arma::drowvec dW = pow(dt, 0.5) * arma::randn<arma::drowvec>(arma::uword(lst_gen - int_gen) * ptn_num);
  
  // declare the mutant allele frequency trajectory
  arma::drowvec ale_frq_pth(arma::uword(lst_gen - int_gen) * ptn_num + 1);
  
  // initialise the mutant allele frequency in generation 0
  ale_frq_pth(0) = int_frq;
  
  for (arma::uword t = 1; t < arma::uword(lst_gen - int_gen) * ptn_num + 1; t++) {
    // calculate the drift coefficient
    double mu = scl_sel_cof * ale_frq_pth(t - 1) * (1 - ale_frq_pth(t - 1)) * ((1 - dom_par) - (1 - 2 * dom_par) * ale_frq_pth(t - 1));
    
    // calculate the diffusion coefficient
    double sigma = pow(ale_frq_pth(t - 1) * (1 - ale_frq_pth(t - 1)), 0.5);
    
    // proceed the Euler-Maruyama scheme
    ale_frq_pth(t) = ale_frq_pth(t - 1) + mu * dt + sigma * dW(t - 1);
    
    // remove the noise from the numerical techniques
    if (ale_frq_pth(t) < 0) {
      ale_frq_pth(t) = 0;
    }
    if (ale_frq_pth(t) > 1) {
      ale_frq_pth(t) = 1;
    }
  }
  
  return ale_frq_pth;
}
/*************************/


/********** BPF **********/
// Calculate the emission probabilities
// [[Rcpp::export]]
double calculateEmissionProb_arma(const arma::icolvec& smp_cnt, const arma::icolvec& ptl_cnt, const int& smp_siz, const double& pop_frq, const double& phi) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  double emn_prob = R::dbinom(ptl_cnt(0), smp_siz, pop_frq, false);
  emn_prob = emn_prob * R::dbinom(smp_cnt(0), ptl_cnt(0), 1 - phi, false);
  emn_prob = emn_prob * R::dbinom(smp_cnt(1), ptl_cnt(1), 1 - phi, false);
  
  return emn_prob;
}

// Run the bootstrap particle filter with the one-locus Wright-Fisher diffusion with selection
// [[Rcpp::export]]
List runBPF_arma(const double& sel_cof, const double& dom_par, const int& pop_siz, const arma::irowvec& smp_gen, const arma::irowvec& smp_siz, const arma::imat& smp_cnt, const arma::uword& ptn_num, const arma::uword& pcl_num) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  double lik = 1;
  
  arma::dmat wght = arma::zeros<arma::dmat>(pcl_num, smp_gen.n_elem);
  arma::dmat part_pre = arma::zeros<arma::dmat>(pcl_num, smp_gen.n_elem);
  arma::dmat part_pst = arma::zeros<arma::dmat>(pcl_num, smp_gen.n_elem);
  
  arma::dcolvec wght_tmp = arma::zeros<arma::dcolvec>(pcl_num);
  arma::dcolvec part_tmp = arma::randu<arma::dcolvec>(pcl_num);
  
  double phi = 1.0 - (double) arma::accu(smp_cnt) / (double) arma::accu(smp_siz);
  
  if (phi == 0) {
    // initialise the particles
    cout << "generation: " << smp_gen(0) << endl;
    for (arma::uword i = 0; i < pcl_num; i++) {
      wght_tmp(i) = R::dbinom(smp_cnt(0, 0), smp_siz(0), part_tmp(i), false);
    }
    
    if (arma::sum(wght_tmp) > 0) {
      arma::dcolvec prob = arma::normalise(wght_tmp, 1);
      arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
      arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, true, prob);
      
      lik = lik * arma::mean(wght_tmp);
      wght.col(0) = wght_tmp;
      part_pre.col(0) = part_tmp;
      part_pst.col(0) = part_tmp.elem(indx);
    } else {
      lik = 0;
      wght.shed_cols(0, smp_gen.n_elem - 1);
      part_pre.shed_cols(0, smp_gen.n_elem - 1);
      part_pst.shed_cols(0, smp_gen.n_elem - 1);
      
      return List::create(Named("lik", lik), 
                          Named("wght", wght), 
                          Named("part_pre_resmp", part_pre), 
                          Named("part_pst_resmp", part_pst));
    }
    
    // run the bootstrap particle filter
    for (arma::uword k = 1; k < smp_gen.n_elem; k++) {
      cout << "generation: " << smp_gen(k) << endl;
      wght_tmp = arma::zeros<arma::dcolvec>(pcl_num);
      for (arma::uword i = 0; i < pcl_num; i++) {
        arma::drowvec path = simulateOLWFDS_arma(sel_cof, dom_par, pop_siz, part_pst(i, k - 1), smp_gen(k - 1), smp_gen(k), ptn_num);
        part_tmp(i) = arma::as_scalar(path.tail(1));
        wght_tmp(i) = R::dbinom(smp_cnt(0, k), smp_siz(k), part_tmp(i), false);
      }
      
      if (arma::sum(wght_tmp) > 0) {
        arma::dcolvec prob = arma::normalise(wght_tmp, 1);
        arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
        arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, true, prob);
        
        lik = lik * arma::mean(wght_tmp);
        wght.col(k) = wght_tmp;
        part_pre.col(k) = part_tmp;
        part_pst.col(k) = part_tmp.elem(indx);
      } else {
        lik = 0;
        wght.shed_cols(k, smp_gen.n_elem - 1);
        part_pre.shed_cols(k, smp_gen.n_elem - 1);
        part_pst.shed_cols(k, smp_gen.n_elem - 1);
        break;
      }
    }
  } else {
    // initialise the particles
    cout << "generation: " << smp_gen(0) << endl;
    arma::imat smp_ale_cnt = arma::zeros<arma::imat>(2, smp_siz(0) - smp_cnt(0, 0) - smp_cnt(1, 0) + 1);
    smp_ale_cnt.row(0) = arma::linspace<arma::irowvec>(smp_cnt(0, 0), smp_siz(0) - smp_cnt(1, 0), smp_siz(0) - smp_cnt(0, 0) - smp_cnt(1, 0) + 1);
    smp_ale_cnt.row(1) = arma::linspace<arma::irowvec>(smp_siz(0) - smp_cnt(0, 0), smp_cnt(1, 0), smp_siz(0) - smp_cnt(0, 0) - smp_cnt(1, 0) + 1);
    for (arma::uword i = 0; i < pcl_num; i++) {
      for (arma::uword j = 0; j < smp_ale_cnt.n_cols; j++) {
        wght_tmp(i) = wght_tmp(i) + calculateEmissionProb_arma(smp_cnt.col(0), smp_ale_cnt.col(j), smp_siz(0), part_tmp(i), phi);
      }
    }
    
    if (arma::sum(wght_tmp) > 0) {
      arma::dcolvec prob = arma::normalise(wght_tmp, 1);
      arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
      arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, true, prob);
      
      lik = lik * arma::mean(wght_tmp);
      wght.col(0) = wght_tmp;
      part_pre.col(0) = part_tmp;
      part_pst.col(0) = part_tmp.elem(indx);
    } else {
      lik = 0;
      wght.shed_cols(0, smp_gen.n_elem - 1);
      part_pre.shed_cols(0, smp_gen.n_elem - 1);
      part_pst.shed_cols(0, smp_gen.n_elem - 1);
      
      return List::create(Named("lik", lik), 
                          Named("wght", wght), 
                          Named("part_pre_resmp", part_pre), 
                          Named("part_pst_resmp", part_pst));
    }
    
    // run the bootstrap particle filter
    for (arma::uword k = 1; k < smp_gen.n_elem; k++) {
      cout << "generation: " << smp_gen(k) << endl;
      wght_tmp = arma::zeros<arma::dcolvec>(pcl_num);
      arma::imat smp_ale_cnt = arma::zeros<arma::imat>(2, smp_siz(k) - smp_cnt(0, k) - smp_cnt(1, k) + 1);
      smp_ale_cnt.row(0) = arma::linspace<arma::irowvec>(smp_cnt(0, k), smp_siz(k) - smp_cnt(1, k), smp_siz(k) - smp_cnt(0, k) - smp_cnt(1, k) + 1);
      smp_ale_cnt.row(1) = arma::linspace<arma::irowvec>(smp_siz(k) - smp_cnt(0, k), smp_cnt(1, k), smp_siz(k) - smp_cnt(0, k) - smp_cnt(1, k) + 1);
      for (arma::uword i = 0; i < pcl_num; i++) {
        arma::drowvec path = simulateOLWFDS_arma(sel_cof, dom_par, pop_siz, part_pst(i, k - 1), smp_gen(k - 1), smp_gen(k), ptn_num);
        part_tmp(i) = arma::as_scalar(path.tail(1));
        for (arma::uword j = 0; j < smp_ale_cnt.n_cols; j++) {
          wght_tmp(i) = wght_tmp(i) + calculateEmissionProb_arma(smp_cnt.col(k), smp_ale_cnt.col(j), smp_siz(k), part_tmp(i), phi);
        }
      }
      
      if (arma::sum(wght_tmp) > 0) {
        arma::dcolvec prob = arma::normalise(wght_tmp, 1);
        arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
        arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, true, prob);
        
        lik = lik * arma::mean(wght_tmp);
        wght.col(k) = wght_tmp;
        part_pre.col(k) = part_tmp;
        part_pst.col(k) = part_tmp.elem(indx);
      } else {
        lik = 0;
        wght.shed_cols(k, smp_gen.n_elem - 1);
        part_pre.shed_cols(k, smp_gen.n_elem - 1);
        part_pst.shed_cols(k, smp_gen.n_elem - 1);
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
double calculateLogLikelihood_arma(const double& sel_cof, const double& dom_par, const int& pop_siz, const arma::irowvec& smp_gen, const arma::irowvec& smp_siz, const arma::imat& smp_cnt, const arma::field<arma::imat>& ptl_cnt, const double& phi, const arma::uword& ptn_num, const arma::uword& pcl_num) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  double log_lik = 0;
  
  arma::dcolvec wght = arma::zeros<arma::dcolvec>(pcl_num);
  arma::dcolvec part_pre = arma::zeros<arma::dcolvec>(pcl_num);
  arma::dcolvec part_pst = arma::zeros<arma::dcolvec>(pcl_num);
  
  if (phi == 0) {
    // initialise the particles
    part_pre = arma::randu<arma::dcolvec>(pcl_num);
    for (arma::uword i = 0; i < pcl_num; i++) {
      wght(i) = R::dbinom(smp_cnt(0, 0), smp_siz(0), part_pre(i), false);
    }
    
    if (arma::sum(wght) > 0) {
      log_lik = log_lik + log(arma::mean(wght));
      arma::dcolvec prob = arma::normalise(wght, 1);
      arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
      arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, true, prob);
      part_pst = part_pre.elem(indx);
    } else {
      log_lik = -(arma::datum::inf);
      return log_lik;
    }
    
    // run the bootstrap particle filter
    for (arma::uword k = 1; k < smp_gen.n_elem; k++) {
      wght = arma::zeros<arma::dcolvec>(pcl_num);
      for (arma::uword i = 0; i < pcl_num; i++) {
        arma::drowvec path = simulateOLWFDS_arma(sel_cof, dom_par, pop_siz, part_pst(i), smp_gen(k - 1), smp_gen(k), ptn_num);
        part_pre(i) = arma::as_scalar(path.tail(1));
        wght(i) = R::dbinom(smp_cnt(0, k), smp_siz(k), part_pre(i), false);
      }
      
      if (arma::sum(wght) > 0) {
        log_lik = log_lik + (log(arma::sum(wght)) - log(pcl_num));
        arma::dcolvec prob = arma::normalise(wght, 1);
        arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
        arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, true, prob);
        part_pst = part_pre.elem(indx);
      } else {
        log_lik = -(arma::datum::inf);
        return log_lik;
      }
    }
  } else {
    // initialise the particles
    part_pre = arma::randu<arma::dcolvec>(pcl_num);
    arma::imat smp_ale_cnt = ptl_cnt(0);
    for (arma::uword i = 0; i < pcl_num; i++) {
      for (arma::uword j = 0; j < smp_ale_cnt.n_cols; j++) {
        wght(i) = wght(i) + calculateEmissionProb_arma(smp_cnt.col(0), smp_ale_cnt.col(j), smp_siz(0), part_pre(i), phi);
      }
    }
    
    if (arma::sum(wght) > 0) {
      log_lik = log_lik + log(arma::mean(wght));
      arma::dcolvec prob = arma::normalise(wght, 1);
      arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
      arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, true, prob);
      part_pst = part_pre.elem(indx);
    } else {
      log_lik = -(arma::datum::inf);
      return log_lik;
    }
    
    // run the bootstrap particle filter
    for (arma::uword k = 1; k < smp_gen.n_elem; k++) {
      wght = arma::zeros<arma::dcolvec>(pcl_num);
      arma::imat smp_ale_cnt = ptl_cnt(k);
      for (arma::uword i = 0; i < pcl_num; i++) {
        arma::drowvec path = simulateOLWFDS_arma(sel_cof, dom_par, pop_siz, part_pst(i), smp_gen(k - 1), smp_gen(k), ptn_num);
        part_pre(i) = arma::as_scalar(path.tail(1));
        for (arma::uword j = 0; j < smp_ale_cnt.n_cols; j++) {
          wght(i) = wght(i) + calculateEmissionProb_arma(smp_cnt.col(k), smp_ale_cnt.col(j), smp_siz(k), part_pre(i), phi);
        }    
      }

      if (arma::sum(wght) > 0) {
        log_lik = log_lik + (log(arma::sum(wght)) - log(pcl_num));
        arma::dcolvec prob = arma::normalise(wght, 1);
        arma::ucolvec elem = arma::linspace<arma::ucolvec>(0, pcl_num - 1, pcl_num);
        arma::ucolvec indx = RcppArmadillo::sample(elem, pcl_num, true, prob);
        part_pst = part_pre.elem(indx);
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
List calculateOptimalParticleNum_arma(const double& sel_cof, const double& dom_par, const int& pop_siz, const arma::irowvec& smp_gen, const arma::irowvec& smp_siz, const arma::imat& smp_cnt, const arma::uword& ptn_num, const arma::uword& pcl_num, const arma::uword& gap_num) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  double phi = 1.0 - (double) arma::accu(smp_cnt) / (double) arma::accu(smp_siz);
  
  arma::field<arma::imat> ptl_cnt(smp_gen.n_elem);
  for (arma::uword k = 0; k < smp_gen.n_elem; k++) {
    arma::imat smp_ale_cnt = arma::zeros<arma::imat>(2, smp_siz(k) - smp_cnt(0, k) - smp_cnt(1, k) + 1);
    smp_ale_cnt.row(0) = arma::linspace<arma::irowvec>(smp_cnt(0, k), smp_siz(k) - smp_cnt(1, k), smp_siz(k) - smp_cnt(0, k) - smp_cnt(1, k) + 1);
    smp_ale_cnt.row(1) = arma::linspace<arma::irowvec>(smp_siz(k) - smp_cnt(0, k), smp_cnt(1, k), smp_siz(k) - smp_cnt(0, k) - smp_cnt(1, k) + 1);
    ptl_cnt(k) = smp_ale_cnt;
  }
  
  arma::drowvec log_lik(300);
  for (arma::uword i = 0; i < 300; i++) {
    log_lik(i) = calculateLogLikelihood_arma(sel_cof, dom_par, pop_siz, smp_gen, smp_siz, smp_cnt, ptl_cnt, phi, ptn_num, pcl_num);
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
        log_lik(i) = calculateLogLikelihood_arma(sel_cof, dom_par, pop_siz, smp_gen, smp_siz, smp_cnt, ptl_cnt, phi, ptn_num, opt_pcl_num(0));
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
        log_lik(i) = calculateLogLikelihood_arma(sel_cof, dom_par, pop_siz, smp_gen, smp_siz, smp_cnt, ptl_cnt, phi, ptn_num, opt_pcl_num(0));
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
        log_lik(i) = calculateLogLikelihood_arma(sel_cof, dom_par, pop_siz, smp_gen, smp_siz, smp_cnt, ptl_cnt, phi, ptn_num, opt_pcl_num(0));
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
        log_lik(i) = calculateLogLikelihood_arma(sel_cof, dom_par, pop_siz, smp_gen, smp_siz, smp_cnt, ptl_cnt, phi, ptn_num, opt_pcl_num(0));
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
List runPMMH_arma(const double& sel_cof, const double& dom_par, const int& pop_siz, const arma::irowvec& smp_gen, const arma::irowvec& smp_siz, const arma::imat& smp_cnt, const arma::uword& ptn_num, const arma::uword& pcl_num, const arma::uword& itn_num) {
  // ensure RNG gets set/reset
  RNGScope scope;
  
  double phi = 1.0 - (double) arma::accu(smp_cnt) / (double) arma::accu(smp_siz);
  
  arma::field<arma::imat> ptl_cnt(smp_gen.n_elem);
  for (arma::uword k = 0; k < smp_gen.n_elem; k++) {
    arma::imat smp_ale_cnt = arma::zeros<arma::imat>(2, smp_siz(k) - smp_cnt(0, k) - smp_cnt(1, k) + 1);
    smp_ale_cnt.row(0) = arma::linspace<arma::irowvec>(smp_cnt(0, k), smp_siz(k) - smp_cnt(1, k), smp_siz(k) - smp_cnt(0, k) - smp_cnt(1, k) + 1);
    smp_ale_cnt.row(1) = arma::linspace<arma::irowvec>(smp_siz(k) - smp_cnt(0, k), smp_cnt(1, k), smp_siz(k) - smp_cnt(0, k) - smp_cnt(1, k) + 1);
    ptl_cnt(k) = smp_ale_cnt;
  }
  
  arma::drowvec sel_cof_chn = arma::zeros<arma::drowvec>(itn_num);
  
  //arma::drowvec log_pri_chn = arma::zeros<arma::drowvec>(itn_num);
  arma::drowvec log_lik_chn = arma::zeros<arma::drowvec>(itn_num);
  
  double sel_cof_sd = 8e-03;
  
  // initialise the population genetic parameters
  cout << "iteration: " << 1 << endl;
  // take the uniform prior and fix the initial value of the selection coefficient to zero
  // or take the beta prior with alpha = 1 and beta = 3
  sel_cof_chn(0) = sel_cof;
  
  log_lik_chn(0) = calculateLogLikelihood_arma(sel_cof_chn(0), dom_par, pop_siz, smp_gen, smp_siz, smp_cnt, ptl_cnt, phi, ptn_num, pcl_num);
  
  double apt_rto = 0;
  for (arma::uword i = 1; i < itn_num; i++) {
    cout << "iteration: " << i + 1 << endl;
    
    // draw the candidate of the selection coefficient from the random walk proposal
    sel_cof_chn(i) = sel_cof_chn(i - 1) + sel_cof_sd * arma::randn();
    
    if (sel_cof_chn(i) > 1) {
      sel_cof_chn(i) = sel_cof_chn(i - 1);
      log_lik_chn(i) = log_lik_chn(i - 1);
    } else {
      // calculate the proposal
      //double log_psl_old_new = log(arma::normpdf(sel_cof_chn(i - 1), sel_cof_chn(i), sel_cof_sd));
      //double log_psl_new_old = log(arma::normpdf(sel_cof_chn(i), sel_cof_chn(i - 1), sel_cof_sd));
      
      // calculate the likelihood
      log_lik_chn(i) = calculateLogLikelihood_arma(sel_cof_chn(i), dom_par, pop_siz, smp_gen, smp_siz, smp_cnt, ptl_cnt, phi, ptn_num, pcl_num);
      
      // calculate the acceptance ratio
      apt_rto = exp(log_lik_chn(i) - log_lik_chn(i - 1));
      //apt_rto = exp((log_pri_chn(i) + log_lik_chn(i) + log_psl_old_new) - (log_pri_chn(i - 1) + log_lik_chn(i - 1) + log_psl_new_old));
      
      if (arma::randu() > apt_rto) {
        sel_cof_chn(i) = sel_cof_chn(i - 1);
        log_lik_chn(i) = log_lik_chn(i - 1);
      }
    }
  }
  
  return List::create(Named("sel_cof_chn", sel_cof_chn));
}
/*************************/
