/* Copyright (C) 2021 Janik Schuettler
 *
 * Three-point return correlations for single-file diffusion.
 *
 * This file is an EXTENSION of BetheSF by Alessio Lapolla and Aljaz Godec
 * (https://data.mendeley.com/datasets/3bs74vf72n/2), whose files carry their
 * own copyright and live alongside this one. Everything in THIS file is mine;
 * nothing in it is theirs.
 *
 * BetheSF diagonalises the Fokker-Planck operator of the full N-particle
 * single file via a coordinate Bethe ansatz and projects onto one tagged
 * particle, giving TWO-point propagators Q(x,t|x0). What is added here is the
 * THREE-point return correlation
 *
 *     C(x,tau,y,t) = Q(x,t+tau ; y,t | x) - Q(x,t+tau|x) Q(y,t|x)
 *
 * and the machinery needed to make it computable on a grid.
 *
 * These are member functions of SingleFile, so their declarations necessarily
 * live in SingleFile.hpp (in a block marked as mine); C++ permits a class's
 * members to be defined across several translation units, which is what keeps
 * the definitions separable like this.
 *
 * You may use, distribute and modify this code under the terms of the MIT
 * license, the same terms as the code it extends.
 */
#include <algorithm>
#include <cmath>
#include <vector>

#include "SingleFile.hpp"

namespace SingleFile {

std::vector<double> SingleFile::joint2dens(const double x, const std::vector<double> t, const double x0)
{
  int n_t = t.size();
  std::vector<double> g2p(n_t, 0.0);

  for(auto& eigenfunction : eigenfunction_store){// all these elements will have the same manybody eigenvalue
    double V0k=V0k_element(eigenfunction, x);
    double Vk0=Vk0_element(eigenfunction, x0);
    double lambda=lambda_fp(eigenfunction);
    for(int i_t = 0; i_t < n_t; i_t++){
      g2p[i_t] += V0k*Vk0*std::exp(-lambda*t[i_t]);
    }
  }
  return g2p;
}

double SingleFile::joint3dens(const double x2, const double t2, const double x1, const double t1, const double x0)
{
  double j3d=0.0;
  for(auto k : eigenfunction_store){
    for(auto l : eigenfunction_store){
      double V0k=V0k_element(k, x2);
      double Vkl=Vkl_element(k, l, x1);
      double Vl0=Vk0_element(l, x0);
      double lambda_k=lambda_fp(k);
      double lambda_l=lambda_fp(l);
      j3d+=V0k*Vkl*Vl0*std::exp(-lambda_k*(t2-t1)-lambda_l*t1);
    }
  }
  return j3d;
}

double SingleFile::correlation(const double x2, const double t2, const double x1, const double t1, const double x0)
{
  double corr=0.0;
  std::vector<int> zeros(N);
  double V00=V0k_element(zeros, x0);
  for(auto& k : eigenfunction_store){
    double V0k=V0k_element(k, x2);
    double Vk0=Vk0_element(k, x0);
    double lambda_k=lambda_fp(k);
    for(auto& l : eigenfunction_store){
      double Vkl=Vkl_element(k, l, x1);
      double Vl0=Vk0_element(l, x0);
      double V0l=V0k_element(l, x1);
      double lambda_l=lambda_fp(l);
      corr+=V0k*Vl0*( Vkl*V00 - V0l*Vk0*std::exp(-lambda_k*t1) )*std::exp(-lambda_k*(t2-t1)-lambda_l*t1);
    }
  }
  return corr / std::pow(V00, 2);
}

std::vector<std::vector<double>> SingleFile::correlation(const double x2, const std::vector<double> tau, const double x1, const std::vector<double> t1, const double x0)
{
  int n1 = t1.size(), n2 = tau.size();
  std::vector<std::vector<double>> corr(n1 , std::vector<double> (n2, 0.0));
  
  std::vector<int> zeros(N);
  double V00=V0k_element(zeros, x0);

  for(auto& k : eigenfunction_store){
    for(auto& l : eigenfunction_store){
      double V0k=V0k_element(k, x2);
      double Vk0=Vk0_element(k, x0);
      double Vkl=Vkl_element(k, l, x1);
      double Vl0=Vk0_element(l, x0);
      double V0l=V0k_element(l, x1);
      double lambda_k=lambda_fp(k);
      double lambda_l=lambda_fp(l);
      // for(auto& t1_ : t1){
      //   for(auto& t2_ : t2){
      for(int i1 = 0; i1 < n1; i1++){
        for(int i2 = 0; i2 < n2; i2++){
          double val = V0k*Vl0*( Vkl*V00 - V0l*Vk0*std::exp(-lambda_k*t1[i1]) )*std::exp(-lambda_k*tau[i2]-lambda_l*t1[i1]) / std::pow(V00, 2);
          corr[i1][i2] = corr[i1][i2] + val;
        }
      }
    }
  }
  return corr;
}

std::vector<std::vector<std::vector<std::vector<double>>>> SingleFile::correlation(const std::vector<double> x2, const std::vector<double> tau, 
  const double x1, const std::vector<double> t1, const std::vector<double> x0)
{
  int n_x0  = x0.size();
  int n_x2  = x2.size();
  int n_t1  = t1.size();
  int n_tau = tau.size();
  //                                                                                                                                                           (n_t1, std::vector<double>(n_tau, 0.0));
  //                                                                                                                   (n_x0, std::vector<std::vector<double>> (n_t1, std::vector<double>(n_tau, 0.0)));
  std::vector<std::vector<std::vector<std::vector<double>>>> corr(n_x0, std::vector<std::vector<std::vector<double>>> (n_x2, std::vector<std::vector<double>> (n_t1, std::vector<double>(n_tau, 0.0))));

  std::vector<int> zeros(N);

  for(auto& k : eigenfunction_store){
    for(auto& l : eigenfunction_store){
      double lambda_k=lambda_fp(k);
      double lambda_l=lambda_fp(l);

      double Vkl=Vkl_element(k, l, x1);
      double V0l=V0k_element(l, x1);

      for(int i_x0 = 0; i_x0 < n_x0; i_x0++){
        double V00=V0k_element(zeros, x0[i_x0]);
        double Vl0=Vk0_element(l, x0[i_x0]);
        double Vk0=Vk0_element(k, x0[i_x0]);

        for(int i_x2 = 0; i_x2 < n_x2; i_x2++){
          double V0k=V0k_element(k, x2[i_x2]);

          for(int i_t1 = 0; i_t1 < n_t1; i_t1++){
            for(int i_tau = 0; i_tau < n_tau; i_tau++){
              double val = V0k*Vl0*( Vkl*V00 - V0l*Vk0*std::exp(-lambda_k*t1[i_t1]) )*std::exp(-lambda_k*tau[i_tau]-lambda_l*t1[i_t1]) / std::pow(V00, 2);
              corr[i_x0][i_x2][i_t1][i_tau] = corr[i_x0][i_x2][i_t1][i_tau] + val;
            }
          }
        }
      }
    }
  }
  return corr;
}

std::vector<std::vector<double>> SingleFile::precompute_Vkl(const double x1, const unsigned n_threads,
 const std::vector<std::vector<int>>& eig_to_exclude) 
{
  const unsigned n_eigenfunctions = eigenfunction_store.size();
  std::vector<std::vector<double>> Vkl_precomputed(n_eigenfunctions, std::vector<double> (n_eigenfunctions));

#ifndef _OPENMP
  (void)n_threads;   // built without OpenMP: the pragma below is ignored
#endif

#pragma omp parallel for num_threads(n_threads) 
  for(unsigned k = 0; k < n_eigenfunctions; k++){
    for(unsigned l = 0; l < n_eigenfunctions; l++){
      auto eig_k = eigenfunction_store[k], eig_l = eigenfunction_store[l];

      if( std::find(eig_to_exclude.begin(), eig_to_exclude.end(), eig_k) == eig_to_exclude.end() || 
        std::find(eig_to_exclude.begin(), eig_to_exclude.end(), eig_l) == eig_to_exclude.end() ) {

        if (k == l){
          Vkl_precomputed[k][l] = Vkl_element(eig_k, eig_k, x1);
        }
        else {
          Vkl_precomputed[k][l] = Vkl_element(eig_k, eig_l, x1);
        }
      }
    }
  }

  return Vkl_precomputed;
}

/* Computes the three-point return correlation C, Eq. (17) of the thesis:
     C(x,tau,y,t) = Q(x,t+tau; y,t | x) - Q(x,t+tau|x) Q(y,t|x)
   over a grid of return points x, waiting times t1 and lag times tau, for a
   fixed intermediate point x1 (=y).

   The expensive overlap elements Vkl(x1) do not depend on x, t1 or tau, so they
   are precomputed once (in parallel) and reused across the whole grid; the
   remaining elements are hoisted as far out of the loop nest as they allow.

   eig_to_exclude drops all (k,l) pairs for which BOTH indices are in the given
   set. Running with a low-eigenvalue set excluded therefore yields the
   incremental contribution of the higher eigenvalue shell, so a convergence
   study can be split across independent jobs and summed.

   The innermost vector holds 5 components:
     [0] C = connected - disconnected   [1] connected   [2] disconnected
     [3] Q(y,t|x)                       [4] Q(x,t+tau|x)                     */
std::vector<std::vector<std::vector<std::vector<double>>>> SingleFile::correlation_looped(const std::vector<double>& x,
  const std::vector<double>& tau, const double x1, const std::vector<double>& t1, const unsigned n_threads,
  const std::vector<std::vector<int>>& eig_to_exclude)
{
  int n_x = x.size(), n_t1 = t1.size(), n_tau = tau.size();

  std::vector<std::vector<std::vector<std::vector<double>>>> corr(n_x, std::vector<std::vector<std::vector<double>>> (n_t1, std::vector<std::vector<double>>(n_tau, std::vector<double> (5, 0.0))));
  std::vector<int> zeros(N);
  int m = multiplicity(zeros);

  auto Vkl_precomputed = precompute_Vkl(x1, n_threads, eig_to_exclude);

  for(unsigned i_k = 0; i_k < eigenfunction_store.size(); i_k++){
    for(unsigned i_l = 0; i_l < eigenfunction_store.size(); i_l++){

      auto eig_k = eigenfunction_store[i_k], eig_l = eigenfunction_store[i_l];

      if( std::find(eig_to_exclude.begin(), eig_to_exclude.end(), eig_k) == eig_to_exclude.end() ||
        std::find(eig_to_exclude.begin(), eig_to_exclude.end(), eig_l) == eig_to_exclude.end() )
      {
        double lambda_k=lambda_fp(eig_k);
        double lambda_l=lambda_fp(eig_l);

        double Vkl=Vkl_precomputed[i_k][i_l];
        double V0l=V0k_element(eig_l, x1);

        for(int i_x = 0; i_x < n_x; i_x++){
          double V00=V0k_element(zeros, x[i_x]);
          double Vl0=Vk0_element(eig_l, x[i_x]);
          double Vk0=Vk0_element(eig_k, x[i_x]);

          double V0k=Vk0 * multiplicity(eig_k) / m; // cheaper than V0k_element(eig_k, x[i_x])

          for(int i_t1 = 0; i_t1 < n_t1; i_t1++){
            for(int i_tau = 0; i_tau < n_tau; i_tau++){

              double g3p  = V0k*Vl0*Vkl * std::exp(-lambda_k*tau[i_tau]-lambda_l*t1[i_t1]) / V00;
              double g2p1 = V0l*Vl0 * std::exp(-lambda_l*t1[i_t1]) / V00;
              double g2p2 = V0k*Vk0 * std::exp(-lambda_k*(tau[i_tau] + t1[i_t1])) / V00;

              corr[i_x][i_t1][i_tau][0] += g3p - g2p1 * g2p2;
              corr[i_x][i_t1][i_tau][1] += g3p;
              corr[i_x][i_t1][i_tau][2] += g2p1 * g2p2;
              corr[i_x][i_t1][i_tau][3] += g2p1;
              corr[i_x][i_t1][i_tau][4] += g2p2;
            }
          }
        }
      }
    }
  }

  return corr;
}

} // namespace SingleFile
