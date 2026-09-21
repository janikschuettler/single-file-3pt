/* Minimal example: three-point return correlations from the spectral solution.
 *
 * Computes  C(x,tau,y,t) = Q(x,t+tau; y,t | x) - Q(x,t+tau|x) Q(y,t|x)
 * for a single file of N particles on [0,1], for one return point x and one
 * intermediate point y, over a log-spaced grid of waiting times t and lag
 * times tau. Writes correlation.csv.
 *
 * Usage:  ./example [N] [tag] [x] [y] [nbins_t] [max_many_eig] [n_threads]
 * Default:        ./example 5 3 0.5 0.6 40 64 4
 */

#include "BetheSF/SingleFileDerived.hpp"
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>

static std::vector<double> log_spaced(int n, double lo, double hi)
{
  std::vector<double> v;
  if (n < 2) { v.push_back(lo); return v; }
  for (int i = 0; i < n; ++i)
    v.push_back(std::exp(std::log(lo) + i * (std::log(hi) - std::log(lo)) / (n - 1)));
  return v;
}

int main(int argc, char* argv[])
{
  const int      N       = argc > 1 ? std::atoi(argv[1]) : 5;
  const int      tag     = argc > 2 ? std::atoi(argv[2]) : 3;
  const double   x       = argc > 3 ? std::atof(argv[3]) : 0.5;
  const double   y       = argc > 4 ? std::atof(argv[4]) : 0.6;
  const int      nbins_t = argc > 5 ? std::atoi(argv[5]) : 40;
  const int      maxeig  = argc > 6 ? std::atoi(argv[6]) : 64;
  const unsigned threads = argc > 7 ? std::atoi(argv[7]) : 4;

  const double D = 1.0;

  // Timescales, for context. The spectral sum is truncated at max_many_eig,
  // which sets the shortest time that can be resolved.
  const double t_collision = 1.0 / D / N / N;
  const double t_relax     = 1.0 / D / std::pow(SingleFile::pi, 2);
  const double t_resolvable = 1.0 / D / std::pow(SingleFile::pi * maxeig, 2);

  const double t_min = 1e-3, t_max = 2.0;
  auto t1  = log_spaced(nbins_t, t_min, t_max);
  auto tau = log_spaced(nbins_t, t_min, t_max);
  std::vector<double> xs{x};

  SingleFile::SingleFileFlat sf(N, tag, maxeig, D);

  std::cout << "N = " << N << ", tag = " << tag << ", x = " << x << ", y = " << y << "\n"
            << "max_many_eig = " << maxeig << " -> " << sf.eigenfunctions().size()
            << " eigenfunctions, " << threads << " thread(s)\n"
            << "t_collision = " << t_collision << ", t_relax = " << t_relax
            << ", shortest resolvable t ~ " << t_resolvable << "\n"
            << "computing " << nbins_t << " x " << nbins_t << " grid ..." << std::endl;

  auto start = std::chrono::high_resolution_clock::now();
  auto c = sf.correlation_looped(xs, tau, y, t1, threads, {});
  auto stop  = std::chrono::high_resolution_clock::now();
  std::cout << "done in "
            << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count() / 1000.0
            << " s" << std::endl;

  std::ofstream out("correlation.csv");
  out << "# N=" << N << " tag=" << tag << " x=" << x << " y=" << y
      << " max_many_eig=" << maxeig << "\n";
  out << "t,tau,C,connected,disconnected\n";
  for (int i = 0; i < nbins_t; ++i)
    for (int j = 0; j < nbins_t; ++j)
      out << t1[i] << "," << tau[j] << "," << c[0][i][j][0] << ","
          << c[0][i][j][1] << "," << c[0][i][j][2] << "\n";
  out.close();

  std::cout << "wrote correlation.csv" << std::endl;
  return 0;
}
