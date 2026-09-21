/* Copyright (C) 2020 Alessio Lapolla and Aljaz Godec - All Rights Reserved
 * You may use, distribute and modify this code under the
 * terms of the MIT license. 
 *
 * You should have received a copy of the MIT license with
 * this file. If not, please write to: alessio.lapolla@mpibpc.mpg.de
 */
#include <cmath>
#include <vector>
#include <set>
#include <utility>
#include <string>
#include <algorithm>
#include <iostream>
#include "combinations.hpp"
#include "SingleFileDerived.hpp"


namespace SingleFile{
//##################################################
// Single File Flat
//##################################################

double SingleFileFlat::eq_prob(const double x) const
{
  return factorial(N)/(factorial(nl)*factorial(nr))*std::pow(x,nl)*std::pow(1.0-x,nr);
}

//Eqs 17a, 17b and 17c
inline double SingleFileFlat::tagged(const int lambda, const double x) const
{
  return (lambda==0) ?   1.0 : std::sqrt(2)*std::cos(lambda*pi*x);
}


inline  double SingleFileFlat::lefttagged(const int lambda, const double x) const
{
  return (lambda==0) ? x : std::sqrt(2)*std::sin(lambda*pi*x)/(lambda*pi);
}

inline double SingleFileFlat::righttagged(const int lambda, const double x) const
{
  return (lambda==0) ? 1.0-x : -std::sqrt(2)*std::sin(lambda*pi*x)/(lambda*pi);
}

//Eqs 17a, 17b and 17c for the general case
double SingleFileFlat::tagged(const int lambda_k, const int lambda_l, const double x) const
{
  if(lambda_k==0)
    return tagged(lambda_l, x);//takes care of the lambda_k==0 && lambda_l==0 as well
  else if(lambda_l==0)
    return tagged(lambda_k, x);
  else
    return 2*std::cos(lambda_l*pi*x)*std::cos(lambda_k*pi*x);
}

double SingleFileFlat::lefttagged(const int lambda_k, const int lambda_l,const double x) const
{
  if(lambda_k==0)
    return lefttagged(lambda_l, x);//takes care of the lambda_k==0 && lambda_l==0 as well
  else if(lambda_l==0)
    return lefttagged(lambda_k, x);
  else if(lambda_k==lambda_l)
    return x+std::sin(2*lambda_k*pi*x)/(2*lambda_k*pi);
  else
    return 2*(lambda_k*std::cos(lambda_l*pi*x)*std::sin(lambda_k*pi*x)-lambda_l*std::cos(lambda_k*pi*x)*std::sin(lambda_l*pi*x))/(pi*(lambda_k*lambda_k-lambda_l*lambda_l));
}

double SingleFileFlat::righttagged(const int lambda_k, const int lambda_l,const double x) const
{
  if(lambda_k==0)
    return righttagged(lambda_l, x);//takes care of the lambda_k==0 && lambda_l==0 as well
  else if(lambda_l==0)
    return righttagged(lambda_k, x);
  else if(lambda_k==lambda_l)
    return 1-x-std::sin(2*lambda_k*pi*x)/(2*lambda_k*pi);
  else
    return 2*(-lambda_k*std::cos(lambda_l*pi*x)*std::sin(lambda_k*pi*x)+lambda_l*std::cos(lambda_k*pi*x)*std::sin(lambda_l*pi*x))/(pi*(lambda_k*lambda_k-lambda_l*lambda_l));
}

//####################################################################################################
// Single file on slope
// For this system V_{k0} and V_{0k} are different
// The implementations of Eqs.13 must be different, here they follow
//####################################################################################################

double SingleFileOnSlope::eq_prob(const double x) const 
{
  double li=(1-std::exp(-g*x/D))/(1-std::exp(-g/D));
  double ri=(std::exp(-g*x/D)-std::exp(-g/D))/(1-std::exp(-g/D));
  double t=g*std::exp(-g*x/D)/(D*(1-std::exp(-g/D)));
  return factorial(N)*t*std::pow(li, nl)*std::pow(ri, nr)/(factorial(nl)*factorial(nr));
}

//Eqs 17a, 17b and 17c for V_k0
double SingleFileOnSlope::tagged(const int lambda, const double x) const
{
  if(lambda==0)
    return g/D*std::exp(-g*x/D)/(1.0-std::exp(-g/D));
  else
    return std::exp(-0.5*g*x/D)*(std::sin(lambda*pi*x)-2*D/g*lambda*pi*std::cos(lambda*pi*x))/std::sqrt(0.5+2*std::pow(D*lambda*pi/g,2));
}

double SingleFileOnSlope::lefttagged(const int lambda, const double x) const
{
  if(lambda==0)
    return (1-std::exp(-g*x/D))/(1-std::exp(-g/D));
  else 
    return -2*D/g*std::exp(-0.5*g*x/D)*std::sin(lambda*pi*x)/std::sqrt(0.5+2*D*std::pow(lambda*pi/g,2));
}

double SingleFileOnSlope::righttagged(const int lambda, const double x) const
{
  if(lambda==0)
    return (std::exp(-g*x/D)-std::exp(-g/D))/(1-std::exp(-g/D));
  else 
    return 2*D/g*std::exp(-0.5*g*x/D)*std::sin(lambda*pi*x)/std::sqrt(0.5+2*D*std::pow(lambda*pi/g,2));
}

//Eqs 17a, 17b and 17c for the general case for V_k0
double SingleFileOnSlope::tagged(const int lambda_k, const int lambda_l, const double x) const
{
  if(lambda_k==0)
    return tagged(lambda_l, x);//takes care of the lambda_k==0 && lambda_l==0 as well
  else if(lambda_l==0)
    return tagged(lambda_k, x);
  else
    return tagged(lambda_k,x)*tagged(lambda_l,x);
}

double SingleFileOnSlope::lefttagged(const int lambda_k, const int lambda_l,const double x) const
{
  if(lambda_k==0)
    return lefttagged(lambda_l, x);//takes care of the lambda_k==0 && lambda_l==0 as well
  else if(lambda_l==0)
    return lefttagged(lambda_k, x);
  else if(lambda_k==lambda_l)
    return 2*(-D*g+0.5*g*g*x+2*std::pow(D*lambda_k*pi,2)*x+D*g*std::cos(2*lambda_k*pi*x)+(D*D*lambda_k*pi-g*g/(4*lambda_k*pi))*std::sin(2*lambda_k*pi*x))/(g*g+std::pow(2*D*lambda_k*pi,2));
  else
    return 2*((g*g+std::pow(2*D*lambda_k*pi,2))*(lambda_l*std::cos(lambda_l*pi*x)*std::sin(lambda_k*pi*x)-lambda_k*std::cos(lambda_k*pi*x)*std::sin(lambda_l*pi*x))-2*D*g*(lambda_k*lambda_k-lambda_l*lambda_l)*pi*std::sin(lambda_k*pi*x)*std::sin(lambda_l*pi*x))/(g*g*pi*(lambda_k*lambda_k-lambda_l*lambda_l)*(1+std::sqrt(1+std::pow(2*D*lambda_k*pi/g,2)))*(1+std::sqrt(1+std::pow(2*D*lambda_l*pi/g,2))));
}

double SingleFileOnSlope::righttagged(const int lambda_k, const int lambda_l,const double x) const
{
  if(lambda_k==0)
    return righttagged(lambda_l, x);//takes care of the lambda_k==0 && lambda_l==0 as well
  else if(lambda_l==0)
    return righttagged(lambda_k, x);
  else if(lambda_k==lambda_l)
    return (2*lambda_k*pi*(g*g+std::pow(2*D*lambda_k*pi,2))*(1-x)+4*D*g*lambda_k*pi*(1-std::cos(2*lambda_k*pi*x))+(g*g-std::pow(2*D*lambda_k*pi,2))*std::sin(2*lambda_k*pi*x))/(2*lambda_k*pi*(g*g+std::pow(2*D*lambda_k*pi,2)));
  else
    return 2*(g*g*lambda_k*std::cos(lambda_k*pi*x)*std::sin(lambda_l*pi*x)-(g*g+std::pow(2*D*lambda_k*pi,2))*std::cos(lambda_l*pi*x)*std::sin(lambda_k*pi*x)+std::pow(2*D*lambda_l*pi,2)*lambda_k*std::cos(lambda_k*pi*x)*std::sin(lambda_l*pi*x)+2*D*g*(lambda_k*lambda_k-lambda_l*lambda_l)*pi*std::sin(lambda_k*pi*x)*std::sin(lambda_l*pi*x))/(g*g*pi*(lambda_k*lambda_k-lambda_l*lambda_l)*(1+std::sqrt(1+std::pow(2*D*lambda_k*pi/g,2)))*(1+std::sqrt(1+std::pow(2*D*lambda_l*pi/g,2))));
}

//Eqs 17a, 17b and 17c for V_0k
inline double SingleFileOnSlope::tagged_k0(const int n, const double x) const
{
  if(n==0)
    return g/D*std::exp(-g*x/D)/(1.0-std::exp(-g/D));
  else
    return g/D*std::exp(-0.5*g*x/D)*(std::sin(n*pi*x)-2*D/g*n*pi*std::cos(n*pi*x))/(std::sqrt(0.5+2*std::pow(n*pi/g,2))*(1-std::exp(-g/D)));
}

inline double SingleFileOnSlope::lefttagged_k0(const int n, const double x) const
{
  if(n==0)
    return (1-std::exp(-g*x/D))/(1-std::exp(-g/D));
  else
    return -2*std::exp(-0.5*g*x/D)*std::sin(n*pi*x)/(std::sqrt(0.5+2*D*std::pow(n*pi/g,2))*(1-std::exp(-g/D)));
}

inline double SingleFileOnSlope::righttagged_k0(const int n, const double x) const
{
  if(n==0)
    return (std::exp(-g*x/D)-std::exp(-g/D))/(1-std::exp(-g/D));
  else
    return 2*std::exp(-0.5*g*x/D)*std::sin(n*pi*x)/(std::sqrt(0.5+2*D*std::pow(n*pi/g,2))*(1-std::exp(-g/D)));
}

//products in EQ.16 for V_k0
double SingleFileOnSlope::left_integral_k0(const std::vector<int>& left_k, const double x) const 
{
  if(left_k.empty())
    return 1.0;
  else{
    double left_contr=1.0;
    for(const auto k : left_k)
      left_contr*=lefttagged_k0( k, x);
    return left_contr;
  }
}

double SingleFileOnSlope::right_integral_k0(const std::vector<int>& right_k, const double x) const 
{
  if(right_k.empty())
    return 1.0;
  else{
    double right_contr=1.0;
    for(const auto k : right_k)
      right_contr*=righttagged_k0(k, x);
    return right_contr;
  }
}

//called by the friend class UCombinations see its file!!!
double SingleFileOnSlope::combinations_contributionk0(const std::vector<int>& remaining, const std::vector<int>& combination, const int unique, const double x) const
{
  std::vector<int> l_comb;
  auto s_comb=combination;
  std::sort(s_comb.begin(), s_comb.end());
  std::set_difference(remaining.cbegin(), remaining.cend(),
    s_comb.cbegin(), s_comb.cend(),
    std::back_inserter(l_comb));
  double n_perm_short=number_permutations(s_comb);
  double n_perm_long=number_permutations(s_comb);
  double t=tagged_k0(unique, x);
  double l, r;
  if(nl<=nr){
    l=left_integral_k0(s_comb, x);
    r=right_integral_k0(l_comb, x);
  }
  else{
    l=left_integral_k0(l_comb, x);
    r=right_integral_k0(s_comb, x);
  }
      
  return n_perm_short*n_perm_long*t*l*r;
}

// //called by the friend class UCombinations see its file!!!
// double SingleFileOnSlope::combinations_contributionk0(const std::vector<int>& remaining, const std::vector<int>& combination, const int unique, const double x) const
// {
//   std::vector<int> l_comb;
//   std::set_difference(remaining.cbegin(), remaining.cend(),
//     combination.cbegin(), combination.cend(),
//     std::inserter(l_comb , l_comb.begin()));
//   double n_perm_short=number_permutations(combination);
//   double n_perm_long=number_permutations(combination);
//   double t=tagged_k0(unique, x);
//   double l, r;
//   if(nl<=nr){
//     l=left_integral_k0(combination, x);
//     r=right_integral_k0(l_comb, x);
//   }
//   else{
//     l=left_integral_k0(l_comb, x);
//     r=right_integral_k0(combination, x);
//   }
      
//   return n_perm_short*n_perm_long*t*l*r;
// }

double SingleFileOnSlope::Vk0_element(std::vector<int>& k_vec, const double x) const
{
  double intres=0.0;
  std::set<int> unique_k_vec_numbers(k_vec.cbegin(), k_vec.cend());
  for(const auto unique : unique_k_vec_numbers){
    std::array<int,1> chosen={unique};
    std::vector<int> all_but_chosen_eigennumbers;
    std::set_difference(k_vec.cbegin(), k_vec.cend(),
        chosen.cbegin(), chosen.cend(),
        std::inserter(all_but_chosen_eigennumbers ,all_but_chosen_eigennumbers.begin()));
    if(tag == 1 || tag == N){
      double n_perm_long=number_permutations(all_but_chosen_eigennumbers);
      double t=tagged_k0(unique, x);
      double rem;
      if(tag == 1)
        rem=right_integral_k0(all_but_chosen_eigennumbers, x);
      else
        rem=left_integral_k0(all_but_chosen_eigennumbers, x);
      intres+=n_perm_long*t*rem;
    }
    else{
      UCombinations<int> com_gen;
      intres+=com_gen.unique_combinations(all_but_chosen_eigennumbers, min, unique, x, *this);
    }//end else clause
  }//end for over unique
  return factorial(N)*intres/(factorial(nl)*factorial(nr));
}

//####################################################################################################
// Single File Harmonic
//####################################################################################################

// //Eqs 17a, 17b and 17c
// double SingleFileHarmonic::tagged(const int lambda, const double x) const
// {
//   if(lambda==0)
//     return std::sqrt(gamma/(2*pi*D))*std::exp(-0.5*gamma*x*x/D);
//   else
//     return tagged(0,x)*std::hermite(lambda,x*std::sqrt(0.5*gamma/D))/(std::pow(2,lambda)*factorial(lambda));
// }

// double SingleFileHarmonic::righttagged(const int lambda, const double y) const
// {
//   const double coeff=1.0/std::sqrt(pi*std::pow(2,lambda)*factorial(lambda));
//   const double x=y*std::sqrt(0.5*gamma/D);
//   if(lambda==0)
//     return coeff*0.5*std::sqrt(pi)*std::erfc(x);
//   else if(lambda%2==0)
//     return coeff*std::exp(-x*x)*std::hermite(lambda-1,x);      
//   else{
//     auto a=std::hermite(lambda-1,0)-std::exp(-x*x)*std::hermite(lambda-1,x);
//     const int m=(lambda-1)/2;
//     auto b=std::pow(-1,m)*std::pow(2,2*m+1)*std::tgamma(m+1.5)/(std::sqrt(pi)*(1+2*m));//*pochhammer(0.5,n)*pochhammer(1.5,n);
//     return coeff*(b-a); 
//   }
// }
    
// double SingleFileHarmonic::lefttagged(const int lambda, const double y) const
// {
//   const double coeff=1.0/std::sqrt(pi*std::pow(2,lambda)*factorial(lambda));
//   if(lambda==0)
//     return coeff*0.5*std::sqrt(pi)*(1+std::erf(y*std::sqrt(0.5*gamma/D)));
//   else if(lambda%2==0)
//     return righttagged(lambda,-y);
//   else
//     return -righttagged(lambda,-y);
// }

// double SingleFileHarmonic::tagged(const int lambda_k, const int lambda_l, const double x) const
// {
//   if(lambda_k==0)
//     return tagged(lambda_l, x);//takes care of the lambda_k==0 && lambda_l==0 as well
//   else if(lambda_l==0)
//     return tagged(lambda_k, x);
//   else
//     throw NotImplementedException();
// }

// double SingleFileHarmonic::lefttagged(const int lambda_k, const int lambda_l, const double x) const
// {
//   if(lambda_k==0)
//     return lefttagged(lambda_l, x);//takes care of the lambda_k==0 && lambda_l==0 as well
//   else if(lambda_l==0)
//     return lefttagged(lambda_k, x);
//   else if(lambda_k==lambda_l)
//     throw NotImplementedException();
//   else
//     throw NotImplementedException();
// }

// double SingleFileHarmonic::righttagged(const int lambda_k, const int lambda_l, const double x) const
// {
//   if(lambda_k==0)
//     return righttagged(lambda_l, x);//takes care of the lambda_k==0 && lambda_l==0 as well
//   else if(lambda_l==0)
//     return righttagged(lambda_k, x);
//   else if(lambda_k==lambda_l)
//     throw NotImplementedException();
//   else
//     throw NotImplementedException();
// }

}//end namespace SingleFIle
