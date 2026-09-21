/* Copyright (C) 2020 Alessio Lapolla and Aljaz Godec - All Rights Reserved
 * You may use, distribute and modify this code under the
 * terms of the MIT license. 
 *
 * You should have received a copy of the MIT license with
 * this file. If not, please write to: alessio.lapolla@mpibpc.mpg.de
 */
#ifndef SINGLEFILE_HPP
#define SINGLEFILE_HPP
#include <vector>
#include <array>
#include <set>
#include <string>
#include <utility>
#include <cmath>
#include <algorithm>
#include <iostream>

namespace SingleFile{

constexpr double pi=3.14159265359;

double factorial(const int n);//helper function

class NotImplementedException : public std::logic_error
{
public:
  NotImplementedException () : std::logic_error{"Function not yet implemented."} {}
};

class NotAllowedParameters : public std::logic_error
{
public:
  NotAllowedParameters () : std::logic_error{"The parameters chosen are not allowed."} {}
};


class IntegerPartitions;
template<typename T>
class UCombinations;

class SingleFile
{
  template<typename T>
  friend class UCombinations;//see file combinations.hpp
  friend class IntegerPartitions;//see file integerpartitions.hpp 
  friend std::ostream& operator<<(std::ostream& os, const SingleFile& sf);
  //prints members of the class
    
public:
  virtual double eq_prob(const double x) const;
  //calculates the equilibrium probability density function P(x)
  double joint2dens(const double  x, const double t, const double x0); 
  /*calculates the two joint density P(x,t;x_0,0)
    cannot be const for possible usage of permutations of the eigennumbers*/
  double green_function(const double x, const double t, const double x0)
    { return joint2dens(x,t,x0)/eq_prob(x0);};
  /*calculates the Green's function G(x,t|x_0,0)
    cannot be const for possible usage of permutations of the eigennumbers */
//=============================================================================
// Extension by Janik Schuettler (2021), MIT licensed -- three-point return
// correlations. These are members of SingleFile, so they must be declared here
// inside the class, but every one of them is DEFINED in
// SingleFileCorrelations.cpp, which contains nothing else. Nothing below this
// marker, down to the matching end marker, is part of upstream BetheSF.
//=============================================================================
  std::vector<double> joint2dens(const double x, const std::vector<double> t, const double x0);
  double joint3dens(const double x2, const double t2, const double x1, const double t1, const double x0);

  double correlation(const double x2, const double t2, const double x1, const double t1, const double x0);
  std::vector<std::vector<double>> correlation(const double x2, const std::vector<double> tau, const double x1, 
              const std::vector<double> t1, const double x0);
  std::vector<std::vector<std::vector<std::vector<double>>>> correlation(const std::vector<double> x2, const std::vector<double> tau, const double x1, 
              const std::vector<double> t1, const std::vector<double> x0);
  std::vector<std::vector<double>> precompute_Vkl(const double x1, const unsigned n_threads, const std::vector<std::vector<int>>& eig_to_exclude);
  std::vector<std::vector<std::vector<std::vector<double>>>> correlation_looped(const std::vector<double>& x, const std::vector<double>& tau, 
              const double x1, const std::vector<double>& t1, const unsigned n_threads, const std::vector<std::vector<int>>& eig_to_exclude);
//===================== end of Janik Schuettler's extension ===================

  //Naive functions, they do not work for asymmetric matrix elements, they can be modified to allow it
  double joint2dens_naive(const double  x, const double t, const double x0);
  //cannot be const for possible usage of permutations of the eigennumbers
  double green_function_naive(const double x, const double t, const double x0)
    { return joint2dens_naive(x,t,x0)/eq_prob(x0);};
//##############################################################################################################
  /*movable to private; here only for benchmark purpouses.
    They compute Eq.16 using the naive or the fast implementation respectively
  */
  double Vkl_element_naive(std::vector<int>& k_vec, std::vector<int>& l_vec, const double x) const;
  //Slavish implementation of Eq.16
  double Vkl_element(std::vector<int>& k_vec, std::vector<int>& l_vec, const double x) const;
  /* input: the multisets k and l and the position of the tagged particle x.
     It chooses the multiset with less permutations. Then it calls permuting_short_string.
     It returns the result of Eq.16       
   */
//##############################################################################################################
  SingleFile(const int N_, const int tag_, const int max_many_eig_, const double D_):
  N(N_), tag(tag_), max_many_eig(max_many_eig_), nl(tag_-1), nr(N_- tag_), min(std::min(nl,nr)), D(D_){
    if(tag<1 || tag>N) throw NotAllowedParameters();
  }
  SingleFile(const int N_, const int tag_):
  N(N_), tag(tag_), max_many_eig(0),  nl(tag_-1), nr(N_- tag_), min(std::min(nl,nr)), D(1){
    if(tag<1 || tag>N) throw NotAllowedParameters();
  }
  SingleFile(const SingleFile& other)=default;
  SingleFile& operator=(const SingleFile& other)=default;
  SingleFile(SingleFile&& other) noexcept = default;
  SingleFile& operator=(SingleFile&& other) noexcept = default ;      
  virtual ~SingleFile()=default;
  //utilities
  void get_eigenfunctions(const SingleFile& sf);
  //copies the eigenfuntion store from another SingleFile
  void get_eigenfunctions(const std::vector<std::vector<int>>& vec);
  //copies the eigenfunction store from the input
  void change_tag(const int i);//changes the tagged particle to i
  void change_DiffCoeff(const double d){
    if(d>0) D=d; else throw NotAllowedParameters();
  }
  void change_MaximumEigenvalue(const int max_many_eig);
  //changes member max_many_eig,resets eigenfunction_store as well,
  unsigned print_numbereigenfunctions(){
    return eigenfunction_store.size();
  }

  // --- added by Janik Schuettler (2021): accessors used by the extension ---
  std::vector<std::vector<int>> eigenfunctions() {
    return eigenfunction_store;
  }

  void print_eigenfunctions() {
    std::cout << "number of eigenfunctions: " << eigenfunction_store.size() << std::endl;
  }
  // --- end ---
    
private:
  std::vector<std::vector<int>> eigenfunction_store;//container of all allowed multisets

  
  //naive implementation
  double left_integral(const std::vector<int>& left_k, const std::vector<int>& left_l, const double x) const;
  /*calculates \prod_{j=1}^{i-1} L_j(x) in Eq.16, input:
    the vector of eigennumbers to the left of the tagged particle for the multiset k (left_k) and l (left_l)
    the position of the tagged particle x*/
  double right_integral(const std::vector<int>& right_k, const std::vector<int>& right_l, const double x) const;
  /*calculates \prod_{j=i+1}^N R_j(x) in Eq.16
    the vector of eigennumbers to the right of the tagged particle for the multiset k (right_k) and l (right_l)
    the position of the tagged particle x*/

  //fast code
  double permuting_short_string(std::vector<int>& short_string,
        std::vector<int>& long_string, const double x) const;
  /*This function is called by Vkl_element.
    Input: long_string is the multiset with more permutations, short_string the multiset with less,
    x is the tagged particle position.
    It combines the two multisets in a multiset of pairs.
    It finds all unique pairs.
    For each unique pair computes the multiset of remaining pairs.
    It uses the class UCombinations to compute the combinations of the remaining multiset of pairs.
    It repeats everything for each permutation of short_string.
    It returns Eq.16 without the prefactor \frac{m_\mathbf{l}}{N_L!N_R!}.
   */
  double combinations_contribution(const std::vector<std::pair<int,int>>& remaining_pairs,
           const std::vector<std::pair<int,int>>& short_combination,
           const std::pair<int, int>& unique_pair, const double x) const;
  /*This function is called by the class UCombinations.
    Input: the multiset of remaining_pairs, a single combination short_combination of the multiset of pairs to 
    the left or to the right of the tagged particle index (it picks the smallest multiset), the unique_pair is 
    the pair with the eigennumbers of the tagged particle, x is the position of the tagged particle.
    The function computes the combination completing the set by difference between remaining_pairs and short_combinations.
    It computes how many identical terms there are in Eq.16.
    Finally, it uses left_integral and right_integral to compute the products in Eq.16.
    It returns the total contribution of this product.
   */
  double left_integral(const std::vector<std::pair<int,int>>& left, const double x) const;
  /*calculates \prod_{j=1}^{i-1} L_j(z) in Eq.16, input:
    the vector of pairs of eigennumbers to the left of the tagged particle 
    the position of the tagged particle x*/
  double right_integral(const std::vector<std::pair<int,int>>& right, const double x) const;
  /*calculates \prod_{j=i+1}^N R_j(z) in Eq.16
    the vector of pairs of eigennumbers to the right of the tagged particle
    the position of the tagged particle x*/
  //##############################
  virtual int eigenfunction_condition(const int i) const=0;
  //function selecting only the allowed multisets
  virtual double lambda_single(const int n) const=0;
  //calculates the single particle eigenvalues \lambda
  double lambda_fp(const std::vector<int>& eig_num) const{
    double s=0.0; 
    for(const auto n : eig_num) 
      s+=lambda_single(n); 
    return s;
  }
  //calculates the many-body eigenvalues \Lambda
    

  /*lambda, lambda_k and lambda_l input parameters in the next 6 functions
     always label a single particle eigennumber, x is the position of the tagged particle.*/    
  /*The next 3 functions implement the contribution of the single particle eigenfunctions 
    integral tagged is Eqs 17a, lefttagged is Eq.17b and righttagged is Eq.17c.*/
  virtual double tagged(const int lambda_k, const int lambda_l, const double x) const=0;
  virtual double lefttagged(const int lambda_k, const int lambda_l, const double x) const=0;
  virtual double righttagged(const int lambda_k, const int lambda_l, const double x) const=0;
  //the next 3 functions assume that one single particle eigenvalues is 0 and are called by the previous 3
  virtual double tagged(const int lambda, const double x) const=0;
  virtual double lefttagged(const int lambda, const double x) const=0;
  virtual double righttagged(const int lambda, const double x) const=0;
  //Specialization of Vkl_element, Eq.13. They assume that one of the multisets is the ground state one
  virtual double V0k_element(std::vector<int>& k_vec, const double x) const;
  virtual double Vk0_element(std::vector<int>& k_vec, const double x) const;

protected:
  const int N;//size of the single file
  int tag;//tagged particle index
  int max_many_eig;//parameter \Lambda_M
  int nl, nr;//number of particle to the left and to the right of the tagged one
  int min;//min between nl and nr
  double D;//Diffusion coefficent
  void eigenfunction_store_init();
  //This function initializes eigenfunction_store. It is called by the constructor of the derived class
  // denominator Eq. A.1
  template<typename T>
  double multiplicity(const std::vector<T>& eigenfunction_numbers) const;
  //Eq A.1
  template<typename T>
  double number_permutations(const std::vector<T>& l_s_combin) const;
};

template<typename T> // denominator Eq. A.1
double SingleFile::multiplicity(const std::vector<T>& eigenfunction_numbers) const
{
  double tot_mult=1.0;
  std::set<T> unique(eigenfunction_numbers.cbegin(), eigenfunction_numbers.cend());
  for(const auto el : unique){
    int sing_mult=std::count(eigenfunction_numbers.cbegin(), eigenfunction_numbers.cend(), el);
    tot_mult*=factorial(sing_mult);
  }
  return tot_mult;
}
    
template<typename T> //Eq A.1
double SingleFile::number_permutations(const std::vector<T>& l_s_combin) const 
{
  double tot_mult=multiplicity(l_s_combin);
  return factorial(l_s_combin.size())/tot_mult;
}    

//prints the members of the class
std::ostream& operator<<(std::ostream& os, const SingleFile& sf);
    
}//end namespace SingleFile
#endif
 
