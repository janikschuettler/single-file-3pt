/* Copyright (C) 2020 Alessio Lapolla and Aljaz Godec - All Rights Reserved
 * You may use, distribute and modify this code under the
 * terms of the MIT license. 
 *
 * You should have received a copy of the MIT license with
 * this file. If not, please write to: alessio.lapolla@mpibpc.mpg.de
 */
#ifndef SINGLEFILEDERIVED_HPP
#define SINGLEFILEDERIVED_HPP
#include "SingleFile.hpp"

namespace SingleFile{
/* all these classes implement the pure virtual functions their base class. 
   See file SingleFile.hpp. */

class SingleFileFlat : public SingleFile
{
friend class IntegerPartitions;
public:
  SingleFileFlat(const int N_, const int tag_, const int max_many_eig_,  const double D_=1.0) :
  SingleFile(N_, tag_, max_many_eig_, D_){ eigenfunction_store_init(); }
  SingleFileFlat(const int N, const int tag):
  SingleFile(N,tag){};
  SingleFileFlat(const SingleFileFlat& other)=default;
  SingleFileFlat& operator=(const SingleFileFlat& other)=default;
  SingleFileFlat(SingleFileFlat&& other) noexcept = default;
  SingleFileFlat& operator=(SingleFileFlat&& other) noexcept = default;
  ~SingleFileFlat()=default;
  
  double eq_prob(const double x) const override;
private:
  int eigenfunction_condition(const int i) const override{
    return i*i;
  }//qudratic constraint
  double lambda_single(const int n) const override {
    return D*std::pow(pi*n,2);
  }
  double tagged(const int lambda, const double x) const override;
  double lefttagged(const int lambda, const double x) const override;
  double righttagged(const int lambda, const double x) const override;
  double tagged(const int lambda_k, const int lambda_l, const double x) const override;
  double lefttagged(const int lambda_k, const int lambda_l, const double x) const override;   
  double righttagged(const int lambda_k, const int lambda_l, const double x) const override;
};

class SingleFileOnSlope : public SingleFile
{
  friend class IntegerPartitions;
  template<typename T>
  friend class UCombinations;
public:
  SingleFileOnSlope(const int N_, const int tag_, const int max_many_eig_, const double D_=1.0, const double g_=1.0):
  SingleFile(N_, tag_, max_many_eig_, D_), g(g_){ eigenfunction_store_init(); }
  SingleFileOnSlope(const int N, const int tag) :
  SingleFile(N, tag), g(1) {}
  SingleFileOnSlope(const SingleFileOnSlope& other)=default;
  SingleFileOnSlope& operator=(const SingleFileOnSlope& other)=default;
  SingleFileOnSlope(SingleFileOnSlope&& other) noexcept = default;
  SingleFileOnSlope& operator=(SingleFileOnSlope&& other) noexcept = default ;        
  ~SingleFileOnSlope()=default;
  
  double eq_prob(const double x) const override;

    
private:
  double g;
  int eigenfunction_condition(const int i) const override{
  return i*i;}//quadratic constraint
  double lambda_single(const int n) const override {
  return (n==0) ? 0 : D*std::pow(pi*n,2)+0.25*g*g/D;}
  double tagged(const int lambda, const double x) const override;
  double lefttagged(const int lambda, const double x) const override;
  double righttagged(const int lambda, const double x) const override;
  double tagged(const int lambda_k, const int lambda_l, const double x) const override;
  double lefttagged(const int lambda_k, const int lambda_l, const double x) const override;   
  double righttagged(const int lambda_k, const int lambda_l, const double x) const override;
  //extra functions to take in account of the asymmetry between V_k0 and V_0k
  double tagged_k0(const int lambda, const double x) const;             
  double lefttagged_k0(const int lambda, const double x) const;     
  double righttagged_k0(const int lambda, const double x) const;
  double left_integral_k0(const std::vector<int>& left_k, const double x) const;
  double right_integral_k0(const std::vector<int>& right_k, const double x) const;
  double combinations_contributionk0(const std::vector<int>& remaining, const std::vector<int>& combination, const int unique, const double x) const;
public:
  double Vk0_element(std::vector<int>& k_vec, const double x) const override;  
};

// class SingleFileHarmonic: public SingleFile
// {
//   friend class IntegerPartitions;
// public:
//   SingleFileHarmonic(const int N_, const int tag_, const int max_many_eig_,  const double D_=1.0, const double gamma=1.0):
//   SingleFile(N_, tag_, max_many_eig_, D_), gamma(gamma) { eigenfunction_store_init(); }
//   SingleFileHarmonic(const int N, const int tag) :
//   SingleFile(N, tag), gamma(1) {}
//   SingleFileHarmonic(const SingleFileHarmonic& other)=default;
//   SingleFileHarmonic& operator=(const SingleFileHarmonic& other)=default;
//   SingleFileHarmonic(SingleFileHarmonic&& other) noexcept = default;
//   SingleFileHarmonic& operator=(SingleFileHarmonic&& other) noexcept = default ;      
//   ~SingleFileHarmonic()=default;

// private:
//   double gamma;
//   int eigenfunction_condition(const int i) const override{
//   return i;}//linear constraint
//   double lambda_single(const int n) const override {
//   return gamma*n;}
//   double tagged(const int lambda, const double x) const override;
//   double lefttagged(const int lambda, const double x) const override;
//   double righttagged(const int lambda, const double x) const override;
//   double tagged(const int lambda_k, const int lambda_l, const double x) const override;
//   double lefttagged(const int lambda_k, const int lambda_l, const double x) const override;   
//   double righttagged(const int lambda_k, const int lambda_l, const double x) const override;
// };



}//end namespace SingleFile

#endif
