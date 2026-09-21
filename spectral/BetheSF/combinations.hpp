/* Copyright (C) 2020 Alessio Lapolla and Aljaz Godec - All Rights Reserved
 * You may use, distribute and modify this code under the
 * terms of the MIT license. 
 *
 * You should have received a copy of the MIT license with
 * this file. If not, please write to: alessio.lapolla@mpibpc.mpg.de
 */
#ifndef COMBINATIONS_HPP
#define COMBINATIONS_HPP
#include <vector>
#include <iostream>
#include <set>
#include "SingleFile.hpp"
#include "SingleFileDerived.hpp"
/*this code is built and adapted for the SingleFile class,
It computes all the unique combinations of a multiset.
Let's suppose our vector is {0, 0, 0, 0, 0, 0, 0, 0, 0, 1} and we want the 2-combinations, 
the result must be
{ {0,0}, {0,1} } and not
{{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 
  1}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 
  1}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 1}, {0, 
  0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 1}, {0, 0}, {0, 0}, {0, 
  0}, {0, 0}, {0, 1}, {0, 0}, {0, 0}, {0, 0}, {0, 1}, {0, 0}, {0, 
  0}, {0, 1}, {0, 0}, {0, 1}, {0, 1}}
*/
namespace SingleFile{

template<typename T>
class UCombinations{
public:
  UCombinations()=default;
  double intres;
  double unique_combinations(const std::vector<T>& remaining, const int k, const T unique, 
                             const double x, const SingleFile& sf); //calls the homonym
  /*
    Input
    unique : the eigenumber(s) of the tagged particle
    remaining: the multiset of (pairs) of eigennumbers of the remaining particles
    The vector remaining must be sorted from the smallest to the larger element, duplicates are welcome!
    k: the size of the k-combination (it is the min between nl and nr)
    x: the point we are interested in
    sf: the single file we are interested in.
   */
 
  void unique_combinations(const std::vector<T>& remaining,  const int k, const T unique, const double x,
                           std::vector<T>& combination, const SingleFile& sf, const int index=0);
  /* This function computes the k-combination (the vector combination) of the vector remaining.
     Then the vector combination, remaining, unique, x and sf are passed to the function
     SingleFile::combinations_contribution (see file SingleFile.hpp)
     index is a parameter to recursively iterate through all the combinations.
   */
  
  /* The following two function are needed for the specialization of SingleFileOnSlope::Vk0_element (see file SingleFileDerived.hpp.
   They are equivalent  to the previous two, the only difference is that 
   they call SingleFileOnSlope::combinations_contributionk0.*/
  double unique_combinations(const std::vector<T>& remaining, const int k, const T unique, const double x, 
                             const SingleFileOnSlope& sf);
  void unique_combinations(const std::vector<T>& remaining,  const int k, const T unique, const double x,
                           std::vector<T>& combination, const SingleFileOnSlope& sf, const int index=0);
};

template<typename T>
void UCombinations<T>::unique_combinations(const std::vector<T>& remaining,  const int k, const T unique,
                                           const double x, std::vector<T>& combination, const SingleFile& sf, 
                                           const int index)
{

  int n=remaining.size();
  if(k==n){ //takes care of the case in which there is only one k-combination
    intres+=sf.combinations_contribution(remaining, combination, unique, x);
    return;
  }
  if(k==0){ //leaf of the tree reached, new combination reached
    intres+=sf.combinations_contribution(remaining, combination, unique, x);
    return;
  }

  for(int j=index; j<n; ++j){ //generates the tree going through all the combinations
    combination.at(k-1)=remaining[j];
    unique_combinations(remaining, k-1, unique, x, combination, sf, j+1);
    while(j<n-1 && remaining[j]==remaining[j+1]) //handles duplicates skipping them
      j++;
  }
}

template<typename T>
double UCombinations<T>::unique_combinations(const std::vector<T>& remaining, const int k, const T unique, 
                                             const double x, const SingleFile& sf)
{
  intres=0.0;
  std::vector<T> combination(k);
  unique_combinations(remaining, k, unique, x, combination, sf);
  return intres;
} 

//Specialized version for the SingleFileOnSlope class.
template<typename T>
void UCombinations<T>::unique_combinations(const std::vector<T>& remaining,  const int k, const T unique, 
                                           const double x, std::vector<T>& combination, const SingleFileOnSlope& sf, 
                                           const int index)
{

  int n=remaining.size();
  if(k==n){
    intres+=sf.combinations_contributionk0(remaining, combination, unique, x);
    return;
  }
  if(k==0){
    intres+=sf.combinations_contributionk0(remaining, combination, unique, x);
    return;
  }

  for(int j=index; j<n; ++j){
    combination.at(k-1)=remaining[j];
    unique_combinations(remaining, k-1, unique, x, combination, sf, j+1);
    while(j<n-1 && remaining[j]==remaining[j+1])
      j++;
  }
}

template<typename T>
double UCombinations<T>::unique_combinations(const std::vector<T>& remaining, const int k, const T unique, 
                                             const double x, const SingleFileOnSlope& sf)
{
  intres=0.0;
  std::vector<T> combination(k);
  unique_combinations(remaining, k, unique, x, combination, sf);
  return intres;
} 

}//end namespace SingleFile
#endif
