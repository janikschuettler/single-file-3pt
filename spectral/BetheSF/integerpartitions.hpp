/* Copyright (C) 2020 Alessio Lapolla and Aljaz Godec - All Rights Reserved
 * You may use, distribute and modify this code under the
 * terms of the MIT license. 
 *
 * You should have received a copy of the MIT license with
 * this file. If not, please write to: alessio.lapolla@mpibpc.mpg.de
 */
#ifndef INTEGERPARTITIONS_HPP
#define INTEGERPARTITIONS_HPP
#include <vector>
#include <iostream>
#include "SingleFile.hpp"
//this code is built and adapted for the SingleFile class
/*this class generates all the integer partitions of a number n in m
  parts, a part can be 0. Code adapted from
  the art of computer science vol4 pag 392*/
namespace SingleFile{
    
class IntegerPartitions{
public:
  static void printer(const std::vector<int>& parts, const SingleFile& sf, std::vector<std::vector<int>>& result);
  static std::vector<std::vector<int>> integer_partitions( const int n, const int m, const SingleFile& sf);
        //n= number to partition, m=size of the partition
};

void IntegerPartitions::printer(const std::vector<int>& parts, const SingleFile& sf, std::vector<std::vector<int>>& result)
{
  auto sum_condition=[&](){
    int s=0;
    for(auto i=0; i<int(parts.size())-1; ++i){
       s+=sf.eigenfunction_condition(parts[i]);}
    return s;
  };
    
  if(sum_condition()<=static_cast<int>(sf.max_many_eig)){
    std::vector<int> res(parts.crbegin()+1, parts.crend());//parts in this way is sorted
    result.push_back(res);
  }
}

std::vector<std::vector<int>> IntegerPartitions::integer_partitions( const int n, const int m, const SingleFile& sf)
{ //generates the values in reverse sorted fashion
  std::vector<std::vector<int>> result;
  if(n==0){
    std::vector<int> parts(m);
    result.push_back(parts);
  }
  else{
    std::vector<int> parts(m+1);
    parts[0]=n;
    parts[m]=-1;//modification to allow 0 in the partition
    while(true)
    {
      printer(parts, sf, result);
      while(parts[1]<parts[0]-1)//tweak h3
      {
        parts[0]-=1;
        parts[1]+=1;
        printer(parts, sf, result);
      }
      int j=2;//step h4
      int s=parts[0]+parts[1]-1;
      while(parts[j]>=parts[0]-1)
      {
        //cout<<"h4"<<endl;
        s+=parts[j];
        ++j;
      }
      if(j>m-1)//exit statement
        break;
      int x=parts[j]+1;//step h5
      parts[j]=x;
      --j;
      while(j>0)//tweak h6
      {
        parts[j]=x;
        s-=x;
        --j;
        parts[0]=s;
      }
    }
  }
  return result;
}

}//end namespace SingleFile
#endif
