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
#include <map>
#include <utility>
#include <string>
#include <algorithm>
#include <iostream>
#include "combinations.hpp"
#include "integerpartitions.hpp"
#include "SingleFile.hpp"
#include <stdexcept>


#include <chrono>


using std::cout; using std::endl;
namespace SingleFile{

    //utility function
// double factorial(const int n)// IT IS NOT THREAD SAFE!!!
// {
//   static std::map<int, double> fact_store={{0,1.0}};
//   if(fact_store.find(n)!=fact_store.end()){
//     return fact_store[n];
//   }
//   else{
//     double p=std::tgamma(n+1);
//     fact_store[n]=p;
//     return p;
//   }
// }

double factorial(const int n)
{
    constexpr std::array<double,171> fact={1.,1.,2.,6.,24.,120.,720.,5040.,40320.,362880.,3.6288e6,3.99168e7,
    4.790016e8,6.2270208e9,8.71782912e10,1.307674368e12,2.0922789888e13,
    3.55687428096e14,6.402373705728e15,1.21645100408832e17,
    2.43290200817664e18,5.109094217170944e19,1.1240007277776077e21,
    2.585201673888498e22,6.204484017332394e23,1.5511210043330986e25,
    4.0329146112660565e26,1.0888869450418352e28,3.0488834461171387e29,
    8.841761993739702e30,2.6525285981219107e32,8.222838654177922e33,
    2.631308369336935e35,8.683317618811886e36,2.9523279903960416e38,
    1.0333147966386145e40,3.7199332678990125e41,1.3763753091226346e43,
    5.230226174666011e44,2.0397882081197444e46,8.159152832478977e47,
    3.345252661316381e49,1.40500611775288e51,6.041526306337383e52,
    2.658271574788449e54,1.1962222086548019e56,5.502622159812089e57,
    2.5862324151116818e59,1.2413915592536073e61,6.082818640342675e62,
    3.0414093201713376e64,1.5511187532873822e66,8.065817517094388e67,
    4.2748832840600255e69,2.308436973392414e71,1.2696403353658276e73,
    7.109985878048635e74,4.0526919504877214e76,2.3505613312828785e78,
    1.3868311854568984e80,8.32098711274139e81,5.075802138772248e83,
    3.146997326038794e85,1.98260831540444e87,1.2688693218588417e89,
    8.247650592082472e90,5.443449390774431e92,3.647111091818868e94,
    2.4800355424368305e96,1.711224524281413e98,1.1978571669969892e100,
    8.504785885678623e101,6.1234458376886085e103,4.4701154615126844e105,
    3.307885441519386e107,2.48091408113954e109,1.8854947016660504e111,
    1.4518309202828587e113,1.1324281178206297e115,8.946182130782976e116,
    7.156945704626381e118,5.797126020747368e120,4.753643337012842e122,
    3.945523969720659e124,3.314240134565353e126,2.81710411438055e128,
    2.4227095383672734e130,2.107757298379528e132,1.8548264225739844e134,
    1.650795516090846e136,1.4857159644817615e138,1.352001527678403e140,
    1.2438414054641308e142,1.1567725070816416e144,1.087366156656743e146,
    1.032997848823906e148,9.916779348709496e149,9.619275968248212e151,
    9.426890448883248e153,9.332621544394415e155,9.332621544394415e157,
    9.42594775983836e159,9.614466715035127e161,9.90290071648618e163,
    1.0299016745145628e166,1.081396758240291e168,1.1462805637347084e170,
    1.226520203196138e172,1.324641819451829e174,1.4438595832024937e176,
    1.588245541522743e178,1.7629525510902446e180,1.974506857221074e182,
    2.2311927486598138e184,2.5435597334721877e186,2.925093693493016e188,
    3.393108684451898e190,3.969937160808721e192,4.684525849754291e194,
    5.574585761207606e196,6.689502913449127e198,8.094298525273444e200,
    9.875044200833601e202,1.214630436702533e205,1.506141741511141e207,
    1.882677176888926e209,2.372173242880047e211,3.0126600184576594e213,
    3.856204823625804e215,4.974504222477287e217,6.466855489220474e219,
    8.47158069087882e221,1.1182486511960043e224,1.4872707060906857e226,
    1.9929427461615188e228,2.6904727073180504e230,3.659042881952549e232,
    5.012888748274992e234,6.917786472619489e236,9.615723196941089e238,
    1.3462012475717526e241,1.898143759076171e243,2.695364137888163e245,
    3.854370717180073e247,5.5502938327393044e249,8.047926057471992e251,
    1.1749972043909107e254,1.727245890454639e256,2.5563239178728654e258,
    3.80892263763057e260,5.713383956445855e262,8.62720977423324e264,
    1.3113358856834524e267,2.0063439050956823e269,3.0897696138473508e271,
    4.789142901463394e273,7.471062926282894e275,1.1729568794264145e278,
    1.853271869493735e280,2.9467022724950384e282,4.7147236359920616e284,
    7.590705053947219e286,1.2296942187394494e289,2.0044015765453026e291,
    3.287218585534296e293,5.423910666131589e295,9.003691705778438e297,
    1.503616514864999e300,2.5260757449731984e302,4.269068009004705e304,
    7.257415615307999e306};
    return fact[n];
}

    //utility function
std::ostream& operator<<(std::ostream& os, const SingleFile& sf)
{
  os<<"Single File of "<<sf.N<<" elements, the tagged particle is the number "<<sf.tag;
  os<<", Nl "<<sf.nl<<", Nr is "<<sf.nr<<", and the min is "<<sf.min;   
  os<<", the eigenfunction considered are: "<<sf.eigenfunction_store.size()<<", ";
  os<<"the maximum eigenvalue is "<<sf.max_many_eig<<".\n";
  // for(auto& v : sf.eigenfunction_store){
  //  for( auto e : v)
  //      os<<e<<" ";
  //  os<<std::endl;
  // }    
  return os;
}


//####################################################################################################
//   Abstract Base Class SingleFile
//####################################################################################################

//####################################################################################################
//   public functions
//####################################################################################################

    
//it saves all the single particles eigennumbers such that they satify eigenfunction_condition and that are
//smaller that \Lambda_M, see file integerpartitions.hpp.
void SingleFile::eigenfunction_store_init()
{
  for(int n=0; n<=max_many_eig; ++n){
    std::vector<std::vector<int>> tmp=IntegerPartitions::integer_partitions(n, N, *this);
    for(const auto& i : tmp)
      eigenfunction_store.push_back(i);
  } 
}

//Special case of Eq.16 for the invariant measure
double SingleFile::eq_prob(const double x) const 
{
  static std::vector<int> zeros(N);
  return V0k_element(zeros, x);
}


//Numerator Eq.10, it uses the fast algorithm
double SingleFile::joint2dens(const double x, const double t, const double x0)
{
  double  j2d=0.0;
  for(auto& eigenfunction : eigenfunction_store){// all these elements will have the same manybody eigenvalue
    double V0k=V0k_element(eigenfunction, x);
    double Vk0=Vk0_element(eigenfunction, x0);
    double lambda=lambda_fp(eigenfunction);
    j2d+=V0k*Vk0*std::exp(-lambda*t);
  }
  return j2d;
}

/* The three-point return correlation functions that used to sit here are an
 * extension by Janik Schuettler and now live in SingleFileCorrelations.cpp,
 * so that this file stays as close to upstream BetheSF as possible. */

//Numerator Eq.10, it uses the naive algorithm
double SingleFile::joint2dens_naive(const double  x, const double t, const double x0)
{
  std::vector<int> zeros(N);
  double  j2d=0.0;
  for(auto& eigenfunction : eigenfunction_store){// all these elements will have the same many-body eigenvalue
    double V0k=Vkl_element_naive(zeros, eigenfunction, x);
    double Vk0=Vkl_element_naive(eigenfunction, zeros, x0);
    double lambda=lambda_fp(eigenfunction);
    j2d+=V0k*Vk0*std::exp(-lambda*t);
  }
  return j2d;
}

//It allows copying eigenfunction_store from another SingleFile
void SingleFile::get_eigenfunctions(const SingleFile& sf)
{
  if(N==sf.N)
    eigenfunction_store=sf.eigenfunction_store;
  else
    throw NotAllowedParameters();
}

//It allows copying eigenfunction_store from the input container
void SingleFile::get_eigenfunctions(const std::vector<std::vector<int>>& vec)
{
  if(N==static_cast<int>(vec[0].size()))
    eigenfunction_store=vec;
  else
    throw NotAllowedParameters();
}

//it changes the tagged-particle
void SingleFile::change_tag(const int i)
{
  if(i>0 && i<=N){
    tag=i;
    nl=tag-1; nr=N-tag;
    min=std::min(nl,nr);
  }
  else
    throw NotAllowedParameters();
}

//Changes the \Lambda_M. The container eigenfunction_store is cleared and refilled
void SingleFile::change_MaximumEigenvalue(const int mme)
{
  if(mme>=0){
    max_many_eig=mme;
    eigenfunction_store.clear();
    for(int n=0; n<=max_many_eig; ++n){
      std::vector<std::vector<int>> tmp=IntegerPartitions::integer_partitions(n, N, *this);
      for(const auto& i : tmp)
        eigenfunction_store.push_back(i);
    }
  }
  else
    throw NotAllowedParameters();
}
    
//####################################################################################################
//   VKl_elemenst_naive
//####################################################################################################

//First product Eq. 16
double SingleFile::left_integral(const std::vector<int>& left_k, const std::vector<int>& left_l, const double x) const 
{
  if(left_k.empty())
    return 1.0;
  else{
    double left_contr=1.0;
    for(unsigned i=0; i<left_k.size(); ++i)
      left_contr*=lefttagged( left_k[i], left_l[i], x);
    return left_contr/factorial(left_k.size());
  }
}
    
//Second product Eq.16
double SingleFile::right_integral(const std::vector<int>& right_k, const std::vector<int>& right_l, const double x) const 
{
  if(right_k.empty())
    return 1.0;
  else{
    double right_contr=1.0;
    for(unsigned i=0; i<right_k.size(); ++i)
      right_contr*=righttagged(right_k[i], right_l.at(i), x);
    return right_contr/factorial(right_k.size());
  }
}

//Slavish implementation of Eq.16
double SingleFile::Vkl_element_naive(std::vector<int>& k_vec, std::vector<int>& l_vec, const double x) const
{
  double intres=0.0;
  do{
    std::vector<int> left_k(k_vec.begin(), k_vec.begin()+tag-1);
    int lam_k=k_vec[tag-1];
    std::vector<int> right_k(k_vec.begin()+tag, k_vec.end());
    do{
        std::vector<int> left_l(l_vec.begin(), l_vec.begin()+tag-1);
        int lam_l=l_vec[tag-1];
        std::vector<int> right_l(l_vec.begin()+tag, l_vec.end());
        double t=tagged(lam_k, lam_l, x);
        double l=left_integral(left_k, left_l, x);
        double r=right_integral(right_k, right_l, x);
        intres+=t*l*r;
    }while(std::next_permutation(l_vec.begin(), l_vec.end()));
  }while(std::next_permutation(k_vec.begin(), k_vec.end()));
  return intres;
  // return multiplicity(l_vec)*intres;
}


//####################################################################################################
//   Vkl_elements fast
//####################################################################################################
//####################################################################################################
//   private functions
//####################################################################################################

//First product Eq.16
double SingleFile::left_integral(const std::vector<std::pair<int,int>>& left, const double x) const
{
  if(left.empty()){
    return 1.0;
  }
  else{
    double res=1.0;
    for(const auto p : left)
      res*=lefttagged(p.first, p.second, x);
    return res;
  }
}

//Second product Eq.16
double SingleFile::right_integral(const std::vector<std::pair<int,int>>& right, const double x) const
{
  if(right.empty()){
    return 1.0;
  }
  else{
    double res=1.0;
    for(const auto p : right)
      res*=righttagged(p.first, p.second, x);
    return res;
  }
}

// syntactic sugar for the function: permuting_short_string
std::vector<std::pair<int,int>> multvec(std::vector<int>& v, std::vector<int>& u)// u and v have the same size
{
  std::vector<std::pair<int,int>> res;
  res.reserve(v.size());
  for(unsigned i=0; i<v.size(); ++i){
    res.push_back(std::make_pair(v[i], u[i]));
  }
  return res;
}

//called by the friend class UCombinations, see its file
double SingleFile::combinations_contribution(const std::vector<std::pair<int,int>>& remaining_pairs,
                         const std::vector<std::pair<int,int>>& short_combination, const std::pair<int, int>& unique_pair, const double x) const
{    
    std::vector<std::pair<int,int>> l_comb;
    auto s_comb=short_combination;
    std::sort(s_comb.begin(),s_comb.end());
    //needs sorted elements, remaining_pairs is sorted
    std::set_difference(remaining_pairs.cbegin(), remaining_pairs.cend(),
                s_comb.cbegin(), s_comb.cend(),
                std::back_inserter(l_comb));
    double n_perm_short=number_permutations(s_comb);//number of equivalent contributions of the particles to the left
    double n_perm_long=number_permutations(l_comb);//number of equivalent contributions of the particles to the right
    double t=tagged(unique_pair.first, unique_pair.second, x); //contribution of the tagged particle
    double l, r;
    if(nl<=nr){ 
        l=left_integral(s_comb, x); //contribution of the particles to the left
        r=right_integral(l_comb, x); //contribution of the particles to the right
    }
    else{
        l=left_integral(l_comb, x);
        r=right_integral(s_comb, x);
    }
    return n_perm_short*n_perm_long*t*l*r;
}

// //called by the friend class UCombinations, see its file
// double SingleFile::combinations_contribution(const std::vector<std::pair<int,int>>& remaining_pairs,
//                const std::vector<std::pair<int,int>>& short_combination, const std::pair<int, int>& unique_pair, const double x) const
// { 
//   std::vector<std::pair<int,int>> l_comb;
//   std::set_difference(remaining_pairs.cbegin(), remaining_pairs.cend(),
//           short_combination.cbegin(), short_combination.cend(),
//           std::inserter(l_comb , l_comb.begin()));
//   double n_perm_short=number_permutations(short_combination);//number of equivalent contributions of the particles to the left
//   double n_perm_long=number_permutations(l_comb);//number of equivalent contributions of the particles to the right
//   double t=tagged(unique_pair.first, unique_pair.second, x); //contribution of the tagged particle
//   double l, r;
//   if(nl<=nr){ 
//     l=left_integral(short_combination, x); //contribution of the particles to the left
//     r=right_integral(l_comb, x); //contribution of the particles to the right
//   }
//   else{
//     l=left_integral(l_comb, x);
//     r=right_integral(short_combination, x);
//   }
      
//   return n_perm_short*n_perm_long*t*l*r;
// }

double SingleFile::permuting_short_string(std::vector<int>& short_string, std::vector<int>& long_string, const double x) const
{
  double intres=0.0;
  do{
    std::vector<std::pair<int,int>> pairs=multvec(long_string,short_string);
    if(std::is_sorted(pairs.cbegin(), pairs.cend())){//Takes care I am not overcounting equivalent pairs
      //std::next_permutation does not distinguish identical pairs
      std::set<std::pair<int,int>> unique_pairs(pairs.cbegin(), pairs.cend());
      for(auto& unique_pair : unique_pairs){
        const std::array<std::pair<int, int>,1> chosen={unique_pair};
        std::vector<std::pair<int,int>> remaining_pairs;
        std::set_difference(pairs.cbegin(), pairs.cend(),
                            chosen.cbegin(), chosen.cend(),
                            std::inserter(remaining_pairs, remaining_pairs.begin()) );
        if(tag == 1 || tag == N){ //faster route if the first or last particle are tagged
          double n_perm_long=number_permutations(remaining_pairs);
          double t=tagged(unique_pair.first, unique_pair.second, x);
          double rem;
          if(tag == 1)
            rem=right_integral(remaining_pairs, x);
          else
            rem=left_integral(remaining_pairs, x);
          intres+=n_perm_long*t*rem;
        }
        else{
          UCombinations<std::pair<int,int>> com_gen;
          intres+=com_gen.unique_combinations(remaining_pairs, min, unique_pair, x, *this);
        }//end else clause if the particle is not the first or the last
      }// end for loop on the unique pairs
    }//end is_sorted
  } while(std::next_permutation(short_string.begin(), short_string.end()));
  return intres;
}

//Final implementation of Eq.16
double SingleFile::Vkl_element(std::vector<int>& k_vec, std::vector<int>& l_vec, const double x) const
{
  double intres;
  double n_perm_k=number_permutations(k_vec), n_perm_l=number_permutations(l_vec);
  if(n_perm_l <= n_perm_k){
    intres=permuting_short_string(l_vec, k_vec, x); 
  }
  else{
    intres=permuting_short_string(k_vec, l_vec, x);
  }
  return multiplicity(l_vec)*intres/(factorial(nl)*factorial(nr));
}

//Eq. 13
double SingleFile::V0k_element(std::vector<int>& k_vec, const double x) const 
{
  std::vector<int> zeros(N);
  return Vkl_element(zeros, k_vec, x);
}
//Eq. 13     
double SingleFile::Vk0_element(std::vector<int>& k_vec, const double x) const 
{
  std::vector<int> zeros(N);
  return Vkl_element(k_vec, zeros, x);
}


}//end namespace SingleFile
