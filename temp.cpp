using namespace std;
#include <fstream>
#include <set>
#include <string>
#include "common.hpp"
#include "graph.hpp"
#include "net_ops.hpp"
#include "net_algo.hpp"
#include "node.hpp"

int main(){
  
  /*pair<int,double> K_tilde(double lambda, int mm, int avg_k);
  
  for(double lambda=2.1;lambda<3;lambda+=0.1){
    pair<int,double> K_delta=K_tilde(lambda,3,4);
    cout<<lambda<<' '<<K_delta.first<<' '<<K_delta.second<<endl;
    }*/
  
  pair<double,double> p_inf_p_c(double avg_k1,double avg_k2);

  for(int k=1;k<=1;k++){
    pair<double,double> pinf_pc=p_inf_p_c(k,k);
    cout<<k<<' '<<pinf_pc.first<<' '<<pinf_pc.second<<endl;
  }
}

pair<int,double> K_tilde(double lambda, int mm, int avg_k){
  int lower=mm, upper=1000000;
  int K=upper;
  double delta=1;
  
  while(upper-lower>1){
    if(delta>0){
      upper=K;
      K=(lower+upper+1)/2;
    }else if(delta<0){
      lower=K;
      K=(lower+upper+1)/2;
    }else if(delta==0){
      break;
    }
    double avg_k_2=0;
    for(int k=mm;k<=K;++k)
      avg_k_2+=double(k)/(pow(K,1-lambda)-pow(mm,1-lambda))*(pow(k+1.,1-lambda)-pow(k,1-lambda));
    delta=avg_k_2-avg_k;
  }
  return make_pair(K,delta);
}


pair<double,double> p_inf_p_c(double avg_k1,double avg_k2){
  double upper=1, lower=0;
  double p_inf=upper;
  double delta=-1;
  while(abs(delta)>0.0001){
    if(delta<0){
      upper=p_inf;
      p_inf=(lower+upper)/2;
    }else {
      lower=p_inf;
      p_inf=(lower+upper)/2;
    }
    delta=p_inf-(exp(avg_k1*p_inf)/avg_k1/avg_k2*(exp(avg_k2-avg_k2*exp(-avg_k1*p_inf))-1));
    cout<<p_inf<<' '<<upper<<' '<<lower<<' '<<delta<<endl;
  }
  double p_c=exp(avg_k1*p_inf)/avg_k1/avg_k2*exp(avg_k2-avg_k2*exp(-avg_k1*p_inf));
  return make_pair(p_inf,p_c);
}
