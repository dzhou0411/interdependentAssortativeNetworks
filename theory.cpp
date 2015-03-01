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
  initsrand(1);

  double ER_g(double p, double k);//1
  void ER_step(int k,double p_0,double p_protect1,double p_protect2);//2
  //void step_si(int num_nodes, int avg_degree,int realization,int type,
  //	       double q1,double q2,double p);//3
  pair<double,double> fc_pc(int k, double q);//4
  void SF_intential_attack_step(double lamda, double p_0);
  const int num_nodes=100000;
  
  const int avg_degree=4;
  const double lamda=2.8;
  const int realization=1;
  int type=1; //ER:1, SF:2
  //---protect portion.
  double q1_start=0.0;
  double q1_end=q1_start+1;
  double q1_delta=0.01;

  double q2_start=0.15;
  double q2_end=q2_start+0.03;
  double q2_delta=0.04;
  //---attack portion
  double p_start=0.10;
  double p_end=p_start+0.04;
  double p_delta=0.05;
  //----------------
  ofstream ofile;
  string file;
  //cout<<P_inf(lamda,p_start)<<endl;
  double K;
  for(double p=9.99;p>0.5;p-=0.01){
    K=2*pow(1-p,1./(1-lamda));
    //cout<<K<<endl;
    cout<<p<<' '<<SF_g_B(lamda,K,2,p)<<endl;
  }
  /*
  //SF_intential_attack_step(lamda,  p_start);
  file="../data/temp.dat";
  ofile.open(file.c_str());
  double SF_g_B(double lamda,double K,double m,double p);
  for(double p=0.99;p>0.9;p-=0.01){
    double pp=1-pow(1-p,(2.-lamda)/(1-lamda));
    double K=2*pow(1-p,1./(1.-lamda));
    cout<<pp<<' '<<K<<endl;
    ofile<<p<<' '<<p*SF_g_B(lamda,K,2,p)<<endl;
    }*/
  //ER_step_th(avg_degree,p_start ,q1_start,q2_start);
  //step_si(num_nodes,avg_degree,realization,type,q1_start,q2_start,p_start);
  //pair<double,double> a=fc_pc(avg_degree, q1_start);
  //cout<<q1_start<<' '<<a.first<<' '<<a.second<<endl;

  /*
  file="../data/pc_muinf/pc_muinf_fc_k"+d_c(avg_degree)+".dat";
  ofile.open(file.c_str());
  pair<double,double> a;
  for(double q=q1_start;q<=q1_end;q+=q1_delta){
    a=fc_pc(avg_degree, q);
    double muinf=q+(1-q)*a.second*(1-a.first);
    ofile<<q<<' '<<a.second<<' '<<muinf<<' '<<a.first<<endl;
    if(a.second==0) break;
  }
  ofile.close();//*/
}


//-----------//
//----1------//
double ER_g(double p, double k){
  double upper=1.,lower=0.;
  double f=0.5;
  double delta=0.1;
  while(delta>0.00001){
    delta=exp(k*p*(f-1.))-f;
    if(delta>0){
      lower=f;
      f=(f+upper)/2.;
    }
    else{
      upper=f;    
      f=(f+lower)/2.;
    }
    delta=abs(delta);
  }
  return (1.-f);
}
//----2------//
void ER_step_th(int k,double p_0,double p_protect1,double p_protect2){
  double q1=p_protect1,q2=p_protect2;
  string file="../data/step/ER_randprotect_Th.dat";
  ofstream ofile(file.c_str());
  vector<double> p(7,0.);
  p[0]=p_0;
  p[1]=q1+(1-q1)*p[0];
  p[2]=p[1]*ER_g(p[1],double(k));
  p[3]=q2+(1-q2)*p[0]*ER_g(p[1],double(k));
  p[4]=p[3]*ER_g(p[3],double(k));
  ///*
  p[5]=q1+(1-q1)*p[0]*ER_g(p[3],double(k));
  p[6]=p[5]*ER_g(p[5],double(k));
  int i=p.size();
  while(abs(p[i-1]-p[i-5])>0.0001){
    if(i%4==1) {
      p.push_back(q1+(1-q1)*p[0]*ER_g(p[i-2],double(k)));
      i=p.size();
      p.push_back(p[i-1]*ER_g(p[i-1],double(k)));
      i=p.size();
    } else if(i%4==3){
      p.push_back(q2+(1-q2)*p[0]*ER_g(p[i-2],double(k)));
      i=p.size();
      p.push_back(p[i-1]*ER_g(p[i-1],double(k)));
      i=p.size();
    }
    cout<<i<<' '<<p[i-1]<<endl;
  }
  //*/
  for(int i=2;i<p.size();i+=4){
    if((p[i]!=0)&&(p[i+2]!=0))
      ofile<<i/4+1<<' '<<p[i]<<' '<<p[i+2]<<endl;
  }
  ofile.close();
}

//-----3------//
pair<double,double> fc_pc(int k, double q){
  double upper=100,lower=0.;
  double f=0.5;
  double delta=0.1;
  while(delta>0.00000001){
    delta=exp((1-k*q*f)*(f-1)/2/f-k*q*(1-f))-f;
    if(delta<0){
      lower=f;
      f=(f+upper)/2.;
    }
    else{
      upper=f;    
      f=(f+lower)/2.;
    }
    delta=abs(delta);
  }
  double p=(1./k-q*f)/2./(1-q)/(1-f)/f;
  if(p<0||p>1) p=0;
  return(make_pair(f,p));
}

//-----4------//
/*void step_si(int num_nodes, int avg_degree,int realization,
	     int type,double q1,double q2,double p){
  double q1_start=q1, q2_start=q2, p_start=p;
  double q1_end=q1+0.03,q2_end=q2+0.03, p_end=p+0.03;
  double q1_delta=0.04,q2_delta=0.04, p_delta=0.04;
  random_protect(num_nodes,avg_degree,realization,type,q1_start,q1_end,
		 q1_delta,q2_start,q2_end,q2_delta,p_start,p_end,p_delta);
  vector<vector<double> > v(realization);
  string file;
  ifstream ifile;
  ofstream ofile;
  for(int i=0;i<realization;++i){
    file="../data/temp"+d_c(i)+".dat";
    ifile.open(file.c_str());
    while(!ifile.eof()){
      if(!ifile.eof()){
	ifile>>p;
	v[i].push_back(p);
      }
    }
    ifile.close();
  }
  int k;
  file="../data/step/ER_randprotect_Si.dat";
  ofile.open(file.c_str());
  for(int i=0;i<v[0].size();++i){
    p=0.;k=0;
    for(int j=0;j<realization;++j)
      if(i<v[j].size()){
	p+=v[j][i];
	k++;
      }
    if(i%2==0) ofile<<i/2+1<<' '<<p/k<<' ';
    if(i%2==1) ofile<<p/k<<endl;
  }
}
*/
