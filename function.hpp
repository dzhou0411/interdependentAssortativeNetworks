#ifndef FUNCTION_HPP
#define FUNCTION_HPP
#include <cmath>
#include <queue>
#include <stack>
#include <string>
#include "common.hpp"
//-------------------------------//
//1. percolation
//-------------------------------//
void percolation(int num_nodes,int avg_degree,double lamda, int realization,int type,int attack,double p_start,double p_end,double p_delta){
  string file,nu;
  ofstream ofile;
  double p,p_inf;
  Graph net_1_prime(num_nodes),net_2_prime(num_nodes);
  Graph net_1(num_nodes),net_2(num_nodes);
  map<double,double> map_p;
  int index;
  vector<pair<int,int> > links;
  for(int i=0;i<realization;++i){
    cout<<i<<endl;
    net_1.clear();net_2.clear();links.clear();
    if(type==1){
      lt_ER_algo2(num_nodes,net_1,avg_degree,links);
      lt_ER_algo(num_nodes,net_2,avg_degree);
      nu="ER"+d_c(avg_degree);
    }
    else if(type==2){
      lt_SF_algo2(num_nodes,net_1,lamda,links);
      lt_SF_algo(num_nodes,net_2,lamda);
      nu="SF"+d_c(int(lamda*10.0001));
    }
    if(attack==1) {nu+="_random";}
    else if(attack==2) {nu+="_maxdeg";}
    else if(attack==3) {nu+="_intentional";}
    for(p=p_start;p<=p_end;p+=p_delta){
      p_inf=0.;
      net_1_prime.clear();net_2_prime.clear();
      net_1_prime=net_1;net_2_prime=net_2;
      if(attack==1) {rand_init_attack(p, net_1_prime);}
      else if(attack==2) {max_degree_init_attack(p, net_1_prime);}
      else if(attack==3) {
	intentional_attack(p,net_1_prime,links);
      }
      p_inf=cascade_failure1(net_1_prime,net_2_prime,i);
      map_p[p]+=p_inf;
    }
  }
  file="../data/percolation/"+nu+".dat";
  cout<<file<<endl;
  ofile.open(file.c_str());
  for(p=p_start;p<=p_end;p+=p_delta)
    ofile<<p<<' '<<map_p[p]/double(realization)<<' '<<map_p[p]/double(realization)/p<<endl;
  ofile.close();
}

//---------------------------//
//2.step simulation
// read temp.dat files and average
//---------------------------//
void step_intential_attack_simu(double p_0,int num_nodes,int type,double lamda,int avg_degree,int realization){
  string file;
  if(type==2) 
    file="../data/step/SF_maxdeg_Sim"+d_c(int(lamda*10.001))+".dat";
  else if(type==1)   
    file="../data/step/ER_maxdeg_Sim"+d_c(int(p_0*100))+"_k"+d_c(avg_degree)+".dat";
  ofstream ofile(file.c_str());
  vector<double> v(90);
  ifstream ifile;
  double a;
  for(int i=0;i<realization;++i){
    file="../data/temp"+d_c(i)+".dat";
    ifile.open(file.c_str());
    int j=0;
    while(!ifile.eof()){
      ifile>>a;
      v[j++]+=a/realization;
    }
    while(j<50)
      v[j++]+=a/realization;
    ifile.close();
  }
  for(int i=0;i<50;++i)
    ofile<<v[i]<<endl;
  ofile.close();
}
/*
//--------------------------//
//3.1 step theory for intentional attack
//--------------------------//
//------random attack g
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
//---intentional attack g_A;
//parameter p_0 initial survive ratio
//---------------------//
double ER_g_A(double p,double p_0, double k){
  double f=log(p_0)/k+1;
  double p_tild=f*exp(k*(f-1));
  double f_A=0.5;
  double upper=1.,lower=0.;
  double delta=0.1;
  while(delta>0.000001){
    delta=1./p_0*exp(k*(f-1+p_tild/p_0*f*(p*f_A-p)))-f_A;
    if(delta>0){
      lower=f_A;
      f_A=(f_A+upper)/2.;
    }else{
      upper=f_A;
      f_A=(f_A+lower)/2.;
    }
    delta=abs(delta);
  }
  return 1-f_A;
}
void step_ER_theo(double p_0,double k){
  string file="../data/step/ER_intentional_Th"+d_c(int(p_0*100))+"_k"+d_c(k)+".dat";
  ofstream ofile(file.c_str());
  vector<double> p(7,0.);
  p[0]=p_0;
  p[1]=p[0];
  p[2]=p[1]*ER_g_A(p[1],p_0,double(k));
  p[3]=p[0]*ER_g_A(p[1],p_0,double(k));
  p[4]=p[3]*ER_g(p[3],double(k));
  ///*
  p[5]=p[0]*ER_g(p[3],double(k));
  p[6]=p[5]*ER_g_A(p[5],p_0,double(k));
  int i=p.size();
  while(abs(p[i-1]-p[i-5])>0.0001){
    if(i%4==1) {
      p.push_back(p[0]*ER_g(p[i-2],double(k)));
      i=p.size();
      p.push_back(p[i-1]*ER_g_A(p[i-1],p_0,double(k)));
      i=p.size();
    } else if(i%4==3){
      p.push_back(p[0]*ER_g_A(p[i-2],p_0,double(k)));
      i=p.size();
      p.push_back(p[i-1]*ER_g(p[i-1],double(k)));
      i=p.size();
    }
  }
  for(int i=2;i<p.size();i+=2){
    if((p[i]!=0)&&(p[i+2]!=0))
      ofile<<i/2-1<<' '<<p[i]<<endl;
  }
  ofile.close();
} 
//--------------------------//
//3.2 step theory for max degree attack in SF
//--------------------------//
double SF_g_A_maxdeg(double p,double p_0, double lamda){
  double K_tild=2*pow(1-p_0,1./(1-lamda));
  //double p_tild=1-pow(1-p_0,(2-lamda)/(1-lamda));
  double a=0,b=0;
  for(int i=2;i<90000;++i){
    a+=i*(pow(i,1-lamda)-pow(i+1,1-lamda))/pow(2.,1-lamda);
  }
  for(int i=2;i<K_tild+1;++i){
    b+=i*(pow(i,1-lamda)-pow(i+1,1-lamda))/pow(2,1-lamda);
  }
  double p_tild=b/a;
  double c=(1-lamda)/(0-pow(2,1-lamda));
  double f_A=0.5;
  double upper=1.,lower=0.;
  double delta=0.1;
  double G1,k=0;
  for(int i=2;i<K_tild+1;++i)
    k+=i*(pow(i,1-lamda)-pow(i+1,1-lamda))/pow(2,1-lamda)/p_0;
  while(delta>0.000001){
    G1=0;
    for(int i=2;i<K_tild+1;++i)
      G1+=i*(pow(i,1-lamda)-pow(i+1,1-lamda))/pow(2,1-lamda)/p_0*pow(1+p_tild/p_0*(f_A*p-p),i-1);
    G1/=k;
    delta=G1-f_A;
    if(delta>0){
      lower=f_A;
      f_A=(f_A+upper)/2.;
    }else{
      upper=f_A;
      f_A=(f_A+lower)/2.;
    }
    delta=abs(delta);
  }
  double G0=0;
  for(int i=2;i<K_tild+1;++i)
    G0+=(pow(i,1-lamda)-pow(i+1,1-lamda))/pow(2,1-lamda)/p_0*pow(1+p_tild/p_0*(f_A*p-p),i);
  return 1-G0;
}
void step_SF_theo(double p_0,double lamda){
  string file="../data/step/SF_maxdeg_Th"+d_c(int(p_0*100))+"_lamda"+d_c(lamda*10.001)+".dat";
  ofstream ofile(file.c_str());
  vector<double> p(7,0.);
  p[0]=p_0;
  p[1]=p[0];
  p[2]=p[1]*SF_g_A_maxdeg(p[1],p_0,lamda);
  p[3]=p[0]*SF_g_A_maxdeg(p[1],p_0,lamda);
  p[4]=p[3]*SF_g_A(p[3],lamda);
  ///*
  p[5]=p[0]*SF_g_A(p[3],lamda);
  p[6]=p[5]*SF_g_A_maxdeg(p[5],p_0,lamda);
  int i=p.size();
  while(abs(p[i-1]-p[i-5])>0.0001){
    if(i%4==1) {
      p.push_back(p[0]*SF_g_A(p[i-2],lamda));
      i=p.size();
      p.push_back(p[i-1]*SF_g_A_maxdeg(p[i-1],p_0,lamda));
      i=p.size();
    } else if(i%4==3){
      p.push_back(p[0]*SF_g_A_maxdeg(p[i-2],p_0,lamda));
      i=p.size();
      p.push_back(p[i-1]*SF_g_A(p[i-1],lamda));
      i=p.size();
    }
  }
  for(int i=2;i<p.size();i+=2){
    if((p[i]!=0)&&(p[i+2]!=0))
      ofile<<i/2-1<<' '<<p[i]<<endl;
  }
  ofile.close();
} 
//----------------------------//
//4 looking for ER critical point
//----------------------------//
//4.1 pg_A(pg_B(x))
void function1(double p_0,int avg_deg){
  double f_B,f_A,y;
  double lower,upper,delta;
  double f=log(p_0)/avg_deg+1;
  double p_tild=f*exp(avg_deg*(f-1));
  string file;ofstream ofile;
  file="../data/pg_A(pg_B(x)).dat";
  ofile.open(file.c_str());
  for(double x=0.;x<=1.;x+=0.02){
    cout<<x<<endl;
    //find f_B
    upper=1.;lower=0;delta=0.1;
    f_B=0.5;
    while(delta>0.000001){
      delta=exp(avg_deg*x*(f_B-1))-f_B;
      if(delta>0){
	lower=f_B;
	f_B=(f_B+upper)/2.;
      }else{
	upper=f_B;
	f_B=(f_B+lower)/2.;
      }
      delta=abs(delta);
    }
    y=p_0*(1-f_B);    
    upper=1.;lower=0;delta=0.1;
    f_A=0.5;
    while(delta>0.000001){
      delta=1./p_0*exp(avg_deg*(p_tild/p_0*f*(f_A-1)*y+f-1))-f_A;
      if(delta>0){
	lower=f_A;
	f_A=(f_A+upper)/2.;
      }else{
	upper=f_A;
	f_A=(f_A+lower)/2.;
      }
      delta=abs(delta);
    }
    ofile<<x<<' '<<p_0*(1-f_A)<<endl;
  }
}*/
#endif
