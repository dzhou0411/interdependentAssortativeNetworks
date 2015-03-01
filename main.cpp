using namespace std;
#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include <string>
#include "common.hpp"
#include "graph.hpp"
#include "net_ops.hpp"
#include "net_algo.hpp"
//#include "net_theory.hpp"
#include "node.hpp"
#include "./statool/srand.hpp"

int main(){
  //-------------------//
  //   net_algo.hpp    //
  //-------------------//
  //a.
  void lt_ER_algo(int num_nodes,Graph & net, const int avg_degree);
  void lt_ER_algo2(int num_nodes,Graph & net, const int avg_degree, vector<pair<int,int> >& links);
  void lt_SF_algo(int num_nodes,Graph & net, const double lembda,int mm);
  void lt_SF_algo2(int num_nodes,Graph & net, const double lembda, vector<pair<int,int> >& links, int mm);
  //-------------------//
  //    net_ops.hpp    //
  //-------------------//
  //b.1.
  void find_biggest_cluster(Graph & net, vector<int> & biggest_cluster);
  //***b.2.1 record size of cluster in each step in temp.dat; output ulti-size of cluster
  double cascade_failure1(Graph &net_1, Graph &net_2, int realization);
  //b.2.2 record the ids of those nodes which die in each step in v1 & v2
  double cascade_failure2(Graph &net_1, Graph &net_2, int real, vector<int> &v1, vector<int> &v2);
  //b.2.3 cascade failure by arbitrary dependance defined between two nets
  double cascade_failure3(Graph &net_1, Graph &net_2,vector<int> &dep1,vector<int> &dep2);
  
  //***b.3.1 initial attack
  void rand_init_attack(double p, Graph &net);
  //b.3.2 
  void max_degree_init_attack(double  p, Graph &net);
  //b.3.3 a node has k*p(k)/<k> possibility being removed
  void intentional_attack(double p, Graph &net,vector<pair<int,int> > &links);
  //b.3.4 a node is removed proportional to k^alpha
  void intentional_attack_alpha(double p,Graph &net,double alpha);

  //b.4.1 protection of nodes randomly
  void random_protect(Graph &net_1,Graph &net_2,double q1,double q2);
  //b.4.2 intentional protection of net_2;
  //    select= 1.degree; 2.btwness; 3.aved
  //    type= 0.rank by absolute value; 1.difference
  void max_protect(int select,int type,Graph &net_1,Graph &net_2,double q);
  
  //b.5 cascade failure simulation
  //b.5.1 without protection
  void percolation_simulation1(double alpha,int num_nodes,int avg_degree,double lembda,int realization,int type,int attack,double p_start,double p_end,double p_delta,int mm);
  //b.5.2 double attack
  void percolation_simulation2(double alpha,int num_nodes,int avg_degree,double lembda,int realization,int type,int attack,double p_start,double p_end,double p_delta,int mm);
  //b.5.3 step of cascades failure by simulation
  void step_simu(double alpha,double p_0,int num_nodes,int type,int attack,double lembda,int avg_degree,int realization,int mm);

  //b.6 pdf of k after removed nodes but not edges
  void pdf_k_intent_alpha(double p_0,int num_nodes,double avg_k,double lembda,double alpha,int M);
  
  
  //--------------------//
  //  net_theory.hpp    //
  //--------------------//
  //c.1.1 ER g(p)
  double ER_g_A(double p, double avg_deg);
  //c.1.2 SF with K=inf, m=2
  double SF_g_A(double p, double lembda,int mm);
  //c.1.3 SF with K>k>m
  double SF_g_B(double lembda,double K,double m,double p);
  //c.1.4 g(p) for intentional attack ER; parameter p_0 initial survive ratio
  double ER_g_A_intent(double p,double p_0, double k);
  //c.1.5 g(p) for intentional attack SF; parameter p_0 initial survive ratio
  double SF_g_A_intent(double p,double p_0,double lembda,int mm);
  //c.1.6 g(p) for maxdeg attack ER
  double ER_g_A_maxdeg(double p,double p_0,double k);
  //c.1.7 g(p) for maxdeg attack SF
  double SF_g_A_maxdeg(double p,double p_0,double k);
  //c.1.8 g(p) for (alpha) intentional attack ER
  double ER_g_A_intent_alpha(double p,double p_0,double avg_k,double alpha);
  double SF_g_A_intent_alpha(double p,double p_0,double lembda,double alpha,int mm);
  //c.1.9 f=G^-1(r) for ER net
  double ER_f_alpha(double alpha,double p_0,double avg_k);
  long double SF_f_alpha(double alpha,double p_0,double lembda,int mm);
  //c.1.10 p_tild for (alpha) ER
  double ER_p_tild_alpha(double alpha,double p_0,double avg_k);


  //c.2.1 normal situation
  void ER_step(int p_0,double k);
  //c.2.2
  void SF_step(double p_0, double lembda,int mm);
  //c.2.3 ER intentional attack situation
  void ER_step_theo_intent(double p_0,double k);
  //c.2.4 SF intentional attack situation
  void SF_step_theo_intent(double p_0,double lembda,int mm);
  //c.2.5 ER maxdeg attack situation
  void ER_step_theo_maxdeg(double p_0,double k);
  //c.2.6 SF maxdeg attack situation
  void SF_step_theo_maxdeg(double p_0,double lembda);
  //c.2.7 ER intentional attack for alpha
  void ER_step_theo_intent_alpha(double p_0,double k,double alpha,int mm);



  //c.3.1 ER-ER,print pg_A(pg_B(x)) for intentional attack for different x
  void function1(double p_0,int avg_deg);
  //c.3.2 ER-ER,pg_A(pg_B(x)) for intentional attack
  double pg_Apg_B(double p_0,double x, double avg_deg,double avg_deg2);
  //c.3.3 ER-ER,find critical point p_c for intentional attack
  double p_c(double avg_deg,double avg_deg2);
  //c.3.4 SF-SF,(3.2)
  //c.3.5 ER_pgApg_B for intent (alpha) attack
  double ERER_pg_Apg_B(double x,double p_0,double alpha,double avg_k,double avg_k2);
  double SFSF_pg_Apg_B(double x,double p_0,double alpha,double lembda,double lembda2,int mm);
  //c.3.6 ER-ER,find p_c for intent (alpha) attack
  double ERER_p_c(double avg_k,double avg_k2,double alpha);
  pair<double,double> SFSF_p_c(double lembda,double lembda2,double alpha,int mm);
  //c.3.7 P_inf for alpha
  double ERER_p_inf_alpha(double p, double avg_k1,double avg_k2,double alpha);
  double SFSF_p_inf_alpha(double p, double lembda1,double lembda2,double alpha,int mm);
  //c.3.8 Double attack for alpha=1
  double ERER_pc_doubleattack(double avg_k);
  double ERER_p_inf_doubleattack(double p,double avg_k);
 //these series of functions work in this way: they try to find the upper intersection of gA(gB(x)) and x. When this intersection position suddenly goes to zero,     then it's the critical point of xc.

 // 031511 ERER pc     
  double graph_xc_ERER(double k);
 // 031511 SFSF pc
  double graph_xc_SFSF(double lamda, int mm);
 // 031511 ER-SF pc 
  double graph_xc_ERSF(double k,double lembda, int mm);
 // 031511 SF-ER pc
  double graph_xc_SFER(double k,double lembda, int mm);
  //c.4.1 
  void ER_step_simu();


  //------------------------------------------------------//
  //------------------------------------------------------//
  initsrand(1);
  const int avg_k=4,avg_k2=4;
  const int num_nodes=10000;
  const double alpha=-1;
  const int realization=20;
  const double lembda=2.5,lembda2=2.5;//
  const int mm=2;
  int type=1;//1:ER, 2:SF
  int attack=1;//1:random,2:max_deg,3:intentional,4:intent alpha
  //---protect portion.
  double q1_start=0.;
  double q1_end=0.0;
  double q1_delta=0.02;
  double q2_start=0.;
  double q2_end=0.0;
  double q2_delta=0.04;
  double p_start=0.613;
  double p_end=0.89;
  double p_delta=0.01;
  //----------------
  
  ofstream ofile,ofile2,ofile3;ifstream ifile;string file;
  vector<pair<int,int> > links;
  

  /* Code of 022111 starts. Aim is to plot pinfinity vs. p for coupled ER-SF network */

  //Graph ERnet2(num_nodes),SFnet2(num_nodes),net1(num_nodes),net2(num_nodes);

  // Graph ERnet1(num_nodes),ERnet2(num_nodes);
  //lt_ER_algo(num_nodes,ERnet1, avg_k);  //E-R network1
  //lt_ER_algo(num_nodes,ERnet2, avg_k); //E-R network2
  //lt_SF_algo(num_nodes,SFnet, 2.5, mm);  //SF network1
  //lt_SF_algo(num_nodes,SFnet2, lembda, mm); //SF network2



//******--------here start the file go through-------******//
for (int     dig1=-10; dig1<=+10;    dig1=dig1+1){
for (int     dig2=85;   dig2<=99;  dig2++){
    
    Graph newnet1(num_nodes),newnet2(num_nodes);    // projgao
    Graph net1(num_nodes),net2(num_nodes);      
     
    string filein,fileout;
    stringstream s1,s2;
    ifstream netdata;
       
    if(dig1==0) s1<<"+00";
    else if(dig1==10) s1<<"+10";
    else if(dig1==-10) s1<<"-10";
    else if(dig1>0) s1<<"+0"<<dig1;
    else if(dig1<0) s1<<"-0"<<(-dig1);   

    if(dig2==0) s2<<"000";
    else if(dig2<10) s2<<"00"<<dig2;
    else if(dig2<100) s2<<"0"<<dig2;
    
    filein="./networkdata/netN10000beta"+s1.str()+"."+s2.str()+".txt";    
    fileout=       "./dataoutput1/PinfVp"+s1.str()+"."+s2.str()+".dat";
/*
    s1<<dig1num;
    s2<<dig2num;

    if (dig1num==-10){
     if (dig2num==0){
     filein="./networkdata/netN10000beta"+s1.str()+".000.txt";    
     fileout=       "./dataoutput/PinfVp"+s1.str()+".000.dat";}
     else if (dig2num<10){
     filein="./networkdata/netN10000beta"+s1.str()+".00"+s2.str()+".txt";    
     fileout=       "./dataoutput/PinfVp"+s1.str()+".00"+s2.str()+".txt";}
     else if (dig2num<100{
     filein="./networkdata/netN10000beta"+s1.str()+".0"+s2.str()+".txt";    
     fileout=       "./dataoutput/PinfVp"+s1.str()+".0"+s2.str()+".txt";}
                    }

    else if (dig1num==+10){
     if (dig2num==0){
     filein="./networkdata/netN10000beta+"+s1.str()+".000.txt";    
     fileout=       "./dataoutput/PinfVp+"+s1.str()+".000.dat";}
     else if (dig2num<10){
     filein="./networkdata/netN10000beta+"+s1.str()+".00"+s2.str()+".txt";    
     fileout=       "./dataoutput/PinfVp+"+s1.str()+".00"+s2.str()+".txt";}
     else if (dig2num<100{
     filein="./networkdata/netN10000beta+"+s1.str()+".0"+s2.str()+".txt";    
     fileout=       "./dataoutput/PinfVp+"+s1.str()+".0"+s2.str()+".txt";}
                    } 



    else if (dig2num==0){
     if (dig1num<0){
     filein="./networkdata/netN10000beta"+s1.str()+".000.txt";    
     fileout=       "./dataoutput/PinfVp"+s1.str()+".000.dat";}
     if (dig1num>=0){
     filein="./networkdata/netN10000beta+"+s1.str()+".000.txt";    
     fileout=       "./dataoutput/PinfVp+"+s1.str()+".000.dat";}    
                    }
    else if (dig2num<10){
    if (dig1num<0){
    filein="./networkdata/netN10000beta"+s1.str()+".00"+s2.str()+".txt";    
    fileout=       "./dataoutput/PinfVp"+s1.str()+".00"+s2.str()+".dat";}
    if (dig1num>=0){
    filein="./networkdata/netN10000beta+"+s1.str()+".00"+s2.str()+".txt";    
    fileout=       "./dataoutput/PinfVp+"+s1.str()+".00"+s2.str()+".dat";}
                    }
    else if (dig2num<100){
    if (dig1num<0){
    filein="./networkdata/netN10000beta"+s1.str()+".0"+s2.str()+".txt";    
    fileout=       "./dataoutput/PinfVp"+s1.str()+".0"+s2.str()+".dat";}
    if (dig1num>=0){
    filein="./networkdata/netN10000beta+"+s1.str()+".0"+s2.str()+".txt";    
    fileout=       "./dataoutput/PinfVp+"+s1.str()+".0"+s2.str()+".dat";}
                    }
*/
    netdata.open(filein.c_str(),ios::in);

    cout<<"doing "<<filein<<endl;

    int a,b;
    int numnodes,numlinks;
    netdata>>numnodes>>numlinks;
    while (netdata>>a>>b){//cout<<a<<" "<<b<<endl;
    newnet1.insert_connection(a,b);}
    newnet2=newnet1;   

// generate dependence vectors dep1 = dep2
vector<int> dep1(num_nodes);
vector<int> dep2(num_nodes);
vector<int> dep_check(num_nodes);
for (int dci=0;dci<num_nodes;dci++) dep_check[dci]=0;
double q=1;
int dep_id;


for (int i=0; i<num_nodes ; i++) { dep1[i]=-1;}


for (int nq=0; nq<num_nodes*q ; nq++) {
	do {dep_id=srand()*num_nodes;} 
	while (dep_check[dep_id]!=0);
        dep1[nq]=dep_id;                //if[dep_id], then same=same; if[nq], then diff=diff 
        dep2[dep_id]=nq;
        dep_check[dep_id]=1;
                                       }

//dep2=dep1;    

//cout<<"??"<<endl;
/*
for (int nq=0; nq<num_nodes ; nq++) {
cout<< nq << "  " << dep1[nq] << "  " <<dep2[dep1[nq]]<<endl;
                                      }
*/

  double p;
  double pinfinity;  
  double pout[101];
  int i;
  for (i=0;i<101;i++) pout[i]=0.;
   

 ofstream fout;
 fout.open(fileout.c_str());



// cout << graph_xc_ERER(avg_k) << endl;
//cout << graph_xc_SFSF(lembda, mm) << endl;
 // cout << graph_xc_ERSF(avg_k,lembda, mm) << endl;
 // cout << graph_xc_SFER(avg_k,lembda, mm) << endl;
  
int realiz;
for (realiz=0;realiz<realization;realiz++){i=0;cout << "projgaotest real = " << realiz << endl;
 for (p=0;p<=1.00000000001;p=p+0.01)
 {//cout<<"doing "<<p<<endl;
  net1.clear();net2.clear();
  net1=newnet1;net2=newnet2;

  rand_init_attack(p,net1);
  pinfinity=cascade_failure3(net1, net2, dep1 ,dep2);
  pout[i]=pout[i]+pinfinity;
  i=i+1;
  //fout << p <<"  "<< pinfinity << endl;
  }
                                              }
 i=0;
 for (p=0;p<=1.00000000001;p=p+0.01){
 fout << p << "  " << pout[i]/realization << endl;
 i=i+1;
                                }

}//end 1st digit
}//end 2nd digit

 
}
