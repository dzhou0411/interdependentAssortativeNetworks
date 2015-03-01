#ifndef NET_OPS_HPP
#define NET_OPS_HPP
#include <cmath>
#include <queue>
#include <stack>
#include <string>
#include "common.hpp"
#include "net_algo.hpp"
//------------------//
//1.
//------------------//
inline void find_biggest_cluster(Graph & net, vector<int> & biggest_cluster) {
  biggest_cluster.clear();
  int num_nodes=net.get_num_vertices();
  bool * marked;
  marked=new bool[num_nodes];
  for (int i=0; i < num_nodes; ++i) marked[i]=false;
  stack<int> search_stack;
  int idx_search;
  Graph::node_neighbor_iterator idx_neighbor;
  vector<int> temp_cluster;

  for (int idx_node=0; idx_node < num_nodes; ++idx_node) {
    if (!marked[idx_node]) {
      temp_cluster.clear();
      temp_cluster.push_back(idx_node);
      search_stack.push(idx_node);
      marked[idx_node]=true;
      
      while (!search_stack.empty()) {
	idx_search=search_stack.top();
	search_stack.pop();
	for (idx_neighbor=net.vertex_neighbor_begin(idx_search);
	     idx_neighbor!=net.vertex_neighbor_end(idx_search); ++idx_neighbor)
	  if (!marked[*idx_neighbor]) {
	    search_stack.push(*idx_neighbor);
	    temp_cluster.push_back(*idx_neighbor);
	    marked[*idx_neighbor]=true;
	  }
      }
      if (temp_cluster.size() > biggest_cluster.size())
	biggest_cluster=temp_cluster;
    }
  }
  delete []marked;
}


//---------------------------------------------------//
inline void find_biggest_and_second_cluster(Graph & net, vector<int> & biggest_cluster, vector<int> & second_biggest_cluster) {
  biggest_cluster.clear();
  second_biggest_cluster.clear();
  int num_nodes=net.get_num_vertices();
  bool * marked;
  marked=new bool[num_nodes];
  for (int i=0; i < num_nodes; ++i) marked[i]=false;
  stack<int> search_stack;
  int idx_search;
  Graph::node_neighbor_iterator idx_neighbor;
  vector<int> temp_cluster;

  for (int idx_node=0; idx_node < num_nodes; ++idx_node) {
    if (!marked[idx_node]) {
      temp_cluster.clear();
      temp_cluster.push_back(idx_node);
      search_stack.push(idx_node);
      marked[idx_node]=true;
      
      while (!search_stack.empty()) {
	idx_search=search_stack.top();
	search_stack.pop();
	for (idx_neighbor=net.vertex_neighbor_begin(idx_search);
	     idx_neighbor!=net.vertex_neighbor_end(idx_search); ++idx_neighbor)
	  if (!marked[*idx_neighbor]) {
	    search_stack.push(*idx_neighbor);
	    temp_cluster.push_back(*idx_neighbor);
	    marked[*idx_neighbor]=true;
	  } // iterate each neighbour end
      } //  while search_stack,empty() end
      if (temp_cluster.size() > biggest_cluster.size())
	{second_biggest_cluster=biggest_cluster;biggest_cluster=temp_cluster;}
         
    }// if !marked[idx_node] end
  } // for each node end
  delete []marked;

  

  



} // function end

//--------------------------------------------//
//2.Cascade Failures
//Given P, output the P_inf
//net_1 was attacked; net_2 is complete
//--------------------------------------------//
//2.1 record size of cluster in each step in temp.dat; 
//    output ulti-size of cluster
//    consider protection.
//--------------------------------------------//
inline double cascade_failure1(Graph &net_1, Graph &net_2, int real){
  double p_inf;
  int num_nodes=net_1.get_num_vertices();
  //cout<<num_nodes<<endl;
  vector<int> biggest_cluster1,biggest_cluster2;
  set<int> set_1,set_2;
  for(int i=0;i<num_nodes;++i){
    biggest_cluster2.push_back(i);
  }
  ofstream ofile;
  string file="../data/temp"+d_c(real)+".dat";
  ofile.open(file.c_str());
  //Start cascade process
  int id;int size;
  while(1){
    size=biggest_cluster1.size();
    biggest_cluster1.clear();set_1.clear();
    find_biggest_cluster(net_1,biggest_cluster1);
    ofile<<double(biggest_cluster1.size())/(num_nodes)<<endl;
    if(biggest_cluster1.size()==size){
      p_inf=double(biggest_cluster1.size())/double(num_nodes);
      return p_inf;
    }
    for(int i=0;i<biggest_cluster1.size();++i){
      set_1.insert(biggest_cluster1[i]);
    }
    for(int i=0;i<biggest_cluster2.size();++i){
      id=biggest_cluster2[i];
      if(!set_1.count(id)){
	net_1.rm_a_node(id);
	if(!net_2.protect_check(id))
	  net_2.rm_a_node(id);
      }
    }
    size=biggest_cluster2.size();
    biggest_cluster2.clear();set_2.clear(); 
    find_biggest_cluster(net_2,biggest_cluster2);
    ofile<<double(biggest_cluster2.size())/(num_nodes)<<endl;
    if(biggest_cluster2.size()==size){
      p_inf=double(biggest_cluster1.size())/double(num_nodes);
      return p_inf;
    }
    for(int i=0;i<biggest_cluster2.size();++i){
      set_2.insert(biggest_cluster2[i]);
    }
    for(int i=0;i<biggest_cluster1.size();++i){
      id=biggest_cluster1[i];
      if(!set_2.count(id)){
	net_2.rm_a_node(id);
	if(!net_1.protect_check(id))
	  net_1.rm_a_node(id);
      }
    }
  }
  ofile.close();
}
//--------------------------------//
//2.2 record the ids of those nodes which die in each step in v1 & v2
//    used to find the properties of nodes that dies in each step
//--------------------------------//
inline double cascade_failure2(Graph &net_1, Graph &net_2, int real, vector<int> &v1, vector<int> &v2){
  double p_inf;
  int num_nodes=net_1.get_num_vertices();
  vector<int> biggest_cluster1,biggest_cluster2,cluster1,cluster2;
  set<int> set_1,set_2;
  for(int i=0;i<num_nodes;++i){
    biggest_cluster2.push_back(i);
  }
  ofstream ofile;
  string file="../data/temp"+d_c(real)+".dat";
  ofile.open(file.c_str());
  //Start cascade process
  int id;int size;
  while(1){
    size=biggest_cluster1.size();
    cluster1.clear();cluster1=biggest_cluster1;
    biggest_cluster1.clear();set_1.clear();
    find_biggest_cluster(net_1,biggest_cluster1);
    ofile<<double(biggest_cluster1.size())/(num_nodes)<<endl;
    if(biggest_cluster1.size()==size){
      p_inf=double(biggest_cluster1.size())/double(num_nodes);
      return p_inf;
    }
    for(int i=0;i<biggest_cluster1.size();++i){
      set_1.insert(biggest_cluster1[i]);
    }
    for(int i=0;i<cluster1.size();++i)
      if(!set_1.count(cluster1[i]))
    	v1.push_back(cluster1[i]);
    v1.push_back(-1);
    for(int i=0;i<biggest_cluster2.size();++i){
      id=biggest_cluster2[i];
      if(!set_1.count(id)){
	net_1.rm_a_node(id);
	if(!net_2.protect_check(id))
	  net_2.rm_a_node(id);
      }
    }
    size=biggest_cluster2.size();
    cluster2.clear();cluster2=biggest_cluster2;
    biggest_cluster2.clear();set_2.clear(); 
    find_biggest_cluster(net_2,biggest_cluster2);
    ofile<<double(biggest_cluster2.size())/(num_nodes)<<endl;
    if(biggest_cluster2.size()==size){
      p_inf=double(biggest_cluster1.size())/double(num_nodes);
      return p_inf;
    }
    for(int i=0;i<biggest_cluster2.size();++i){
      set_2.insert(biggest_cluster2[i]);
    }
    for(int i=0;i<cluster2.size();++i)
      if(!set_2.count(cluster2[i]))
	v2.push_back(cluster2[i]);
    v2.push_back(-1);
    for(int i=0;i<biggest_cluster1.size();++i){
      id=biggest_cluster1[i];
      if(!set_2.count(id)){
	net_2.rm_a_node(id);
	if(!net_1.protect_check(id))
	  net_1.rm_a_node(id);
      }
    }
  }
  ofile.close();
}
//-----------------------------//
//2.3 cascade failure by arbitrary dependance defined between two nets
//-----------------------------//
 double cascade_failure3(Graph &net_1, Graph &net_2,vector<int> &dep1,vector<int> &dep2){
  double p_inf;
  int num_nodes=net_1.get_num_vertices();
  vector<int> biggest_cluster1,biggest_cluster2;
  set<int> set_1,set_2;
  for(int i=0;i<num_nodes;++i){
    biggest_cluster2.push_back(i);
  }
  ofstream ofile;
  string file="../data/temp.dat";
  ofile.open(file.c_str());
  //Start cascade process
  int id;int size;
  while(1){
    size=biggest_cluster1.size();
    biggest_cluster1.clear();set_1.clear();
    find_biggest_cluster(net_1,biggest_cluster1);
    //ofile<<double(biggest_cluster1.size())<<endl;
    if(biggest_cluster1.size()==size){
      p_inf=double(biggest_cluster1.size())/double(num_nodes);
      return p_inf;
    }
    for(int i=0;i<biggest_cluster1.size();++i){
      set_1.insert(biggest_cluster1[i]);
    }
    for(int i=0;i<biggest_cluster2.size();++i){
      id=biggest_cluster2[i]; if (dep2[id]>=0) {
      if(!set_1.count(dep2[id])){
	net_1.rm_a_node(dep2[id]);
	if(!net_2.protect_check(id))
	  net_2.rm_a_node(id);
      }                                        }
    }
    size=biggest_cluster2.size();
    biggest_cluster2.clear();set_2.clear(); 
    find_biggest_cluster(net_2,biggest_cluster2);
    //ofile<<double(biggest_cluster2.size())<<endl;
    if(biggest_cluster2.size()==size){
      p_inf=double(biggest_cluster1.size())/double(num_nodes);
      return p_inf;
    }
    for(int i=0;i<biggest_cluster2.size();++i){
      set_2.insert(biggest_cluster2[i]);
    }
    for(int i=0;i<biggest_cluster1.size();++i){
      id=biggest_cluster1[i]; if (dep1[id]>=0)    {
      if(!set_2.count(dep1[id])){
	if(!net_1.protect_check(id))
	  net_1.rm_a_node(id);
	net_2.rm_a_node(dep1[id]);
      }                                           }
    }
  }
  ofile.close();
}//*/
//------------Cascade of Failure with Number Of Iteration (NOI) recorded------------------
double cascade_failure3_NOI(Graph &net_1, Graph &net_2,vector<int> &dep1,vector<int> &dep2, int & NOI){
  NOI=0;
  double p_inf;
  int num_nodes=net_1.get_num_vertices();
  vector<int> biggest_cluster1,biggest_cluster2;
  set<int> set_1,set_2;
  for(int i=0;i<num_nodes;++i){
    biggest_cluster2.push_back(i);
  }
  ofstream ofile;
  string file="../data/temp.dat";
  ofile.open(file.c_str());
  //Start cascade process
  int id;int size;
  while(1){
    size=biggest_cluster1.size();
    biggest_cluster1.clear();set_1.clear();
    find_biggest_cluster(net_1,biggest_cluster1);
    //ofile<<double(biggest_cluster1.size())<<endl;
    if(biggest_cluster1.size()==size){
      p_inf=double(biggest_cluster1.size())/double(num_nodes);
      return p_inf;
    }
   NOI=NOI+1;
    for(int i=0;i<biggest_cluster1.size();++i){
      set_1.insert(biggest_cluster1[i]);
    }
    for(int i=0;i<biggest_cluster2.size();++i){
      id=biggest_cluster2[i]; if (dep2[id]>=0) {
      if(!set_1.count(dep2[id])){
	net_1.rm_a_node(dep2[id]);
	if(!net_2.protect_check(id))
	  net_2.rm_a_node(id);
      }                                        }
    }
    size=biggest_cluster2.size();
    biggest_cluster2.clear();set_2.clear(); 
    find_biggest_cluster(net_2,biggest_cluster2);
    //ofile<<double(biggest_cluster2.size())<<endl;
    if(biggest_cluster2.size()==size){
      p_inf=double(biggest_cluster1.size())/double(num_nodes);
      return p_inf;
    }
    NOI=NOI+1;
    for(int i=0;i<biggest_cluster2.size();++i){
      set_2.insert(biggest_cluster2[i]);
    }
    for(int i=0;i<biggest_cluster1.size();++i){
      id=biggest_cluster1[i]; if (dep1[id]>=0)    {
      if(!set_2.count(dep1[id])){
	if(!net_1.protect_check(id))
	  net_1.rm_a_node(id);
	net_2.rm_a_node(dep1[id]);
      }                                           }
    }
  }
  ofile.close();
}//*/

//----------------cascade failure3, with dep links, NOI, and second largest cluster recorded
double cascade_failure3_NOI_sec(Graph &net_1, Graph &net_2,vector<int> &dep1,vector<int> &dep2, int & NOI, double & second_p){
  second_p=0.;
  NOI=0;
  double p_inf;
  int num_nodes=net_1.get_num_vertices();
  vector<int> biggest_cluster1,biggest_cluster2;
  vector<int> second_cluster;
  set<int> set_1,set_2;
  for(int i=0;i<num_nodes;++i){
    biggest_cluster2.push_back(i);
  }
  ofstream ofile;
  string file="../data/temp.dat";
  ofile.open(file.c_str());
  //Start cascade process
  int id;int size;
  while(1){
    size=biggest_cluster1.size();
    biggest_cluster1.clear();set_1.clear();
    second_cluster.clear();
    find_biggest_and_second_cluster(net_1,biggest_cluster1, second_cluster);
    //ofile<<double(biggest_cluster1.size())<<endl;
    if(biggest_cluster1.size()==size){
      p_inf=double(biggest_cluster1.size())/double(num_nodes);
      second_p=double(second_cluster.size())/double(num_nodes);
      return p_inf;
    }
   NOI=NOI+1;
    for(int i=0;i<biggest_cluster1.size();++i){
      set_1.insert(biggest_cluster1[i]);
    }
    for(int i=0;i<biggest_cluster2.size();++i){
      id=biggest_cluster2[i]; if (dep2[id]>=0) {
      if(!set_1.count(dep2[id])){
	net_1.rm_a_node(dep2[id]);
	if(!net_2.protect_check(id))
	  net_2.rm_a_node(id);
      }                                        }
    }
    size=biggest_cluster2.size();
    biggest_cluster2.clear();set_2.clear();
    second_cluster.clear(); 
    find_biggest_and_second_cluster(net_2,biggest_cluster2, second_cluster);
    //ofile<<double(biggest_cluster2.size())<<endl;
    if(biggest_cluster2.size()==size){
      p_inf=double(biggest_cluster1.size())/double(num_nodes);
      second_p=double(second_cluster.size())/double(num_nodes);
      return p_inf;
    }
    NOI=NOI+1;
    for(int i=0;i<biggest_cluster2.size();++i){
      set_2.insert(biggest_cluster2[i]);
    }
    for(int i=0;i<biggest_cluster1.size();++i){
      id=biggest_cluster1[i]; if (dep1[id]>=0)    {
      if(!set_2.count(dep1[id])){
	if(!net_1.protect_check(id))
	  net_1.rm_a_node(id);
	net_2.rm_a_node(dep1[id]);
      }                                           }
    }
  }
  ofile.close();
}//*/


//-------------------------------//
// 3. Initial attacks
// p: percentage alive
//-------------------------------//
inline void rand_init_attack(double p, Graph &net){ 
  int num_nodes=net.get_num_vertices();
  set<int> set_1;
  int init_attack=int((1.-p)*num_nodes);
  set_1.clear();
  while(set_1.size()<init_attack){
    int node_idx=int(srand()*num_nodes);//index of node removed
    if(net.get_deg_vertex(node_idx)&&!net.protect_check(node_idx))
      net.rm_a_node(node_idx);
    set_1.insert(node_idx);
  }
  set_1.clear();
}
inline void max_degree_init_attack(double  p, Graph &net){
  int num_nodes=net.get_num_vertices();
  int init_attack=int((1.-p)*num_nodes);
  vector<pair<int,double> > v;
  for(int i=0;i<net.size();++i)
    v.push_back(make_pair(i,net.get_deg_vertex(i)));
  bubblesort2(v);
  for(int i=0;i<init_attack;i++){
    net.rm_a_node(v[i].first);
  }
}
//attack nodes proportional to degree
inline void intentional_attack(double p, Graph &net,vector<pair<int,int> > &links){
  set<int> set_1;
  int num_nodes=net.get_num_vertices();
  int init_attack=int((1.-p)*num_nodes);
  int node_idx;
  while(set_1.size()<init_attack){
    node_idx=(srand()>0.5)?links[int(srand()*links.size())].first:links[int(srand()*links.size())].second;
    if(net.get_deg_vertex(node_idx))
      net.rm_a_node(node_idx);
    set_1.insert(node_idx);
  }
}

//attack nodes proportional to k^alpha

//------------------------------//
//good for alpha very large
//------------------------------//
/*
void intentional_attack_alpha(double p_0,Graph &net,double alpha){
  list<pair<int,double> > l1;
  int size=net.size();
  double range=0,temp;
  double m;
  for(int i=0;i<size;++i){
    temp=pow(net.get_deg_vertex(i),alpha)/pow(4,alpha);
    l1.push_back(make_pair(i,temp));
    range+=temp;
  }
  cout<<':'<<l1.size()<<' '<<range<<endl;
  int num=0;
  int NUM=(1-p_0)*net.size();
  list<pair<int,double> >::iterator it;
  while(num<NUM){
    //for(list<pair<int,double> >::iterator it2=l1.begin();it2!=l1.end();++it2)
    //cout<<it2->first<<' '<<it2->second<<endl;
    m=srand()*range;
    it=l1.begin();
    temp=it->second;    
    while(temp<m){
      it++;
      temp+=it->second;
    }
    //cout<<num<<' '<<it->first<<' '<<it->second<<endl;
    range-=it->second;    
    net.rm_a_node(it->first);
    l1.erase(it);
    num++;
    
  }
}//*/

//------------------------------//
//good for alpha small or negative
//------------------------------//
//*
void intentional_attack_alpha(double p_0,Graph &net,double alpha){
  vector<pair<double,int> > v(net.size()+1);
  vector<int> v1(net.size()+1,1);
  v[0].first=0;
  for(int i=0;i<net.size();++i){
    v[i].second=i-1;
    if(net.get_deg_vertex(i)==0) v[i+1].first=v[i].first;
    else v[i+1].first=v[i].first+pow(net.get_deg_vertex(i),alpha);
  }
  int num=0;
  double m;int upper,lower;
  int b,c;
  int NUM=(1-p_0)*net.size();
  while(num<NUM){
    m=srand()*v[net.size()-num].first;
    b=int((net.size()-num)/2);
    upper=net.size()-num;lower=0;
    while(abs(upper-lower)>1){
      if(m>v[b].first){
	lower=b;
	b=int((b+upper+1)/2);
      }else if(m<v[b].first){
	upper=b;
	b=int((b+lower+1)/2);
      }else if(m==v[b].first){
	break;
      }
    }
    net.rm_a_node(v[b].second);
    num++;
    double delta=v[b].first-v[b-1].first;
    for(int i=0;i<net.size()-num-b;++i){
      v[b+i].first=v[b+i+1].first-delta;
      v[b+i].second=v[b+i+1].second;
    }
  }

}
//*/

/*
void intentional_attack_alpha(double p_0,Graph &net,double alpha){
  vector<double> v(net.size()+1);
  vector<int> v1(net.size()+1,1);
  v[0]=0;//pow(net.get_deg_vertex(0),alpha);
  for(int i=0;i<net.size();++i){
    if(net.get_deg_vertex(i)==0) v[i+1]=v[i];
    else v[i+1]=v[i]+pow(net.get_deg_vertex(i),alpha);
    //cout<<v[i]<<endl;
  }
  int num=0;
  double m,upper,lower;
  int b,c;
  while(num<(1-p_0)*net.size()){
    m=srand()*v[net.size()];
    b=int(net.size()/2);
    upper=net.size();lower=0;
    while(abs(upper-lower)>1){
      if(m>v[b]){
	lower=b;
	b=int((b+upper+1)/2);
      }else if(m<v[b]){
	upper=b;
	b=int((b+lower+1)/2);
      }else if(m==v[b]){
	break;
      }
    }
    if(v1[b]==1){
      num++;
      v1[b]=0;
      net.rm_a_node(b-1);
    }
  }
}
//*/

//-------------------------------//
//4. Protect Strategies
//q1,q2: fraction of nodes protected
//-------------------------------//
//4.1
//-------------------------------//
void random_protect(Graph &net_1,Graph &net_2,double q1,double q2){
  int index;
  net_1.unprotect_allnodes();
  int num_nodes=net_1.size();
  int num_protect=int(q1*num_nodes);
  for(int j=0;j<num_protect;){
    index=int(srand()*num_nodes);
    if(!net_1.protect_check(index)){
      net_1.protect_node(index);
      ++j;
    }
  }
  net_2.unprotect_allnodes();
  num_nodes=net_2.size();
  num_protect=int(q2*num_nodes);
  for(int j=0;j<num_protect;){
    index=int(srand()*num_nodes);
    if(!net_2.protect_check(index)){
      net_2.protect_node(index);
      ++j;
    }
  }
}
//--------------------------------//
//4.2 intentional protection 
// only protect net_2.
//select= 1.degree; 2.btwness; 3.aved
//type= 0.rank by absolute value; 1.difference
//--------------------------------//
void max_protect(int select,int type,Graph &net_1,Graph &net_2,double q){
  vector<vector<pair<int,double> > > v(2);//v1 is absolute value
  if(select==1){//degree
    for(int i=0;i<net_2.size();++i){
      v[0].push_back(make_pair(i,net_2.get_deg_vertex(i)));
      v[1].push_back(make_pair(i,net_2.get_deg_vertex(i)-net_1.get_deg_vertex(i)));
      }
  }
  else if(select==2){//max btw
    net_2.set_btwness();net_1.set_btwness();
    for(int i=0;i<net_2.size();++i){
      v[0].push_back(make_pair(i,net_2.btw(i)));
      v[1].push_back(make_pair(i,net_2.btw(i)-net_1.btw(i)));
    }
  }
  else if(select==3){//min aved
    net_2.set_aved();net_1.set_aved();
    for(int i=0;i<net_2.size();++i){
      v[0].push_back(make_pair(i,100-net_2.aved(i)));
      v[1].push_back(make_pair(i,net_1.aved(i)-net_2.aved(i)));
    }
  }
  bubblesort2(v[type]);
  cout<<"bubblesort finished with size: "<<v.size()<<endl;
  int num_nodes=net_2.size();
  int num_protect=int(q*num_nodes);
  net_2.unprotect_allnodes();
  for(int i=0;i<num_protect;i++) 
    net_2.protect_node(v[type][i].first);
}
//------------------------------------------//
//5. cascade failure simulation
//------------------------------------------//
//5.1 without protection
//------------------------------------------//
void percolation_simulation1(double alpha,int num_nodes,int avg_degree,double lamda,int realization,int type,int attack,double p_start,double p_end,double p_delta,int mm){
  string file,nu,nu1;
  ofstream ofile;
  Graph net_1(num_nodes),net_2(num_nodes),net_1_prime(num_nodes),net_2_prime(num_nodes);
  map<double,double> map_p;
  int index;
  vector<pair<int,int> > links;
  for(int i=0;i<realization;++i){
    cout<<i<<endl;
    net_1.clear();net_2.clear();
    if(type==1){
      lt_ER_algo2(num_nodes,net_1,avg_degree,links);
      lt_ER_algo(num_nodes,net_2,avg_degree);
      nu="ER";
      nu1="k_"+d_c(avg_degree);
    }
    else if(type==2){
      lt_SF_algo2(num_nodes,net_1,lamda,links,mm);
      lt_SF_algo(num_nodes,net_2,lamda,mm);
      nu="SF";
      nu1="lamda_"+d_c(int(lamda*10.0001));
    }
    for(double p=p_start;p<=p_end;p+=p_delta){
      double p_inf=0.;
      net_1_prime.clear();net_2_prime.clear();
      net_1_prime=net_1;net_2_prime=net_2;
      if(attack==1) rand_init_attack(p, net_1_prime);
      else if(attack==2) max_degree_init_attack(p, net_1_prime);
      else if(attack==3) intentional_attack(p,net_1_prime,links);
      else if(attack==4) intentional_attack_alpha(p,net_1_prime,alpha);
      cout<<i<<' '<<p<<endl;
      p_inf=cascade_failure1(net_1_prime,net_2_prime,i);
      cout<<p_inf<<endl;
      map_p[p]+=p_inf;
    }
  }
  file="../data/percolation/"+nu+"_noprotection_"+nu1+"_"+d_c(int(alpha))+"mm_"+d_c(mm)+".dat";
  ofile.open(file.c_str());
  for(double p=p_start;p<=p_end;p+=p_delta){
    ofile<<p<<' '<<map_p[p]/double(realization)<<' '<<map_p[p]/double(realization)/p<<endl;
  }
  ofile.close();
}
//------------------------------------------//
//5.2 double attack
//------------------------------------------//
void percolation_simulation2(double alpha,int num_nodes,int avg_degree,double lamda,int realization,int type,int attack,double p_start,double p_end,double p_delta,int mm){
  string file,nu,nu1;
  ofstream ofile;
  Graph net_1(num_nodes),net_2(num_nodes),net_1_prime(num_nodes),net_2_prime(num_nodes);
  map<double,double> map_p;
  int index;
  vector<pair<int,int> > links,links2;
  for(int i=0;i<realization;++i){
    cout<<i<<endl;
    net_1.clear();net_2.clear();
    if(type==1){
      lt_ER_algo2(num_nodes,net_1,avg_degree,links);
      lt_ER_algo2(num_nodes,net_2,avg_degree,links2);
      nu="ER_doubleattack_";
      nu1="k_"+d_c(avg_degree);
    }
    else if(type==2){
      lt_SF_algo2(num_nodes,net_1,lamda,links,mm);
      lt_SF_algo(num_nodes,net_2,lamda,mm);
      nu="SF_doubleattack_";
      nu1="lamda_"+d_c(int(lamda*10.001));
    }
    for(double p=p_start;p<=p_end;p+=p_delta){
      double p_inf=0.;
      net_1_prime.clear();net_2_prime.clear();
      net_1_prime=net_1;net_2_prime=net_2;
      if(attack==1) rand_init_attack(p, net_1_prime);
      else if(attack==2) max_degree_init_attack(p, net_1_prime);
      else if(attack==3) {
	intentional_attack(p,net_1_prime,links);
	intentional_attack(p,net_2_prime,links2);
      }
      else if(attack==4) intentional_attack_alpha(p,net_1_prime,alpha);
      cout<<i<<' '<<p<<endl;
      p_inf=cascade_failure1(net_1_prime,net_2_prime,i);
      cout<<p_inf<<endl;
      map_p[p]+=p_inf;
    }
  }
  file="../data/percolation/"+nu+"_noprotection_"+nu1+"_alpha_"+d_c(int(alpha*1.001))+".dat";
  ofile.open(file.c_str());
  for(double p=p_start;p<=p_end;p+=p_delta){
    ofile<<p<<' '<<map_p[p]/double(realization)<<' '<<map_p[p]/double(realization)/p<<endl;
  }
  ofile.close(); 
}
//---------------------------//
//5.3.step simulation
//---------------------------//
void step_simu(double alpha,double p_0,int num_nodes,int type,int attack,double lamda,int avg_degree,int realization,int mm){
  string file,nu;
  
  switch (attack){
  case 1:
    nu="rand";
    break;
  case 2:
    nu="maxdeg";
    break;
  case 3:
    nu="intent";
    break;
  case 4:
    nu="intent_alpha";
  }
  switch (type){
  case(1):
    file="../data/step/ER_k"+d_c(avg_degree)+"_p0"+d_c(int(p_0*100.001))+"_sim_"+nu+".dat";
    break;
  case(2):
    file="../data/step/SF_lamda"+d_c(lamda*10.001)+"_p0"+d_c(int(p_0*100.001))+"_sim_"+nu+".dat";
    break;
  }
  percolation_simulation1(alpha,num_nodes,avg_degree,lamda,realization,type,attack,p_0,p_0+0.00001,0.1,mm);
  vector<double> v(100,0);
  ofstream ofile(file.c_str());
  ifstream ifile;
  double a;
  for(int i=0;i<realization;++i){
    file="../data/temp"+d_c(i)+".dat";
    ifile.open(file.c_str());
    int j=0;
    while(!ifile.eof()){
      ifile>>a;
      v[j++]+=a/realization;
      //cout<<i<<' '<<v[j-1]<<endl;
    }
    while(j<100)
      v[j++]+=a/realization;
    ifile.close();
  }
  for(int i=0;i<100;++i){
    if(v[i-1]-v[i]>0.0001)
      ofile<<i<<' '<<v[i]<<endl;
  }
  ofile.close();
}
//---------------------------//
//6. pdf of k after removed nodes but not edges
//---------------------------//
void pdf_k_intent_alpha(double p_0,int num_nodes,double avg_k,double lamda,double alpha,int M){
  Graph net(num_nodes);
  net.clear();
  //lt_ER_algo(num_nodes,net,avg_k);
  lt_SF_algo(num_nodes,net,lamda,M);
  vector<double> v(num_nodes+1);
  vector<int> v1(num_nodes+1,1);
  v[0]=0;//pow(net.get_deg_vertex(0),alpha);
  for(int i=0;i<num_nodes;++i){
    v[i+1]=v[i]+pow(net.get_deg_vertex(i),alpha);
  }
  int num=0;
  double m,upper,lower;
  int b,c;
  //cout<<(1-p_0)*num_nodes<<endl;
  while(num<(1-p_0)*num_nodes){
    m=srand()*v[num_nodes];
    b=int(num_nodes/2);
    upper=num_nodes;lower=0;
    while(abs(upper-lower)>1){
      if(m>v[b]){
	lower=b;
	b=int((b+upper+1)/2);
      }else if(m<v[b]){
	upper=b;
	b=int((b+lower+1)/2);
      }else if(m==v[b]){
	break;
      }
    }
    if(v1[b]==1){
      num++;
      v1[b]=0;
    }
  }
  
  map<double,double> a;
  for(int i=1;i<num_nodes+1;++i)
    if(v1[i]==1)
      a[net.get_deg_vertex(i-1)]++;
  string file;ofstream ofile;
  file="../data/temp.dat";
  ofile.open(file.c_str());
  int w=0,w1=0;
  for(map<double,double>::iterator it=a.begin();it!=a.end();++it){
    ofile<<it->first<<' '<<it->second/double(p_0*num_nodes)<<endl;
    w+=pow(it->first,alpha)*it->second;
    w1+=it->first*it->second;
  }
  //cout<<"Simulation Avg_deg: "<<w/double(p_0*num_nodes)<<endl;
  //cout<<w1/double(p_0*num_nodes)<<endl;
  ofile.close();
}

//7. nodes' degree uncorrelated after intentional attack
/*
void uncorrelated_proof(int num_nodes,double alpha,int avg_degree){
  
  map<int,double> N;
  map<pair<int,int>,double> N_link;
  vector<int> v;
  double t_n_links=0;
  Graph net1(num_nodes);
  
  
  double p=0.5;
  lt_ER_algo(num_nodes,net1,avg_degree);
    cout<<"good"<<endl;
  for(int i=0;i<num_nodes;++i){
    cout<<i<<' '<<net1.get_deg_vertex(i)<<endl;
    v.clear();
    N[net1.get_deg_vertex(i)]+=1;
    
    net1.show_neighbors(i,v);
    //cout<<"bad"<<endl;
    t_n_links+=v.size();
    cout<<"neighbors: ";
    for(int j=0;j<v.size();++j){
      cout<<v[j]<<' ';
      int a=min(net1.get_deg_vertex(v[j]),net1.get_deg_vertex(i));
      int b=max(net1.get_deg_vertex(v[j]),net1.get_deg_vertex(i));
      N_link[make_pair(a,b)]+=1;
    }
    cout<<endl;
  }
  cout<<"numnodes and links: "<<num_nodes<<' '<< t_n_links<<endl;
  cout<<"good"<<endl;
  for(map<int,double>::iterator it=N.begin();it!=N.end();++it){
    cout<<it->first<<' '<<it->second<<endl;
    it->second/=double(num_nodes);
  }
  cout<<"good"<<endl;
  for(map<pair<int,int>,double>::iterator it=N_link.begin();it!=N_link.end();++it){
    //it->second/=t_n_links;
    cout<<it->first.first<<' '<<N[it->first.first]<<' '
	<<' '<<it->first.second<<' '<<N[it->first.second]<<' '
	<<' '<<it->second<<endl;
    it->second-=N[it->first.first]*N[it->first.second];
    //cout<<it->first.first<<' '<<it->first.second<<' '<<it->second<<endl;
  }
  
}
*/

void deg_deg_correlation(int num_nodes, int avg_degree,double alpha){
  Graph net(num_nodes),net1(num_nodes);
  vector<pair<int,int> > links;
  vector<pair<double,double> > v;
  vector<int> vv;
  int R=100;

  double before=0;
  vector<double> after(9,0);
  for(int realization=0;realization<R;++realization){
  cout<<realization<<endl;
    net.clear();net.resize(num_nodes);
    cout<<"111"<<endl;
    net1.clear();net1.resize(num_nodes);
    cout<<"111"<<endl;
    links.clear();
    cout<<"111"<<endl;
    v.clear();
    
    lt_ER_algo2(num_nodes, net, avg_degree, links);
    cout<<"222"<<endl;
    for(int i=0;i<links.size();++i)
      v.push_back(make_pair(net.get_deg_vertex(links[i].first),net.get_deg_vertex(links[i].second)));  
    before+=correlation(v);
    
    cout<<"good"<<endl;
    double p=0.1;
    for(int j=0;j<9;j++){
      p=0.1+j*0.1;
      net1=net;
      intentional_attack_alpha(p,net1,alpha);
      double t_n_links=0;
      v.clear();
      for(int i=0;i<num_nodes;++i){   
	vv.clear();
	net1.show_neighbors(i,vv);
	t_n_links+=vv.size();
	for(int j=0;j<vv.size();++j){
	  int a=net1.get_deg_vertex(vv[j]);
	  int b=net1.get_deg_vertex(i);
	  v.push_back(make_pair(a,b));
	}
      }   
      after[j]+=correlation(v);
    }
    cout<<"one real finish"<<endl;
  }
  string file="../data/deg_deg_correlation_ER_"+d_c(int(alpha))+".dat";
  ofstream ofile(file.c_str());
  cout<<"good"<<endl;
  for(int i=0;i<9;i++)
    ofile<<after[i]/double(R)<<endl;
  ofile<<before/double(R)<<endl;
}

void deg_deg_correlation2(int num_nodes, double lambda,double alpha, int M){
  Graph net(num_nodes),net1(num_nodes);
  vector<pair<int,int> > links;
  vector<pair<double,double> > v;
  vector<int> vv;
  int R=100;

  double before=0;
  vector<double> after(9,0);
  for(int realization=0;realization<R;++realization){
  cout<<realization<<endl;
    net.clear();net.resize(num_nodes);
    cout<<"111"<<endl;
    net1.clear();net1.resize(num_nodes);
    cout<<"111"<<endl;
    links.clear();
    cout<<"111"<<endl;
    v.clear();
    
    lt_SF_algo2(num_nodes, net, lambda, links, M);
    cout<<"222"<<endl;
    for(int i=0;i<links.size();++i)
      v.push_back(make_pair(net.get_deg_vertex(links[i].first),net.get_deg_vertex(links[i].second)));  
    before+=correlation(v);
    
    cout<<"good"<<endl;
    double p=0.1;
    for(int j=0;j<9;j++){
      p=0.1+j*0.1;
      net1=net;
      intentional_attack_alpha(p,net1,alpha);
      double t_n_links=0;
      v.clear();
      for(int i=0;i<num_nodes;++i){   
	vv.clear();
	net1.show_neighbors(i,vv);
	t_n_links+=vv.size();
	for(int j=0;j<vv.size();++j){
	  int a=net1.get_deg_vertex(vv[j]);
	  int b=net1.get_deg_vertex(i);
	  v.push_back(make_pair(a,b));
	}
      }   
      after[j]+=correlation(v);
    }
    cout<<"one real finish"<<endl;
  }
  string file="../data/deg_deg_correlation_SF_"+d_c(int(lambda*10.0001))+"_"+d_c(int(alpha))+".dat";
  ofstream ofile(file.c_str());
  cout<<"good"<<endl;
  for(int i=0;i<9;i++)
    ofile<<after[i]/double(R)<<endl;
  ofile<<before/double(R)<<endl;
}

/*
void random_protect(int num_nodes,int avg_degree,double lamda, int realization,int type,int attack,double q1_start,double q1_end,double q1_delta,double q2_start,double q2_end,double q2_delta,double p_start,double p_end,double p_delta){
  string file,nu;
  ofstream ofile;
  double p,p_inf;
  Graph net_1_prime(num_nodes),net_2_prime(num_nodes);
  Graph net_1(num_nodes),net_2(num_nodes);
  map<pair<pair<double,double>,double>,double> map_p;
  //map<double,double> map_p;
  int index;
  for(int i=0;i<realization;++i){
    net_1.clear();net_2.clear();
    if(type==1){
      lt_ER_algo(num_nodes,net_1,avg_degree);
      lt_ER_algo(num_nodes,net_2,avg_degree);
      nu="ER";
    }
    else if(type==2){
      lt_SF_algo(num_nodes,net_1,lamda);
      lt_SF_algo(num_nodes,net_2,lamda);
      nu="SF";
    }
    for(double f1=q1_start;f1<=q1_end;f1+=q1_delta){
      //cout<<f1<<endl;
      net_1.unprotect_allnodes();
      int num_protect=int(f1*net_1.size());
      for(int j=0;j<num_protect;){
	index=int(srand()*num_nodes);
	if(!net_1.protect_check(index)){
	  net_1.protect_node(index);
	  ++j;
	}
      }
      for(double f=q2_start;f<=q2_end;f+=q2_delta){
	net_2.unprotect_allnodes();
	num_protect=int(f*net_2.size());
	for(int j=0;j<num_protect;){
	  index=int(srand()*num_nodes);
	  if(!net_2.protect_check(index)){
	    net_2.protect_node(index);
	    ++j;
	  }
	}
	for(p=p_start;p<=p_end;p+=p_delta){
	  p_inf=0.;
	  net_1_prime.clear();net_2_prime.clear();
	  net_1_prime=net_1;net_2_prime=net_2;
	  if(attack==1) rand_init_attack(p, net_1_prime);
	  else if(attack==2) max_degree_init_attack(p, net_1_prime);
	  p_inf=cascade_failure1(net_1_prime,net_2_prime,i);
	  map_p[make_pair(make_pair(f1,f),p)]+=p_inf;
	}
      }
    }
  }
  
  for(double f1=q1_start;f1<=q1_end;f1+=q1_delta){
    for(double f=q2_start;f<=q2_end;f+=q2_delta){
      //double f=f1;  
      file="../data/percolation/"+nu+"_rand_q1_"+d_c(int(f1*1000.001))+"_q2_"+d_c(int(f*1000.001))+".dat";
      ofile.open(file.c_str());
      for(p=p_start;p<=p_end;p+=p_delta){
	ofile<<p<<' '<<map_p[make_pair(make_pair(f1,f),p)]/double(realization)<<' '<<map_p[make_pair(make_pair(f1,f),p)]/double(realization)/p<<endl;
      }
      ofile.close();
    }
  }//
}*/
/*
//-----------------------------//
//max protection:
//1.degree; 2.btwness; 3.aved
//-----------------------------//
void max_protect(int select,int num_nodes,int avg_degree,int realization,int type,double f_start,double f_end,double f_delta,double p_start,double p_end,double p_delta){
  string file,nu;
  ofstream ofile;
  double p,p_inf;
  Graph net_1_prime(num_nodes),net_2_prime(num_nodes);
  Graph net_1(num_nodes),net_2(num_nodes);
  map<pair<double,double>,double> map_p,map_p2;
  vector<pair<int,double> > v,v2;
  for(int i=0;i<realization;++i){
    cout<<"realization: "<<i<<endl;
    net_1.clear();net_2.clear();
    if(type==1){
      lt_ER_algo(num_nodes,net_1,avg_degree);
      lt_ER_algo(num_nodes,net_2,avg_degree);
      nu="ER";
    }
    v.clear();
    if(select==1){//degree
      for(int i=0;i<net_2.size();++i){
	v.push_back(make_pair(i,net_2.get_deg_vertex(i)));
	v2.push_back(make_pair(i,net_2.get_deg_vertex(i)-net_1.get_deg_vertex(i)));
      }
      nu+="_deg";
    }
    else if(select==2){//btw
      net_2.set_btwness();net_1.set_btwness();
      for(int i=0;i<net_2.size();++i){
	v.push_back(make_pair(i,net_2.btw(i)));
	v2.push_back(make_pair(i,net_2.btw(i)-net_1.btw(i)));
      }
      nu+="_btw";
    }
    else if(select==3){//aved
      net_2.set_aved();net_1.set_aved();
      for(int i=0;i<net_2.size();++i){
	v.push_back(make_pair(i,100-net_2.aved(i)));
	v2.push_back(make_pair(i,net_1.aved(i)-net_2.aved(i)));
      }
      nu+="_aved";
    }
    bubblesort2(v);
    bubblesort2(v2);
    cout<<"bubblesort good "<<v.size()<<endl;
    for(double f=f_start;f<=f_end;f+=f_delta){
      int num_protect=int(f*net_2.size());
      net_2.unprotect_allnodes();
      for(int i=0;i<num_protect;++i)
	net_2.protect_node(v[i].first);
      for(p=p_start;p<=p_end;p+=p_delta){
	p_inf=0.;
	net_1_prime.clear();net_2_prime.clear();
	net_1_prime=net_1;net_2_prime=net_2;
	rand_init_attack(p, net_1_prime);
	p_inf=cascade_failure1(net_1_prime,net_2_prime,1);
	map_p[make_pair(f,p)]+=p_inf;
      }
    }
    for(double f=f_start;f<=f_end;f+=f_delta){
      int num_protect=int(f*net_2.size());
      net_2.unprotect_allnodes();
      for(int i=0;i<num_protect;++i)
	net_2.protect_node(v2[i].first);
      for(p=p_start;p<=p_end;p+=p_delta){
	p_inf=0.;
	net_1_prime.clear();net_2_prime.clear();
	net_1_prime=net_1;net_2_prime=net_2;
	rand_init_attack(p, net_1_prime);
	p_inf=cascade_failure1(net_1_prime,net_2_prime,1);
	map_p2[make_pair(f,p)]+=p_inf;
      }
    }
  }
  for(double f=f_start;f<=f_end;f+=f_delta){
    file="../data/percolation/"+nu+"_"+d_c(int(f*1000.001))+".dat";
    ofile.open(file.c_str());
    for(p=p_start;p<=p_end;p+=p_delta){
      ofile<<p<<' '<<map_p[make_pair(f,p)]/double(realization)<<' '<<map_p[make_pair(f,p)]/double(realization)/p<<endl;
    }
    ofile.close();
  }
  for(double f=f_start;f<=f_end;f+=f_delta){
    file="../data/percolation/"+nu+"diff_"+d_c(int(f*1000.001))+".dat";
    ofile.open(file.c_str());
    for(p=p_start;p<=p_end;p+=p_delta){
      ofile<<p<<' '<<map_p2[make_pair(f,p)]/double(realization)<<' '<<map_p2[make_pair(f,p)]/double(realization)/p<<endl;
    }
    ofile.close();
  }
}
//rank by degree difference between two nodes. (k_2-k_1)
/*void max_degdiff_protect(int num_nodes,int avg_degree,int realization,int type,double f_start,double f_end,double f_delta,double p_start,double p_end,double p_delta){
  string file,nu;
  ofstream ofile;
  double p,p_inf;
  Graph net_1_prime(num_nodes),net_2_prime(num_nodes);
  Graph net_1(num_nodes),net_2(num_nodes);
  map<pair<double,double>,double> map_p;
  vector<pair<int,double> > v;
  vector<int> v1,v2;
  for(int i=0;i<realization;++i){
    net_1.clear();net_2.clear();
    if(type==1){
      lt_ER_algo(num_nodes,net_1,avg_degree);
      lt_ER_algo(num_nodes,net_2,avg_degree);
      nu="ER";
    }
    v.clear();
    dependant_relation(net_1,net_2,v1,v2);
    for(int i=0;i<net_2.size();++i){
      v.push_back(make_pair(i,net_2.get_deg_vertex(i)-net_1.get_deg_vertex(i)));
    }
    bubblesort2(v);
    for(double f=f_start;f<=f_end;f+=f_delta){
      int num_protect=int(f*net_2.size());
      for(int i=0;i<num_protect;++i){
	net_2.protect_node(v[i].first);
      }
      for(p=p_start;p<=p_end;p+=p_delta){
	p_inf=0.;
	net_1_prime.clear();net_2_prime.clear();
	net_1_prime=net_1;net_2_prime=net_2;
	rand_init_attack(p, net_1_prime);
	p_inf=cascade_failure3(net_1_prime,net_2_prime,v1,v2);
	map_p[make_pair(f,p)]+=p_inf;
      }
    }
  }
  for(double f=f_start;f<=f_end;f+=f_delta){
    file="../data/percolation/"+nu+"_degdiff_"+d_c(int(f*1000))+".dat";
    ofile.open(file.c_str());
    for(p=p_start;p<=p_end;p+=p_delta){
      ofile<<p<<' '<<map_p[make_pair(f,p)]/double(realization)<<' '<<map_p[make_pair(f,p)]/double(realization)/p<<endl;
    }
    ofile.close();
  }
}
void max_degnear_protect(int num_nodes,int avg_degree,int realization,int type,double f_start,double f_end,double f_delta,double p_start,double p_end,double p_delta){
  string file,nu;
  ofstream ofile;
  double p,p_inf;
  Graph net_1_prime(num_nodes),net_2_prime(num_nodes);
  Graph net_1(num_nodes),net_2(num_nodes);
  map<pair<double,double>,double> map_p;
  vector<pair<int,double> > v;
  vector<int> v1,v2;
  for(int i=0;i<realization;++i){
    net_1.clear();net_2.clear();
    if(type==1){
      lt_ER_algo(num_nodes,net_1,avg_degree);
      lt_ER_algo(num_nodes,net_2,avg_degree);
      nu="ER";
    }
    v.clear(); 
    dependant_relation(net_1,net_2,v1,v2);
    list<int>::iterator it;
    for(int i=0;i<net_2.size();++i){
      double a=0;
      for(it=net_2.vertex_neighbor_begin(i);it!=net_2.vertex_neighbor_end(i);++it){
	a+=net_2.get_deg_vertex(*it);
      }
      v.push_back(make_pair(i,a));
    }
    bubblesort2(v);
    //for(int i=0;i<v.size()/2;++i)
    //swap(v[i],v[v.size()-i]);
    for(double f=f_start;f<=f_end;f+=f_delta){
      int num_protect=int(f*net_2.size());
      for(int i=0;i<num_protect;++i){
	net_2.protect_node(v[i].first);
      }
      for(p=p_start;p<=p_end;p+=p_delta){
	p_inf=0.;
	net_1_prime.clear();net_2_prime.clear();
	net_1_prime=net_1;net_2_prime=net_2;
	rand_init_attack(p, net_1_prime);
	p_inf=cascade_failure3(net_1_prime,net_2_prime,v1,v2);
	map_p[make_pair(f,p)]+=p_inf;
      }
    }
  }
  for(double f=f_start;f<=f_end;f+=f_delta){
    file="../data/percolation/"+nu+"_degnear_"+d_c(int(f*1000))+".dat";
    ofile.open(file.c_str());
    for(p=p_start;p<=p_end;p+=p_delta){
      ofile<<p<<' '<<map_p[make_pair(f,p)]/double(realization)<<' '<<map_p[make_pair(f,p)]/double(realization)/p<<endl;
    }
    ofile.close();
  }
  }*/

/*
void max_aved_protect(int num_nodes,int avg_degree,int realization,int type,double f_start,double f_end,double f_delta,double p_start,double p_end,double p_delta){
  for(int i=0;i<net_2.size();++i)
    net_2.unprotect_node(i);
  double p,p_inf;
  Graph net_1_prime(net_1.size()),net_2_prime(net_2.size());
  vector<pair<int,double> > v1,v2;
  net_1.set_aved();net_2.set_aved();
  for(int i=0;i<net_2.size();++i){
    v1.push_back(make_pair(i,net_2.aved(i)));
    v2.push_back(make_pair(i,net_2.aved(i)-net_1.aved(i)));
  }
  bubblesort2(v1);
  bubblesort2(v2);
  cout<<"sort finished!"<<endl;
  string file;
  ofstream ofile;
  for(double f=start;f<end;f+=delta){
    int num_protect=int(f*net_2.size());
    //---------------aved protect
    for(int i=0;i<num_protect;++i){
      net_2.protect_node(v1[i].first);
    }
    for(int i=num_protect;i<v1.size();++i)
      net_2.unprotect_node(v1[i].first);
    file="../data/percolation/avedprotect_"+d_c(int(f*100))+"_"+nu+".dat";
    ofile.open(file.c_str());
    for(p=0.1;p<0.9;p+=0.05){
      p_inf=0.;
      for(int i=0;i<realization;++i){
	net_1_prime.clear();net_2_prime.clear();
	net_1_prime=net_1;net_2_prime=net_2;
	rand_init_attack(p, net_1_prime);
	p_inf+=cascade_failure3(net_1_prime,net_2_prime,v3,v4);
      }
      ofile<<p<<' '<<p_inf/double(realization)<<' '<<p_inf/double(realization)/p<<endl;
    }
    ofile.close();
    //--------------aved difference protect
    for(int i=0;i<num_protect;++i){
      net_2.protect_node(v2[i].first);
    }
    for(int i=num_protect;i<v2.size();++i)
      net_2.unprotect_node(v2[i].first);
    file="../data/percolation/aveddiffprotect_"+d_c(int(f*100))+"_"+nu+".dat";
    ofile.open(file.c_str());
    for(p=0.1;p<0.9;p+=0.05){
      p_inf=0.;
      for(int i=0;i<realization;++i){
	net_1_prime.clear();net_2_prime.clear();
	net_1_prime=net_1;net_2_prime=net_2;
	rand_init_attack(p, net_1_prime);
	p_inf+=cascade_failure(net_1_prime,net_2_prime);
      }
      ofile<<p<<' '<<p_inf/double(realization)<<' '<<p_inf/double(realization)/p<<endl;
    }
    ofile.close();
  }
}


















  //----------------
/*  for(int i=0;i<net_2.size();++i)
    net_2.unprotect_node(i);
  double p,p_inf;
  Graph net_1_prime(net_1.size()),net_2_prime(net_2.size());
  vector<pair<int,double> > v1,v2;
  net_1.set_btwness();net_2.set_btwness();
  for(int i=0;i<net_2.size();++i){
    v1.push_back(make_pair(i,net_2.btw(i)));
    //v2.push_back(make_pair(i,net_2.btw(i)-net_1.btw(i)));
  }
  bubblesort2(v1);
  //bubblesort2(v2);
  cout<<"sort finished!"<<endl;
  string file;
  ofstream ofile;
  for(double f=0.00;f<0.2;f+=0.04){
    int num_protect=int(f*net_2.size());
    //---------------btw protect
    for(int i=0;i<num_protect;++i){
      net_2.protect_node(v1[i].first);
    }
    for(int i=num_protect;i<v1.size();++i)
      net_2.unprotect_node(v1[i].first);
    file="../data/percolation/ER_btw_"+d_c(int(f*100))+"_"+nu+".dat";
    ofile.open(file.c_str());
    for(p=0.1;p<0.9;p+=0.05){
      p_inf=0.;
      for(int i=0;i<realization;++i){
	net_1_prime.clear();net_2_prime.clear();
	net_1_prime=net_1;net_2_prime=net_2;
	rand_init_attack(p, net_1_prime);
	p_inf+=cascade_failure(net_1_prime,net_2_prime);
      }
      ofile<<p<<' '<<p_inf/double(realization)<<' '<<p_inf/double(realization)/p<<endl;
    }
    ofile.close();
    //--------------aved difference protect
    for(int i=0;i<num_protect;++i){
      net_2.protect_node(v2[i].first);
    }
    for(int i=num_protect;i<v2.size();++i)
      net_2.unprotect_node(v2[i].first);
    file="../data/percolation/btwdiffprotect_"+d_c(int(f*100))+"_"+nu+".dat";
    ofile.open(file.c_str());
    for(p=0.1;p<0.9;p+=0.05){
      p_inf=0.;
      for(int i=0;i<realization;++i){
	net_1_prime.clear();net_2_prime.clear();
	net_1_prime=net_1;net_2_prime=net_2;
	rand_init_attack(p, net_1_prime);
	p_inf+=cascade_failure(net_1_prime,net_2_prime);
      }
      ofile<<p<<' '<<p_inf/double(realization)<<' '<<p_inf/double(realization)/p<<endl;
    }
    ofile.close();
  }
}


*/



/*
//--------------------------------------//
//Theoretical cascade step by step
//--------------------------------------//
void ER_step(int k,double p_0){
  string file="../data/step/ER_theo_"+d_c(k)+".dat";
  ofstream ofile(file.c_str());
  double p1,p2,p3=2.;
  p1=p_0*ER_g_A(p_0,double(k));
  p2=p1;
  while(abs(p1-p3)>0.000001){
    p3=p1;
    p1=p2*ER_g_A(p2,double(k));
    p2=p_0*ER_g_A(p2,double(k));
    ofile<<p3/p_0<<endl;
    cout<<p3<<endl;
  }
  ofile.close();
}
void SF_step(double lamda, double p_0){
  string file="../data/step/SF_theo_"+d_c(int(lamda*10.001))+".dat";
  ofstream ofile(file.c_str());
  double p1,p2,p3=2.;
  p1=p_0*SF_g_A(lamda,p_0);
  p2=p1;
  while(abs(p1-p3)>0.0001){
    p3=p1;
    p1=p2*SF_g_A(lamda,p2);
    p2=p_0*SF_g_A(lamda,p2);
    ofile<<p3<<endl;
    cout<<p1<<' '<<p2<<' '<<p3<<endl;
  }
  ofile.close();
}
void SF_maxdeg_attack_step(double lamda, double p_0){
  string file="../data/step/SF_maxdeg_theo_"+d_c(int(lamda*10.001))+"_"+d_c(int(p_0*100.001))+".dat";
  ofstream ofile(file.c_str());
  double m=2;
  double K=m*pow(p_0,1./(1-lamda));
  double pp=1-pow(1-p_0,(2.-lamda)/(1-lamda));//equivalent initial attack
  double p1,p2,p3,p4,p5=2.;

  p1=pp*SF_g_B(lamda,K,m,pp);
  p2=p1;
  while(abs(p1-p5)>0.0001){
    p5=p1;
    p3=p2*SF_g_A(lamda,p2);
    p4=pp*SF_g_A(lamda,p2);
    p1=p2*SF_g_B(lamda,K,m,p4);
    p2=pp*SF_g_B(lamda,K,m,p4);
    ofile<<p3<<endl;
    cout<<p1<<' '<<p2<<' '<<p3<<endl;
  }
  ofile.close();
}
//*/
/*
void step_simulation(int type,double lamda,int avg_degree,int realization){
  string file;
  if(type==1) file="../data/step/ER_simu_maxdeg_"+d_c(avg_degree)+".dat";
  else if(type==2) file="../data/step/SF_simu_maxdeg_"+d_c(int(lamda*10.001))+".dat";
  ofstream ofile(file.c_str());
  vector<double> v(50);
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
}//*/


/*
//------------------------------//
// mu_inf=p[1-exp(-k*mu_inf)]^2
//------------------------------//
double mu_inf(double avg_deg,double p){
  double mu_max=1.0,mu_min=0.;
  double mu=0.5;
  double f=p*pow((1-exp(0-avg_deg*mu)),2.)-mu;
  while(abs(f)>0.0001){
    if(f>0)
      mu_max=mu;
    else
      mu_min=mu;
    mu=(mu_max+mu_min)/2.;
    f=p*pow((1-exp(0-avg_deg*mu)),2.)-mu;
  }
  return mu;
}
//-------------------------------//
double P(double lamda,double k){
  return pow((k/2.),1-lamda)-pow((k+1.)/2.,1-lamda);
}
double pp(double lamda, double p){
  1-pow(p,(2.-lamda)/(1-lamda));
}
double PP(double lamda, double K, double p, double k){
  double a=0.,b,c,d=1;
  for(double i=k;i>0;i--)
    d*=i;
  for( double k_0=k;k_0<=K;++k_0){
    b=1.;c=1.;
    for(double i=k_0;i>0;i--)
      b*=i;
    for(double i=k_0-k;i>0;i--)
      c*=i;
    a+=P(lamda,k_0)*b/c/d*pow(1-pp(lamda,p),k)*pow(pp(lamda,p),k_0-k);
  }
  return a;
}
double u(double lamda){
  double avg_k=0;
  double a,b;
  for(int k=2;k<90000;++k){
    avg_k+=k*P(lamda,k);
  }
  double up=2,down=0.01;
  a=(up+down)/2.;
  double delta=2.;
  while(abs(delta)>0.00001){
    cout<<delta<<endl;
    cout<<a<<endl;
    b=0;
    for(int k=2;k<100000;++k)
      b+=k*P(lamda,k)*pow(a,k-1);
    delta=b-avg_k*a;
    if(delta>0){
      up=a;
      a=(a+down)/2.;
    }else if(delta<0){
      down=b;
      a=(a+up)/2.;
    }
  }
  cout<<"a"<<a<<endl;
  return a;
}
    
double P_inf(double lamda, double p){
  double K=2*pow(p,1./(1.-lamda));
  double uu=1;//=u(lamda);
  cout<<"uu: "<<uu<<endl;
  double a=0;
  for(int k=2;k<1000000;++k){
    a+=PP(lamda,K,p,k)*pow(uu,k-1);
  }
  return 1-a;
}


double ff(double lamda,double p){
  double avg_k=0;
  double K=2*pow(p,1./(1.-lamda));
  for(int k=2;k<1000000;++k){
    avg_k+=k*P(lamda,k);
  }
  string file="../data/temp.dat";
  ofstream ofile(file.c_str());
  for(double u=0;u<1.01;u+=0.02){
    double a=0;
    for(int k=2;k<1000000;++k){
      a+=k*PP(lamda,K,p,k)*pow(u,k);
    }
    a=a/avg_k;
    cout<<u<<' '<<a<<endl;
    ofile<<u<<' '<<a<<endl;
  }
}
double f(double lamda, double p){
  double K=2*pow(p,1./(1.-lamda));
  double a,b=0;
  for(int k=2;k<K;++k)
    b+=k*(pow(k+1,1-lamda)-pow(k,1-lamda));
  for(double x=0;x<1.01;x+=0.02){
    a=0;
    for(int k=2;k<K;++k)
      a+=k*(pow(k+1,1-lamda)-pow(k,1-lamda))*pow(x,k-1);
    a=a/b;
    cout<<x<<' '<<a<<endl;
  }
}
*/
#endif
