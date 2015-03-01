//2 D
#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <vector>
#include <list>
#include <queue>
#include "./statool/srand.hpp"
#include "node.hpp"


class Graph {
  
protected:
  //total num < 2*10^9
  vector<Node> vertices;
  
public:
  int get_deg_vertex(int idx_vertex);
  int get_num_vertices();
  int size();
  int get_num_validvertices();
  bool protect_check(int node_idx);//check if node is protected or not
  bool connection_check(int node1,int node2);//check if node1 and node2 are connected or not
  double distance(int node_idx);
  double aved(int node_idx);
  double btw(int node_idx);
  void insert_connection(int n1,int n2);
  void add_a_dependant(int n1,int n2){vertices[n1].add_a_dependant(n2);}//*
  void rm_a_node(int node_idx);//**
  void show_neighbors(int node_idx,vector<int> &v);
  void clear() {vertices.clear();}
  void protect_node(int node_idx);
  void unprotect_node(int node_idx);
  void unprotect_allnodes();
  void set_aved();
  void set_btwness();
  void resize(int num_nodes);
  double average_degree_of_network();

  typedef list<int>::iterator node_neighbor_iterator;
  node_neighbor_iterator vertex_neighbor_begin(size_t index) { 
    return vertices[index].neighbor_begin(); }
  node_neighbor_iterator vertex_neighbor_end(size_t index) { 
    return vertices[index].neighbor_end(); }
  Graph(int num_nd) {
    vertices.resize(num_nd);
  }
};

//---------------------------------------------
void Graph::insert_connection(int n1,int n2) {
  vertices[n1].add_a_neighbor(n2);
  //vertices[n2].add_a_neighbor(n1);         //this is no need in reading file of projectgao
}
//---------------------------------------------
void Graph::rm_a_node(int node_idx){
  node_neighbor_iterator it;
  for(it=vertex_neighbor_begin(node_idx);it!=vertex_neighbor_end(node_idx);++it){
    vertices[*it].rm_a_neighbor(node_idx);
  }
  vertices[node_idx].rm_all_neighbors();
  return;
}
//--------------------------------------------
void Graph::show_neighbors(int node_idx,vector<int> &v){
  //vector<int> v;
  node_neighbor_iterator it; 
  for(it=vertex_neighbor_begin(node_idx);it!=vertex_neighbor_end(node_idx);++it){
    {//cout<<node_idx<<' '<<*it<<endl;
      v.push_back(*it);
    }
  }
}
//--------------------------------------------
int Graph::get_deg_vertex(int idx_vertex) {
  return vertices[idx_vertex].get_degree();
}
//--------------------------------------------
int Graph::get_num_vertices(){
  return vertices.size();
} 
int Graph::size() {
  return vertices.size();
}
//--------------------------------------------
int Graph::get_num_validvertices(){
  int a=0;
  for(int i=0;i<vertices.size();++i){
    if(get_deg_vertex(i)){
      a++;//cout<<i<<endl;
    }
  }
  return a;
}
//-------------------------------------------
double Graph::distance(int node_idx){return vertices[node_idx].distance();}
double Graph::aved(int node_idx){return vertices[node_idx].aved();}
double Graph::btw(int node_idx){return vertices[node_idx].tbtw();}
//-------------------------------------------
void Graph::protect_node(int node_idx){
  vertices[node_idx].protect_node();
}
void Graph::unprotect_node(int node_idx){
  vertices[node_idx].unprotect_node();
}
void Graph::unprotect_allnodes(){
  for(int i=0;i<vertices.size();++i)
    vertices[i].unprotect_node();
}
bool Graph::protect_check(int node_idx){
  return vertices[node_idx].protect();
}
//----------------------------------------
bool Graph::connection_check(int node1,int node2){
  node_neighbor_iterator it;
  for(it=vertex_neighbor_begin(node1);it!=vertex_neighbor_end(node1);++it)
    if(*it==node2)
      return true;
  return false;
}
void Graph::resize(int num_nodes){
  vertices.resize(num_nodes);
}
//----------------------------------------
void Graph::set_aved(){
  list<int> s;
  list<int>::iterator it;
  int count,k;
  vector<int> distance(vertices.size(),-1);
  for(int node_idx=0;node_idx<vertices.size();++node_idx){
    //cout<<node_idx<<endl;
    for(int j=0;j<vertices.size();++j)
      vertices[j].set_distance(-1);
    vertices[node_idx].set_distance(0);
    s.clear();
    s.push_back(node_idx);
    count=1;
    while(s.size()){
      k=s.front();s.pop_front();
      if(vertices[k].get_degree())
	for(it=vertices[k].neighbor_begin();it!=vertices[k].neighbor_end();++it){
	  if(vertices[*it].distance()==-1){
	    vertices[*it].set_distance(vertices[k].distance()+1);
	    vertices[node_idx].set_aved(vertices[node_idx].aved()+vertices[*it].distance());
	    s.push_back(*it);
	    ++count;
	  }
	}
    }
    if(count<(vertices.size()*0.8))
      vertices[node_idx].set_aved(0);
    else
      vertices[node_idx].set_aved(vertices[node_idx].aved()/double(vertices.size()-1));
  }
}
//------------------------------------
void Graph::set_btwness(){
  list<int> q;
  list<int>::iterator it;
  int j,k,max_d=0;

  for(int i=0;i<vertices.size();++i){
    for(j=0;j<vertices.size();++j)
      vertices[i].set_distance(0);
    q.push_back(i);
    vertices[i].set_distance(1);
    max_d=0;
    //set distance;
    while(!q.empty()){
      k=q.front();q.pop_front();
      for(it=vertices[k].neighbor_begin();it!=vertices[k].neighbor_end();++it)
	if(vertices[*it].distance()==0){
	  vertices[*it].set_distance(vertices[k].distance()+1);
	  max_d=vertices[*it].distance();
	  q.push_back(*it);
	}
    }
    q.clear();
    //end setting distance
    for(j=0;j<vertices.size();++j){
      vertices[j].set_btw(1.);
      if(vertices[j].distance()==max_d) q.push_back(j);
    }
    while(!q.empty()){
      k=q.front();q.pop_front();
      j=0;
      for(it=vertices[k].neighbor_begin();it!=vertices[k].neighbor_end();++it){
	if((vertices[k].distance()-vertices[*it].distance()==1)&&(find(q.begin(),q.end(),*it)==q.end()))
	  ++j;
      }
      for(it=vertices[k].neighbor_begin();it!=vertices[k].neighbor_end();++it)
	if((vertices[k].distance()-vertices[*it].distance()==1)&&(find(q.begin(),q.end(),*it)==q.end())){
	  vertices[*it].set_btw(vertices[*it].btw()+vertices[k].btw()/double(j));
	  //	  vertices[*it].btw+=vertices[k].btw/double(j);
	  q.push_back(*it);
	}
    }
    for(j=0;j<vertices.size();++j)
      vertices[j].set_tbtw(vertices[j].tbtw()+vertices[j].btw());
  }
}

//-------------------------------------
double Graph::average_degree_of_network(){
int i;
double totaldegree=0;

for(i=0;i<vertices.size();i++) {totaldegree=totaldegree+get_deg_vertex(i);}

totaldegree=totaldegree/vertices.size();

return totaldegree;}


#endif
    
