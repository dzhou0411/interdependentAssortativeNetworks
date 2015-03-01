#ifndef NODE_HPP
#define NODE_HPP

class Node {

protected:
  list<int> neighbors;
  list<int> dependants;
  int status;//0:free,1:protected
  double btwness;
  double tbtwness;
  double d;
  double ave_d;

public:
  void add_a_neighbor(int other_node_index);
  void rm_a_neighbor(int node_index);
  void rm_all_neighbors();
  int get_degree();
  void protect_node();
  void unprotect_node();
  bool protect();
  double distance();
  double aved();
  double btw();
  double tbtw();
  void set_distance(double a);
  void set_aved(double a);
  void set_btw(double a);
  void set_tbtw(double a);
  void rm_a_dependant(int node_index);
  void rm_all_dependants(){dependants.clear();}
  void add_a_dependant(int node_index){dependants.push_back(node_index);}

  typedef list<int>::iterator neighbor_iterator;
  neighbor_iterator neighbor_begin() { return neighbors.begin(); }
  neighbor_iterator neighbor_end() { return neighbors.end(); }
  Node() {
    status=0;
    btwness=0;
    tbtwness=0;
    d=0;
    ave_d=0;
  }
};
//---------------------------------------
void Node::rm_a_neighbor(int node_index){
  neighbor_iterator it;
  for(it=neighbor_begin();it!=neighbor_end();++it){
    if(*it==node_index){
      break;
    }
  }
  if(it!=neighbor_end()){
    it=neighbors.erase (it);
  }
  else{
    cerr<<"don't have such neighbor!"<<endl;
    exit(1);
  }
}
void Node::rm_a_dependant(int node_index){
  neighbor_iterator it;
  for(it=dependants.begin();it!=dependants.end();++it){
    if(*it==node_index){
      break;
    }
  }
  if(it!=dependants.end()){
    it=dependants.erase (it);
  }
  else{
    cerr<<"don't have such dependant!"<<endl;
    exit(1);
  }
}
//---------------------------------------
void Node::rm_all_neighbors(){
  neighbors.clear();
  return;
}
//--------------------------------------
void Node::add_a_neighbor(int other_node_index) {
  neighbors.push_back(other_node_index);
}
//---------------------------------------
int Node::get_degree() { return neighbors.size(); }
//---------------------------------------
void Node::protect_node(){status=1;}
void Node::unprotect_node(){status=0;}
bool Node::protect(){
  if(status) return true;
  else return false;}
//---------------------------------------
double Node::distance(){  return d;}
double Node::aved(){return ave_d;}
double Node::btw(){return btwness;}
double Node::tbtw(){return tbtwness;}
void Node::set_distance(double a){d=a;}
void Node::set_aved(double a){ave_d=a;}
void Node::set_btw(double a){btwness=a;}
void Node::set_tbtw(double a){tbtwness=a;}
#endif
