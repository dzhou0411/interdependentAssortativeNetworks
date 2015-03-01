#ifndef NET_ALGO_HPP
#define NET_ALGO_HPP


//-----------------------------//
// Build Random Regular network//
//-----------------------------//
void lt_RR_algo(int num_nodes,Graph & lt_net, const int avg_degree){
 lt_net.clear();
  lt_net.resize(num_nodes);

  int a,b;
  for(int i=0;i<num_nodes*avg_degree/2;i++){//cout << "done link no. " << "  " << i << endl;
       do {a=int(srand()*num_nodes);} while (lt_net.get_deg_vertex(a)==4);
       do {do {b=int(srand()*num_nodes);} while (b==a);} while (lt_net.get_deg_vertex(b)==4);
  lt_net.insert_connection(a,b); 
                                      }
                                                                        }     

//---------------------------//
// Build ER network
//---------------------------//
inline void lt_ER_algo(int num_nodes,Graph & lt_net, const int avg_degree){
  lt_net.clear();
  lt_net.resize(num_nodes);
  int a,b;
  for(int i=0;i<num_nodes*avg_degree/2;){
    a=int(srand()*num_nodes);
    b=int(srand()*num_nodes);
    if((a!=b)&&(!lt_net.connection_check(a,b))){
      ++i;
      lt_net.insert_connection(a,b);
    }
  }
}
inline void lt_ER_algo2(int num_nodes,Graph & lt_net, const int avg_degree, vector<pair<int,int> >& links){
  lt_net.clear();
  lt_net.resize(num_nodes);
  int a,b;
  for(int i=0;i<num_nodes*avg_degree/2;){
    a=int(srand()*num_nodes);
    b=int(srand()*num_nodes);
    if((a!=b)&&(!lt_net.connection_check(a,b))){
      ++i;
      lt_net.insert_connection(a,b);
      links.push_back(make_pair(a,b));
    }
  }
}
//----------------------------//
// Build SF network
//----------------------------//
inline void lt_SF_algo(int num_nodes,Graph &net, const double lamda, int M){
  int k,h; //degree of a node
  vector<int> v;
  net.clear();net.resize(num_nodes);
  for(int j=0;j<num_nodes;++j){
    k=int(M*pow(1-srand(),1./(1-lamda)));
    for(int i=0;i<k;++i)
      v.push_back(j);
  }
  cout<<"size at the beginning: "<<v.size()<<endl;
  int mm=0;
  int n=0,size=v.size();
  bool mark=true;
  while(mark){
    k=0;h=0;
    while((v[k]==v[h])||net.connection_check(v[h],v[k])){
      if(++n>100000){
	mark=false;break;
      }
      k=int(srand()*size);
      h=int(srand()*size);
    }
    if(!mark) {break;}
    net.insert_connection(v[k],v[h]);
    if(k<(size-2)&&h<(size-2)) {
      v[k]=v[size-1];
      v[h]=v[size-2];
    }
    else if(k>=size-2 && h>=size-2)
      ;
    else if(h>=size-2){
      n=(h==size-2)?size-1:size-2;
      v[k]=v[n];
    }
    else if(k>=size-2){
      n=(k==size-2)?size-1:size-2;
      v[h]=v[n];
    }
    size-=2;n=0;    
  }
  cout<<"size at the end: "<<size<<endl;
  return; 
}
inline void lt_SF_algo2(int num_nodes,Graph &net, const double lamda, vector<pair<int,int> >& links,int M){
  int k,h; //degree of a node
  vector<int> v;
  net.clear();net.resize(num_nodes);
  for(int j=0;j<num_nodes;++j){
    k=int(M*pow(1-srand(),1./(1-lamda)));
    for(int i=0;i<k;++i)
      v.push_back(j);
  }
  cout<<"size at the beginning: "<<v.size()<<endl;
  int mm=0;
  int n=0,size=v.size();
  bool mark=true;
  while(mark){
    k=0;h=0;
    while((v[k]==v[h])||net.connection_check(v[h],v[k])){
      if(++n>100000){
	mark=false;break;
      }
      k=int(srand()*size);
      h=int(srand()*size);
    }
    if(!mark) {break;}
    net.insert_connection(v[k],v[h]);
    links.push_back(make_pair(v[k],v[h]));
    if(k<(size-2)&&h<(size-2)) {
      v[k]=v[size-1];
      v[h]=v[size-2];
    }
    else if(k>=size-2 && h>=size-2)
      ;
    else if(h>=size-2){
      n=(h==size-2)?size-1:size-2;
      v[k]=v[n];
    }
    else if(k>=size-2){
      n=(k==size-2)?size-1:size-2;
      v[h]=v[n];
    }
    size-=2;n=0;    
  }
  cout<<"size at the end: "<<size<<endl;
  return; 
  }
//---------------
//realation between two nets
//---------------
void dependant_relation2(Graph &net_1,Graph &net_2,vector<int> &dep1,vector<int> &dep2){
  dep1.clear();dep2.clear();
  dep1.resize(net_1.size());dep2.resize(net_2.size());
  /*
  for(int i=0;i<net_1.size();++i){//random connect
    dep1[i]=i;
    dep2[i]=i;
  }//*/
  ///*
  vector<pair<int,double> > v1,v2;
  for(int i=0;i<net_1.size();++i){
    v1.push_back(make_pair(i,net_1.get_deg_vertex(i)));
    v2.push_back(make_pair(i,net_2.get_deg_vertex(i)));
  }
  bubblesort2(v1);
  bubblesort2(v2);
  /*
  for(int i=0;i<v1.size();++i){//nested
    dep1[v1[i].first]=v2[v2.size()-1-i].first;
    dep2[v2[i].first]=v1[v1.size()-1-i].first;
  }
  //*/
  for(int i=0;i<v1.size();++i){//anti-nested
    dep1[v1[i].first]=v2[i].first;
    dep2[v2[i].first]=v1[i].first;
  }//*/
}
#endif
