#ifndef COMMON_HPP
#define COMMON_HPP
#include <string>
#include <vector>
#include <set>
#include <map>
#include <cmath>
#include <iostream>
#include <algorithm>
//-------------------------------//
//1.1 Convert digital num to char//
// num<999,999
//-------------------------------//
string d_c(int num){
  string nu,c="y";
  int mark=0;
  if(num<0) {mark=1;num=-num;}
  c[0]=num%10+48;
  nu=c+nu;
  if(num>9){c[0]=(num%100)/10+48;nu=c+nu;}
  if(num>99){c[0]=(num%1000)/100+48;nu=c+nu;}
  if(num>999){c[0]=(num%10000)/1000+48;nu=c+nu;}
  if(num>9999){c[0]=(num%100000)/10000+48;nu=c+nu;}
  if(num>99999){c[0]=(num%1000000)/100000+48;nu=c+nu;}
  
  if(mark) {c="-";nu=c+nu;}
  return nu;
}

//-------------------------------//
//1.2 Average and std
// double data stored in a vector
//-------------------------------//

pair<double,double> average_std(vector<double> a){
  pair<double,double> avg_std;
  long double sum=0.;
  for(int i=0;i<a.size();++i){
    sum+=a[i];
  }
  avg_std.first=sum/double(a.size());
  sum=0.;
  for(int i=0;i<a.size();++i){
    sum+=(a[i]-avg_std.first)*(a[i]-avg_std.first);
  }
  if(a.size()==1) sum=0.;
  else sum=sum/double(a.size()-1);
  avg_std.second=sqrt(sum);
  return avg_std;
}
pair<double,double> average_std2(vector<pair<int,double> > a){
  pair<double,double> avg_std;
  long double sum=0.;
  for(int i=0;i<a.size();++i){
    sum+=a[i].second;
  }
  avg_std.first=sum/double(a.size());
  sum=0.;
  for(int i=0;i<a.size();++i){
    sum+=(a[i].second-avg_std.first)*(a[i].second-avg_std.first);
  }
  if(a.size()==1) sum=0.;
  else sum=sum/double(a.size()-1);
  avg_std.second=sqrt(sum);
  return avg_std;
}

double correlation(vector<pair<double,double> > &v){
  vector<double> v1,v2;
  double Exy=0;
  for(int i=0;i<v.size();++i){
    //cout<<v[i].first<<' '<<v[i].second<<endl;
    v1.push_back(v[i].first);
    v2.push_back(v[i].second);
    Exy+=v[i].first*v[i].second;
  }
  Exy/=double(v.size());
  pair<double,double> a1,a2;
  a1=average_std(v1);
  a2=average_std(v2);
  double cor=(Exy-a1.first*a2.first)/a1.second/a2.second;
  return cor;
}

//-------------------------//
//1.3 Bubble descendent sort //
//-------------------------//
void bubblesort(vector<double> &a)
{
  int i,j,flag=1;
  double temp;
  int arraylength=a.size();

  for(i=1;(i<=arraylength)&&flag;i++){
    flag=0;
    for(j=0;j<(arraylength-1);j++){
      if(a[j+1]>a[j]){
	temp=a[j];
	a[j]=a[j+1];
	a[j+1]=temp;
	flag=1;
      }
    }
  }
  return;
}

void bubblesort2(vector<pair<int,double> > &a){
  int i,j,flag=1;
  pair<int,double> temp;
  int arraylength=a.size();

  for(i=1;(i<=arraylength)&&flag;i++){
    flag=0;
    for(j=0;j<(arraylength-1);j++){
      if(a[j+1].second>a[j].second){
	temp=a[j];
	a[j]=a[j+1];
	a[j+1]=temp;
	flag=1;
      }
    }
  }
  return;
}
//-----------------------------------------//
//1.4 Correlation                           //
//-----------------------------------------//
double correlation(vector<double> a,vector<double> b){
  double ave=0.,ave1=0.,ave2=0.,std1=0.,std2=0.;
  double corr;
  int arraylength;
  if(a.size()!=b.size()) cout<<"Correlation ERREO!"<<endl;
  else arraylength=a.size();

  for(int i=0;i<arraylength;i++){
    ave1=ave1+a[i];
    ave2=ave2+b[i];
    ave=ave+a[i]*b[i];
  }
  ave1=ave1/arraylength;
  ave2=ave2/arraylength;
  ave=ave/arraylength;
  for(int i=0;i<arraylength;i++){
    std1=std1+(a[i]-ave1)*(a[i]-ave1);
    std2=std2+(b[i]-ave2)*(b[i]-ave2);
  }
  std1=sqrt(std1/arraylength);
  std2=sqrt(std2/arraylength);
  corr=(ave-ave1*ave2)/std1/std2;
  return corr;
}
//-------------------------------------------//
//1.5 CDF and PDF
//num is the number of points what to have
//-------------------------------------------//
vector<pair<double,double> > cdf(vector<double> a,int num){
  double min=*min_element(a.begin(),a.end());
  double max=*max_element(a.begin(),a.end());
  double delta=(max-min)/double(num)+(max-min)/1000.;
  vector<double> pdf(num,0.);
  vector<pair<double,double> > b(num);
  
  for(int i=0;i<a.size();++i){
    ++pdf[int((a[i]-min)/delta)];
  }
  for(int i=num-2;i>=0;--i)
    pdf[i]+=pdf[i+1];
  for(int i=0;i<num;++i){
    b[i].first=min+i*delta;
    b[i].second=pdf[i]/a.size()/delta;
  }
  return b;
}
vector<pair<double,double> > pdf(vector<double> a,int num){
  double min=*min_element(a.begin(),a.end());
  double max=*max_element(a.begin(),a.end());
  double delta=(max-min)/double(num)+(max-min)/1000.;
  vector<double> pdf(num,0.);
  vector<pair<double,double> > b(num);
 
  for(int i=0;i<a.size();++i){
    ++pdf[int((a[i]-min)/delta)];
  }
  for(int i=0;i<num;++i){
    b[i].first=min+i*delta;
    b[i].second=pdf[i]/a.size()/delta;
  }
  return b;
}
void logpdf(vector<double> &a,int num,vector<pair<double,double> > &b){
  double min=*min_element(a.begin(),a.end());
  double max=*max_element(a.begin(),a.end());
  double delta=(log10(max)-log10(min))/double(num)+(log10(max)-log10(min))/1000.;
  vector<double> pdf(num,0.);
  b.clear();b.resize(num);
  for(int i=0;i<a.size();++i){
    ++pdf[int((log10(a[i])-log10(min))/delta)];
  }
  for(int i=0;i<num;++i){
    b[i].first=pow(10,log10(min)+i*delta);
    b[i].second=pdf[i]/a.size()/(pow(10,log10(min)+i*delta)-pow(10,log10(min)+(i-1)*delta));
  }
}
//--------------------------------------//
//1.6
//--------------------------------------//
int factorial (int num)
{
 if (num==1)
   return 1;
 return factorial(num-1)*num;
}

#endif
