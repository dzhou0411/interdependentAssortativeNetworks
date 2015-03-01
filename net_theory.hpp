#ifndef NET_THEORY_HPP
#define NET_THEORY_HPP

double delta(double x, double k, double p)
{return exp(k*p*(x-1.))-x;}

//---------------------//
//1.g(p) 
//---------------------//
//1.1 g(p) for ER normal situation
//---------------------//
double ER_g_A(double p, double k){
  double upper=1.,lower=0.;
  double f=0.5;
  double delta=0.1;
  while(delta>0.000000001){
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


double ER_g_Adi(double p, double k){
   double upper=1.,lower=0.;
   double f;
   double deltaf=1.;
   while (deltaf>0.0000000001){
      f=(upper+lower)/2.;
      deltaf=abs(delta(f,k,p));
      if (delta(f,k,p)*delta(lower,k,p) < 0) { upper=f;}
      else                           { lower=f;}
      }
      
   return (1.-f);                                 }






//-----------------//
//1.2 g(p) for SF without upper cutoff but k>2
//-----------------//
double SF_G_A_1(double lembda,double z,int mm){
  double a=0,b=0,c=0;
  for(int k=mm;k<50000;++k){
    b+=k*(pow(double(mm)/k,lembda-1.)-pow(double(mm)/(k+1),lembda-1))*pow(z,k-1);
    c+=k*(pow(double(mm)/k,lembda-1)-pow(double(mm)/(k+1),lembda-1));
  }
  a=b/c;
  return a;
}
double SF_G_A_0(double lembda,double z,int mm){
  double a=0;
  for(int k=mm;k<50000;++k){
    a+=(pow(double(mm)/k,lembda-1)-pow(double(mm)/(k+1),lembda-1))*pow(z,k);
  }
  return a;
}
double SF_g_A(double p,double lembda,int mm){
  double upper=1.,lower=0.;
  double f=0;
  double delta=0.1;
  while(abs(delta)>0.00000001){
    if(delta>0){
      lower=f;
      f=(f+upper)/2.;
    }
    else{
      upper=f;    
      f=(f+lower)/2.;
    }
    delta=SF_G_A_1(lembda,1-p*(1-f),mm)-f;
  }
  return (1.-SF_G_A_0(lembda,1-p*(1-f),mm));
}
/*//------------------//
//1.3 g(p) for SF with upper and lower cutoff K>k>m
//------------------//
double SF_G_B_1(double lembda,double K,double m,double z){
  double a=0,b=0,c=0;
  for(int k=m;k<K;++k){
    b+=(pow(k+1,1-lembda)-pow(k,1-lembda))*k*pow(z,k-1);
    c+=(pow(k+1,1-lembda)-pow(k,1-lembda))*k;
  }
  a=b/c;
  return a;
}
double SF_G_B_0(double lembda,double K,double m,double z){
  double a=0;
  for(int k=m;k<K;++k){
    a+=(pow(k+1,1-lembda)-pow(k,1-lembda))*pow(z,k);
  }
  a/=(pow(K,1-lembda)-pow(m,1-lembda));
  return a;
}
double SF_g_B(double lembda,double K,double m,double p){//p is actually pp
  double upper=1.,lower=0.;
  double f=0.5;
  double delta=0.1;
  while(delta>0.00001){
    delta=SF_G_B_1(lembda,K,m,1-p*(1-f))-f;
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
  return (1.-SF_G_B_0(lembda,K,m,1-p*(1-f)));
}
*/

//------------------------------//
//1.4 g(p) for intentional attack ER
//parameter p_0 initial survive ratio
//---------------------//
double ER_g_A_intent(double p,double p_0, double k){
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
//-----------------------------//
//1.5 g(p) for intentional attack SF
//parameter p_0 initial survive ratio
//-----------------------------//
double SF_g_A_intent(double p,double p_0,double lembda,int mm){

}
//-----------------------------//
//1.6 g(p) for maxdeg attack ER
//-----------------------------//
double ER_g_A_maxdeg(double p,double p_0,double k){
  
}
//-----------------------------//
//1.7 g(p) for maxdeg attack SF
//-----------------------------//
double SF_g_A_maxdeg(double p,double p_0,double lembda){
  double K_tild=2*pow(1-p_0,1./(1-lembda));
  double p_tild=1-pow(1-p_0,(2-lembda)/(1-lembda));
  double c=(1-lembda)/(0-pow(2,1-lembda));
  double f_A=0.5;
  double upper=1.,lower=0.;
  double delta=0.1;
  double G1,k=0;
  for(int i=2;i<K_tild+1;++i)
    k+=i*(pow(i,1-lembda)-pow(i+1,1-lembda))/pow(2,1-lembda)/p_0;
  while(delta>0.000001){
    G1=0;
    for(int i=2;i<K_tild+1;++i)
      G1+=i*(pow(i,1-lembda)-pow(i+1,1-lembda))/pow(2,1-lembda)/p_0*pow(1+p_tild/p_0*(f_A*p-p),i-1);
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
    G0+=(pow(i,1-lembda)-pow(i+1,1-lembda))/pow(2,1-lembda)/p_0*pow(1+p_tild/p_0*(f_A*p-p),i);
  return 1-G0;
}
//c.1.8 g(p) for (alpha) intentional attack ER
double ER_g_A_intent_alpha(double p,double p_0,double avg_k,double alpha){
  long double ER_f_alpha(double alpha,double p_0,double avg_k);
  long double ER_p_tild_alpha(double alpha,double p_0,double avg_k);
  long double f_alpha=ER_f_alpha(alpha,p_0,avg_k);
  long double p_tild=ER_p_tild_alpha(alpha,p_0,avg_k);
  long double c=0;
  long double temp1;  
  for(int i=1;i<999;i++){
    temp1=1;
    for(int j=1;j<=i;++j)
      temp1*=avg_k/j;
    if(temp1<=0) break;
    c+=temp1*exp(-avg_k)*i*pow(f_alpha,pow(i,alpha));
  }
  double f_A=0.5;
  double upper=1.,lower=0.,delta=0.1;
  double G1;
  while(abs(delta)>0.00001){
    //cout<<"Upper:"<<upper<<" Lower:"<<lower<<' '<<delta<<endl;
    if(upper==0.||lower==1.) {cout<<"fail!"<<endl;break;}
    G1=0;
    for(int i=1;i<999;++i){
      temp1=1;
      for(int j=1;j<=i;++j)
	temp1*=avg_k/j;
      if(temp1<=0) break;
      G1+=temp1*exp(-avg_k)*i*pow(f_alpha,pow(i,alpha))*pow(1+p_tild/p_0*p*(f_A-1),i-1);
    }
    G1/=c;
    delta=G1-f_A;
    if(delta>0){
      lower=f_A;
      f_A=(f_A+upper)/2.;
    }else{
      upper=f_A;
      f_A=(f_A+lower)/2.;
    }
  }
  double G0=0;
  for(int i=1;i<999;++i){
    temp1=1;
    for(int j=1;j<=i;++j)
      temp1*=avg_k/j;
    if(temp1<=0) break;
    G0+=temp1*exp(-avg_k)*pow(f_alpha,pow(i,alpha))*pow(1+p_tild/p_0*p*(f_A-1),i);
  }
  G0+=exp(-avg_k);
  G0=G0/p_0;
  return 1-G0;
}
double SF_g_A_intent_alpha(double p,double p_0,double lembda,double alpha,int mm){
  long double M=mm;
  long double SF_f_alpha(double alpha,double p_0,double lembda,int mm);
  long double SF_p_tild_alpha(double alpha,double p_0,double lembda,int mm);
  long double f_alpha=SF_f_alpha(alpha,p_0,lembda,mm);
  long double p_tild=SF_p_tild_alpha(alpha,p_0,lembda,mm);
  long double c=0;// sum of bottom in f_A equation
  long double temp1;
  for(int i=M;i<100000;i++)
    c+=(pow(i,1-lembda)-pow(i+1,1-lembda))/(pow(M,1-lembda))*i*pow(f_alpha,pow(i,alpha));
  long double f_A=0.5;
  long double upper=1.,lower=0.,delta=0.1;
  long double G1;
  while(abs(delta)>0.00001){
    //cout<<f_A<<" Upper:"<<upper<<" Lower:"<<lower<<' '<<delta<<endl;
    if(upper==0.||lower==1.) {cout<<"fail!"<<endl;break;}
    G1=0;
    for(long double i=M;i<10000;i+=1.){
      G1+=(pow(i,1-lembda)-pow(i+1,1-lembda))/(pow(M,1-lembda))*i*pow(f_alpha,pow(i,alpha))*pow(1+p_tild/p_0*p*(f_A-1),i-1);
    }
    G1/=c;
    delta=G1-f_A;
    if(delta>0){
      lower=f_A;
      f_A=(f_A+upper)/2.;
    }else{
      upper=f_A;
      f_A=(f_A+lower)/2.;
    }
  }
  double G0=0;
  for(int i=M;i<10000;++i){
    G0+=(pow(i,1-lembda)-pow(i+1,1-lembda))/(pow(M,1-lembda))*pow(f_alpha,pow(i,alpha))*pow(1+p_tild/p_0*p*(f_A-1),i);
  }
  G0=G0/p_0;
  return 1-G0;
}

//-----------------------------------------------
//1.9 reversed function f=G^-1(r) for ER network
//----------------------------------------------
long double ER_f_alpha(double alpha,double p,double avg_k){
  long double f=1,upper=1,lower=0.000;
  long double delta=1;
  long double G,temp1,a=10,b=-4949;
  long double c=pow(a,b);
  while(abs(delta)>0.00001){
    if(delta>0){
      upper=f;
      f=(f+lower)/2;
    }else{
      lower=f;
      f=(f+upper)/2;
    }
    if((upper-lower)<c) break;
    G=0;
    for(int i=1;i<1000;++i){
      temp1=1;
      for(int j=1;j<=i;++j){
	if(j==0) temp1*=avg_k;
	else temp1*=avg_k/j;
      }
      if(temp1<=0) break;
      G+=temp1*exp(0-avg_k)*pow(f,pow(i,alpha));
    }
    G+=exp(0-avg_k);// sum from k=0
    delta=G-p;
  }
  return f;
}
long double SF_f_alpha(double alpha,double p,double lembda,int mm){
  int M=mm;
  long double f=1,upper=100,lower=0.0;
  long double delta=1;
  long double G,temp1;
  long double a=10.,b=-4949;
  long double c=pow(a,b);
  while(abs(delta)>0.0001){
    if(delta>0){
      upper=f;
      f=(f+lower)/2;
    }else{
      lower=f;
      f=(f+upper)/2;
    }
    if((upper-lower)<c) break;
    G=0;
    for(int i=M;i<1000;++i){
      G+=(pow(i,1-lembda)-pow(i+1,1-lembda))/(pow(M,1-lembda))*pow(f,pow(i,alpha));
    }
    delta=G-p;
  }
  return f;
}
//--------------------------------------------
//1.10 p_tild for (alpha) ER
//-------------------------------------------
double ER_p_tild_alpha(double alpha,double p,double avg_k){
  double f_alpha=ER_f_alpha(alpha,p,avg_k);
  double m=0,temp=1;
  for(int i=1;i<500;++i){
    temp=1;
    for(int j=1;j<=i;++j){//**pay attention! j<=i
      temp*=avg_k/j;
    }
    if(temp<=0) break;
    m+=temp*exp(0-avg_k)*i*pow(f_alpha,pow(i,alpha));
  }
  double p_tild=m/avg_k;
  return p_tild;
}
double SF_p_tild_alpha(double alpha,double p,double lembda,int mm){
  int M=mm;
  double f_alpha=SF_f_alpha(alpha,p,lembda,mm);
  double m=0,avg_k=0;
  for(int i=M;i<1000;++i){
    m+=(pow(i,1-lembda)-pow(i+1,1-lembda))/pow(M,1-lembda)*i*pow(f_alpha,pow(i,alpha));
    avg_k+=(pow(i,1-lembda)-pow(i+1,1-lembda))/(pow(M,1-lembda))*i;
  }
  double p_tild=m/avg_k;
  return p_tild;
}
//-----------------------------------//
//2. Theoretical cascade step by step
//-----------------------------------//
//2.1 ER normal situation
//--------------------------------------//
void ER_step_theo(int p_0,double k){
  string file="../data/step/ER_k"+d_c(k)+"_theo.dat";
  ofstream ofile(file.c_str());
  double p1,p2,p3=2.;
  p1=p_0*ER_g_A(p_0,double(k));
  p2=p1;
  while(abs(p1-p3)>0.000001){
    p3=p1;
    p1=p2*ER_g_A(p2,double(k));
    p2=p_0*ER_g_A(p2,double(k));
    ofile<<p3/p_0<<endl;
  }
  ofile.close();
}
//2.2 SF normal situation
void SF_step_theo(double p_0, double lembda,int mm){
  string file="../data/step/SF_lembda"+d_c(int(lembda*10.001))+"_theo.dat";
  ofstream ofile(file.c_str());
  double p1,p2,p3=2.;
  p1=p_0*SF_g_A(lembda,p_0,mm);
  p2=p1;
  while(abs(p1-p3)>0.0001){
    p3=p1;
    p1=p2*SF_g_A(lembda,p2,mm);
    p2=p_0*SF_g_A(lembda,p2,mm);
    ofile<<p3<<endl;
  }
  ofile.close();
}

//2.3 ER intentional attack situation
void ER_step_theo_intent(double p_0,double k){
  string file="../data/step/ER_k"+d_c(k)+"_p0"+d_c(int(p_0*100.01))+"_theo_intent.dat";
  ofstream ofile(file.c_str());
  vector<double> p(7,0.);
  p[0]=p_0;
  p[1]=p[0];
  p[2]=p[1]*ER_g_A_intent(p[1],p_0,double(k));
  p[3]=p[0]*ER_g_A_intent(p[1],p_0,double(k));
  p[4]=p[3]*ER_g_A(p[3],double(k));
  ///*
  p[5]=p[0]*ER_g_A(p[3],double(k));
  p[6]=p[5]*ER_g_A_intent(p[5],p_0,double(k));
  int i=p.size();
  while(abs(p[i-1]-p[i-5])>0.0001){
    if(i%4==1) {
      p.push_back(p[0]*ER_g_A(p[i-2],double(k)));
      i=p.size();
      p.push_back(p[i-1]*ER_g_A_intent(p[i-1],p_0,double(k)));
      i=p.size();
    } else if(i%4==3){
      p.push_back(p[0]*ER_g_A_intent(p[i-2],p_0,double(k)));
      i=p.size();
      p.push_back(p[i-1]*ER_g_A(p[i-1],double(k)));
      i=p.size();
    }
  }
  for(int i=2;i<p.size();i+=2){
    if((p[i]!=0)&&(p[i+2]!=0))
      ofile<<i/2-1<<' '<<p[i]<<endl;
  }
  ofile.close();
} //*/
//2.4 SF intentional attack situation
void SF_step_theo_intent(double p_0,double lembda){
}
//2.5 ER maxdeg attack situation
void ER_step_theo_maxdeg(double p_0,double k){
}
//2.6 SF maxdeg attack situation
void SF_step_theo_maxdeg(double p_0,double lembda,int mm){
  string file="../data/step/SF_lembda"+d_c(lembda*10.001)+"_p0"+d_c(int(p_0*100.01))+"_theo_maxdeg.dat";
  ofstream ofile(file.c_str());
  vector<double> p(7,0.);
  p[0]=p_0;
  p[1]=p[0];
  p[2]=p[1]*SF_g_A_maxdeg(p[1],p_0,lembda);
  p[3]=p[0]*SF_g_A_maxdeg(p[1],p_0,lembda);
  p[4]=p[3]*SF_g_A(p[3],lembda,mm);
  p[5]=p[0]*SF_g_A(p[3],lembda,mm);
  p[6]=p[5]*SF_g_A_maxdeg(p[5],p_0,lembda);
  int i=p.size();
  while(abs(p[i-1]-p[i-5])>0.0001){
    if(i%4==1) {
      p.push_back(p[0]*SF_g_A(p[i-2],lembda,mm));
      i=p.size();
      p.push_back(p[i-1]*SF_g_A_maxdeg(p[i-1],p_0,lembda));
      i=p.size();
    } else if(i%4==3){
      p.push_back(p[0]*SF_g_A_maxdeg(p[i-2],p_0,lembda));
      i=p.size();
      p.push_back(p[i-1]*SF_g_A(p[i-1],lembda,mm));
      i=p.size();
    }
  }
  for(int i=2;i<p.size();i+=2){
    if((p[i]!=0)&&(p[i+2]!=0))
      ofile<<i/2-1<<' '<<p[i]<<endl;
  }
  ofile.close();
}
//2.7 ER intentional attack situation
void ER_step_theo_intent_alpha(double p_0,double k,double alpha){
  string file="../data/step/ER_k"+d_c(k)+"_p0"+d_c(int(p_0*100.01))+"_theo_intent_alpha.dat";
  ofstream ofile(file.c_str());
  vector<double> p(7,0.);
  p[0]=p_0;
  p[1]=p[0];
  p[2]=p[1]*ER_g_A_intent_alpha(p[1],p_0,double(k),alpha);
  p[3]=p[0]*ER_g_A_intent_alpha(p[1],p_0,double(k),alpha);
  p[4]=p[3]*ER_g_A(p[3],double(k));
  ///*
  p[5]=p[0]*ER_g_A(p[3],double(k));
  p[6]=p[5]*ER_g_A_intent_alpha(p[5],p_0,double(k),alpha);
  int i=p.size();
  while(abs(p[i-1]-p[i-5])>0.0001){
    if(i%4==1) {
      p.push_back(p[0]*ER_g_A(p[i-2],double(k)));
      i=p.size();
      p.push_back(p[i-1]*ER_g_A_intent_alpha(p[i-1],p_0,double(k),alpha));
      i=p.size();
    } else if(i%4==3){
      p.push_back(p[0]*ER_g_A_intent_alpha(p[i-2],p_0,double(k),alpha));
      i=p.size();
      p.push_back(p[i-1]*ER_g_A(p[i-1],double(k)));
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
//3 looking for ER critical point
//----------------------------//
//3.1 print pg_A(pg_B(x)) for intentional attack for different x
void function1(double p_0,int avg_deg){
  double f_B,f_A,y;
  double lower,upper,delta;
  double f=log(p_0)/avg_deg+1;
  double p_tild=f*exp(avg_deg*(f-1));
  string file;ofstream ofile;
  file="../data/pg_A(pg_B(x))_intent.dat";
  ofile.open(file.c_str());
  for(double x=0.;x<=1.;x+=0.02){
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
}
//3.2 pg_A(pg_B(x)) for intentional attack
double pg_Apg_B(double p_0,double x, double avg_deg,double avg_deg2){
  double f_B,f_A,y;
  double lower,upper,delta;
  double f=log(p_0)/avg_deg+1;
  double p_tild=f*p_0;//f*exp(avg_deg*(f-1));
  string file;ofstream ofile;
  upper=1.;lower=0;delta=0.1;
  f_B=0.5;
  while(delta>0.000001){
    delta=exp(avg_deg2*x*(f_B-1))-f_B;
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
    delta=exp(avg_deg*f*f*(f_A-1)*y)-f_A;
    if(delta>0){
      lower=f_A;
      f_A=(f_A+upper)/2.;
    }else{
      upper=f_A;
      f_A=(f_A+lower)/2.;
    }
    delta=abs(delta);
  }
   return p_0*(1-f_A);
}
//3.3 find critical point p_c for intentional attack
// 
double p_c(double avg_deg,double avg_deg2){
  double L=0.,U=1.,D=1.;
  double p_0,x,y,y_x,temp;
  p_0=0.5;
  while(abs(D)>0.005){
    //--find x
    double lower=0.01;
    double upper=1.,delta1=1.;
    y=0;
    while(y<0.0001){
      lower=lower+0.01;
      //cout<<lower<<endl;
      y=pg_Apg_B(p_0,lower,avg_deg,avg_deg2);
      if(lower>=1.){
	goto ooo;
      }
    }
    x=0.999;
    y=pg_Apg_B(p_0,x,avg_deg,avg_deg2);
    temp=y-x;
    upper=x;
    x=(x+lower)/2.;
    while(abs(delta1)>0.00001){
      y=pg_Apg_B(p_0,x,avg_deg,avg_deg2);
      y_x=y-x;
      delta1=y_x-temp;
      temp=y_x;
      if(delta1>0){
	upper=x;
	x=(x+lower)/2.;
      }else{
	lower=x;
	x=(x+upper)/2.;
      }
    }
    //---------
    D=y_x;
    if(D>0){
      U=p_0;
      p_0=(p_0+L)/2.;
    }else{
    ooo:
      L=p_0;
      p_0=(p_0+U)/2.;
    }
  }
  return p_0;
}


//3.4 SF-SF, pg_A(pg_B(x)) for intentional attack
//3.5
double ERER_pg_Apg_B(double x,double p_0,double alpha,double avg_k,double avg_k2){
  double y=0,z=0;
  y=p_0*ER_g_A(x,avg_k2);
  if(y<0.01) return 0;
  z=p_0*ER_g_A_intent_alpha(y,p_0,avg_k,alpha);
  return z;
}
long double SFSF_pg_Apg_B(double x,double p_0,double alpha,double lembda,double lembda2,int mm){
  long double y,z;
  y=p_0*SF_g_A(x,lembda2,mm);
  //cout<<"y:"<<y<<endl;
  if(y<0.001) return 0;
  z=p_0*SF_g_A_intent_alpha(y,p_0,lembda,alpha,mm);
  //cout<<"z:"<<z<<endl;
  return z; 
}
//3.6
/*double ERER_p_c(double avg_deg,double avg_deg2,double alpha){
  double L=0.,U=1.,D=1.;
  double p_0,x,y,y_x,temp;
  p_0=1;
  while(abs(D)>0.0001){//0.005
    if(D>0){
      U=p_0;
      p_0=(p_0+L)/2.;
    }else{
    ooo:
      L=p_0;
      p_0=(p_0+U)/2.;
    }    
    //--find x
    double lower=0.01;
    double upper=1.,delta1=1.;
    
    if(U-L<0.000001) {break;}
    y=0;
    while(y<0.0001){//0.0001
      lower=lower+0.01;
      y=ERER_pg_Apg_B(lower,p_0,alpha,avg_deg,avg_deg2);
      if(lower>=1.){
	goto ooo;cout<<"jump"<<endl;
      }
    }
    x=0.99999;
    y=ERER_pg_Apg_B(x,p_0,alpha,avg_deg,avg_deg2);
    temp=y-x;
    upper=x;
    x=(x+lower)/2.;
    while(abs(delta1)>0.00001){//0.001
      //cout<<x<<' '<<upper<<' '<<lower<<' '<<delta1<<endl;
      y=ERER_pg_Apg_B(x,p_0,alpha,avg_deg,avg_deg2);
      y_x=y-x;
      delta1=y_x-temp;
      temp=y_x;
      if(delta1>0){
	upper=x;
	x=(x+lower)/2.;
      }else{
	lower=x;
	x=(x+upper)/2.;
      }
    }
    //---------
    D=y_x;
    //cout<<p_0<<' '<<L<<' '<<U<<' '<<D<<endl;
  }
  return p_0;
  }//*/
//-------------------------
double ERER_p_c(double avgk,double avgk2,double alpha){
  long double L=0.,U=1.,D=1.;
  long double p_0,x,x1,y,y_plus,y_x,temp;
  p_0=1;
  while(abs(D)>0.0001){
    if(U-L<0.0001) {break;}
    if(D>0){
      U=p_0;
      p_0=(p_0+L)/2.;
    }else{
    ooo:
      L=p_0;
      p_0=(p_0+U)/2.;
    }
    //--find x
    long double lower=0.01;
    long double upper=1.,delta1=1.;
    
    y=0;
    while(y<0.001){//0.0001
      lower=lower+0.01;
      y=ERER_pg_Apg_B(lower,p_0,alpha,avgk,avgk2);
      if(lower>=1.){
	cout<<"jump"<<endl;
	goto ooo;
      }
    }
    cout<<"lower: "<<lower<<endl;
    x=1;
    y=1;y_plus=0;
    delta1=1;
    while(upper-lower>0.00001){
      if(y>y_plus){
	upper=x;
	x=(x+lower)/2.;
      }else if(y==y_plus){
	break;
      }else{
	lower=x;
	x=(x+upper)/2.;
      }
      y=ERER_pg_Apg_B(x,p_0,alpha,avgk,avgk2)-x;
      x1=x+0.001;
      y_plus=ERER_pg_Apg_B(x1,p_0,alpha,avgk,avgk2)-x1;
      cout<<"x,y,y+,upper,lower: "<<' '<<x<<' '<<y<<' '<<y_plus<<' '<<upper<<' '<<lower<<endl;
    }
    //---------
    D=y;
    cout<<"p0,L,U,D:"<<p_0<<' '<<L<<' '<<U<<' '<<D<<endl;
  }
  return p_0;
}
//-------------------------
pair<double,double> SFSF_p_c(double lembda,double lembda2,double alpha,int mm){
  long double L=0.,U=1.,D=1.;
  long double p_0,x,x1,y,y_plus,y_x,temp;
  long double const1,const2;
  const1=pow(10.,-6);const2=pow(10.,-10);
  p_0=1;
  while(abs(D)>0.0001){
    if(U-L<const1) {break;}
    if(D>0){
      U=p_0;
      p_0=(p_0+L)/2.;
    }else{
    ooo:
      L=p_0;
      p_0=(p_0+U)/2.;
    }
    //--find x
    long double lower=0.01;
    long double upper=1.,delta1=1.;
    
    y=0;
    while(y<0.001){//0.0001
      lower=lower+0.01;
      y=SFSF_pg_Apg_B(lower,p_0,alpha,lembda,lembda2,mm);
      if(lower>=1.){
	cout<<"jump"<<endl;
	goto ooo;
      }
    }
    cout<<"lower: "<<lower<<endl;
    x=1;
    y=1;y_plus=0;
    delta1=1;
    while(upper-lower>const2){
      if(y>y_plus){
	upper=x;
	x=(x+lower)/2.;
      }else if(y==y_plus){
	break;
      }else{
	lower=x;
	x=(x+upper)/2.;
      }
      y=SFSF_pg_Apg_B(x,p_0,alpha,lembda,lembda2,mm)-x;
      x1=x+(upper-lower)/100.;
      y_plus=SFSF_pg_Apg_B(x1,p_0,alpha,lembda,lembda2,mm)-x1;
      //cout<<"x,y,y+,y-y+,upper,lower: "<<' '<<x<<' '<<y<<' '<<y_plus<<' '<<y-y_plus<<' '<<upper<<' '<<lower<<endl;
    }
    //---------
    D=y;
    cout<<"p0,L,U,D:"<<p_0<<' '<<L<<' '<<U<<' '<<D<<endl;
  }
  return make_pair(x,p_0);
}//*/
//3.7. p_inf for alpha attack
double ERER_p_inf_alpha(double p_0, double k1,double k2,double alpha){
  vector<double> p(7,0.);
  p[0]=p_0;
  p[1]=p[0];
  p[2]=p[1]*ER_g_A_intent_alpha(p[1],p_0,double(k1),alpha);
  p[3]=p[0]*ER_g_A_intent_alpha(p[1],p_0,double(k1),alpha);
  p[4]=p[3]*ER_g_A(p[3],double(k2));
  ///*
  p[5]=p[0]*ER_g_A(p[3],double(k2));
  p[6]=p[5]*ER_g_A_intent_alpha(p[5],p_0,double(k1),alpha);
  int i=p.size();
  while(abs(p[i-1]-p[i-5])>0.0001){
    if(i%4==1) {
      p.push_back(p[0]*ER_g_A(p[i-2],double(k2)));
      i=p.size();
      p.push_back(p[i-1]*ER_g_A_intent_alpha(p[i-1],p_0,double(k1),alpha));
      i=p.size();
    } else if(i%4==3){
      p.push_back(p[0]*ER_g_A_intent_alpha(p[i-2],p_0,double(k1),alpha));
      i=p.size();
      p.push_back(p[i-1]*ER_g_A(p[i-1],double(k2)));
      i=p.size();
    }
  }
  double a=p.size();a-=1;
  return p[a];
} 
double SFSF_p_inf_alpha(double p_0, double lembda1,double lembda2,double alpha,int mm){
  vector<double> p(7,0.);
  p[0]=p_0;
  p[1]=p[0];
  p[2]=p[1]*SF_g_A_intent_alpha(p[1],p_0,double(lembda1),alpha,mm);
  p[3]=p[0]*SF_g_A_intent_alpha(p[1],p_0,double(lembda1),alpha,mm);
  p[4]=p[3]*SF_g_A(p[3],double(lembda2),mm);
  ///*
  p[5]=p[0]*SF_g_A(p[3],double(lembda2),mm);
  p[6]=p[5]*SF_g_A_intent_alpha(p[5],p_0,double(lembda1),alpha,mm);
  int i=p.size();
  while(abs(p[i-1]-p[i-5])>0.002){
    if(i%4==1) {
      p.push_back(p[0]*SF_g_A(p[i-2],double(lembda2),mm));
      i=p.size();
      p.push_back(p[i-1]*SF_g_A_intent_alpha(p[i-1],p_0,double(lembda1),alpha,mm));
      i=p.size();
    } else if(i%4==3){
      p.push_back(p[0]*SF_g_A_intent_alpha(p[i-2],p_0,double(lembda1),alpha,mm));
      i=p.size();
      p.push_back(p[i-1]*SF_g_A(p[i-1],double(lembda2),mm));
      i=p.size();
    }
    //cout<<p[i-1]<<' '<<p[i-5]<<' '<<p[i-1]-p[i-5]<<endl;
  }
  double a=p.size();a-=1;
  return p[a];
}



// 3.8. Double attack
//P_c double attack
double ERER_pc_doubleattack(double avg_k){
  double pc=1,upper=1,lower=0,delta=1;
  while(abs(delta)>0.0001){
    if(delta>0){
      upper=pc;
      pc=(pc+lower)/2;
    }else{
      lower=pc;
      pc=(pc+upper)/2;
    }
    delta=avg_k*pow(pc,2)*pow(log(pc)/avg_k+1,2)-2.4554;
    cout<<pc<<' '<<upper<<' '<<lower<<' '<<delta<<endl;
  }
  return pc;
}

// p_inf double attack
double ERER_p_inf_doubleattack(double p,double avg_k){
  double f=log(p)/avg_k+1;
  double upper=1,lower=0.000001,delta=1;
  double p_inf=0.5;
  while(abs(delta)>0.000001){
    delta=p_inf-pow(p,2)*pow(1-exp(-avg_k*f*f*p_inf),2);
    //cout<<p_inf<<' '<<lower<<' '<<upper<<' '<<delta<<endl;
    if(delta<0){
      lower=p_inf;
      p_inf=(p_inf+upper)/2;
    }else{
      upper=p_inf;
      p_inf=(p_inf+lower)/2;
    }
  }
  return p_inf;
}

// single network P_c for SF

double pc_SF_single(double lembda,int mm,double alpha){
  double delta=1,upper=1000000,lower=0,qc=0;
  double temp;
  double b=0;//avegrage degree
  for(int i=mm;i<100000;++i){
    b+=pow(i,2-lembda)-pow(i+1,1-lembda)*i;
  }
  b/=pow(mm,1-lembda);
  cout<<b<<endl;
  while(abs(delta)>0.000001){
    if(delta<0){
      upper=qc;
      qc=(qc+lower)/2;
    }else{
      lower=qc;
      qc=(qc+upper)/2;
    }
    temp=0;
    for(int i=mm;i<100000;++i)
      temp+=(pow(double(i)/double(mm),1-lembda)-pow(double(i+1)/double(mm),1-lembda))*i*(i-1)*exp(-qc*pow(i,alpha));
    temp/=b;
    delta=temp-1;
    cout<<qc<<' '<<upper<<' '<<lower<<' '<<delta<<endl;
  }
  b=0;
  cout<<"qc: "<<qc<<endl;
  for(int i=mm;i<100000;++i){
    b+=(pow(double(i)/double(mm),1-lembda)-pow(double(i+1)/double(mm),1-lembda))*exp(-qc*pow(i,alpha));
  }
  //  b/=pow(mm,1-lembda);
  return b;
}
      
      
      
/*
//c.1.8 g(p) for (alpha) intentional attack ER with approximation method

double ER_g_A_intent_alpha(double p,double p_0,double avg_k,double alpha){
  double k_alpha=0;//<k^alpha>
  double k_alpha1=0;//<k^(alpha+1)>
  double c=0;// sum(k*p_tild(k))
  double temp1,temp2;
  for(int i=1;i<500;++i){
    temp1=1;
    for(int j=1;j<=i;++j){
      temp1*=avg_k/j;
    }
    if(temp1<=0) {break;}
    k_alpha+=temp1*exp(0.-avg_k)*pow(i,alpha);
    //k_alpha+=pow(avg_k,i)*exp(0-avg_k)/factorial(i)*pow(i,alpha);
    k_alpha1+=temp1*exp(0-avg_k)*pow(i,alpha+1);
  }
  for(int i=1;i<500;++i){
    temp1=1;
    for(int j=1;j<=i;++j){
      temp1*=avg_k/j;
    }
    if(temp1<=0) {break;}
    c+=i*(1-(1-p_0)*pow(i,alpha)/k_alpha)/p_0*temp1*exp(0-avg_k);
  }
  double p_tild=1-(1-p_0)*k_alpha1/k_alpha/avg_k;
  cout<<"c:"<<c<<endl;

  double f_A=0.5;
  double upper=1.,lower=0.,delta=0.1;
  double G1;
  while(delta>0.00001){
    //cout<<"Upper:"<<upper<<" Lower:"<<lower<<' '<<delta<<endl;
    if(upper==lower) break;
    G1=0;
    for(int i=1;i<999999;++i){
      temp1=1;
      for(int j=1;j<=i;++j){
	temp1*=avg_k/j;
      }
      if(temp1<=0) {break;}
      G1+=i*(1-(1-p_0)*pow(i,alpha)/k_alpha)/p_0*temp1*exp(0-avg_k)*pow(1+p_tild*(f_A-1),i-1);
    }
    G1/=c;
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
  for(int i=0;i<999999;++i){
    temp1=1;
    for(int j=1;j<=i;++j){
      temp1*=avg_k/j;
    }
    if(temp1<=0) {break;}
    G0+=(1-(1-p_0)*pow(i,alpha)/k_alpha)/p_0*temp1*exp(0-avg_k)*pow(1+p_tild*(f_A-1),i);
  }
  return 1-G0;
  }*/

//031511 calculate the xc value when x=y and x=pgA(pgB(x)) are tangential
//these series of functions work in this way: they try to find the upper intersection of gA(gB(x)) and x. When this intersection position suddenly goes to zero,     then it's the critical point of xc. 
//

double graph_xc_ERER(double k)
{ double xupper=0.99,xlower=0.1;
  double xcupper=1., xclower=0.;
  double delta=0.1, deltaderi=0.1,deltaxc=1.;
  double p=0.8;
  double x,y,deri;
  int counter=0;
  double deltap=0.00001;
  double pinfi=0;
 
  while (xcupper>0.0001)
  {
   pinfi=xcupper*ER_g_A(xcupper,k);
   p=p-deltap;
   cout << "p = " << p << endl;
        x=xupper; 	
	delta=0.1; counter=0;
        while (delta>0.0000001){
   	y=p*ER_g_A(ER_g_A(x,k)*p,k);
        delta=x-y;
        x=y;
   	counter=counter+1;
   	}
        xcupper=x; 
       cout << "xcupper" << xcupper << endl;      
   
      
 }
cout<< "P-infinity = " << pinfi << endl;
cout<< "pc = " << p+deltap << endl;
return (p+deltap);
 }




double graph_xc_SFSF(double lembda, int mm)
{ double xupper=0.99;
  double xcupper=0.;
  double delta=0.1, deltaderi=0.1,deltaxc;
  double p=0.8;
  double x,y,deri;
  int counter=0;
  double deltap=0.00001;
  double pinfi=0;
 
   while (1)
  {pinfi=xcupper*SF_g_A(xcupper,lembda,mm);
   p=p-deltap;
   cout << "p = " << p << endl;
        x=xupper; 	
	delta=0.1; 
        while (delta>0.0000001){
        y=p*SF_g_A(SF_g_A(x,lembda,mm)*p,lembda,mm);
        delta=x-y;
        x=y;
   	counter=counter+1;
   	}
        xcupper=x; 
        cout << "xcupper" << xcupper << endl;      
 
        if (xcupper<=0.0001) break;
  
 }
cout<< "p-infinity = " << pinfi << endl; 
return (p+deltap);
 }


double graph_xc_ERSF(double k,double lembda, int mm)
{ double xupper=0.99,xlower=0.1;
  double xcupper=0., xclower=0.;
  double delta=0.1, deltaderi=0.1,deltaxc=1.;
  double p=0.8;
  double x,y,deri;
  int counter=0;
  double deltap=0.001;
  double pinfi;
 
  while (1)
  {p=p-deltap;
   cout << "p = " << p << endl;
        x=xupper; 	
	delta=0.1; counter=0;
        pinfi=xcupper*SF_g_A(xcupper,lembda,mm);
        while (delta>0.0000001){
   	y=p*ER_g_A(SF_g_A(x,lembda,mm)*p,k);
        delta=x-y;
        x=y;
   	counter=counter+1;
   	}
        xcupper=x; 
       cout << "xcupper" << xcupper << endl;      
        if (xcupper<=0.0001) break;
      
 }
cout<< "p-infinity = " << pinfi << endl; 
return (p+deltap);
 }

double graph_xc_SFER(double k,double lembda, int mm)
{ double xupper=0.99,xlower=0.1;
  double xcupper=0., xclower=0.;
  double delta=0.1, deltaderi=0.1,deltaxc=1.;
  double p=0.8;
  double x,y,deri;
  int counter=0;
  double deltap=0.00001;
  double pinfi;
 
  while (1)
  {p=p-deltap;
   cout << "p = " << p << endl;
        x=xupper; 	
	delta=0.1; counter=0;
        pinfi=xcupper*ER_g_A(xcupper,k);
        while (delta>0.0000001){
   	y=p*SF_g_A(ER_g_A(x,k)*p,lembda,mm);
        delta=x-y;
        x=y;
   	counter=counter+1;
   	}
        xcupper=x; 
       cout << "xcupper" << xcupper << endl;      
        if (xcupper<=0.0001) break;
      
 }
cout<< "p-infinity = " << pinfi << endl;  
return (p+deltap);
 } 


void plot_graph_xc_ERER(double k)
{ double xupper=0.99,xlower=0.1;
  double xcupper=1., xclower=0.;
  double delta=0.1, deltaderi=0.1,deltaxc=1.;
  double p=0.9;
  double x,y,deri;
  int counter=0;
  double deltap=0.00001;
  double pinfi=0.;
  
  ofstream fout;
 fout.open("./theoERER.dat"); 


  while (1)
  {
   //pinfi=xcupper*ER_g_Adi(xcupper,k);
   p=p-deltap;
   cout << "p = " << p << endl;
        x=xupper; 	
	delta=0.1; counter=0;
        while (1){
   	y=p*ER_g_A(ER_g_A(x,k)*p,k);
        delta=x-y;
        if (delta<0.000000001) break;
        x=y;
   	counter=counter+1;
   	}
        xcupper=x; 
        fout << p << "  " << x*ER_g_A(x,k) <<endl;   
        if (xcupper<=0.0001) break;
      
 }

 }

void plot_graph_xc_SFSF(double lembda, int mm)
{ double xupper=0.99;
  double xcupper=0.;
  double delta=0.1, deltaderi=0.1,deltaxc;
  double p=0.8;
  double x,y,deri;
  int counter=0;
  double pinfi=0;
  double deltap=0.001;
 
   ofstream fout;
 fout.open("./theoSFSF.dat"); 

   while (1)
  {pinfi=xcupper*SF_g_A(xcupper,lembda,mm);
   p=p-deltap;
   cout << "p = " << p << endl;
        x=xupper; 	
	delta=0.1; 
        while (delta>0.0000001){
        y=p*SF_g_A(SF_g_A(x,lembda,mm)*p,lembda,mm);
        delta=x-y;
        x=y;
   	counter=counter+1;
   	}
        xcupper=x; 
        fout << p << "  " << x*SF_g_A(x,lembda,mm)<<endl;    
        if (xcupper<=0.0001) break;
  
 }

 }

void plot_graph_xc_ERSF(double k,double lembda, int mm)
{ double xupper=0.99,xlower=0.1;
  double xcupper=0., xclower=0.;
  double delta=0.1, deltaderi=0.1,deltaxc=1.;
  double p=0.8;
  double x,y,deri;
  int counter=0;
  double deltap=0.001;
  double pinfi;
 
  ofstream fout;
 fout.open("./theoERSF.dat"); 

  while (1)
  {p=p-deltap;
   cout << "p = " << p << endl;
        x=xupper; 	
	delta=0.1; counter=0;
        pinfi=xcupper*SF_g_A(xcupper,lembda,mm);
        while (delta>0.0000001){
   	y=p*ER_g_A(SF_g_A(x,lembda,mm)*p,k);
        delta=x-y;
        x=y;
   	counter=counter+1;
   	}
        xcupper=x; 
        fout << p << "  " << x*SF_g_A(x,lembda,mm)<<endl;     
        if (xcupper<=0.0001) break;
      
 }
 }


void plot_graph_xc_SFER(double k,double lembda, int mm)
{ double xupper=0.99,xlower=0.1;
  double xcupper=0., xclower=0.;
  double delta=0.1, deltaderi=0.1,deltaxc=1.;
  double p=0.8;
  double x,y,deri;
  int counter=0;
  double deltap=0.001;
  double pinfi;
 
  ofstream fout;
 fout.open("./theoSFER.dat"); 

  while (1)
  {p=p-deltap;
   cout << "p = " << p << endl;
        x=xupper; 	
	delta=0.1; counter=0;
        pinfi=xcupper*ER_g_A(xcupper,k);
        while (delta>0.0000001){
   	y=p*SF_g_A(ER_g_A(x,k)*p,lembda,mm);
        delta=x-y;
        x=y;
   	counter=counter+1;
   	}
        xcupper=x; 
        fout << p << "  " << x*ER_g_A(x,k)<<endl;     
        if (xcupper<=0.0001) break;
      
 }
 }

void plot_betavsp_graph_xc_ERER(double k)
{ double xupper=0.99,xlower=0.1;
  double xcupper=1., xclower=0.;
  double delta=0.1, deltaderi=0.1,deltaxc=1.;
  double p=0.61395,pc=0.;
  double x,y,deri;
  int counter=0;
  double deltap=0.000001;
  double pinfi=0;
  double beta=0.;
  
 ofstream fout;
 fout.open("./theoBETAvspERER.dat"); 


  while (1)
  {
   pinfi=xcupper*ER_g_A(xcupper,k);
   p=p-deltap;
   cout << "p = " << p << endl;
        x=xupper; 	
	delta=0.1; counter=0;
        while (1){
   	y=p*ER_g_A(ER_g_A(x,k)*p,k);
        delta=x-y;
        if (delta<0.000000001) break;
        x=y;
   	counter=counter+1;
   	}
        xcupper=x; 
        //beta=log(x*ER_g_A(x,k) - 0.314497)/log(p-0.61385);
        //fout << p << "  " << beta <<endl;   
        fout << p-0.613852 << "  " << x*ER_g_A(x,k) - 0.314505 << endl;
        if (xcupper<=0.0001) break;
      
 }

   //pinfi=xcupper*ER_g_A(xcupper,k);
   pc=p+deltap;

cout << "pinfi" << pinfi << endl;
cout << "pc" << pc << endl;

 }



void plot_betavsp_graph_xc_SFSF(double lembda, int mm)
{ double xupper=0.99,xlower=0.1;
  double xcupper=1., xclower=0.;
  double delta=0.1, deltaderi=0.1,deltaxc=1.;
  double p=0.8,pc=0.;
  double x,y,deri;
  int counter=0;
  double deltap=0.001;
  double pinfi=0;
  double beta=0.;
  
 ofstream fout;
 fout.open("./theoBETAvspSFSF.dat"); 


  while (1)
  {pinfi=xcupper*SF_g_A(xcupper,lembda,mm);
   p=p-deltap;
   cout << "p = " << p << endl;
        x=xupper; 	
	delta=0.1; counter=0;
        while (delta>0.0000001){
	y=p*SF_g_A(SF_g_A(x,lembda,mm)*p,lembda,mm);
        delta=x-y;
        x=y;
   	counter=counter+1;
   	}
        xcupper=x; 
        beta=log(x*SF_g_A(x,lembda,mm) - 0.407091)/log(p-0.753);
        fout << p << "  " << beta <<endl;   
        if (xcupper<=0.0001) break;
      
 }

   //pinfi=xcupper*ER_g_A(xcupper,k);
   pc=p+deltap;

cout << "pinfi" << pinfi << endl;
cout << "pc" << pc << endl;

 }


void plot_betavsp_graph_xc_ERSF(double k, double lembda, int mm)
{ double xupper=0.99,xlower=0.1;
  double xcupper=1., xclower=0.;
  double delta=0.1, deltaderi=0.1,deltaxc=1.;
  double p=0.8,pc=0.;
  double x,y,deri;
  int counter=0;
  double deltap=0.001;
  double pinfi=0;
  double beta=0.;
  
 ofstream fout;
 fout.open("./theoBETAvspERSF.dat"); 


  while (1)
  {pinfi=xcupper*SF_g_A(xcupper,lembda,mm);
   p=p-deltap;
   cout << "p = " << p << endl;
        x=xupper; 	
	delta=0.1; counter=0;
        while (delta>0.0000001){
        y=p*ER_g_A(SF_g_A(x,lembda,mm)*p,k);
        delta=x-y;
        x=y;
   	counter=counter+1;
   	}
        xcupper=x; 
        beta=log(x*SF_g_A(x,lembda,mm) - 0.346486)/log(p-0.681);
        fout << p << "  " << beta <<endl;   
        if (xcupper<=0.0001) break;
      
 }

   //pinfi=xcupper*ER_g_A(xcupper,k);
   pc=p+deltap;

cout << "pinfi" << pinfi << endl;
cout << "pc" << pc << endl;

 }



void plot_betavsp_graph_xc_SFER(double k, double lembda, int mm)
{ double xupper=0.99,xlower=0.1;
  double xcupper=1., xclower=0.;
  double delta=0.1, deltaderi=0.1,deltaxc=1.;
  double p=0.8,pc=0.;
  double x,y,deri;
  int counter=0;
  double deltap=0.001;
  double pinfi=0;
  double beta=0.;
  
 ofstream fout;
 fout.open("./theoBETAvspSFER.dat"); 


  while (1)
  {pinfi=xcupper*ER_g_A(xcupper,k);
   p=p-deltap;
   cout << "p = " << p << endl;
        x=xupper; 	
	delta=0.1; counter=0;
        while (delta>0.0000001){
        y=p*SF_g_A(ER_g_A(x,k)*p,lembda,mm);
        delta=x-y;
        x=y;
   	counter=counter+1;
   	}
        xcupper=x; 
        beta=log(x*ER_g_A(x,k) - 0.346486)/log(p-0.681);
        fout << p << "  " << beta <<endl;   
        if (xcupper<=0.0001) break;
      
 }

   //pinfi=xcupper*ER_g_A(xcupper,k);
   pc=p+deltap;

cout << "pinfi" << pinfi << endl;
cout << "pc" << pc << endl;

 }















  

#endif
