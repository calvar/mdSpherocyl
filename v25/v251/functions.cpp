#include "functions.hpp"
#include <cstdio>
#include <cstring>
#include <sys/stat.h>
//extern "C"{
#include "jacobi.h"
//}

bool r_input(int& Nsteps, int& pre_steps, int& iprint, int& signal, part_num& NP, Dynvec<double>& length, 
	     Dynvec<double>& diameter, Dynvec<double>& mass, Dynvec<Matrix<double> >& dip, double& dens, double& temp, 
	     double& Dt, double& sr_cut, int& K_max, double& alpha, double& bar_w, double& ArLx, double& ArLy, 
	     int& corrf, double& Q, double& gamma, int& tcfd, int& tcfv, double B[3], int& cl_fr,
	     double& min_s){
  std::string line;
  ifstream In("input.txt");
  if(! In){
    ofstream Err("Error.dat");
    Err << "Couldn't open " << "input.txt" << endl;
    Err.close();
    In.close();
    return false;
  }
  int n_comp;
  Dynvec<int> n_part;
  In >> line >> line >> line >> Nsteps;
  In >> line >> line >> line >> pre_steps;
  In >> line >> signal;
  In >> line >> line >> line >> n_comp;
  n_part.set_size(n_comp);
  In >> line >> line >> line;
  for(int i = 0; i < n_comp; ++i)
    In >> n_part[i];
  length.set_size(n_comp);
  In >> line;
  for(int i = 0; i < n_comp; ++i)
    In >> length[i];
  diameter.set_size(n_comp);
  In >> line;
  for(int i = 0; i < n_comp; ++i)
    In >> diameter[i];
  mass.set_size(n_comp);
  In >> line;
  for(int i = 0; i < n_comp; ++i)
    In >> mass[i];
  In >> line >> line;
  for(int i = 0; i < 3; ++i)
    In >> B[i];    
  In >> line >> dens;
  In >> line >> temp;
  In >> line >> Q;
  In >> line >> gamma;
  In >> line >> line >> line >> line >> line;
  In >> ArLx >> ArLy; 
  In >> line >> line >> line >> Dt;
  In >> line >> line >> sr_cut;
  In >> line >> line >> line >> K_max;
  In >> line >> line >> alpha;
  In >> line >> line >> line >> bar_w;
  In >> line >> line >> line >> corrf;
  In >> line >> line >> tcfd;
  In >> line >> line >> tcfv;
  In >> line >> line >> cl_fr;
  In >> line >> line >> min_s;
  In >> line >> line >> iprint;
  //read dipole configurations
  dip.set_size(n_comp);
  In >> line;
  In >> line >> line >> line >> line >> line;
  int ndp;
  for(int i = 0; i < n_comp; ++i){
    In >> line >> ndp;
    dip[i].set_size(ndp, 5);
    for(int j = 0; j < ndp; ++j){
      In >> dip[i][j][0];
      In >> line;
      for(int a = 1; a < 4; ++a)
	In >> dip[i][j][a];
      In >> line;
      In >> dip[i][j][4];
    }
  }
  In.close();
  if(! NP.init(n_comp, n_part))
    return false;
  return true;
}

double spherocyl_dist(const Spherocyl& A, const Spherocyl& B, const double L[3], const double rAB[3], double dAB[3],
		      double st[2]){
  double oA[3], oB[3], ldA[2], ldB[2], doto, dotrA, dotrB;
  A.get_ori(oA); A.get_shape(ldA);
  B.get_ori(oB); B.get_shape(ldB);
  dot(oA, oB, doto);
  dot(rAB, oA, dotrA);
  dot(rAB, oB, dotrB);
//   //
//   cout<<"r("<<rAB[0]<<","<<rAB[1]<<","<<rAB[2]
//       <<") par("<<oA[0]*dotrA<<","<<oA[1]*dotrA<<","<<oA[2]*dotrA<<") per("
//       <<rAB[0]-oA[0]*dotrA<<","<<rAB[1]-oA[1]*dotrA<<","<<rAB[2]-oA[2]*dotrA<<")\n";
//   //
  
  double a = ldA[0] * ldA[0];
  double b = doto * ldA[0] * ldB[0];
  double c = ldB[0] * ldB[0];
  double d = dotrA*ldA[0] - a/2 + b/2;
  double e = dotrB*ldB[0] + c/2 - b/2;
  
  double det = a*c - b*b;
  double s = b*e - c*d;
  double t = a*e - b*d;
//   //
//   cout<<"doto="<<doto<<" dotrA="<<dotrA<<" dotrB="<<dotrB<<endl;
//   cout<<"a="<<a<<" b="<<b<<" c="<<c<<" d="<<d<<" e="<<e<<endl;
//   //
  double dis[3], dsq;
  if(det > 1.e-12){ // Not parallel
    double s0, t0, dsq0;
    s0 = s/det; t0 = t/det; 
    for(int i = 0; i < 3; ++i)
      dis[i] = rAB[i] + (s0-0.5)*ldA[0]*oA[i] - (t0-0.5)*ldB[0]*oB[i];
    dot(dis, dis, dsq0);
//     //
//     cout<<"s "<<s0<<" t "<<t0<<" det "<<det<<" "<<dsq0<<endl;
//     //

    if(dsq0 < 1.e-12){ //if lines intersect (distance between line segment and end point)
      double rA[3], rB[3];
      A.get_pos(rA); B.get_pos(rB);
      double pp[4][3];
      for(int a = 0; a < 3; ++a){
	pp[0][a] = rA[a] + ldA[0]*oA[a]/2; //s=1
	pp[1][a] = rA[a] - ldA[0]*oA[a]/2; //s=0
	pp[2][a] = rB[a] + ldB[0]*oB[a]/2; //t=1
	pp[3][a] = rB[a] - ldB[0]*oB[a]/2; //t=0
      }
      double dd[4], d[4][3];
      for(int a = 0; a < 3; ++a){
	d[0][a] = pp[0][a] - pp[2][a]; d[0][a] -= L[a] * floor(d[0][a] / L[a] + 0.5); 
	d[1][a] = pp[0][a] - pp[3][a]; d[1][a] -= L[a] * floor(d[1][a] / L[a] + 0.5);
	d[2][a] = pp[1][a] - pp[2][a]; d[2][a] -= L[a] * floor(d[2][a] / L[a] + 0.5);
	d[3][a] = pp[1][a] - pp[3][a]; d[3][a] -= L[a] * floor(d[3][a] / L[a] + 0.5);
	//rAB[a] - (ldA[0]*oA[a] - ldB[0]*oB[a]) / 2;
      }
      dot(d[0], d[0], dd[0]);
      dot(d[1], d[1], dd[1]);
      dot(d[2], d[2], dd[2]);
      dot(d[3], d[3], dd[3]);
      int im = 0;
      double mind = 1.e6; 
      for(int b = 0 ; b < 4; ++b){
	if(dd[b] < mind){ im = b; mind = dd[b]; }
      }
      double vec[3];
      switch(im){
	case 0:
	  t = 1;
	  vec[0] = pp[1][0]-pp[2][0]; vec[0] -= L[0] * floor(vec[0] / L[0] + 0.5); 
	  vec[1] = pp[1][1]-pp[2][1]; vec[1] -= L[1] * floor(vec[1] / L[1] + 0.5);
	  vec[2] = pp[1][2]-pp[2][2]; vec[2] -= L[2] * floor(vec[2] / L[2] + 0.5);
	  dot(vec, oA, s);
	  s /= -ldA[0];
	  if(s >= 0 && s <= 1){
	    for(int i = 0; i < 3; ++i)
	      dAB[i] = rAB[i] + (s-0.5)*ldA[0]*oA[i] - (t-0.5)*ldB[0]*oB[i];
	    dot(dAB, dAB, dsq);
	  }else{
	    s = 1;
	    vec[0] = pp[3][0]-pp[0][0]; vec[0] -= L[0] * floor(vec[0] / L[0] + 0.5); 
	    vec[1] = pp[3][1]-pp[0][1]; vec[1] -= L[1] * floor(vec[1] / L[1] + 0.5);
	    vec[2] = pp[3][2]-pp[0][2]; vec[2] -= L[2] * floor(vec[2] / L[2] + 0.5);
	    dot(vec, oB, t);
	    t /= -ldB[0];
	    if(t >= 0 && t <= 1){
	      for(int i = 0; i < 3; ++i)
	      dAB[i] = rAB[i] + (s-0.5)*ldA[0]*oA[i] - (t-0.5)*ldB[0]*oB[i];
	      dot(dAB, dAB, dsq);
	    }else{
	      t = 1;
	      for(int i = 0; i < 3; ++i)
		dAB[i] = d[0][i];
	      dot(dAB, dAB, dsq);
	    }
	  }                                  //cout<<s<<" "<<t<<endl;
	  break;
	case 1:
	  t = 0;
	  vec[0] = pp[1][0]-pp[3][0]; vec[0] -= L[0] * floor(vec[0] / L[0] + 0.5);
	  vec[1] = pp[1][1]-pp[3][1]; vec[1] -= L[1] * floor(vec[1] / L[1] + 0.5);
	  vec[2] = pp[1][2]-pp[3][2]; vec[2] -= L[2] * floor(vec[2] / L[2] + 0.5);
	  dot(vec, oA, s);
	  s /= -ldA[0];
	  if(s >= 0 && s <= 1){
	    for(int i = 0; i < 3; ++i)
	      dAB[i] = rAB[i] + (s-0.5)*ldA[0]*oA[i] - (t-0.5)*ldB[0]*oB[i];
	    dot(dAB, dAB, dsq);
	  }else{
	    s = 1;
	    vec[0] = pp[3][0]-pp[0][0]; vec[0] -= L[0] * floor(vec[0] / L[0] + 0.5);
	    vec[1] = pp[3][1]-pp[0][1]; vec[1] -= L[1] * floor(vec[1] / L[1] + 0.5);
	    vec[2] = pp[3][2]-pp[0][2]; vec[2] -= L[2] * floor(vec[2] / L[2] + 0.5);
	    dot(vec, oB, t);
	    t /= -ldB[0];
	    if(t >= 0 && t <= 1){
	      for(int i = 0; i < 3; ++i)
	      dAB[i] = rAB[i] + (s-0.5)*ldA[0]*oA[i] - (t-0.5)*ldB[0]*oB[i];
	      dot(dAB, dAB, dsq);
	    }else{
	      t = 0;
	      for(int i = 0; i < 3; ++i)
		dAB[i] = d[1][i];
	      dot(dAB, dAB, dsq);
	    }
	  }
	  break;
	case 2:
	  t = 1;
	  vec[0] = pp[1][0]-pp[2][0]; vec[0] -= L[0] * floor(vec[0] / L[0] + 0.5);
	  vec[1] = pp[1][1]-pp[2][1]; vec[1] -= L[1] * floor(vec[1] / L[1] + 0.5);
	  vec[2] = pp[1][2]-pp[2][2]; vec[2] -= L[2] * floor(vec[2] / L[2] + 0.5);
	  dot(vec, oA, s);
	  s /= -ldA[0];
	  if(s >= 0 && s <= 1){
	    for(int i = 0; i < 3; ++i)
	      dAB[i] = rAB[i] + (s-0.5)*ldA[0]*oA[i] - (t-0.5)*ldB[0]*oB[i];
	    dot(dAB, dAB, dsq);
	  }else{
	    s = 0;
	    vec[0] = pp[3][0]-pp[1][0]; vec[0] -= L[0] * floor(vec[0] / L[0] + 0.5);
	    vec[1] = pp[3][1]-pp[1][1]; vec[1] -= L[1] * floor(vec[1] / L[1] + 0.5);
	    vec[2] = pp[3][2]-pp[1][2]; vec[2] -= L[2] * floor(vec[2] / L[2] + 0.5);
	    dot(vec, oB, t);
	    t /= -ldB[0];
	    if(t >= 0 && t <= 1){
	      for(int i = 0; i < 3; ++i)
	      dAB[i] = rAB[i] + (s-0.5)*ldA[0]*oA[i] - (t-0.5)*ldB[0]*oB[i];
	      dot(dAB, dAB, dsq);
	    }else{
	      t = 1;
	      for(int i = 0; i < 3; ++i)
		dAB[i] = d[2][i];
	      dot(dAB, dAB, dsq);
	    }
	  }
	  break;
	case 3:
	  t = 0;
	  vec[0] = pp[1][0]-pp[3][0]; vec[0] -= L[0] * floor(vec[0] / L[0] + 0.5);
	  vec[1] = pp[1][1]-pp[3][1]; vec[1] -= L[1] * floor(vec[1] / L[1] + 0.5);
	  vec[2] = pp[1][2]-pp[3][2]; vec[2] -= L[2] * floor(vec[2] / L[2] + 0.5);
	  dot(vec, oA, s);
	  s /= -ldA[0];
	  if(s >= 0 && s <= 1){
	    for(int i = 0; i < 3; ++i)
	      dAB[i] = rAB[i] + (s-0.5)*ldA[0]*oA[i] - (t-0.5)*ldB[0]*oB[i];
	    dot(dAB, dAB, dsq);
	  }else{
	    s = 0;
	    vec[0] = pp[3][0]-pp[1][0]; vec[0] -= L[0] * floor(vec[0] / L[0] + 0.5);
	    vec[1] = pp[3][1]-pp[1][1]; vec[1] -= L[1] * floor(vec[1] / L[1] + 0.5);
	    vec[2] = pp[3][2]-pp[1][2]; vec[2] -= L[2] * floor(vec[2] / L[2] + 0.5);
	    dot(vec, oB, t);
	    t /= -ldB[0];
	    if(t >= 0 && t <= 1){
	      for(int i = 0; i < 3; ++i)
	      dAB[i] = rAB[i] + (s-0.5)*ldA[0]*oA[i] - (t-0.5)*ldB[0]*oB[i];
	      dot(dAB, dAB, dsq);
	    }else{
	      t = 0;
	      for(int i = 0; i < 3; ++i)
		dAB[i] = d[3][i];
	      dot(dAB, dAB, dsq);
	    }
	  }
      }
      
    }else{ //if lines do not intersect
      //s|1|2|3|
      //^|4|5|6|
      // |7|8|9|
      //    > t
      int region = 5;
      if(s >= 0){
	if(s <= det){
	  if(t >= 0){
	    if(t <= det)
	      region = 5;
	    else
	      region = 6;
	  }else{
	    region = 4;
	  }
	}else{
	  if(t >= 0){
	    if(t <= det)
	      region = 2;
	    else
	      region = 3;
	  }else{
	    region = 1;
	  }
	}
      }else{
	if(t >= 0){
	  if(t <= det)
	    region = 8;
	  else
	    region = 9;
	}else{
	  region = 7;
	}
      }
      switch(region){
	case 2:
	  s = 1; t = (e + b) / c;
	  if(t > 1) t = 1;
	  else if(t < 0) t = 0;
	  for(int i = 0; i < 3; ++i)
	    dAB[i] = rAB[i] + (s-0.5)*ldA[0]*oA[i] - (t-0.5)*ldB[0]*oB[i];
	  dot(dAB, dAB, dsq);
	  break;
	case 8:
	  s = 0; t = e / c;
	  if(t > 1) t = 1;
	  else if(t < 0) t = 0;
	  for(int i = 0; i < 3; ++i)
	    dAB[i] = rAB[i] + (s-0.5)*ldA[0]*oA[i] - (t-0.5)*ldB[0]*oB[i];
	  dot(dAB, dAB, dsq);
	  break;
	case 6:
	  t = 1; s = -(d - b) / a;
	  if(s > 1) s = 1;
	  else if(s < 0) s = 0;
	  for(int i = 0; i < 3; ++i)
	    dAB[i] = rAB[i] + (s-0.5)*ldA[0]*oA[i] - (t-0.5)*ldB[0]*oB[i];
	  dot(dAB, dAB, dsq);
	  break;
	case 4:
	  t = 0; s = -d / a;
	  if(s > 1) s = 1;
	  else if(s < 0) s = 0;
	  for(int i = 0; i < 3; ++i)
	    dAB[i] = rAB[i] + (s-0.5)*ldA[0]*oA[i] - (t-0.5)*ldB[0]*oB[i];
	  dot(dAB, dAB, dsq);
	  break;
	case 1:
	  dsq = dsq0 = 1.e6;
	  s = 1; t = (e + b) / c;
	  if(t < 0) t = 0;
	  else if(t > 1) t = 1;
	  for(int i = 0; i < 3; ++i)
	    dAB[i] = rAB[i] + (s-0.5)*ldA[0]*oA[i] - (t-0.5)*ldB[0]*oB[i];
	  dot(dAB, dAB, dsq);
	  
	  t0 = 0; s0 = -d / a;
	  if(s0 < 0) s0 = 0;
	  else if(s0 > 1) s0 = 1;
	  for(int i = 0; i < 3; ++i)
	    dis[i] = rAB[i] + (s0-0.5)*ldA[0]*oA[i] - (t0-0.5)*ldB[0]*oB[i];
	  dot(dis, dis, dsq0);
	  
	  if(dsq0 < dsq){
	    for(int a = 0; a < 3; ++a)
	      dAB[a] = dis[a];
	    s = s0; t = t0; dsq = dsq0; 
	  }
	  break;
	case 3:
	  dsq = dsq0 = 1.e6;
	  s = 1; t = (e + b) / c;
	  if(t < 0) t = 0;
	  else if(t > 1) t = 1;
	  for(int i = 0; i < 3; ++i)
	    dAB[i] = rAB[i] + (s-0.5)*ldA[0]*oA[i] - (t-0.5)*ldB[0]*oB[i];
	  dot(dAB, dAB, dsq);
	 
	  t0 = 1; s0 = -(d - b) / a;
	  if(s0 < 0) s0 = 0;
	  else if(s0 > 1) s0 = 1;
	  for(int i = 0; i < 3; ++i)
	    dis[i] = rAB[i] + (s0-0.5)*ldA[0]*oA[i] - (t0-0.5)*ldB[0]*oB[i];
	  dot(dis, dis, dsq0);
	  
	  if(dsq0 < dsq){
	    for(int a = 0; a < 3; ++a)
	      dAB[a] = dis[a];
	    s = s0; t = t0; dsq = dsq0; 
	  }
	  break;
	case 7:
	  dsq = dsq0 = 1.e6;
	  t = 0; s = -d / a;
	  if(s < 0) s = 0;
	  else if(s > 1) s = 1;
	  for(int i = 0; i < 3; ++i)
	    dAB[i] = rAB[i] + (s-0.5)*ldA[0]*oA[i] - (t-0.5)*ldB[0]*oB[i];
	  dot(dAB, dAB, dsq);
	  
	  s0 = 0; t0 = e / c;
	  if(t0 < 0) t0 = 0;
	  else if(t0 > 1) t0 = 1;
	  for(int i = 0; i < 3; ++i)
	    dis[i] = rAB[i] + (s0-0.5)*ldA[0]*oA[i] - (t0-0.5)*ldB[0]*oB[i];
	  dot(dis, dis, dsq0);
	  
	  if(dsq0 < dsq){
	    for(int a = 0; a < 3; ++a)
	      dAB[a] = dis[a];
	    s = s0; t = t0; dsq = dsq0; 
	  }
	  break;
	case 9:
	  dsq = dsq0 = 1.e6;
	  t = 1; s = -(d - b) / a;
	  if(s < 0) s = 0;
	  else if(s > 1) s = 1;
	  for(int i = 0; i < 3; ++i)
	    dAB[i] = rAB[i] + (s-0.5)*ldA[0]*oA[i] - (t-0.5)*ldB[0]*oB[i];
	  dot(dAB, dAB, dsq);
	  
	  s0 = 0; t0 = e / c;
	  if(t0 < 0) t0 = 0;
	  else if(t0 > 1) t0 = 1;
	  for(int i = 0; i < 3; ++i)
	    dis[i] = rAB[i] + (s0-0.5)*ldA[0]*oA[i] - (t0-0.5)*ldB[0]*oB[i];
	  dot(dis, dis, dsq0);
	  
	  if(dsq0 < dsq){
	    for(int a = 0; a < 3; ++a)
	      dAB[a] = dis[a];
	    s = s0; t = t0; dsq = dsq0; 
	  }
	  break;
	case 5:
	  s /= det; t /= det; 
	  for(int i = 0; i < 3; ++i)
	    dAB[i] = rAB[i] + (s-0.5)*ldA[0]*oA[i] - (t-0.5)*ldB[0]*oB[i];
	  dot(dAB, dAB, dsq);
	  break;
	default:
	  ofstream Err("Error.dat");
	  Err << "distance. Case " << region << "not found." << endl;
	  Err.close();
      }
    }
    
  }else{ // if parallel
    double le = (ldA[0] + ldB[0]) / 2;
    if(fabs(dotrA) < le){ //no unique solution for minimum distance
      for(int i = 0; i < 3; ++i)
	dAB[i] = rAB[i] - dotrA * oA[i]; 
      dot(dAB, dAB, dsq);
      s = t = -1.e+6;
    }else{ //min dist between endpoints
      double d2[3], dsq0;
      if(doto > 0){ //head to tail or tail to head
	for(int i = 0; i < 3; ++i){
	  dAB[i] = rAB[i] + ldA[0]*oA[i]/2 + ldB[0]*oB[i]/2; //headA tailB (+/+)
	  d2[i] = rAB[i] - ldA[0]*oA[i]/2 - ldB[0]*oB[i]/2; //tailA headB (-/-)
	}
	dot(dAB, dAB, dsq);
	dot(d2, d2, dsq0);
	s = t = -1.e+10; //(+/+)
	if(dsq0 < dsq){
	  for(int a = 0; a < 3; ++a)
	    dAB[a] = d2[a];
	  dsq = dsq0;
	  s = t = -2.e+10; //(-/-)
	}
      }else{ //head to head or tail to tail
	for(int i = 0; i < 3; ++i){
	  dAB[i] = rAB[i] + ldA[0]*oA[i]/2 - ldB[0]*oB[i]/2; //headA headB (+/-)
	  d2[i] = rAB[i] - ldA[0]*oA[i]/2 + ldB[0]*oB[i]/2; //tailA tailB (-/+)
	}
	dot(dAB, dAB, dsq);
	dot(d2, d2, dsq0);
	s = t = -3.e+10; //(+/-)
	if(dsq0 < dsq){
	  for(int a = 0; a < 3; ++a)
	    dAB[a] = d2[a];
	  dsq = dsq0;
	  s = t = -4.e+10; //(-/+)
	}
      }
    }
  }
  st[0] = s; st[1] = t;
  return dsq;
}

int ini_pos(Dynvec<Dynvec<Spherocyl> >& part, const double L[3], const part_num& NP){
  int NC = NP.get_Ncomp();                                  
  double len[NC], dia[NC], sh[2];
  for(int n = 0; n < NC; ++n){
    part[n][0].get_shape(sh);
    len[n] = sh[0];
    dia[n] = sh[1];
  }
  
  double r[3];
  Matrix<double> o(3,3,0);
  o[0][0] = o[1][1] = o[2][2] = 1;  
  double sq = sqrt(3.) / 2.;
  double sqt = 1./ (2 * sqrt(3.));
  double sqtt = sqrt(6.) / 3;
  for(int n = 0; n < NC; ++n){
    int NN = NP.get_Npart(n);
    int num1 = static_cast<int>(max(floor(L[0]/(dia[n])), 1.));
    int num2 = static_cast<int>(max(floor(L[1]/(dia[n]*sq)), 1.));
    int num3 = static_cast<int>(max(floor(L[2]/(dia[n]*sqtt+len[n])), 1.));
    int cont = 0;
    while(cont < NN){
      for(int nz = 0; nz < num3; ++nz){
	for(int ny = 0; ny < num2; ++ny){
	  for(int nx = 0; nx < num1; ++nx){
	    if(cont < NN){
	      if(fmod(static_cast<double>(nz),2.) == 0){
		if(fmod(static_cast<double>(ny),2.) == 0){
		  r[0] = dia[n] * static_cast<double>(nx) * 1.0001;
		}else{
		  r[0] = dia[n] * (static_cast<double>(nx) + 0.5) * 1.0001;
		}
		r[1] = dia[n] * static_cast<double>(ny) * sq * 1.0001;
	      }else{
		if(fmod(static_cast<double>(ny),2.) != 0){
		  r[0] = dia[n] * static_cast<double>(nx) * 1.0001;
		}else{
		  r[0] = dia[n] * (static_cast<double>(nx) + 0.5) * 1.0001;
		}
		r[1] = dia[n] * (static_cast<double>(ny) * sq + sqt) * 1.0001;
	      }
	      r[2] = (dia[n] * sqtt + len[n]) * static_cast<double>(nz) * 1.0001;
	      
	      r[0] -= L[0] * floor(r[0] / L[0] + 0.5);
	      r[1] -= L[1] * floor(r[1] / L[1] + 0.5);
	      r[2] -= L[2] * floor(r[2] / L[2] + 0.5);

	      part[n][cont].set_pos(r);
	      part[n][cont].set_omat(o);
	      
	      if(n > 0){
		int con = 0;
		double rj[3], rij[3];
		for(int m = 0; m < n; ++m){
		  int MM = NP.get_Npart(m);
		  double dd = (dia[n] + dia[m]) / 2;
		  double rsq = dd * dd;
		  bool test = false;
		  for(int j = 0; j < MM; ++j){
		    part[m][j].get_pos(rj);
		    
		    for(int a = 0; a < 3; ++a){
		      rij[a] = r[a] - rj[a];
		      rij[a] -= L[a] * floor(rij[a] / L[a] + 0.5);
		    }
		    double st[2], dij[3];
		    double dsq = spherocyl_dist(part[n][cont], part[m][j], L, rij, dij, st);
// 		    //
// 		    cout<<"r2:"<<rsq<<" d2:"<<dsq<<" n("<<nx<<","<<ny<<","<<nz<<")\n";
// 		    //
		    if(dsq < rsq){
		      test = true;
		      break;
		    }
		    
		  }
		  if(test){
		    break;
		  }
		  if(m == n-1){
		    cont++;
		  }
		  con++;
		}
	      }else{
		cont++;
	      }	  
	    }

	  }
	}
      }
    }
  }
 
  //print
  ofstream InPos("iniPos.dat");
  if(! InPos){
    ofstream Err("Error.dat", ios::app);
    Err << "Couldn't open " << "iniPos.dat" << endl;
    Err.close();
    InPos.close();
    return -1;
  }
  for(int n = 0; n < NC; ++n){
    int NN = NP.get_Npart(n);
    for(int i = 0; i < NN; i++){
      part[n][i].get_pos(r);
      InPos << r[0] << " " << r[1] << " " << r[2] 
	    << " 0" << " 0" << " 0" << " ";
      for(int a = 0; a < 3; ++a){
	for(int b = 0; b < 3; ++b)
	  InPos << o[a][b] << " ";
      }
      InPos << " 0" << " 0" << " 0" << " 1" << endl; 
    }
  }
  InPos << L[0] << " " << L[1] << " " << L[2] << endl;
  InPos.close();

  //test for overlaps
  double ri[3], rj[3], rij[3];
  //intra
  for(int n = 0; n < NC; ++n){
    int NN = NP.get_Npart(n);
    double dd = dia[n];
    double rsq = dd * dd;
    for(int i = 0; i < NN-1; ++i){
      part[n][i].get_pos(ri);
      for(int j = i+1; j < NN; ++j){
	part[n][j].get_pos(rj);
	
	for(int a = 0; a < 3; ++a){
	  rij[a] = ri[a] - rj[a];
	  rij[a] -= L[a] * floor(rij[a] / L[a] + 0.5);
	}
	double st[2], dij[3];
	double dsq = spherocyl_dist(part[n][i], part[n][j], L, rij, dij, st);
	
	if(dsq < rsq){
	  ofstream Err("Error.dat", ios::app);
	  Err << "initial conf. overlap\n"
	      << "intra;"<< n << "| " << i << " " << j << " d:" << sqrt(dsq) 
	      << " " << dia[n] << endl;
	  Err.close();
	  return 1;
	}
      }
    }
  }
  //inter
  int con = 0;
  for(int m = 0; m < NC-1; ++m){
    int MM = NP.get_Npart(m);
    for(int n = m+1; n < NC; ++n){
      int NN = NP.get_Npart(n);
      double dd = (dia[m] + dia[n]) / 2;
      double rsq = dd * dd;
      for(int i = 0; i < MM; ++i){
	part[m][i].get_pos(ri);
	for(int j = 0; j < NN; ++j){
	  part[n][j].get_pos(rj);
	  
	  for(int a = 0; a < 3; ++a){
	    rij[a] = ri[a] - rj[a];
	    rij[a] -= L[a] * floor(rij[a] / L[a] + 0.5);
	  }
	  double st[2], dij[3];
	  double dsq = spherocyl_dist(part[m][i], part[n][j], L, rij, dij, st);
	 
	  if(dsq < rsq){
	    ofstream Err("Error.dat", ios::app);
	    Err << "initial cof. overlap\n"
		<< "inter;"<< m << "," << n << "| " << i << " " << j 
		<< " d:" << sqrt(dsq) << " " << dd << endl;
	    Err.close();
	    return 1;
	  }
	}
      }
      con++;
    }
  }
  return 0;
}

void ini_mom(Dynvec<Dynvec<Spherocyl> >& part, const part_num& NP, double temp,
	     const gsl_rng* ran){
  int NC = NP.get_Ncomp();
  int Ntot = NP.get_Tpart();
  
  //for translation velocities************************************************
  Dynvec<Dynvec<double> > rnd[3];
  rnd[0].set_size(NC); rnd[1].set_size(NC); rnd[2].set_size(NC); 
  for(int n = 0; n < NC; ++n){
    for(int j = 0; j < 3; ++j)
      rnd[j][n].set_size(NP.get_Npart(n)+1);
  }
  for(int n = 0; n < NC; ++n){
    int NN = NP.get_Npart(n);
    for(int i = 0; i < NN; ++i){
      for(int j = 0; j < 3; ++j)
	rnd[j][n][i] = gsl_rng_uniform(ran);
    }
  }
  for(int n = 0; n < NC; ++n){
    int NN = NP.get_Npart(n);
    double mtemp = temp * part[n][0].get_mass();
    for(int i = 0; i < NN; i += 2){
      if((fmod(static_cast<float>(NN),2) != 0.) && (i == NN-1)){
	for(int j = 0; j < 3; j++){
	  double rr = gsl_rng_uniform(ran);
	  gauss(mtemp, rr, rnd[j][n][i]);
	}
      }else{
	for(int j = 0; j < 3; j++){
	  gauss(mtemp, rnd[j][n][i], rnd[j][n][i+1]);
	}
      }
    }
  }

  double P[3] = {0., 0., 0.};
  double pp[3][Ntot+1];
  
  int cnt = 0;
  for(int n = 0; n < NC; ++n){
    int NN = NP.get_Npart(n);
    for(int i = 0; i < NN; ++i){
      for(int j = 0; j < 3; j++){
	pp[j][cnt] = rnd[j][n][cnt];
	P[j] += pp[j][cnt];
      }
      cnt++;
    }
  }

  double pr = 0.1;
  
//   //
//   cout<<"begin linear X..\n";
//   //
  
 //  // X part/////////////////////////////////////////////////
  int II = 0;/*, nn=0, ii=0;*/
  long scont = 0;
  while(fabs(P[0]) > pr){ //NOTE: in fabsP > x, if the value of x is very 
                           // small the program will have an error. 
                           //(what is small depends on Ntot) 
                           // for N=256 x=0.1 is fine.
    while(P[0] > pr){
      if(pp[0][II] > 0){ 
	pp[0][II] = -pp[0][II];
// 	rev_II(NP, nn, ii, II);
// 	rnd[0][nn][ii] *= -1;
	P[0] += 2 * pp[0][II];
      }
      II++;
      if(II == Ntot - 1) II = 0;
    }
    while(P[0] < -pr){
      if(pp[0][II] < 0){ 
	pp[0][II] = -pp[0][II];
// 	rev_II(NP, nn, ii, II);
// 	rnd[0][nn][ii] *= -1;
	P[0] += 2 * pp[0][II];
      }
      II++;
      if(II == Ntot - 1) II = 0;
    }
    if(scont > 1e8) break;
    scont++;
  }
  P[0] -= pp[0][Ntot-1];
  pp[0][Ntot-1] = - P[0];
//   rev_II(NP, nn, ii, Ntot-1);
//   rnd[0][nn][ii] = pp[0][Ntot-1];
  P[0] = 0.;
  cnt = 0;
  for(int n = 0; n < NC; ++n){
    int NN = NP.get_Npart(n);
    for(int i = 0; i < NN; ++i){
      rnd[0][n][i] = pp[0][cnt];
      P[0] += rnd[0][n][i];
      cnt++;
    }
  }
  if(fabs(P[0]) > 1e-12){
    ofstream Err("Error.dat", ios::app);
    Err << "Px!=0 !!\n";
    Err.close();
  }
  
//   //
//   cout<<"begin linear Y..\n";
//   //
  
  //  // Y part/////////////////////////////////////////////////
  II = 0;
  scont = 0;
  while(fabs(P[1]) > pr){ //NOTE: in fabsP > x, if the value of x is very 
                           // small the program will have an error. 
                           //(what is small depends on Ntot) 
                           // for N=256 x=0.1 is fine.
    while(P[1] > pr){
      if(pp[1][II] > 0){ 
	pp[1][II] = -pp[1][II];
// 	rev_II(NP, nn, ii, II);
// 	rnd[1][nn][ii] *= -1;
	P[1] += 2 * pp[1][II];
      }
      II++;
      if(II == Ntot - 1) II = 0;
    }
    while(P[1] < -pr){
      if(pp[1][II] < 0){ 
	pp[1][II] = -pp[1][II];
// 	rev_II(NP, nn, ii, II);
// 	rnd[1][nn][ii] *= -1;
	P[1] += 2 * pp[1][II];
      }
      II++;
      if(II == Ntot - 1) II = 0;
    }
    if(scont > 1e8) break;
    scont++;
  }
  P[1] -= pp[1][Ntot-1];
  pp[1][Ntot-1] = - P[1];
//   rev_II(NP, nn, ii, Ntot-1);
//   rnd[1][nn][ii] = pp[1][Ntot-1];
  P[1] = 0.;
  cnt = 0;
  for(int n = 0; n < NC; ++n){
    int NN = NP.get_Npart(n);
    for(int i = 0; i < NN; ++i){
      rnd[1][n][i] = pp[1][cnt];
      P[1] += rnd[1][n][i];
      cnt++;
    }
  }
  if(fabs(P[1]) > 1e-12){
    ofstream Err("Error.dat", ios::app);
    Err << "Py!=0 !!\n";
    Err.close();
  }

//   //
//   cout<<"begin linear Z..\n";
//   //
  
  //  // Z part/////////////////////////////////////////////////
  II = 0;
  scont = 0;
  while(fabs(P[2]) > pr){ //NOTE: in fabsP > x, if the value of x is very 
                           // small the program will have an error. 
                           //(what is small depends on Ntot) 
                           // for N=256 x=0.1 is fine.
    while(P[2] > pr){
      if(pp[2][II] > 0){ 
	pp[2][II] = -pp[2][II];
// 	rev_II(NP, nn, ii, II);
// 	rnd[2][nn][ii] *= -1;
	P[2] += 2 * pp[2][II];
      }
      II++;
      if(II == Ntot - 1) II = 0;
    }
    while(P[2] < -pr){
      if(pp[2][II] < 0){ 
	pp[2][II] = -pp[2][II];
// 	rev_II(NP, nn, ii, II);
// 	rnd[2][nn][ii] *= -1;
	P[2] += 2 * pp[2][II];
      }
      II++;
      if(II == Ntot - 1) II = 0;
    }
    if(scont > 1e8) break;
    scont++;
  }
  P[2] -= pp[2][Ntot-1];
  pp[2][Ntot-1] = - P[2];
//   rev_II(NP, nn, ii, Ntot-1);
//   rnd[2][nn][ii] = pp[2][Ntot-1];
  P[2] = 0;
  cnt = 0;
  for(int n = 0; n < NC; ++n){
    int NN = NP.get_Npart(n);
    for(int i = 0; i < NN; ++i){
      rnd[2][n][i] = pp[2][cnt];
      P[2] += rnd[2][n][i];
      cnt++;
    }
  }
  if(fabs(P[2]) > 1e-12){
    ofstream Err("Error.dat", ios::app);
    Err << "Pz!=0 !!\n";
    Err.close();
  }

  double pi[3];
//   ofstream InVel("iniVel.dat");
//   if(! InVel){
//     ofstream Err("Error.dat", ios::app);
//     Err << "Couldn't open " << "iniVel.dat" << endl;
//     Err.close();
//     InVel.close();
//     return;
//   }
  for(int n = 0; n < NC; ++n){
    int NN = NP.get_Npart(n);
    //double mm = part[n][0].get_mass();
    for(int i = 0; i < NN; i++){
      for(int j = 0; j < 3; j++)
	pi[j] = rnd[j][n][i]; 
      part[n][i].set_mom(pi);
//       InVel << pi[0]/mm << " " << pi[1]/mm << " " << pi[2]/mm <<endl; 
    }
  }
//   InVel.close();
  //***************************************************************************
  
  //for angular velocities*****************************************************
  Dynvec<Dynvec<double> > ond[3];
  ond[0].set_size(NC); ond[1].set_size(NC); ond[2].set_size(NC); 
  for(int n = 0; n < NC; ++n){
    for(int j = 0; j < 3; ++j)
      ond[j][n].set_size(NP.get_Npart(n)+1);
  }
  for(int n = 0; n < NC; ++n){
    int NN = NP.get_Npart(n);
    for(int i = 0; i < NN; ++i){
      for(int j = 0; j < 3; ++j)
	ond[j][n][i] = gsl_rng_uniform(ran);
    }
  }
  for(int n = 0; n < NC; ++n){
    int NN = NP.get_Npart(n);
    double ine[3]; part[n][0].get_inertia(ine);
    double mtemp[3] = {temp*ine[0], temp*ine[1], temp*ine[2]};
    for(int i = 0; i < NN; i += 2){
      if((fmod(static_cast<float>(NN),2) != 0.) && (i == NN-1)){
	for(int j = 0; j < 3; j++){
	  double rr = gsl_rng_uniform(ran);
	  gauss(mtemp[j], rr, ond[j][n][i]);
	}
      }else{
	for(int j = 0; j < 3; j++){
	  gauss(mtemp[j], ond[j][n][i], ond[j][n][i+1]);
	}
      }
    }
  }
  
  double T[3] = {0., 0., 0.};
  double tt[3][Ntot+1];
  
  cnt = 0;
  for(int n = 0; n < NC; ++n){
    int NN = NP.get_Npart(n);
    for(int i = 0; i < NN; ++i){
      for(int j = 0; j < 3; j++){
	tt[j][cnt] = ond[j][n][cnt];
	T[j] += tt[j][cnt];
      }
      cnt++;
    }
  }

  double tr = 0.1;
  
//   //
//   cout<<"begin rot X..\n";
//   //
  
 //  // X part/////////////////////////////////////////////////
  II = 0; /*nn=0; ii=0;*/
  scont = 0;
  while(fabs(T[0]) > tr){ //NOTE: in fabsP > x, if the value of x is very 
                           // small the program will have an error. 
                           //(what is small depends on Ntot) 
                           // for N=256 x=0.1 is fine.
    while(T[0] > tr){
      if(tt[0][II] > 0){ 
	tt[0][II] = -tt[0][II];
// 	rev_II(NP, nn, ii, II);
// 	ond[0][nn][ii] *= -1;
	T[0] += 2 * tt[0][II];
      }
      II++;
      if(II == Ntot - 1) II = 0;
    }
    while(T[0] < -tr){
      if(tt[0][II] < 0){ 
	tt[0][II] = -tt[0][II];
// 	rev_II(NP, nn, ii, II);
// 	ond[0][nn][ii] *= -1;
	T[0] += 2 * tt[0][II];
      }
      II++;
      if(II == Ntot - 1) II = 0;
    }
    if(scont > 1e8) break;
    scont++;
  }
  T[0] -= tt[0][Ntot-1];
  tt[0][Ntot-1] = - T[0];
//   rev_II(NP, nn, ii, Ntot-1);
//   ond[0][nn][ii] = tt[0][Ntot-1];
  T[0] = 0.; 
  cnt = 0;
  for(int n = 0; n < NC; ++n){
    int NN = NP.get_Npart(n);
    for(int i = 0; i < NN; ++i){
      ond[0][n][i] = tt[0][cnt];
      T[0] += ond[0][n][i];
      cnt++;
    }
  }
  if(fabs(T[0]) > 1e-12){
    ofstream Err("Error.dat", ios::app);
    Err << "Tx!=0 !!\n";
    Err.close();
  }
  
//   //
//   cout<<"begin rot Y..\n";
//   //
  
  //  // Y part/////////////////////////////////////////////////
  II = 0;
  scont = 0;
  while(fabs(T[1]) > tr){ //NOTE: in fabsP > x, if the value of x is very 
                           // small the program will have an error. 
                           //(what is small depends on Ntot) 
                           // for N=256 x=0.1 is fine.
    while(T[1] > tr){
      if(tt[1][II] > 0){ 
	tt[1][II] = -tt[1][II];
// 	rev_II(NP, nn, ii, II);
// 	ond[1][nn][ii] *= -1;
	T[1] += 2 * tt[1][II];
      }
      II++;
      if(II == Ntot - 1) II = 0;
    }
    while(T[1] < -tr){
      if(tt[1][II] < 0){ 
	tt[1][II] = -tt[1][II];
// 	rev_II(NP, nn, ii, II);
// 	ond[1][nn][ii] *= -1;
	T[1] += 2 * tt[1][II];
      }
      II++;
      if(II == Ntot - 1) II = 0;
    }
    if(scont > 1e8) break;
    scont++;
  }
  T[1] -= tt[1][Ntot-1];
  tt[1][Ntot-1] = - T[1];
//   rev_II(NP, nn, ii, Ntot-1);
//   ond[1][nn][ii] = tt[1][Ntot-1];
  T[1] = 0.;
  cnt = 0;
  for(int n = 0; n < NC; ++n){
    int NN = NP.get_Npart(n);
    for(int i = 0; i < NN; ++i){
      ond[1][n][i] = tt[1][cnt];
      T[1] += ond[1][n][i];
      cnt++;
    }
  }
  if(fabs(T[1]) > 1e-12){
    ofstream Err("Error.dat", ios::app);
    Err << "Ty!=0 !!\n";
    Err.close();
  }
  
//   //
//   cout<<"begin rot Z..\n";
//   //
  
  //  // Z part/////////////////////////////////////////////////
  II = 0;
  scont = 0;
  while(fabs(T[2]) > tr){ //NOTE: in fabsP > x, if the value of x is very 
                           // small the program will have an error. 
                           //(what is small depends on Ntot) 
                           // for N=256 x=0.1 is fine.
    while(T[2] > tr){
      if(tt[2][II] > 0){ 
	tt[2][II] = -tt[2][II];
// 	rev_II(NP, nn, ii, II);
// 	ond[2][nn][ii] *= -1;
	T[2] += 2 * tt[2][II];
      }
      II++;
      if(II == Ntot - 1) II = 0;
    }
    while(T[2] < -tr){
      if(tt[2][II] < 0){ 
	tt[2][II] = -tt[2][II];
// 	rev_II(NP, nn, ii, II);
// 	ond[2][nn][ii] *= -1;
	T[2] += 2 * tt[2][II];
      }
      II++;
      if(II == Ntot - 1) II = 0;
    }
    if(scont > 1e8) break;
    scont++;
  }
  T[2] -= tt[2][Ntot-1];
  tt[2][Ntot-1] = - T[2];
//   rev_II(NP, nn, ii, Ntot-1);
//   ond[2][nn][ii] = tt[2][Ntot-1];
  T[2] = 0.; 
  cnt = 0;
  for(int n = 0; n < NC; ++n){
    int NN = NP.get_Npart(n);
    for(int i = 0; i < NN; ++i){
      ond[2][n][i] = tt[2][cnt];
      T[2] += ond[2][n][i];
      cnt++;
    }
  }
  if(fabs(T[2]) > 1e-12){
    ofstream Err("Error.dat", ios::app);
    Err << "Tz!=0 !!\n";
    Err.close();
  }
  

  double li[3];
//   ofstream InWel("iniWel.dat");
//   if(! InWel){
//     ofstream Err("Error.dat", ios::app);
//     Err << "Couldn't open " << "iniWel.dat" << endl;
//     Err.close();
//     InWel.close();
//     return;
//   }
  for(int n = 0; n < NC; ++n){
    int NN = NP.get_Npart(n);
    //double in = part[n][0].get_inertia();
    for(int i = 0; i < NN; i++){
      for(int j = 0; j < 3; j++)
	li[j] = ond[j][n][i]; 
      part[n][i].set_angm(li);
//       InWel << li[0]/in << " " << li[1]/in << " " << li[2]/in <<endl; 
    }
  }
//   InWel.close();
}

//Box-Muller algorithm
void gauss(double variance, double& rnd1, double& rnd2){
  double rndA, rndB;
  rndA = sqrt(- 2 * log(rnd1)) * cos(2 * M_PI * rnd2);
  rndB = sqrt(- 2 * log(rnd1)) * sin(2 * M_PI * rnd2);
  double sqrtVar = sqrt(variance);
  rnd1 = rndA * sqrtVar;
  rnd2 = rndB * sqrtVar;
}

int calc_II(const part_num& NP, int nn, int ii){
  int II = 0;
  for(int n = 0; n < nn; ++n)
    II += NP.get_Npart(n);
  II += ii;
  return II;
}

void rev_II(const part_num& NP, int& nn, int& ii, int II){
  int NC = NP.get_Ncomp();
  int Nto=0, Ntn=0;
  for(int n = 0; n < NC; ++n){
    Nto = Ntn;
    int NN = NP.get_Npart(n);
    Ntn += NN; 
    if(II < Ntn){
      nn = n;
      ii = II - Nto;
      break;
    }
  }
}

bool print_conf(const Dynvec<Dynvec<Spherocyl> >& part, Thermostat thermo[2], const part_num& NP, 
		const double L[3], bool fin){
  int NC = NP.get_Ncomp();
  char name[30];
  sprintf(name, "conf.dat");
  std::ofstream Con(name);
  if(! Con){
    Con.close();
    return false;
  }
  double r[3], p[3], l[3];
  Matrix<double> o(3, 3);
  Dipole u;
  Con << L[0] << " " << L[1] << " " << L[2] << endl;
  for(int n = 0; n < NC; ++n){
    int n_dip = part[n][0].get_ndips();
    int NN = NP.get_Npart(n);
    for(int i = 0; i < NN; ++i){
      part[n][i].get_pos(r); part[n][i].get_mom(p);
      part[n][i].get_omat(o); part[n][i].get_angm(l);
      Con << setiosflags(ios::fixed) << setprecision(12);
      Con << r[0] << " " << r[1] << " " << r[2] << " "
	  << p[0] << " " << p[1] << " " << p[2] << " ";
      for(int a = 0; a < 3; ++a){
	for(int b = 0; b < 3; ++b)
	  Con << o[a][b] << " ";
      }
      Con << l[0] << " " << l[1] << " " << l[2] << " "
	  << " | " << part[n][i].get_mass() << " ";
      for(int j = 0; j < n_dip; ++j){
	u = part[n][i].get_dip(j);
	u.get_pos(r); u.get_ori(p); 
	Con << r[0] << " " << r[1] << " " << r[2] << " "
	    << p[0] << " " << p[1] << " " << p[2] << " "
	    << u.get_dip() << " . ";
      }
      Con << part[n][i].get_clust();
      Con << endl;
    }
  }
  for(unsigned i = 0; i < 2; ++i){
    Con << thermo[i].get_coo() << " " << thermo[i].get_strn() << endl;
  }
  Con.close();
  //print state of random generators//
  for(unsigned i = 0; i < 2; ++i){
    if(! thermo[i].print_ran(i, fin)){
      ofstream Err("Error.dat", ios::app);
      Err << "problem printing ran_state" << i << ".dat\n"; 
      Err.close();
      return false;
    }
// 	//
// 	cout<<" t"<<thermo[i].get_coo()<<" ";
// 	//
  }
  return true;
}

bool chrg_conf(Dynvec<Dynvec<Spherocyl> >& part, Thermostat thermo[2], const part_num& NP, 
		double L[3]){
  std::string line;
  int NC = NP.get_Ncomp();
  char name[30];
  sprintf(name, "conf.dat");
  std::ifstream Con(name);
  if(! Con){
    Con.close();
    return false;
  }
  int cl;
  double r[3], p[3], l[3], m;
  Matrix<double> o(3,3);
  Dipole u;
  Con >> L[0] >> L[1] >> L[2];
  for(int n = 0; n < NC; ++n){
    int n_dip = part[n][0].get_ndips();
    int NN = NP.get_Npart(n);
    for(int i = 0; i < NN; ++i){
      Con >> r[0] >> r[1] >> r[2]
	  >> p[0] >> p[1] >> p[2];
      for(int a = 0; a < 3; ++a){
	for(int b = 0; b < 3; ++b)
	  Con >> o[a][b];
      }
      Con >> l[0] >> l[1] >> l[2]
	  >> line >> m; 
      part[n][i].set_pos(r); part[n][i].set_mom(p);
      part[n][i].set_omat(o); part[n][i].set_angm(l);
      part[n][i].set_mass(m);
      part[n][i].set_inertia();
      for(int j = 0; j < n_dip; ++j){
	for(int a = 0; a < 8; ++a) Con >> line;
      }
      Con >> cl;
      part[n][i].set_clust(cl);
    }
  }
  double s, cn;
  for(unsigned i = 0; i < 2; ++i){
    Con >> s >> cn;
    thermo[i].set_coo(s); 
    thermo[i].set_strn(cn);
  }
  Con.close();
  //charge thermo rng state
  for(unsigned i = 0; i < 2; ++i){
    if(! thermo[i].read_ran(i)){
      ofstream Err("Error.dat", ios::app);
      Err << "problem reading ran_state" << i << ".dat\n"; 
      Err.close();
      return false;
    }
  }
  return true;
}

void neigh_list(const Dynvec<Dynvec<Spherocyl> >& part, const part_num& NP, const double L[3], 
		deque<int> nblist[], double r_list){
  int Tpart = NP.get_Tpart();
  double rl_sq = r_list * r_list;
  for(int i = 0; i < Tpart; ++i)
    nblist[i].clear(); 
  
  #pragma omp parallel default(shared)
  {
    int tasks2 = 16 * omp_get_num_threads();
    int chunk = static_cast<int>(ceil(static_cast<float>(Tpart)/tasks2));
    #pragma omp for schedule(dynamic, chunk)
    for(int f = 0; f < Tpart-1; ++f){
      double ri[3], rj[3], rij[3];
      // m and i
      int m=0, i=0;
      rev_II(NP, m, i, f);
      part[m][i].get_pos(ri);
      for(int g = f+1; g < Tpart; ++g){
	// n and j
	int n=0, j=0; 
	rev_II(NP, n, j, g);
	part[n][j].get_pos(rj);
	
	for(int a = 0; a < 3; ++a){
	  rij[a] = ri[a] - rj[a];
	  rij[a] -= L[a] * floor(rij[a] / L[a] + 0.5);
	}
	double st[2], dij[3];
	double dsq = spherocyl_dist(part[m][i], part[n][j], L, rij, dij, st);
	
	if(dsq < rl_sq)
	  nblist[f].push_back(g); //neigbor list
      }
    }
  }//end parallel
}

bool savec_list(const Dynvec<Dynvec<Spherocyl> >& part, const part_num& NP){
  int NC = NP.get_Ncomp();
  double r[3];
  ofstream Conf("compConf.dat");
  if(! Conf){
    Conf.close();
    return false;
  }
  Conf << setiosflags(ios::fixed) << setprecision(12);
  for(int n = 0; n < NC; ++n){
    int NN = NP.get_Npart(n);
    for(int i = 0; i < NN; ++i){
      part[n][i].get_pos(r);
      Conf << r[0] << " " << r[1] << " " << r[2] << endl;
    }
  }
  Conf.close();
  return true;
}

bool check_list(const Dynvec<Dynvec<Spherocyl> >& part, const part_num& NP, const double L[3], 
		double r_list_min, double r_list_max, bool& test){
  int NC = NP.get_Ncomp();
  test = false;
  double r[3], rx, ry, rz;
  double Dx, Dy, Dz;
  ifstream Conf("compConf.dat");
  if(! Conf){
    Conf.close();
    return false;
  }
  double dismax = 0.;
  for(int n = 0; n < NC; ++n){
    int NN = NP.get_Npart(n);
    for(int i = 0; i < NN; ++i){
      part[n][i].get_pos(r);
      Conf >> rx >> ry >> rz;
      Dx = r[0] - rx;
      Dx -= L[0] * floor(Dx / L[0] + 0.5);
      Dy = r[1] - ry;
      Dy -= L[1] * floor(Dy / L[1] + 0.5);
      Dz = r[2] - rz;
      Dz -= L[2] * floor(Dz / L[2] + 0.5);
      dismax = max(fabs(Dx), dismax);
      dismax = max(fabs(Dy), dismax);
      dismax = max(fabs(Dz), dismax);
    }
  }
  Conf.close();
  //test of the list skin crossing
  dismax = 2. * sqrt(3.*dismax*dismax);
  if(dismax > r_list_max-r_list_min){
    test = true;
//     //
//     cout<<dismax<<" "<<r_list-r_cut<<" ";
//     //
  }
  return true;
}


bool force_sc(Dynvec<Dynvec<Spherocyl> >& part, const part_num& NP, const double L[3], 
	   double r_cut, double& U, double& W, double PT[3][3], const deque<int> nblist[]){
  int Tpart = NP.get_Tpart();
  double r_cut_sq = r_cut * r_cut;
  double Uc = 0.;
  double Wc = 0.;
  double P00 = 0., P01 = 0., P02 = 0.;
  double P10 = 0., P11 = 0., P12 = 0.;
  double P20 = 0., P21 = 0., P22 = 0.;
  double srt = pow(2., 1./6);
  double srtsq = srt * srt;
  
  bool abort = false;
  #pragma omp parallel default(shared)
  {
    int tasks2 = 16 * omp_get_num_threads();
    int chunk = static_cast<int>(floor(static_cast<float>(Tpart)/tasks2));
    #pragma omp for schedule(dynamic, chunk) reduction(+: Uc,Wc,P00,P01,P02,P10,P11,P12,P20,P21,P22)
    for(int id = 0; id < Tpart; ++id){
      #pragma omp flush(abort)
      if(!abort){
	double ri[3], rj[3], f[3]; 
	double oi[3], oj[3], t[3];
	int m=0, n=0, i=0, j=0;
	int size = nblist[id].size();
	rev_II(NP, m, i, id);
	part[m][i].get_pos(ri);
	part[m][i].get_ori(oi);
	for(int s = 0; s < size; ++s){
	  int jd = nblist[id][s];
	  rev_II(NP, n, j, jd);
	  part[n][j].get_pos(rj);
	  part[n][j].get_ori(oj);
	  
	  double shi[2], shj[2];
	  part[m][i].get_shape(shi);
	  part[n][j].get_shape(shj);
	  double sigma = (shi[1] + shj[1]) / 2;
	  double sig_sq = sigma * sigma;
	  double rsh = srtsq * sig_sq;
	  
	  double rij[3], rsq=0;
	  for(int a = 0; a < 3; ++a){
	    rij[a] = ri[a] - rj[a];
	    rij[a] -= L[a] * floor(rij[a] / L[a] + 0.5);
	    rsq += rij[a]*rij[a];
	  }
	  double st[2], dij[3];
	  double dsq = spherocyl_dist(part[m][i], part[n][j], L, rij, dij, st);

// 	  //
// 	  if(i==0 && j==1){
// 	    std::ofstream Dis("dist.dat", ios::app);
// 	    Dis<<dij[0]<<" "<<dij[1]<<" "<<dij[2]<< " ; "<<st[0]<<" "<<st[1]<<endl;
// 	    Dis.close();
// 	  }
// 	 
	  if(dsq < 0.49 * sig_sq){ // (0.7*sigma)^2
	    //
	    double rr; dot(rij,rij,rr);
	    ofstream Err("Error.dat", ios::app);
	    Err<<" rij("<<rij[0]<<","<<rij[1]<<","<<rij[2]<<") dij("<<dij[0]<<","<<dij[1]<<","<<dij[2]<<")"<<endl;
	    Err<<dsq<<"<0.7225*"<<sig_sq<<" rsq="<<rr<<" s="<<st[0]<<" t="<<st[1]<<endl;
	    Err<<"m "<<m<<" i "<<i<<endl; Err<<"n "<<n<<" j "<<j<<endl;
	    Err.close();
	    //
	    abort = true;
	    #pragma omp flush(abort) //only safe with bool var
	  }
	  if(dsq > rsq*1.00001){
	    ofstream Err("Error.dat", ios::app);
	    Err <<" rij("<<rij[0]<<","<<rij[1]<<","<<rij[2]<<") dij("<<dij[0]<<","<<dij[1]<<","<<dij[2]<<")"<<endl; 
	    Err <<" rsq="<<rsq<<" > dsq="<<dsq<<"!!\n";
	    Err.close();
	    abort = true;
	    #pragma omp flush(abort) //only safe with bool var
	  }
	  
	  if(dsq < r_cut_sq){ 
	    //Lennard-Jones shifted potential and virial//
	    double SR2s = sig_sq / dsq;
	    double SR6s = SR2s * SR2s * SR2s;
	    double pot = 4 * SR6s * (SR6s - 1.);
	    double UIJ = 0.;
	    double WIJ = 0.;
	    double ddr; dot(rij, dij, ddr);
	    if(dsq < rsh){
	      UIJ = pot + 1.;//shifted potential (epsilon = 1)
	      WIJ = 6 * (pot + 4*SR6s*SR6s) * ddr / dsq;
	    }
	    Uc += UIJ;
	    Wc += WIJ;
	    //Forces and torques are not well defined when the line segments are parallel (infinite points at which 
	    //the distance is minimum), and problems with energy conservation appear (although the probability of such 
	    //problematic configurations in the simulation is very small). As an example run a simulation with two 
	    //particles and the following initial conditions in the conf.dat file:
	    //0 0 0 0 0 0 1 0 0 0 1 0 0 0 1 0.1 0 0 1
	    //1.1 0.25 0 0 0 0 1 0 0 0 1 0 0 0 1 0.1 0 0 1
	    
	    //force//
	    double FIJ = WIJ / ddr;
	    for(int a = 0; a < 3; ++a)
	      f[a] = FIJ * dij[a];
	    part[m][i].add_F(f);
	    //Pressure tensor
	    P00 += rij[0] * f[0]; P01 += rij[0] * f[1]; P02 += rij[0] * f[2];
	    P10 += rij[1] * f[0]; P11 += rij[1] * f[1]; P12 += rij[1] * f[2];
	    P20 += rij[2] * f[0]; P21 += rij[2] * f[1]; P22 += rij[2] * f[2];
	    for(int a = 0; a < 3; ++a)
	      f[a] *= -1;
	    part[n][j].add_F(f);
	    
	    //torque//
	    double icroj[3], icror[3], jcror[3];
	    cros(oi, oj, icroj);
	    cros(oi, rij, icror);
	    cros(oj, rij, jcror);
	    double fact1 = (st[0]-0.5)*shi[0];
	    double fact2 = -(st[1]-0.5)*shj[0];
	    if(st[0] == -1.e+6){ //if parallel (no unique solution for min dist)
	      fact1 = fact2 = 0.;
	    }
	    //if parallel and min dist is between end points
	    if(st[0] == -1.e+10){ //(+/+)
	      fact1 = 0.5*shi[0];
	      fact2 = 0.5*shj[0];
	    }else if(st[0] == -2.e+10){ //(-/-)
	      fact1 = -0.5*shi[0];
	      fact2 = -0.5*shj[0];
	    }else if(st[0] == -3.e+10){ //(+/-)
	      fact1 = 0.5*shi[0];
	      fact2 = -0.5*shj[0];
	    }else if(st[0] == -3.e+10){ //(-/+)
	      fact1 = -0.5*shi[0];
	      fact2 = 0.5*shj[0];
	    }
	    double fmul = fact1 * fact2;
	    for(int a = 0; a < 3; ++a)
	      t[a] = FIJ * (fact1*icror[a] + fmul*icroj[a]);
	    part[m][i].add_T(t);
	    for(int a = 0; a < 3; ++a)
	      t[a] = FIJ * (fact2*jcror[a] - fmul*icroj[a]);
	    part[n][j].add_T(t); 
	  }
	}
      }
    }
  }//end parallel
  if(abort)
    return false;
  U = Uc;
  W = Wc;
  PT[0][0] = P00; PT[0][1] = P01; PT[0][2] = P02;
  PT[1][0] = P10; PT[1][1] = P11; PT[1][2] = P12;
  PT[2][0] = P20; PT[2][1] = P21; PT[2][2] = P22;
  return true;
}

void force_dipR(Dynvec<Dynvec<Spherocyl> >& part, const part_num& NP, const double L[3], 
	        const double Tab[32768][3], double& U, double PT[3][3]){
  int Tpart = NP.get_Tpart();
  double Uc = 0.;
  double P00 = 0., P01 = 0., P02 = 0.;
  double P10 = 0., P11 = 0., P12 = 0.;
  double P20 = 0., P21 = 0., P22 = 0.;
  
  double Lh = 0.5 * (200*1.001) / 32768;
  
  #pragma omp parallel default(shared)
  {
    int tasks2 = 16 * omp_get_num_threads();
    int chunk = static_cast<int>(floor(static_cast<float>(Tpart)/tasks2));
    #pragma omp for schedule(dynamic, chunk) reduction(+: Uc,P00,P01,P02,P10,P11,P12,P20,P21,P22) nowait
    for(int k = 0; k < Tpart-1; ++k){
      double ri[3], rj[3], f[3];
      double oi[3], oj[3], t[3];
      double di[3], dj[3];
      Dipole ui, uj;
      // m and i
      int m=0, i=0;
      rev_II(NP, m, i, k);
      
      for(int l = k+1; l < Tpart; ++l){
	// n and j
	int n=0, j=0; 
	rev_II(NP, n, j, l);
	
	for(int v = 0; v < part[m][i].get_ndips(); ++v){
	  part[m][i].get_lev(v, di);  // lever arm
	  ui = part[m][i].get_dip(v);  
	  ui.get_pos(ri);
	  ui.get_ori(oi);      
	  if(ui.get_dip() != 0.){
	    for(int w = 0; w < part[n][j].get_ndips(); ++w){
	      part[n][j].get_lev(w, dj);  // lever arm
	      uj = part[n][j].get_dip(w);
	      uj.get_pos(rj);
	      uj.get_ori(oj);
	      if(uj.get_dip() != 0.){
	  
		double rij[3];
		for(int a = 0; a < 3; ++a){
		  rij[a] = ri[a] - rj[a];
		  rij[a] -= L[a] * floor(rij[a] / L[a] + 0.5);
		}
		double RIJSQ, UIR, UJR, UIUJ;
		dot(rij, rij, RIJSQ);
		
		double uu = ui.get_dip() * uj.get_dip();
		double SR2 = 1. / RIJSQ;
		double RIJ = sqrt(RIJSQ);
		int ri = static_cast<int>(floor(RIJ / Lh));
		double R_L0 = RIJ - Lh * ri;
		double R_L1 = RIJ - Lh * (ri+1);
	  
		double A = (Tab[ri][0] * R_L0 - Tab[ri-1][0] * R_L1) / Lh;
		double B = (Tab[ri][1] * R_L0 - Tab[ri-1][1] * R_L1) / Lh;
		double C = (Tab[ri][2] * R_L0 - Tab[ri-1][2] * R_L1) / Lh;
		
		dot(oi, rij, UIR);
		dot(oj, rij, UJR);
		double UIJR = UIR * UJR;
		dot(oi, oj, UIUJ);
		double UIJ = 0.;
		
		UIJ = uu * (UIUJ * A - UIJR * B);
		Uc += UIJ;   /*cout<<ui<<" "<<uj<<" "<<UIUJ<<" "<<UIJR<<" "<<UIJ<<endl;*/
	  
		//force//
		double fact = (UIUJ - 5. * UIJR * SR2) * B - UIJR * SR2 * C;
		for(int a = 0; a < 3; ++a)
		  f[a] = uu * (fact * rij[a] + (UIR * oj[a] + UJR * oi[a]) * B);
		part[m][i].add_F(f);
		//Pressure tensor
		P00 += rij[0] * f[0]; P01 += rij[0] * f[1]; P02 += rij[0] * f[2];
		P10 += rij[1] * f[0]; P11 += rij[1] * f[1]; P12 += rij[1] * f[2];
		P20 += rij[2] * f[0]; P21 += rij[2] * f[1]; P22 += rij[2] * f[2];
		for(int a = 0; a < 3; ++a)
		  f[a] *= -1;
		part[n][j].add_F(f);
		
		//torque//
		double icror[3]; cros(oi, rij, icror);
		double dicror[3]; cros(di, rij, dicror);
		double dicroi[3]; cros(di, oi, dicroi);
		
		double jcror[3]; cros(oj, rij, jcror);
		double djcror[3]; cros(dj, rij, djcror);
		double djcroj[3]; cros(dj, oj, djcroj);
		
		double dicroj[3]; cros(di, oj, dicroj);
		double djcroi[3]; cros(dj, oi, djcroi);
		double icroj[3]; cros(oi, oj, icroj);
		for(int a = 0; a < 3; ++a)
		  icroj[a] *= A;
		
		fact = -UIR * UJR * (5*B + C) * SR2 + UIUJ * B; 
		double fact1 = 0.;
		
		for(int a = 0; a < 3; ++a){
		  fact1 = (UJR * (icror[a] + dicroi[a]) + UIR * dicroj[a]) * B;
		  t[a] = uu * (fact * dicror[a] + fact1 - icroj[a]);
		}
		part[m][i].add_T(t);
		for(int a = 0; a < 3; ++a){
		  fact1 = (UIR * (jcror[a] - djcroj[a]) - UJR * djcroi[a]) * B; 
		  t[a] = uu * (-fact * djcror[a] + fact1 + icroj[a]);
		}
		part[n][j].add_T(t);         /*part[n][j].get_T(t);  cout<<t[0]<<","<<t[1]<<","<<t[2]<<endl;*/

	      }
	    }
	  }
	}
      }
    }
  }//end parallel
  U = Uc;
  PT[0][0] = P00; PT[0][1] = P01; PT[0][2] = P02;
  PT[1][0] = P10; PT[1][1] = P11; PT[1][2] = P12;
  PT[2][0] = P20; PT[2][1] = P21; PT[2][2] = P22;
}

void tabul(double Tab[32768][3], double alpha){
  double extL = 200 * 1.001;
  double alsq = alpha * alpha;
  double alft = alsq * alsq;
  double spi = sqrt(M_PI);
  for(int i = 1; i < 32769; i++){
    double r =  0.5 * extL * i / 32768;
    double rsq = r * r;
    double sr1 = 1.0 / r;
    double sr2 = sr1 * sr1;
    double sr3 = sr2 * sr1;
    double ERFC = erfc(alpha * r);
    double aer = 2.0 * alpha * exp(- alsq * rsq) / spi;
    double aa = 2.0 * alsq + 3.0 * sr2;
    Tab[i-1][0] = ERFC * sr3 + aer * sr2;
    Tab[i-1][1] = (3.0 * ERFC * sr3 + aer * aa) * sr2;
    Tab[i-1][2] = 4. * alft * aer; 
  }
  // ofstream Table("TabBC.dat");
//   for(int i = 0; i < 32768; i++){
//     Table << 0.5*extL*(i+1)/32768 << " " << Tab[i][0] << " " << Tab[i][1] 
//           << " " << Tab[i][2] << endl;
//   }
  // Table.close();
}

void force_dipK(Dynvec<Dynvec<Spherocyl> >& part, const part_num& NP, const double L[3], 
	   int K_max, double alpha, const Dynvec<Matrix<double> >& Kvec, double& U, double PT[3][3]){
  int NC = NP.get_Ncomp();
//   double Lmax = max(L[2], max(L[1], L[0]));
  double V = L[0] * L[1] * L[2];
  double P2[3] = {2.*M_PI/L[0], 2.*M_PI/L[1], 2.*M_PI/L[2]};
//   double K_msq = 4*M_PI*M_PI*K_max*K_max/(Lmax*Lmax) * 1.001;
  int K_max2 = K_max * K_max;
  int K_max3 = K_max2 * K_max;
  double a4 = 1. / (4. * alpha);
  double Uc = 0.;
  double P00 = 0., P01 = 0., P02 = 0.;
  double P10 = 0., P11 = 0., P12 = 0.;
  double P20 = 0., P21 = 0., P22 = 0.;
  
  #pragma omp parallel default(shared)
  {
    int chunk = 50;
    #pragma omp for schedule(dynamic, chunk) reduction(+: Uc,P00,P01,P02,P10,P11,P12,P20,P21,P22) nowait
    for(int kn = 0; kn < K_max3; ++kn){
      int nz = static_cast<int>(floor(static_cast<float>(kn)/K_max2));
      float l = kn - K_max2 * nz;
      int ny = static_cast<int>(floor(l/K_max));
      int nx = static_cast<int>(l) - K_max * ny;
      double Kz = P2[2] * nz;
      double Ky = P2[1] * ny;
      double Kx = P2[0] * nx;
      double Ksq = Kx * Kx + Ky * Ky + Kz * Kz;
      
      if(Ksq != 0/* && Ksq < K_msq*/){ 
	double KK = Kvec[nx][ny][nz];
	double val = 2. * KK / V;  //mult by 2 for symmetry
	double pval = 2 * (a4 + 1. / Ksq);
	
	double K[4][3];
	K[0][0] = Kx; K[0][1] = Ky; K[0][2] = Kz;
	K[1][0] = -Kx; K[1][1] = Ky; K[1][2] = Kz;
	K[2][0] = Kx; K[2][1] = -Ky; K[2][2] = Kz;
	K[3][0] = -Kx; K[3][1] = -Ky; K[3][2] = Kz;
	
	//sum array
	complex<double> cc(0., 0.);
	complex<double> sum[4][NC];
	complex<double> uve[4][3]; //not need to separate for different components
	for(int b = 0; b < 4; ++b){
	  for(int k = 0; k < 3; ++k)
	    uve[b][k] = cc;
	  for(int n = 0; n < NC; ++n) 
	    sum[b][n] = cc;
	}
	double self[4][NC]; //extract the self energy corresponding to each term (Kmax finite)
	double uvslf[4][3];
	for(int b = 0; b < 4; ++b){
	  for(int k = 0; k < 3; ++k)
	    uvslf[b][k] = 0.;
	  for(int n = 0; n < NC; ++n)
	    self[b][n] = 0.;
	}
	
	//individual values array FOR EACH DIPOLE (not just particle)!
	Dynvec<Dynvec<Matrix<complex<double> > > > M(4);
	for(int b = 0; b < 4; ++b){
	  M[b].set_size(NC);
	  for(int n = 0; n < NC; ++n){
	    int NN = NP.get_Npart(n);
	    M[b][n].set_size(NN, part[n][0].get_ndips());
	    M[b][n].set_to(cc);
	  }
	}
	
	double ri[3] = {0., 0., 0.}, oi[3] = {0., 0., 0.}; 
	Dipole u;
	
	//generate values
	for(int n = 0; n < NC; ++n){
	  int NN = NP.get_Npart(n);
	  for(int i = 0; i < NN; ++i){
	    complex<double> selfpart[4] = {cc, cc, cc, cc};
	    complex<double> uvslfpart[4][3];
	    for(int b = 0; b < 4; ++b){
	      for(int k = 0; k < 3; ++k)
		uvslfpart[b][k] = cc;
	    }
	    
	    for(int w = 0; w < part[n][i].get_ndips(); ++w){
	      u = part[n][i].get_dip(w);
	      u.get_pos(ri);
	      for(int k = 0; k < 3; ++k)
		ri[k] -= L[k] * floor(ri[k] / L[k] + 0.5);
	      u.get_ori(oi);      /*cout<<"("<<oi[0]<<","<<oi[1]<<","<<oi[2]<<") ";*/
	    
	      if(u.get_dip() != 0.){
		double RIK[4], UIK[4];
		complex<double> z[4];
		for(int b = 0; b < 4; ++b){
		  dot(ri, K[b], RIK[b]); //r.k
		  dot(oi, K[b], UIK[b]); //u.k
		  UIK[b] *= u.get_dip();
		  z[b] = polar(1., RIK[b]); //exp(i r.k)
		}
	      
		for(int b = 0; b < 4; ++b){
		  M[b][n][i][w] = conj(z[b]); //exp(-i r.k)
		  complex<double> ukexp = UIK[b] * M[b][n][i][w];
		  sum[b][n] += ukexp; //sum( (u.k)exp(-i r.k) )
		  selfpart[b] += ukexp; //sum( (u.k)exp(-i r.k) ) just for particle (n, i)
		  for(int k = 0; k < 3; ++k){
		    ukexp = u.get_dip() * oi[k] * M[b][n][i][w];
		    uve[b][k] += ukexp; //sum( u_{a}exp(-i r.k) )
		    uvslfpart[b][k] += ukexp; //sum( u_{a}exp(-i r.k) ) just for particle (n, i)
		  }
		}
	      }
	    }
	    //the interaction between two dipoles (same or different)of the same particle has to be removed!!
	    for(int b = 0; b < 4; ++b){
	      self[b][n] += norm(selfpart[b]); //sum of self energies of particles, including all its internal interactions
	      for(int k = 0; k < 3; ++k)
		uvslf[b][k] += real(conj(selfpart[b]) * uvslfpart[b][k]);
	    }
	  }
	}
	
	//sum over all species
	complex<double> adds[4] = {cc, cc, cc, cc};
	double addslf[4] = {0., 0., 0., 0.};
	for(int b = 0; b < 4; ++b){
	  for(int n = 0; n < NC; ++n){
	    adds[b] += sum[b][n];
	    addslf[b] += self[b][n];
	  }        
	}
	
	//compute energy
	double Upar = 0;
	for(int b = 0; b < 4; ++b){
	  double en = val * 0.5 * (norm(adds[b]) - addslf[b]);
	  Upar += en;
	}       
	if((nx == 0 && ny == 0) || (nx == 0 && nz == 0) || (ny == 0 && nz == 0)){
	  Upar /= 4.;
	}else if((nz == 0 && nx != 0 && ny != 0) || (ny == 0 && nz != 0 && nx != 0) ||
	  (nx == 0 && ny != 0 && nz != 0)){
	  Upar /= 2.;
	}
	Uc += Upar;       
	
	//force and torque
	double f[3], t[3], di[3];
	complex<double> local_adds[4];
	for(int n = 0; n < NC; ++n){
	  int nd = part[n][0].get_ndips();
	  int NN = NP.get_Npart(n);
	  for(int i = 0; i < NN; ++i){
	    //self terms (this was done in the "generate values step". 
	    // Is it better to do it twice like here? or save it (local_adds AND UIK) from before and use more memory?)
	    double UIK[4][nd];
	    for(int b = 0; b < 4; ++b){
	      local_adds[b] = cc;
	      for(int w = 0; w < nd; ++w){
		u = part[n][i].get_dip(w);
		u.get_ori(oi); 
	      
		if(u.get_dip() != 0){
		  dot(oi, K[b], UIK[b][w]);
		  UIK[b][w] *= u.get_dip();
		  local_adds[b] += UIK[b][w] * M[b][n][i][w]; //sum( (u.k)exp(-i r.k) ) just for particle (n, i)
		}
	      }
	    }
	    
	    //compute forces and torques due to each dipole
	    for(int w = 0; w < part[n][i].get_ndips(); ++w){
	      part[n][i].get_lev(w, di);
	      u = part[n][i].get_dip(w);
	      u.get_ori(oi);
	    
	      if(u.get_dip() != 0){
		double icrok[4][3], dicrok[4][3];
		for(int b = 0; b < 4; ++b){
		  cros(di, K[b], dicrok[b]);
		  cros(oi, K[b], icrok[b]);
		  for(int j = 0; j < 3; ++j)
		    icrok[b][j] *= u.get_dip();
		}   
		double factF = 0., factT = 0.;
	      
		//as force and torque are not scalars different values of "a" are not always equivalent, as with the energy
		for(int j = 0; j < 3; ++j){
		  f[j] = 0.;
		  t[j] = 0.;
		}
		for(int b = 0; b < 4; ++b){
		  factF = -val * imag(UIK[b][w] * M[b][n][i][w] * (conj(adds[b]) - conj(local_adds[b])) ); //self force term due to interactions between dipoles in the same particle  
		  for(int j = 0; j < 3; ++j)
		    f[j] += factF * K[b][j];
		
		  factT = -val * real(M[b][n][i][w] * (conj(adds[b]) - conj(local_adds[b])) ); //substract self term
		  for(int j = 0; j < 3; ++j)
		    t[j] += factT * icrok[b][j] + factF * dicrok[b][j];
		}
	      
		if((nx == 0 && ny == 0) || (nx == 0 && nz == 0) || (ny == 0 && nz == 0)){
		  for(int j = 0; j < 3; ++j){
		    f[j] /= 4.;
		    t[j] /= 4.;
		  }
		}else if((nz == 0 && nx != 0 && ny != 0) || (ny == 0 && nz != 0 && nx != 0) ||
		  (nx == 0 && ny != 0 && nz != 0)){
		  for(int j = 0; j < 3; ++j){
		    f[j] /= 2.;
		    t[j] /= 2.;
		  }
		}
		part[n][i].add_F(f);
		part[n][i].add_T(t);
	      }
	    }
	  }
	}
	//pressure tensor
	double P[3][3];
	for(int j = 0; j < 3; ++j){
	  for(int k = 0; k < 3; ++k)
	    P[j][k] = 0.;
	}
	for(int b = 0; b < 4; ++b){
	  for(int j = 0; j < 3; ++j){
	    for(int k = 0; k < 3; ++k)
	      P[j][k] += val * 0.5 * ( (delta(j,k)-pval*K[b][j]*K[b][k])*(norm(adds[b])-addslf[b]) + 2*(real(conj(adds[b])*uve[b][j])-uvslf[b][j])*K[b][k] );
	  }
	}
	if((nx == 0 && ny == 0) || (nx == 0 && nz == 0) || (ny == 0 && nz == 0)){
	  for(int j = 0; j < 3; ++j){
	    for(int k = 0; k < 3; ++k)
	      P[j][k] /= 4.;
	  }
	}else if((nz == 0 && nx != 0 && ny != 0) || (ny == 0 && nz != 0 && nx != 0) ||
	  (nx == 0 && ny != 0 && nz != 0)){
	  for(int j = 0; j < 3; ++j){
	    for(int k = 0; k < 3; ++k)
	      P[j][k] /= 2.;
	  }
	}
	P00 += P[0][0]; P01 += P[0][1]; P02 += P[0][2];
	P10 += P[1][0]; P11 += P[1][1]; P12 += P[1][2];
	P20 += P[2][0]; P21 += P[2][1]; P22 += P[2][2]; 
	
      }
    }
  }//end parallel
  U = Uc;
  PT[0][0] = P00; PT[0][1] = P01; PT[0][2] = P02;
  PT[1][0] = P10; PT[1][1] = P11; PT[1][2] = P12;
  PT[2][0] = P20; PT[2][1] = P21; PT[2][2] = P22;
}

void wave_v(Dynvec<Matrix<double> >& Kvec, const double L[3], double alpha, int K_max){
  double c = 4. * alpha * alpha;
  double P2[3] = {2.*M_PI/L[0], 2.*M_PI/L[1], 2.*M_PI/L[2]};
  for(int nz = 0; nz <= K_max; ++nz){
    double Kz = P2[2] * nz;
    for(int ny = 0; ny <= K_max; ++ny){
      double Ky = P2[1] * ny;
      for(int nx = 0; nx <= K_max; ++nx){
	double Kx = P2[0] * nx;
	double Ksq = Kx * Kx + Ky * Ky + Kz * Kz;
	if(Ksq != 0.){
	  Kvec[nx][ny][nz] = 4. * M_PI * exp(- Ksq / c) / Ksq; 
	}
      }
    }
  }
//   ofstream Vec("Kvec.dat");
//   for(int z = 0; z <= K_max; ++z){
//     for(int y = 0; y <= K_max; ++y){
//       for(int x = 0; x <= K_max; ++x){
// 	Vec << Kvec[x][y][z] << " ";
//       }
//       Vec << "\n";
//     }
//     Vec << "\n";
//   }
//   Vec.close();
}

void ext_field(Dynvec<Dynvec<Spherocyl> >& part, const part_num& NP, const double B[3],
	       double& U){
  int Tpart = NP.get_Tpart();
  
  double Uext = 0.;
  #pragma omp parallel default(shared)
  {
    int tasks2 = 16 * omp_get_num_threads();
    int chunk = static_cast<int>(floor(static_cast<float>(Tpart)/tasks2));
    #pragma omp for schedule(dynamic, chunk) reduction(-: Uext)
    for(int k = 0; k < Tpart; ++k){
       // m and i
      int m=0, i=0;
      rev_II(NP, m, i, k);
      
      double oi[3], t[3], doti;
      Dipole u;
      
      for(int w = 0; w < part[m][i].get_ndips(); ++w){
	u = part[m][i].get_dip(w);
	u.get_ori(oi);
      
	dot(oi, B, doti);
	Uext -= u.get_dip() * doti;
	//Force = 0 for a constant parallel field;
	//torque//
	cros(oi, B, t);
	t[0] *= u.get_dip(); t[1] *= u.get_dip(); t[2] *= u.get_dip();
	part[m][i].add_T(t);
      }
    } 
  }//end parallel
  U = Uext;
}

double self_ener(const Dynvec<Dynvec<Spherocyl> >& part, const part_num& NP){
  int NC = NP.get_Ncomp();
  double ri[3], rj[3];
  double oi[3], oj[3];
  Dipole ui, uj;
  double U = 0.;
  
  for(int n = 0; n < NC; ++n){
    int NN = NP.get_Npart(n);
    double Un = 0.;
    for(int v = 0; v < part[n][0].get_ndips()-1; ++v){
      ui = part[n][0].get_dip(v);  
      ui.get_pos(ri);
      ui.get_ori(oi); 
      for(int w = v+1; w < part[n][0].get_ndips(); ++w){
	uj = part[n][0].get_dip(w);  
	uj.get_pos(rj);
	uj.get_ori(oj); 
	
	double uu = ui.get_dip() * uj.get_dip();
	
	double rij[3];
	for(int a = 0; a < 3; ++a)
	  rij[a] = ri[a] - rj[a];
	  
	double RIJSQ, UIR, UJR, UIUJ;
	dot(rij, rij, RIJSQ);
	double RIJ = sqrt(RIJSQ);
	
	double SR2 = 1. / RIJSQ;
	double SR3 = SR2 / RIJ;
	dot(oi, oj, UIUJ);
	dot(oi, rij, UIR);
	dot(oj, rij, UJR);
	
	Un += uu * (UIUJ * SR3 - 3. * UIR * UJR * SR2 * SR3);
      }
    }
    
    U += NN * Un;
  }
  return U;
}

void moveX(Dynvec<Dynvec<Spherocyl> >& part, const part_num& NP, double Dt){
  int Tpart = NP.get_Tpart();
  
  #pragma omp parallel default(shared)
  {
    int tasks2 = 16 * omp_get_num_threads();
    int chunk = static_cast<int>(floor(static_cast<float>(Tpart)/tasks2));
    #pragma omp for schedule(dynamic, chunk)
    for(int id = 0; id < Tpart; ++id){
      //n and i
      int n=0, i=0;
      rev_II(NP, n, i, id);
      
      //translation
      part[n][i].translate(Dt);
      
      //rotation
      part[n][i].rotate(Dt);
      
    }
  }//end parallel 
}

void moveP(Dynvec<Dynvec<Spherocyl> >& part, Thermostat thermo[2], const part_num& NP, double Dt, 
	   double K[2], double PT[3][3], bool last){
  int Tpart = NP.get_Tpart();
  double Dt2 = Dt / 2;
  double Kc0 = 0., Kc1 = 0.;
  double P00 = 0., P01 = 0., P02 = 0.;
  double P11 = 0., P12 = 0.;
  double P22 = 0.;
  
  #pragma omp parallel default(shared)
  {
    int tasks2 = 16 * omp_get_num_threads();
    int chunk = static_cast<int>(floor(static_cast<float>(Tpart)/tasks2));
    #pragma omp for schedule(dynamic, chunk) reduction(+: Kc0,Kc1,P00,P01,P02,P11,P12,P22)
    for(int id = 0; id < Tpart; ++id){
      //n and i
      int n=0, i=0;
      rev_II(NP, n, i, id);
      
      double m = part[n][i].get_mass();
      double I[3]; part[n][i].get_inertia(I);
      double p[3], f[3];
      part[n][i].get_mom(p);
      part[n][i].get_F(f);
      double l[3], t[3];
      part[n][i].get_angm(l);                  
      part[n][i].get_T(t);
      
      Matrix<double> o(3,3); part[n][i].get_omat(o);
      Matrix<double> I0(3,3,0.), Ib(3,3,0.), ll(3,1,0.), lres(1,1,0.);
      I0[0][0] = 1./I[0]; I0[1][1] = 1./I[1]; I0[2][2] = 1./I[2];
      //transform I^-1 to box frame
      Ib = transpose(o) * I0 * o;
      double psq = 0., lsq = 0.;
      
      if(!last){
	//translational
	double s = thermo[0].get_coo();
	double sfac = 1. + Dt2 * s;       
	for(int j = 0; j < 3; j++){
	  p[j] = (p[j] + Dt2 * f[j]) / sfac;
	  psq += p[j] * p[j];
	}
	part[n][i].set_mom(p);                   
      
	//rotational
	s = thermo[1].get_coo();
	sfac = 1. + Dt2 * s;
	for(int j = 0; j < 3; j++)
	  l[j] = (l[j] + Dt2 * t[j]) / sfac; 
	
	for(int a = 0; a < 3; ++a)
	  ll[a][0] = l[a];
	lres = transpose(ll) * Ib * ll;
	lsq += lres[0][0];
	part[n][i].set_angm(l);                  /*cout<<"b("<<w[0]<<","<<w[1]<<","<<w[2]<<") ";*/
      }else{
	//translational
	double s = thermo[0].get_coo();
	double sfac = 1. - Dt2 * s;
	for(int j = 0; j < 3; j++){
	  p[j] = p[j] * sfac + Dt2 * f[j];
	  psq += p[j] * p[j];
	}
	part[n][i].set_mom(p);
	P00 += p[0] * p[0]/m; P01 += p[0] * p[1]/m; P02 += p[0] * p[2]/m;
			      P11 += p[1] * p[1]/m; P12 += p[1] * p[2]/m;
						    P22 += p[2] * p[2]/m;
	
	//rotational
	s = thermo[1].get_coo();
	sfac = 1. - Dt2 * s;
	for(int j = 0; j < 3; j++)
	  l[j] = l[j] * sfac + Dt2 * t[j];
	for(int a = 0; a < 3; ++a)
	  ll[a][0] = l[a];
	lres = transpose(ll) * Ib * ll;
	lsq += lres[0][0];
	part[n][i].set_angm(l);   
      }
      
      Kc0 += psq / m; 
      Kc1 += lsq;
    }
  }//end parallel 
  K[0] = Kc0 / 2;
  K[1] = Kc1 / 2;
  PT[0][0] = P00; PT[0][1] = P01; PT[0][2] = P02;
  PT[1][0] = P01; PT[1][1] = P11; PT[1][2] = P12;
  PT[2][0] = P02; PT[2][1] = P12; PT[2][2] = P22;
}

void moveTher(Thermostat thermo[2], double temp, double Dt, const double K[2], double Tpart, 
	      int mode){
  double Q, s, g, ex1;
  double rnd[2], rnd1;
  for(int i = 0; i < 2; ++i){
    Q = thermo[i].get_mass();
    switch(mode){
      case 0:
	rnd[0] = thermo[i].get_rand();
	rnd[1] = thermo[i].get_rand();
	gauss(1, rnd[0], rnd[1]);          //cout<<rnd[0]<<" "<<rnd[1]<<endl;
	thermo[i].set_strn(rnd[1]);
	g = thermo[i].get_shr();
	ex1 = exp(-g*Dt/2);
	s = ex1*thermo[i].get_coo() + sqrt(2.*temp*g*Dt/Q)*rnd[0]; //cout<<sqrt(2.*temp*g/Q)*rnd[0]<<" ";
	thermo[i].set_coo(s);         
	break;
      case 1:
	s = thermo[i].get_coo() + Dt*(2.*K[i] - 3.*Tpart*temp)/Q;
	thermo[i].set_coo(s);
	break;
      case 2:
	rnd1 = thermo[i].get_strn();
	g = thermo[i].get_shr();
	ex1 = exp(-g*Dt/2);
	s = ex1*thermo[i].get_coo() + sqrt(2.*temp*g*Dt/Q)*rnd1; //cout<<sqrt(2.*temp*g/Q)*rnd1<<" .";
	thermo[i].set_coo(s);
	break;
      default:
	break;
    }
  }
}

bool save_cont(double Time, long cont, int Ntcf, int Ntcv, int mcont, const accumulator& AcT, 
	       const accumulator& AcTsq, const long Npart_in_clust[2]){
  double PT[3][3];
  char name[30];
  sprintf(name, "counters.dat");
  ofstream Con(name);
  if(! Con){
    Con.close();
    return false;
  }
  Con << setiosflags(ios::fixed) << setprecision(12);
  Con << Time << " " << cont << " " << Ntcf << " " << Ntcv << " "
      << mcont << " ";
  Con << AcT.get_K() << " " << AcTsq.get_K() << " "
      << AcT.get_U() << " " << AcTsq.get_U() << " "
      << AcT.get_P() << " " << AcTsq.get_P() << " "
      << AcT.get_Uss() << " " << AcTsq.get_Uss() << " "
      << AcT.get_Pss() << " " << AcTsq.get_Pss() << " "
      << AcT.get_UdipR() << " " << AcTsq.get_UdipR() << " "
      << AcT.get_PdipR() << " " << AcTsq.get_PdipR() << " "
      << AcT.get_UdipK() << " " << AcTsq.get_UdipK() << " "
      << AcT.get_PdipK() << " " << AcTsq.get_PdipK() << " "
      << AcT.get_Uext() << " " << AcTsq.get_Uext() << " ";
  AcT.get_PT(PT);
  for(int i = 0; i < 3; ++i){
    for(int j = 0; j < 3; ++j)
      Con << PT[i][j] << " ";
  }
  AcTsq.get_PT(PT);
  for(int i = 0; i < 3; ++i){
    for(int j = 0; j < 3; ++j)
      Con << PT[i][j] << " ";
  }
  AcT.get_PTss(PT);
  for(int i = 0; i < 3; ++i){
    for(int j = 0; j < 3; ++j)
      Con << PT[i][j] << " ";
  }
  AcTsq.get_PTss(PT);
  for(int i = 0; i < 3; ++i){
    for(int j = 0; j < 3; ++j)
      Con << PT[i][j] << " ";
  }
  AcT.get_PTdipR(PT);
  for(int i = 0; i < 3; ++i){
    for(int j = 0; j < 3; ++j)
      Con << PT[i][j] << " ";
  }
  AcTsq.get_PTdipR(PT);
  for(int i = 0; i < 3; ++i){
    for(int j = 0; j < 3; ++j)
      Con << PT[i][j] << " ";
  }
  AcT.get_PTdipK(PT);
  for(int i = 0; i < 3; ++i){
    for(int j = 0; j < 3; ++j)
      Con << PT[i][j] << " ";
  }
  AcTsq.get_PTdipK(PT);
  for(int i = 0; i < 3; ++i){
    for(int j = 0; j < 3; ++j)
      Con << PT[i][j] << " ";
  }
  AcT.get_PTkin(PT);
  for(int i = 0; i < 3; ++i){
    for(int j = 0; j < 3; ++j)
      Con << PT[i][j] << " ";
  }
  AcTsq.get_PTkin(PT);
  for(int i = 0; i < 3; ++i){
    for(int j = 0; j < 3; ++j)
      Con << PT[i][j] << " ";
  }
  Con << AcT.get_G1() << " " << AcTsq.get_G1() << " "
      << AcT.get_G2() << " " << AcTsq.get_G2() << " ";
  Con << Npart_in_clust[0] << " " << Npart_in_clust[1] << " ";
  Con << endl;
  Con.close();
  return true;
}

bool chrg_cont(double& Time, long& cont, int& Ntcf, int& Ntcv, int& mcont, accumulator& AcT, 
	       accumulator& AcTsq, long Npart_in_clust[2]){
  double var;
  double PT[3][3];
  char name[30];
  sprintf(name, "counters.dat");
  ifstream Con(name);
  if(! Con){
    Con.close();
    return false;
  }
  Con >> Time >> cont >> Ntcf >> Ntcv >> mcont;
  Con >> var; AcT.set_K(var);
  Con >> var; AcTsq.set_K(var);
  Con >> var; AcT.set_U(var);
  Con >> var; AcTsq.set_U(var);
  Con >> var; AcT.set_P(var);
  Con >> var; AcTsq.set_P(var);
  Con >> var; AcT.set_Uss(var);
  Con >> var; AcTsq.set_Uss(var);
  Con >> var; AcT.set_Pss(var);
  Con >> var; AcTsq.set_Pss(var);
  Con >> var; AcT.set_UdipR(var);
  Con >> var; AcTsq.set_UdipR(var);
  Con >> var; AcT.set_PdipR(var);
  Con >> var; AcTsq.set_PdipR(var);
  Con >> var; AcT.set_UdipK(var);
  Con >> var; AcTsq.set_UdipK(var);
  Con >> var; AcT.set_PdipK(var);
  Con >> var; AcTsq.set_PdipK(var);
  Con >> var; AcT.set_Uext(var);
  Con >> var; AcTsq.set_Uext(var);
  for(int i = 0; i < 3; ++i){
    for(int j = 0; j < 3; ++j){
      Con >> var; PT[i][j] = var;
    }
  }
  AcT.set_PT(PT);
  for(int i = 0; i < 3; ++i){
    for(int j = 0; j < 3; ++j){
      Con >> var; PT[i][j] = var;
    }
  }
  AcTsq.set_PT(PT);
  for(int i = 0; i < 3; ++i){
    for(int j = 0; j < 3; ++j){
      Con >> var; PT[i][j] = var;
    }
  }
  AcT.set_PTss(PT);
  for(int i = 0; i < 3; ++i){
    for(int j = 0; j < 3; ++j){
      Con >> var; PT[i][j] = var;
    }
  }
  AcTsq.set_PTss(PT);
  for(int i = 0; i < 3; ++i){
    for(int j = 0; j < 3; ++j){
      Con >> var; PT[i][j] = var;
    }
  }
  AcT.set_PTdipR(PT);
  for(int i = 0; i < 3; ++i){
    for(int j = 0; j < 3; ++j){
      Con >> var; PT[i][j] = var;
    }
  }
  AcTsq.set_PTdipR(PT);
  for(int i = 0; i < 3; ++i){
    for(int j = 0; j < 3; ++j){
      Con >> var; PT[i][j] = var;
    }
  }
  AcT.set_PTdipK(PT);
  for(int i = 0; i < 3; ++i){
    for(int j = 0; j < 3; ++j){
      Con >> var; PT[i][j] = var;
    }
  }
  AcTsq.set_PTdipK(PT);
  for(int i = 0; i < 3; ++i){
    for(int j = 0; j < 3; ++j){
      Con >> var; PT[i][j] = var;
    }
  }
  AcT.set_PTkin(PT);
  for(int i = 0; i < 3; ++i){
    for(int j = 0; j < 3; ++j){
      Con >> var; PT[i][j] = var;
    }
  }
  AcTsq.set_PTkin(PT);
  Con >> var; AcT.set_G1(var);
  Con >> var; AcTsq.set_G1(var);
  Con >> var; AcT.set_G2(var);
  Con >> var; AcTsq.set_G2(var);
  Con >> Npart_in_clust[0] >> Npart_in_clust[1];
  Con.close();
  return true;
}

void acc_correl(const Dynvec<Dynvec<Spherocyl> >& part, Dynvec<Matrix<double> >& Hista, 
	       Dynvec<Matrix<double> >& Histe, const double L[3], double bar_w, 
		const part_num& NP, const cross Vc[]){
  int NC = NP.get_Ncomp();
  int Tpart = NP.get_Tpart();
  int cst, n_bars;
  Hista[0].get_size(cst, n_bars);
  
  #pragma omp parallel default(shared) //private(P)
  {
    int tasks2 = 16 * omp_get_num_threads();
    int chunk = static_cast<int>(static_cast<float>(Tpart)/tasks2);
    #pragma omp for schedule(dynamic, chunk) nowait
    for(int f = 0; f < Tpart-1; ++f){
      double ri[3], rj[3];
      double oi[3], oj[3];
      // m and i
      int m=0, i=0;
      rev_II(NP, m, i, f);
      part[m][i].get_pos(ri);
      part[m][i].get_ori(oi);
      for(int g = f+1; g < Tpart; ++g){
	// n and j
	int n=0, j=0;
	rev_II(NP, n, j, g);
	part[n][j].get_pos(rj);
	part[n][j].get_ori(oj);
	
	double XIJ = ri[0] - rj[0];
	XIJ -= L[0] * floor(XIJ / L[0] + 0.5);
	double YIJ = ri[1] - rj[1];
	YIJ -= L[1] * floor(YIJ / L[1] + 0.5);
	double ZIJ = ri[2] - rj[2];
	ZIJ -= L[2] * floor(ZIJ / L[2] + 0.5);
	double RIJSQ = XIJ * XIJ + YIJ * YIJ + ZIJ * ZIJ;
	double RIJ = sqrt(RIJSQ);
	
	double dot = oi[0]*oj[0] + oi[1]*oj[1] + oi[2]*oj[2];
	double dotsq = dot * dot; 
	double doti = (oi[0]*XIJ + oi[1]*YIJ + oi[2]*ZIJ) / RIJ;
	double dotj = (oj[0]*XIJ + oj[1]*YIJ + oj[2]*ZIJ) / RIJ;
	
	int bin = static_cast<int>(floor(RIJ / bar_w));
	if(m == n){//intra
	  #pragma omp atomic
	  Hista[n][0][bin] += 2.; // 000
	  #pragma omp atomic
	  Hista[n][1][bin] += 2. * dot; // 110
	  #pragma omp atomic
	  Hista[n][2][bin] += 2. * (3.*doti*dotj - dot); // 112
	  #pragma omp atomic
	  Hista[n][3][bin] += 3.*dotsq - 1.; // 220
	}else{//inter
	  int s;
	  for(s = 0; s < NC-1; ++s){
	    if(Vc[m].P[s] == n) break;
	  }
	  int cc = Vc[m].C[s];
	  #pragma omp atomic
	  Histe[cc][0][bin] += 1.; // 000
	  #pragma omp atomic
	  Histe[cc][1][bin] += dot; // 110
	  #pragma omp atomic
	  Histe[cc][2][bin] += 3.*doti*dotj - dot; // 112
	  #pragma omp atomic
	  Histe[cc][3][bin] += (3.*dotsq - 1.) / 2; // 220
	}
      }
    }
  }//end parallel
}

bool nor_correl(const Dynvec<Matrix<double> >& Hista, const Dynvec<Matrix<double> >& Histe,
		double Lmin, double bar_w, const part_num& NP, double cor_cont){
  int NC = NP.get_Ncomp();
  int cst, n_bars;
  Hista[0].get_size(cst, n_bars);
  char name[30];
  sprintf(name, "correlation.dat");
  ofstream HIST(name);
  if(! HIST){
    HIST.close();
    return false;
  }
  double pd = 4. * M_PI;
  double bar3 = bar_w * bar_w * bar_w / 12.;
  double sq2 = sqrt(2.);
  double hnrma, hnrme;
  
  for(int bin = 0; bin < n_bars; ++bin){
    double pos = (static_cast<double>(bin) + 0.5) * bar_w;
    double norm;
    if(pos >= Lmin / 2. && pos <= Lmin / sq2){
      norm = pd * pos * pos * ( 3 * Lmin / (2. * pos) - 2.) * bar_w;
    }// else if(pos > Lmin / sq2){
      //}
    else
      norm = pd * (pos * pos * bar_w + bar3);
    double nnorm = norm * cor_cont;
    HIST << cor_cont << " " 
	 << setiosflags(ios::fixed) << setprecision(8) << norm << " "
	 << setiosflags(ios::fixed) << setprecision(4) << pos << " & ";
    HIST << setiosflags(ios::fixed) << setprecision(8);
    //intra
    for(int n = 0; n < NC; ++n){
      int NN = NP.get_Npart(n);
      int NN1 = NN - 1;
      if(NN1 > 0){
	hnrma = 1. * Hista[n][0][bin] / (nnorm * NN * NN1);
	HIST << hnrma << " ";
	hnrma = 3. * Hista[n][1][bin] / (nnorm * NN * NN1);
	HIST << hnrma << " ";
	hnrma = (3./2) * Hista[n][2][bin] / (nnorm * NN * NN1);
	HIST << hnrma << " ";
	hnrma = 5. * Hista[n][3][bin] / (nnorm * NN * NN1);
	HIST << hnrma << " ";
      }else{
	HIST << 0 << " ";
      }
      HIST << "| ";
    }
    //inter
    int cont = 0;
    for(int m = 0; m < NC-1; ++m){
      int MM = NP.get_Npart(m);
      for(int n = m+1; n < NC; ++n){
	int NN = NP.get_Npart(n);
	hnrme = 1. * Histe[cont][0][bin] / (nnorm * NN * MM);
	HIST << hnrme << " "; 
	hnrme = 3. * Histe[cont][1][bin] / (nnorm * NN * MM);
	HIST << hnrme << " ";
	hnrme = (3./2) * Histe[cont][2][bin] / (nnorm * NN * MM);
	HIST << hnrme << " ";
	hnrme = 5. * Histe[cont][3][bin] / (nnorm * NN * MM);
	HIST << hnrme << " ";
	cont++;
	HIST << "| ";
      }
    }
    HIST << endl;
  }
  HIST.close();
  return true;
}

bool chrg_correl(double& cor_cont, Dynvec<Matrix<double> >& Hista, Dynvec<Matrix<double> >& Histe,
		 const part_num& NP){
  int NC = NP.get_Ncomp();
  int cst, n_bars;
  Hista[0].get_size(cst, n_bars);
  char name[30];
  sprintf(name, "correlation.dat");
  ifstream HIST(name);
  if(! HIST){
    HIST.close();
    return false;
  }
  string line;
  double norm;
  double pos;
  double hnrma, hnrme;
  
  for(int bin = 0; bin < n_bars; ++bin){
    HIST >> cor_cont >> norm >> pos >> line;
    double nnorm = norm * cor_cont;
    //intra
    for(int n = 0; n < NC; ++n){
      int NN = NP.get_Npart(n);
      int NN1 = NN - 1;
      HIST >> hnrma;
      Hista[n][0][bin] = hnrma * nnorm * NN * NN1 / 1.;
      HIST >> hnrma;
      Hista[n][1][bin] = hnrma * nnorm * NN * NN1 / 3.;
      HIST >> hnrma;
      Hista[n][2][bin] = hnrma * nnorm * NN * NN1 / (3./2);
      HIST >> hnrma;
      Hista[n][3][bin] = hnrma * nnorm * NN * NN1 / 5.;
      HIST >> line;
    }
    //inter
    int cont = 0;
    for(int m = 0; m < NC-1; ++m){
      int MM = NP.get_Npart(m);
      for(int n = m+1; n < NC; ++n){
	int NN = NP.get_Npart(n);
	HIST >> hnrme;
	Histe[cont][0][bin] = hnrme * nnorm * NN * MM / 1.;
	HIST >> hnrme;
	Histe[cont][1][bin] = hnrme * nnorm * NN * MM / 3.;
	HIST >> hnrme;
	Histe[cont][2][bin] = hnrme * nnorm * NN * MM / (3./2);
	HIST >> hnrme;
	Histe[cont][3][bin] = hnrme * nnorm * NN * MM / 5.;
	cont++;
	HIST >> line;
      }
    }
  }
  HIST.close();
  return true;
}
//for computing mean square disp. and velocity autocorrelation
bool pos_vel(const Dynvec<Dynvec<Spherocyl> >& part, const Dynvec<int>& np_tc, /*const double L[3],*/
	     double Dt, int num){
  int NC = np_tc.get_size();
  
  char name[130];
  getcwd(name,130);
  sprintf(name+strlen(name),"/pos_vel");
  mkdir(name, 0777);
  sprintf(name+strlen(name), "/pos_vel_%d.dat", num);
  std::ofstream Con(name);
  if(! Con){
    Con.close();
    return false;
  }
  Con << setiosflags(ios::fixed) << setprecision(12);
  double r[3], p[3], I[3];
  Matrix<double> o(3,3); 
  Matrix<double> I0(3,3,0.), Ib(3,3,0.), ll(3,1,0.), omega(3,1,0.);
  for(int n = 0; n < NC; ++n){
    int NN = np_tc[n];
    double m = part[n][0].get_mass();
    for(int i = 0; i < NN; ++i){
      part[n][i].get_pos(r); part[n][i].get_mom(p);
      Con << r[0] << " " << r[1] << " " << r[2] << " "
	  << p[0]/m << " " << p[1]/m << " " << p[2]/m << " | ";
      part[n][i].get_ori(r); part[n][i].get_angm(p);
      part[n][i].get_omat(o); part[n][i].get_inertia(I);
      I0.set_to(0.);
      I0[0][0] = 1./I[0]; I0[1][1] = 1./I[1]; I0[2][2] = 1./I[2];
      //transform I^-1 to box frame
      Ib = transpose(o) * I0 * o;
      for(int a = 0; a < 3; ++a) ll[a][0] = p[a];
      omega = Ib * ll; //angular velocity
      Con << r[0] << " " << r[1] << " " << r[2] << " "
	  << omega[0][0] << " " << omega[1][0] << " " << omega[2][0] 
	  << endl;
    }
  }
  Con << Dt << endl;
  Con.close();
  return true;
}

bool print_ptens(const Dynvec<Dynvec<Spherocyl> >& part, const part_num& NP, double time, 
		 const double PT[3][3]){
  int NC = NP.get_Ncomp();
  double r[3], p[3];
  //to compute Helfand moment (Einstein relation)
  double HFsum[3][3];
  for(int i = 0; i < 3; ++i){
    for(int j = 0; j < 3; ++j)
      HFsum[i][j] = 0.;
  }
  for(int n = 0; n < NC; ++n){
    int NN = NP.get_Npart(n);
    for(int k = 0; k < NN; ++k){
      part[n][k].get_pos(r); part[n][k].get_mom(p);
      for(int i = 0; i < 3; ++i){
	for(int j = 0; j < 3; ++j)
	  HFsum[i][j] += r[i] * p[j];
      }
    }
  }
  
  std::ofstream Out("ptens.dat", ios::app);
  if(! Out){
    Out.close();
    return false;
  }
  Out << time << " ";
  Out << setiosflags(ios::fixed) << setprecision(12);
  for(int i = 0; i < 3; ++i){
    for(int j = 0; j < 3; ++j)
      Out << HFsum[i][j] << " ";
  }
  for(int i = 0; i < 3; ++i){
    for(int j = 0; j < 3; ++j)
      Out << PT[i][j] << " ";
  }
  Out << endl;
  Out.close();
  return true;
}

bool print_magmom(double time, const double M[3]){
  std::ofstream Out("magmom.dat", ios::app);
  if(! Out){
    Out.close();
    return false;
  }
  Out << time << " ";
  Out << setiosflags(ios::fixed) << setprecision(12);
  for(int i = 0; i < 3; ++i){
    Out << M[i] << " ";
  }
  Out << endl;
  Out.close();
  return true;
}

void oo_params(const Dynvec<Dynvec<Spherocyl> >& part, const part_num& NP, double& G1, double& G2,
	       double M[3]){
  int NC = NP.get_Ncomp();
  int Tpart = NP.get_Tpart();
  
  Dipole dip;
  double o[3];
  
  G1 = 0., G2 = 0.;
  float* Q[4];
  float Q_0[4], Q_1[4], Q_2[4], Q_3[4];
  float* V_[4];
  float V_0[4], V_1[4], V_2[4], V_3[4];
  Q[0] = Q_0; Q[1] = Q_1; Q[2] = Q_2; Q[3] = Q_3; 
  V_[0] = V_0; V_[1] = V_1; V_[2] = V_2; V_[3] = V_3;
  float d[4], Ve[3], Vesq;
  int nrot;
  
  //total magnetic moment
  for(int n = 0; n < NC; ++n){
    int NN = NP.get_Npart(n);
    for(int i = 0; i < NN; ++i){
      int size = part[n][i].get_ndips();
      for(int s = 0; s < size; ++s){
	dip = part[n][i].get_dip(s);
	dip.get_ori(o);
	for(int l = 0; l < 3; ++l)
	  M[l] += dip.get_dip() * o[l];
      }
    }
  }
    
  //ori. order parameters
  Qmat_3d(part, NP, Q);
  jacobi(Q, 3, d, V_, &nrot);
  
  float e1 = fabs(d[1]), e2 = fabs(d[2]), e3 = fabs(d[3]);
  if(e1 >= e2 && e1 >= e3){
    G2 += static_cast<double>(e1);
    Ve[0] = V_[1][1]; Ve[1] = V_[2][1]; Ve[2] = V_[3][1];
  }else if(e2 > e1 && e2 >= e3){
    G2 += static_cast<double>(e2);
    Ve[0] = V_[1][2]; Ve[1] = V_[2][2]; Ve[2] = V_[3][2];
  }else{
    G2 += static_cast<double>(e3);
    Ve[0] = V_[1][3]; Ve[1] = V_[2][3]; Ve[2] = V_[3][3];
  }
  Vesq = sqrt(Ve[0] * Ve[0] + Ve[1] * Ve[1] + Ve[2] * Ve[2]);
  Ve[0] /= Vesq; Ve[1] /= Vesq; Ve[2] /= Vesq;
  
  for(int n = 0; n < NC; ++n){
    int NN = NP.get_Npart(n);
    for(int i = 0; i < NN; ++i){
      part[n][i].get_ori(o);
      G1 += Ve[0] * o[0] + Ve[1] * o[1] + Ve[2] * o[2];
    }
  }
  G1 = fabs(G1) / Tpart;
}

bool print_out(double Time, const accumulator& var, long telap, long tcpu, int iprint, 
	       const double Mom[6]){
  //time difference evaluation
  double t_elap = (telap / iprint) * 1e-6;
  double t_cpu = (static_cast<double>(tcpu) / iprint) * 1e-6;
  //print
  char name[30];
  sprintf(name, "output.dat");
  ofstream Out(name, ios::app);
  if(! Out){
    Out.close();
    return false;
  }
  Out << setiosflags(ios::fixed) << setprecision(12);
  Out << Time << ". " << t_elap << " " << t_cpu << " | "
      << var.get_K() << " " << var.get_U() << " " << var.get_P() << " " 
      << var.get_Uss() << " " << var.get_Pss() << " " 
      << var.get_UdipR() << " " << var.get_PdipR() << " "
      << var.get_UdipK() << " " << var.get_PdipK() << " "
      << Mom[0]*Mom[0]+Mom[1]*Mom[1]+Mom[2]*Mom[2] << " "
      << Mom[3]*Mom[3]+Mom[4]*Mom[4]+Mom[5]*Mom[5] << endl;
  Out.close();
  return true;
}

bool print_final(double Time, double avK, double flK, double avU, double flU, double avP, 
		 double flP, double avUss, double flUss, double avPss, double flPss, 
		 double avUdipR, double flUdipR, double avPdipR, double flPdipR, 
		 double avUdipK, double flUdipK, double avPdipK, double flPdipK, 
		 double avUext, double flUext,
		 double avPT[3][3], double flPT[3][3],
		 double avPTss[3][3], double flPTss[3][3],
		 double avPTdipR[3][3], double flPTdipR[3][3],
		 double avPTdipK[3][3], double flPTdipK[3][3],
		 double avPTkin[3][3], double flPTkin[3][3],
		 double avG1, double flG1, double avG2, double flG2,
		 const Dynvec<double>& length, const Dynvec<double>& diameter, 
		 const Dynvec<double>& mass, const Dynvec<Matrix<double> >& dip, 
		 const double L[3], const part_num& NP, double temp, double Q, 
		 double gamma, const double B[3], double pol_deg){
  int NC = NP.get_Ncomp();
  int TP = NP.get_Tpart();
  //print
  char name[30];
  sprintf(name, "output.dat");
  ofstream Out(name, ios::app);
  if(! Out){
    Out.close();
    return false;
  }
  Out << "NVT ensemble\n";
  Out << "Temperature: " << temp << endl;
  double V = L[0] * L[1] * L[2];
  Out << "Density: " << TP / V << endl;
  Out << "Volume (Lx,Ly,Lz): " 
      << L[0] << " x " << L[1] << " x " << L[2] << " = " 
      << V << endl;
  Out << "Elapsed time: " << Time << endl;
  Out << "Component | No. of particles | mass | L | D | Vol. frac.\n"; 
  for(int n = 0; n < NC; ++n){
    double rad = diameter[n] / 2;
    int NN =NP.get_Npart(n);
    Out << n << " | " << NN << " | " << mass[n] << " | " 
	<< length[n] << " | " << diameter[n] << " | " 
	<< NN * M_PI * (4.*rad/3 + length[n]) * rad * rad / V << endl;
  }
  Out << "Dipoles: (position ; orientation ; dip.moment)\n";
  for(int n = 0; n < NC; ++n){
    int a, b;
    dip[n].get_size(a, b);
    Out << "species " << n << endl;
    for(int j = 0; j < a; ++j)
      Out << dip[n][j][0] << " ; " 
	  << dip[n][j][1] << " " << dip[n][j][2] << " " <<  dip[n][j][3]
	  << " ; " << dip[n][j][4] << endl;
  }
  Out << "External field: (" << B[0] << "," << B[1] << "," << B[2] << ")\n";
  Out << "Thermostat mass: " << Q << " Thermostat shear rate: " << gamma << endl;
  Out << "------------------\n";
  Out << "Av. kinetic energy: " << avK << "+/-" << flK << endl;
  Out << "Av. potential energy: " << avU << "+/-" << flU << endl;
  Out << "Av. soft sphere en.: " << avUss << "+/-" << flUss 
      << " Av. dipolarR en.: " << avUdipR << "+/-" << flUdipR 
      << " Av. dipolarK en.: " << avUdipK << "+/-" << flUdipK << endl;
  Out << "Av. Ext.Field en.: " << avUext << "+/-" << flUext << endl;
  Out << "Av. pressure: " << avP << "+/-" << flP << endl;
  Out << "Av. soft sphere pres.: " << avPss << "+/-" << flPss 
      << " Av. dipolarR pres.: " << avPdipR << "+/-" << flPdipR
      << " Av. dipolarK pres.: " << avPdipK << "+/-" << flPdipK << endl;
  Out << "Av. Pressure tensor: \n";
  Out << "| " << avPT[0][0]<<"+/-"<< flPT[0][0] << " " << avPT[0][1]<<"+/-"<< flPT[0][1] << " " << avPT[0][2]<<"+/-"<< flPT[0][2] << " |\n"
      << "| " << avPT[1][0]<<"+/-"<< flPT[1][0] << " " << avPT[1][1]<<"+/-"<< flPT[1][1] << " " << avPT[1][2]<<"+/-"<< flPT[1][2] << " |\n"
      << "| " << avPT[2][0]<<"+/-"<< flPT[2][0] << " " << avPT[2][1]<<"+/-"<< flPT[2][1] << " " << avPT[2][2]<<"+/-"<< flPT[2][2] << " |\n";
  Out << "kinetic part: \n";
  Out << "| " << avPTkin[0][0]<<"+/-"<< flPTkin[0][0] << " " << avPTkin[0][1]<<"+/-"<< flPTkin[0][1] << " " << avPTkin[0][2]<<"+/-"<< flPTkin[0][2] << " |\n"
      << "| " << avPTkin[1][0]<<"+/-"<< flPTkin[1][0] << " " << avPTkin[1][1]<<"+/-"<< flPTkin[1][1] << " " << avPTkin[1][2]<<"+/-"<< flPTkin[1][2] << " |\n"
      << "| " << avPTkin[2][0]<<"+/-"<< flPTkin[2][0] << " " << avPTkin[2][1]<<"+/-"<< flPTkin[2][1] << " " << avPTkin[2][2]<<"+/-"<< flPTkin[2][2] << " |\n";
  Out << "soft core part: \n";
  Out << "| " << avPTss[0][0]<<"+/-"<< flPTss[0][0] << " " << avPTss[0][1]<<"+/-"<< flPTss[0][1] << " " << avPTss[0][2]<<"+/-"<< flPTss[0][2] << " |\n"
      << "| " << avPTss[1][0]<<"+/-"<< flPTss[1][0] << " " << avPTss[1][1]<<"+/-"<< flPTss[1][1] << " " << avPTss[1][2]<<"+/-"<< flPTss[1][2] << " |\n"
      << "| " << avPTss[2][0]<<"+/-"<< flPTss[2][0] << " " << avPTss[2][1]<<"+/-"<< flPTss[2][1] << " " << avPTss[2][2]<<"+/-"<< flPTss[2][2] << " |\n";
  Out << "real dip part: \n";
  Out << "| " << avPTdipR[0][0]<<"+/-"<< flPTdipR[0][0] << " " << avPTdipR[0][1]<<"+/-"<< flPTdipR[0][1] << " " << avPTdipR[0][2]<<"+/-"<< flPTdipR[0][2] << " |\n"
      << "| " << avPTdipR[1][0]<<"+/-"<< flPTdipR[1][0] << " " << avPTdipR[1][1]<<"+/-"<< flPTdipR[1][1] << " " << avPTdipR[1][2]<<"+/-"<< flPTdipR[1][2] << " |\n"
      << "| " << avPTdipR[2][0]<<"+/-"<< flPTdipR[2][0] << " " << avPTdipR[2][1]<<"+/-"<< flPTdipR[2][1] << " " << avPTdipR[2][2]<<"+/-"<< flPTdipR[2][2] << " |\n";
  Out << "reciprocal dip part: \n";
  Out << "| " << avPTdipK[0][0]<<"+/-"<< flPTdipK[0][0] << " " << avPTdipK[0][1]<<"+/-"<< flPTdipK[0][1] << " " << avPTdipK[0][2]<<"+/-"<< flPTdipK[0][2] << " |\n"
      << "| " << avPTdipK[1][0]<<"+/-"<< flPTdipK[1][0] << " " << avPTdipK[1][1]<<"+/-"<< flPTdipK[1][1] << " " << avPTdipK[1][2]<<"+/-"<< flPTdipK[1][2] << " |\n"
      << "| " << avPTdipK[2][0]<<"+/-"<< flPTdipK[2][0] << " " << avPTdipK[2][1]<<"+/-"<< flPTdipK[2][1] << " " << avPTdipK[2][2]<<"+/-"<< flPTdipK[2][2] << " |\n";  
  Out << "Average G1(P) parameter:\n";
  Out << avG1 << " +/- " << flG1 << "\n";
  Out << "Average G2(S) mayer-saupe parameter:\n";
  Out << avG2 << " +/- " << flG2 << "\n";
  Out << "Polymerization degree: " << pol_deg << endl;
  Out.close();
  return true;
}

//cluster************************************

//Equal probability for all particles
void rnd_part(const std::vector<std::vector<int> >& NI, int& comp, 
              int& I, int& loc, const gsl_rng* ran){
  int Tot = 0;
  for(unsigned n = 0; n < NI.size(); ++n){
    Tot += NI[n].size();
  }
  double rnd = gsl_rng_uniform(ran);
  int ip = static_cast<int>(floor(rnd*Tot));
  int Nto=0, Ntn=0;
  for(unsigned n = 0; n < NI.size(); ++n){
    Nto = Ntn;
    int Np = NI[n].size();
    Ntn += Np;
    if(ip < Ntn){
      comp = n;
      loc = ip - Nto;
      I = NI[n][loc];
      break;
    }
  }
}

// //Part. in a less numerous species are more probable
// void rnd_part(const std::std::vector<std::std::vector<int> >& NI, int& comp, 
// 	      int& I, int& loc, const gsl_rng* ran){
//   double rnd;
//   int size = 0;
//   while(size < 1){
//     rnd = gsl_rng_uniform(ran);
//     comp = static_cast<int>(floor(rnd*NI.size()));
//     size = NI[comp].size();
//   }
//   rnd = gsl_rng_uniform(ran);
//   loc = static_cast<int>(floor(rnd*size));
//   I = NI[comp][loc];
// }

bool dist_test(const Spherocyl& A, const Spherocyl& B, const double L[3], double min_s){
  double ri[3], rj[3]; 
  A.get_pos(ri);
  B.get_pos(rj);
  
  double shi[2], shj[2];
  A.get_shape(shi);
  B.get_shape(shj);
  double sigma = (shi[1] + shj[1]) / 2 + min_s;
  double sig_sq = sigma * sigma;
  
  double rij[3]/*, rsq=0*/;
  for(int a = 0; a < 3; ++a){
    rij[a] = ri[a] - rj[a];
    rij[a] -= L[a] * floor(rij[a] / L[a] + 0.5);
//     rsq += rij[a]*rij[a];
  }
  
  double st[2], dij[3];
  double dsq = spherocyl_dist(A, B, L, rij, dij, st);
  
  if(dsq < sig_sq){
//     //
//   cout<<sig_sq<<" < "<<dsq<<endl;
//   //
    return true;
  }

  return false;
}

int clust_gen(Dynvec<Dynvec<Spherocyl> >& part, 
	      const part_num& NP, std::vector<std::vector<int> >& NI, 
	      Dynvec<bool>& added, std::vector<std::vector<cl_element> >& memb, 
	      const double L[3], double min_s, const gsl_rng* ran, 
	      int n_clust){
  int NC = NP.get_Ncomp();
  int TP = NP.get_Tpart();

  //temporal element container
  cl_element t_el;
  t_el.link.reserve(TP); //to reduce allocation time
  
  //choose seed particle
  int comp=0, I=0, Loc=0;
  rnd_part(NI, comp, I, Loc, ran);
  
  //track of added to this cluster
  Dynvec<bool> locadd(TP);
  for(int k = 0; k < TP; ++k)
    locadd[k] = false;
  
  //set seed added to true
  int II = calc_II(NP, comp, I);
  added[II] = true;
  locadd[II] = true;
  part[comp][I].set_clust(n_clust);
  int cl_size = 1;
  
  //stack. Add seed to the stack
  deque<int> queu[2];       //deque because I need pop_front
  queu[0].clear(); queu[1].clear();
  queu[0].push_back(comp);
  queu[1].push_back(I);
  
  //store members
  std::vector<cl_element> local;
  local.resize(TP);
  local.clear();
  t_el.type = comp; t_el.num = I;
  t_el.gen = 0;
//   t_el.end = false;
  t_el.state = false;
  t_el.link.clear();
  local.push_back(t_el);
  
  //erase added particles from NI
  NI[comp].erase(NI[comp].begin()+Loc);
//    //
//     cout<<NI[comp].size()<<"("<<comp<<","<<I<<")l"<<Loc<<" ";
//     //
  
  // 	// //
  // 	// 	    cout<<"i "<<queu[0].size()<<endl;
  // 	// 	    //
  
  ////construct cluster
  int index = 0;
  while(! queu[0].empty()){ //while the stack is not empty..
    int nn = queu[0][0];
    int ii = queu[1][0];
   
    //check with all non-added interacting partilces 
    // (according to cluster criteria)
    II = 0;
    for(int o = 0; o < NC; ++o){
      int NN = NP.get_Npart(o);
      for(int j_ = 0; j_ < NN; ++j_){
	if(! locadd[II]){
	  //apply criteria (calculate the bond forming probability 
	  // ratio p(i->f)/p(f->i)) 
	  bool test = dist_test(part[nn][ii], part[o][j_], L, min_s);
	  // //
	    // 			cout<<"("<<test<<"a"<<cl_size<<")";
	    // 			//
	  
	  if(test){
	    if(! added[II]){
	      //if criteria accepted add to cluster
	      // and the bond energy remains the same
	      queu[0].push_back(o);
	      queu[1].push_back(j_);
	      t_el.type = o; t_el.num = j_;
	      t_el.gen = local[index].gen + 1; 
	      local.push_back(t_el);
	      local[index].link.push_back(local.size()-1);
	      int loc = 0;
	      for(unsigned i = 0; i < NI[o].size(); ++i){
		if(NI[o][i] == j_){
		  loc = i;
		  break;
		}
	      }
	      //erase added particles from NI
	      NI[o].erase(NI[o].begin()+loc);
	      // //
	      // 		  cout<<"("<<o<<","<<j_<<","<<loc<<") ";
	      // 		  //
	      locadd[II] = true;
	      added[II] = true;
	      
	      part[o][j_].set_clust(n_clust);
	      cl_size++;
	    }else{
	      ofstream Err("Error.dat", ios::app);
	      Err << "particle (" << o << ","
		  << j_ << ") classifies for 2 clusters!" 
		  << endl;
	      Err.close();
	    }
	  }
	}
	II++;
      }
    }
    if(! queu[0].empty()){
      queu[0].pop_front();
      queu[1].pop_front();
    }
//     //NOTE: this is not the real connectivity as particles already in the 
//     // cluster are not checked twice! This is just a particular realization
//     // of the cluster.
//     int connectivity = local[index].link.size() + 1;
//     if(index == 0)
//       connectivity = local[index].link.size();
//     if(connectivity < 2)
//       local[index].end = true;
    index++;
  }
  memb.push_back(local);
  return cl_size;
}

void gb_check(const Dynvec<Dynvec<Spherocyl> >& part, const part_num& NP, std::vector<cl_element>& memb,
	      const double L[3], double min_s, int k){
  //int NC = NP.get_Ncomp();
  int m = memb[k].type;
  int i = memb[k].num;
  
  for(unsigned l = 0; l < memb.size(); ++l){
    bool test = true;
    
    if(memb[k].gen > memb[l].gen){ //if the generation of the other particle is less than than of the current one...
      test = false;
    }
    
    for(unsigned t = 0 ; t < memb[k].link.size(); ++t){
      if(static_cast<int>(l) == memb[k].link[t]) //if there is already a non-ghost bond...
	test = false;
    }
 
    if(test){
      int n = memb[l].type;
      int j = memb[l].num;
      if((m != n) || ((m == n) && (i != j))){
	if(dist_test(part[m][i], part[n][j], L, min_s)){
	  memb[k].link.push_back(l); // add the new bond
// 	  //
// 	  cout<<"add_gba{"/*<<m<<","<<i<<":"<<n<<","<<j<<"|"*/<<k<<"-"<<l<<"} ";
// 	  //
	}
      }
    }
  }
}

bool percolation(const Dynvec<Dynvec<Spherocyl> >& part, const part_num& NP, 
		 std::vector<cl_element>& memb, const double L[3], double dmax, double min_s){
  double Lmin = min(L[2], min(L[1], L[0]));
  double mm = Lmin - (dmax + min_s); //particles below this distance from the initial particle cannot for a bond that closes aroun the PBC(dmax is the largest dimension of a particle in the system)
  int min_n = ceil(mm / (dmax + min_s));
  //minimum generation whose distance from the initial particle is large enough to go through the PBC
  int min_gen = 0;
  if(fmod(static_cast<float>(min_n), 2) == 0){
    min_gen = min_n / 2 - 1;
  }else
    min_gen = (min_n - 1) / 2;
  
//   //
//   cout<<"mm:"<<mm<<" min_n:"<<min_n<<" min_gen:"<<min_gen<<endl;  
//   //
  
  bool percol = false;
  double Dr[3] = {0.0}; 
  std::vector<Walker> wkers;
  Walker wk;
  wk.set_actual(0);
  wk.set_next(memb[wk.get_actual()].link[0]);
  wkers.push_back(wk);
  int w = 0;
  
  while(wkers.size() > 0){
    int k = wkers[w].get_actual();
    wkers[w].get_Dr(Dr);
    for(int i = 0; i < 3; ++i)
      memb[k].dis[i] = Dr[i]; 
    
//     //
//     cout<<"w:"<<w<<" ("<<Dr[0]<<","<<Dr[1]<<","<<Dr[2]<<") "
// 	<<"clus:"<<n_clust<<" part(k:"<<k<<",s:"<<memb[n_clust][k].state<<",g:"<<memb[n_clust][k].gen<<") ";
//     //
    
    if(memb[k].link.size() == 0){ //if no next element
      memb[k].state = true;  //mark as visited
      wkers.erase(wkers.begin()+w); //erase walker from the path
      w = wkers.size()-1; // take anothe walker
//       //
//       cout<<" kill."<<wkers.size()<<endl;
//       //
    }else if(memb[k].link.size() > 0){ //if more elements follow

      if(! memb[k].state){  //if not visited
	for(unsigned l = 1; l < memb[k].link.size(); ++l){ //create new walkers redy to go from the new nodes
	  wk = wkers[w];  //This line...
	  wkers.push_back(wk); //and this line are to add a copy of the element and then...
	  wkers[wkers.size()-1].set_next(memb[k].link[l]); //put a different "next element"
// 	  //
// 	  cout<<" rep("<<k<<") ";
// 	  //
	}
      }
//       //
//       cout<<"[act:"<<wkers[w].get_actual()<<",nxt:"<<wkers[w].get_next()<<"]<="<<memb[n_clust].size()-1<<" ";
//       //
      wkers[w].move(part, L, memb); //measure distance and go to next element
      int l = wkers[w].get_actual();
      //check for ghost bonds at new node (NOTE: if a gb is formed between parts i and j it will be formed again between i and j
      // so there is a posibility for a walker to travel backwards through this bond before dying. This cosumes one additional 
      // distance computation, but does not seems to be a problem beyond that)
      if(memb[l].gen >= min_gen && ! memb[l].state)
	gb_check(part, NP, memb, L, min_s, l); //check for ghost bonds
      if(memb[l].link.size() > 0) //put the links of node l in the path
	wkers[w].set_next(memb[l].link[0]);
    
      //set old node as visited
      memb[k].state = true;
	
      if(memb[wkers[w].get_actual()].state){ //if we passed by this new node already...
	wkers[w].get_Dr(Dr);
// 	//
// 	cout<<"loop-"/*<<"Dr("<<Dr[0]<<","<<Dr[1]<<","<<Dr[2]<<")d("
// 	    <<memb[l].dis[0]<<","<<memb[l].dis[1]<<","<<memb[l].dis[2]<<") "*/;  
// 	//
	double difX = Dr[0] - memb[l].dis[0];
	double difY = Dr[1] - memb[l].dis[1];
	double difZ = Dr[2] - memb[l].dis[2];
	if((fabs(difX) > 1e-8) || (fabs(difY) > 1e-8) || (fabs(difZ) > 1e-8)){ //if this is not a non-spaning closed loop...
	  percol = true;
// 	  //
// 	  cout<<"percol! ";  
// 	  //
	}
	wkers.erase(wkers.begin()+w); //erase this walker
	w = wkers.size()-1; //take another walker
// 	//
// 	cout<<"-kill."<<wkers.size()<<endl;
// 	//
      }
      if(percol)
	break;
    }
  }
  return percol;
}

bool clust_conf(Dynvec<Dynvec<Spherocyl> >& part, const part_num& NP, const double L[3],
		double min_s, double dmax, const gsl_rng* ran, long Npart_in_clust[2], bool p_check){
  int NC = NP.get_Ncomp();
  int TP = NP.get_Tpart();
  std::vector<std::vector<cl_element> > memb;
  memb.reserve(TP);
  memb.clear(); 

  //particle index array
  std::vector<std::vector<int> > NI(NC);
  for(int n = 0; n < NC; ++n){
    int NN = NP.get_Npart(n);
    NI[n].reserve(NN);
    NI[n].clear();
    for(int j = 0; j < NN; ++j)
      NI[n].push_back(j);
  }
  //To take track of added particle, maybe redundant with member m_clust of spheroid
  Dynvec<bool> added(TP);
  for(int k = 0; k < TP; ++k)
    added[k] = false;      
  //generate clusters
  int n_clust = 0;
  int tsize = 1;
  while(tsize > 0){
    clust_gen(part, NP, NI, added, memb, L, min_s, ran, n_clust);
    tsize = 0;
    for(unsigned o = 0; o < NI.size(); ++o){
      tsize += NI[o].size();
    }
    n_clust++;
  }
  
  //degree of polymerization
  Npart_in_clust[0]++;
  for(unsigned n = 0; n < memb.size(); ++n){
    if(memb[n].size() > 1)
      Npart_in_clust[1] += memb[n].size();
  }
  
  //check for percolation
  bool perc = false;
  if(p_check){
    for(unsigned n = 0; n < memb.size(); ++n){
      if(percolation(part, NP, memb[n], L, dmax, min_s))
	perc = true;
    }
  }
  
  return perc;
}

////////////////////////////////

void GCAmix(Dynvec<Dynvec<Spherocyl> >& part, const part_num& NP, const double L[3], 
	    const gsl_rng* ran){
  int NC = NP.get_Ncomp();
  int Tpart = NP.get_Tpart();
  //added particles
  bool added[Tpart];
  for(int i = 0; i < Tpart; ++i)
    added[i] = false;
  //choose pivot
  double pivot[3];
  double rnd = gsl_rng_uniform(ran);
  pivot[0] = (2 * rnd - 1) * L[0] / 2;
  rnd = gsl_rng_uniform(ran);
  pivot[1] = (2 * rnd - 1) * L[1] / 2;
  rnd = gsl_rng_uniform(ran);
  pivot[2] = (2 * rnd - 1) * L[2] / 2;
  //choose seed particle at random
  rnd = gsl_rng_uniform(ran);
  int id = static_cast<int>(floor(rnd*Tpart));
  int nn = 0, ii = 0, Nto = 0;
  for(int n = 0; n < NC; ++n){
    int NN = NP.get_Npart(n);
    Nto += NN;
    if(id < Nto){
      nn = n;
      int loc = Nto - id;
      ii = NN - loc;
      break;
    }
  }
  int II = part[nn][ii].get_id();
  added[II] = true;
  deque<int> queue[2];
  queue[0].push_back(nn);
  queue[1].push_back(ii);
  
  //construct cluster
  double ri[3], rj[3], rij[3];
  while(queue[0].size() > 0){
    nn = queue[0][0];
    ii = queue[1][0];
    //reflect
    part[nn][ii].get_pos(ri);
    for(int a = 0; a < 3; ++a){
      rij[a] = pivot[a] - ri[a];
      rij[a] -= L[a] * floor(rij[a] / L[a] + 0.5);
      
      ri[a] += 2*rij[a];
      ri[a] -= L[a] * floor(ri[a] / L[a] + 0.5);
    }
    part[nn][ii].set_pos(ri);
    
    for(int n = 0; n < NC; ++n){
      double shn[2]; part[n][0].get_shape(shn);
      double shnn[2]; part[nn][0].get_shape(shnn);
      double dd = (shnn[1] + shn[1]) / 2;
      double rsq = dd * dd;
      int NN = NP.get_Npart(n);
      for(int i = 0; i < NN; ++i){
	II = part[n][i].get_id();
	if(! added[II]){
	  part[n][i].get_pos(rj);
	  
	  for(int a = 0; a < 3; ++a){
	    rij[a] = ri[a] - rj[a];
	    rij[a] -= L[a] * floor(rij[a] / L[a] + 0.5);
	  }
	  double st[2], dij[3];
	  double dsq = spherocyl_dist(part[nn][ii], part[n][i], L, rij, dij, st);
	    
	  if(dsq < rsq){
	    queue[0].push_back(n);
	    queue[1].push_back(i);
	    added[II] = true;
	  }
	}
      }
    }
    queue[0].pop_front();
    queue[1].pop_front();
  }
}

void dot(const double a[3], const double b[3], double& ans){
  ans = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

void cros(const double a[3], const double b[3], double ans[3]){
  ans[0] = a[1]*b[2] - a[2]*b[1];
  ans[1] = a[2]*b[0] - a[0]*b[2];
  ans[2] = a[0]*b[1] - a[1]*b[0];
}

int delta(int a, int b){
  if(a == b) return 1;
  else return 0;
}

void Qmat_3d(const Dynvec<Dynvec<Spherocyl> >& part, const part_num& NP, float** Q){
  int NC = NP.get_Ncomp();
  int Tpart = NP.get_Tpart();
  for(int i = 0; i < 4; i++){
    for(int j = 0; j < 4; j++){
      Q[i][j] = 0;
    }
  }
  double o[3];
  for(int n = 0; n < NC; ++n){
    int Npart = NP.get_Npart(n);
    for(int k = 0; k < Npart; ++k){
      part[n][k].get_ori(o);
      for(int i = 1; i < 4; ++i){
	for(int j = 1; j < 4; ++j){
	  Q[i][j] += (static_cast<float>(3 * o[i-1] * o[j-1]) 
		      - static_cast<float>(delta(i, j))) / (2 * Tpart);
	}
      }
    }
  }
  // for(int i = 1; i < 4; ++i){
  //     for(int j = 1; j < 4; ++j){
  //       cout << Q[i][j] << " ";
  //     }
  //     cout << endl;
  //   }
}

///////////////////////////////

//Matrix functions********************************************************
// err = 0: bad size, err = -1: not square, err = 1: matrix singular

Matrix<double> transpose(const Matrix<double>& RESTRICT M){
  Matrix<double> result(M.m_cols, M.m_rows);
  for(int m = 0; m < M.m_cols; ++m){
    for(int n = 0; n < M.m_rows; ++n){
      result[m][n] = M[n][m];
    }
  }
  return result;
}

Matrix<double> LU_dec(const Matrix<double>& RESTRICT M, int* RESTRICT order, 
		      int& sign, bool& reg){
  if(M.m_rows != M.m_cols){
    Matrix<double> err(1, 1, -1);
    return err;
  }
  for(int i = 0; i < M.m_rows; ++i)
    order[i] = i;
  reg = true;
  sign = 1;
  Matrix<double> cpy(M.m_rows, M.m_cols);
  cpy = M;
  Matrix<double> LU(M.m_rows, M.m_cols);
  //separate in lower and upper (LU, LxU=M) triangular matrices and put in LU
  //-----------
  //U00 U01 U02
  //L10 U11 U12
  //L20 L21 U22
  //-----------
  //L00=L11=L22=1
  for(int m = 0; m < M.m_rows; ++m){
    //swap rows if needed
    if(cpy[m][m] == 0){
      bool test = true;
      int i = 1;
      while(test){
	if(cpy[m+i][m] != 0){
	  double temp;
	  for(int n = 0; n < M.m_cols; ++n){
	    temp = cpy[m][n];
	    cpy[m][n] = cpy[m+1][n];
	    cpy[m+1][n] = temp;
	  }
	  int o = order[m];
	  order[m] = order[m+i];
	  order[m+i] = o;
	  sign *= -1;
	  test = false;
	}
      }
      reg = false;
    }
    //Crout's method
    for(int n = 0; n < M.m_cols; ++n){
      if(m <= n){
	LU[m][n] = cpy[m][n];
	for(int i = 0; i < m; ++i)
	  LU[m][n] -= LU[m][i] * LU[i][n];
      }
      else{
	LU[m][n] = cpy[m][n];
	for(int j = 0; j < n; ++j)
	  LU[m][n] -= LU[m][j] * LU[j][n];
	if(LU[m][n] != 0)
	  LU[m][n] /= LU[n][n];
      }
    }
  }
  return LU;
}

Matrix<double> fwd_bkd(const Matrix<double>& RESTRICT M, 
		       const double* RESTRICT b){
  if(M.m_rows != M.m_cols){
    Matrix<double> err(1, 1, -1);
    return err;
  }    
  double y[M.m_rows];
  //Forward substitution
  y[0] = b[0];
  for(int m = 1; m < M.m_rows; ++m){
    y[m] = b[m];
    for(int i = 0; i < m; ++i)
      y[m] -= M[m][i] * y[i]; 
  }
  Matrix<double> x(M.m_rows, 1);
  //Backward substitution
  x[M.m_rows-1][0] = y[M.m_rows-1] / M[M.m_rows-1][M.m_rows-1];
  for(int m = M.m_rows-2; m >= 0; --m){
    x[m][0] = y[m];
    for(int i = m+1; i < M.m_rows; ++i)
      x[m][0] -= M[m][i] * x[i][0];
    if(x[m][0] != 0)
      x[m][0] /= M[m][m];
  }
  return x;
}

Matrix<double> solve(const Matrix<double>& RESTRICT M, 
		const Matrix<double>& RESTRICT b){
  //solves A x = b.
  if(M.m_rows != M.m_cols){
    Matrix<double> err(2, 2, -1);
    return err;
  } 
  if(b.m_cols != 1 || b.m_rows != M.m_rows){
    Matrix<double> err(1, 1, 0);
    return err;
  }
  Matrix<double> LU(M.m_rows, M.m_cols);
  bool reg;
  int sign;
  int order[M.m_rows];
  LU = LU_dec(M, order, sign, reg);
  double cpyb[M.m_rows];
  for(int i = 0; i < M.m_rows; ++i)
    cpyb[i] = b[order[i]][0];
  Matrix<double> x(M.m_rows, 1);
  x = fwd_bkd(LU, cpyb);
  return x;
}

Matrix<double> inverse(const Matrix<double>& RESTRICT M){
  if(M.m_rows != M.m_cols){
    Matrix<double> err(1, 1, -1);
    return err;
  } 
  Matrix<double> inverse(M.m_rows, M.m_cols);
  Matrix<double> LU(M.m_rows, M.m_cols);
  bool reg;
  int sign;
  int order[M.m_rows];
  LU = LU_dec(M, order, sign, reg);
  if(!reg){
    Matrix<double> err(1, 1, 1);
    return err;
  }
  double col[M.m_rows];
  for(int i = 0; i < M.m_rows; ++i)
    col[i] = 0;
  for(int n = 0; n < M.m_cols; ++n){
    if(n != 0)
      col[n-1] = 0;
    col[n] = 1;
    Matrix<double> y(M.m_rows, 1);
    y = fwd_bkd(LU, col);
    for(int i = 0; i < M.m_rows; ++i)
      inverse[i][n] = y[i][0];
  }
  return inverse;
}
