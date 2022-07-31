#include <iostream>
#include <fstream>
#include <deque>
#include <iomanip>
#include <string>
using namespace std;

#include "Matrix.hpp"

bool diffusion(int Nfiles, const deque<int>& np_tc){ //v autocorrel is non-normalized!
  int NC = np_tc.size();
  string line;
  char name[30];
  
  Dynvec<Matrix<double> > r0, v0;
  r0.set_size(NC);
  v0.set_size(NC);
  for(int n = 0; n < NC; ++n){
    int NN = np_tc[n];
    r0[n].set_size(NN, 3);
    v0[n].set_size(NN, 3);
  }
  double rt[3], vt[3], d[3], Dt;
  
  Matrix<int> norm(NC,Nfiles+1, 0); 
  Dynvec<Matrix<double> > msd(NC);
  Dynvec<Matrix<double> > vacf(NC);
  for(int m = 0; m < NC; ++m){
    msd[m].set_size(Nfiles+1,3);
    msd[m].set_to(0.);
    vacf[m].set_size(Nfiles+1,3);
    vacf[m].set_to(0.);
  }
  
  for(int t0 = 0; t0 < Nfiles; ++t0){
    sprintf(name, "pos_vel/pos_vel_%d.dat", t0);
    std::ifstream Con(name);
    if(! Con){
      Con.close();
      return false;
    }
    
    for(int m = 0; m < NC; ++m){
      int MM = np_tc[m];
      for(int i = 0; i < MM; ++i){
	for(int k = 0; k < 3; ++k){
	  Con >> r0[m][i][k];          
	}                              
	for(int k = 0; k < 3; ++k){
	  Con >> v0[m][i][k];         
	}
	for(int k = 0; k < 7; ++k){
	  Con >> line;
	}
      }
    }
    if(t0 == 0)
      Con >> Dt;    
    Con.close();
    
//     for(int m = 0; m < NC; ++m){
//       int MM = np_tc[m];
//       for(int i = 0; i < MM; ++i){
// 	for(int k = 0; k < 3; ++k)
// 	  cout << r0[m][i][k] << " ";
// 	for(int k = 0; k < 3; ++k){
// 	  cout << v0[m][i][k] << " ";
// 	}
// 	for(int k = 0; k < 7; ++k){
// 	  cout << " . ";
// 	}
// 	cout<<endl;
//       }
//       cout << Dt << endl;
//     }
    
    //autocorrelation Dt=0
    for(int n = 0; n < NC; ++n){
      int NN = np_tc[n];
      for(int j = 0; j < NN; ++j){
	for(int k = 0; k < 3; ++k){
	  d[k] = 0.;
	  msd[n][0][k] += 0.;
	  vacf[n][0][k] += v0[n][j][k]*v0[n][j][k];
	}
	norm[n][0]++;
      }
    }
    //correlation with other times
    for(int t = t0+1; t < Nfiles; ++t){
      sprintf(name, "pos_vel/pos_vel_%d.dat", t);
      std::ifstream Con1(name);
      if(! Con1){
	Con1.close();
	return false;
      }
      
      for(int n = 0; n < NC; ++n){
	int NN = np_tc[n];
	for(int j = 0; j < NN; ++j){
	  for(int l = 0; l < 3; ++l)
	    Con1 >> rt[l];
	  for(int l = 0; l < 3; ++l)
	    Con1 >> vt[l];
	  for(int k = 0; k < 7; ++k)
	    Con1 >> line;
	  
	  for(int k = 0; k < 3; ++k){
	    d[k] = rt[k]-r0[n][j][k];
	    msd[n][t-t0][k] += d[k]*d[k];
	    vacf[n][t-t0][k] += v0[n][j][k]*vt[k];
	  }
	  norm[n][t-t0]++;
	}
      }
      Con1.close();
    }
  }
  
  sprintf(name, "difudata.dat");
  std::ofstream Out(name);
  if(! Out){
    Out.close();
    return false;
  }
  Out << setiosflags(ios::fixed) << setprecision(6);
  for(int t = 0; t < Nfiles; ++t){
    Out << t*Dt << " |";
    for(int n = 0; n < NC; ++n){ 
      Out << "|r2 ";
      for(int k = 0; k < 3; ++k)
	Out << msd[n][t][k]/norm[n][t] << " ";
      Out << "|vv ";
      for(int k = 0; k < 3; ++k)
	Out << vacf[n][t][k]/norm[n][t] << " "; 
      Out << "|";
    }
    Out << endl;
  }
  Out.close();
  
  return true;
}

bool shear_visc(int Ndat, double V, double pres){
  Dynvec<double> time(Ndat+1);
  Dynvec<Matrix<double> > HFdat(Ndat+1); //helfand moments are not correct for periodic bc, extra terms have to be added.
					 // (JChemPhys 126, 184512). The extra term (Int Fij Lij) cannot be computed efficiently in reciprocal space!
					 // even though r is not bounded by the box,  mshf seems to be bounded! why?
  Dynvec<Matrix<double> > Pdat(Ndat+1);
  for(int i = 0; i < Ndat+1; ++i){
    HFdat[i].set_size(3,3);
    HFdat[i].set_to(0.);
    Pdat[i].set_size(3,3);
    Pdat[i].set_to(0.);
  }
  std::ifstream In("ptens.dat");
  if(! In){
    In.close();
    return false;
  }
  for(int n = 0; n < Ndat; ++n){
    In >> time[n];
    for(int i = 0; i < 3; ++i){
      for(int j = 0; j < 3; ++j)
	In >> HFdat[n][i][j];
    }
    for(int i = 0; i < 3; ++i){
      for(int j = 0; j < 3; ++j){
	In >> Pdat[n][i][j];
	if(i == j) Pdat[n][i][j] -= pres;
      }
    }
  }
  In.close();
  
  Dynvec<int> norm(Ndat+1);
  Dynvec<Matrix<double> > mshf(Ndat+1);
  Dynvec<Matrix<double> > pacf(Ndat+1);
  for(int i = 0; i < Ndat+1; ++i){
    mshf[i].set_size(3,3);
    mshf[i].set_to(0.);
    pacf[i].set_size(3,3);
    pacf[i].set_to(0.);
  }
  double d[3][3];
  for(int n = 0; n < Ndat; ++n){
    norm[n] = 0;
    for(int i = 0; i < 3; ++i){
      for(int j = 0; j < 3; ++j){
	mshf[n][i][j] = 0.;
	pacf[n][i][j] = 0.;
      }
    }
  }
  for(int t0 = 0; t0 < Ndat; ++t0){
    for(int t = t0; t < Ndat; ++t){
      for(int i = 0; i < 3; ++i){
	for(int j = 0; j < 3; ++j){
	  d[i][j] = (HFdat[t][i][j] - HFdat[t0][i][j]);
	  mshf[t-t0][i][j] += d[i][j] * d[i][j] / V; //take into account that it is divided by V when computing shear visc.!
	  pacf[t-t0][i][j] += Pdat[t0][i][j] * Pdat[t][i][j];
	}
      }
      norm[t-t0]++;
    }
  }
  
  std::ofstream Out("viscodata.dat");
  if(! Out){
    Out.close();
    return false;
  }
  Out << setiosflags(ios::fixed) << setprecision(6);
  for(int t = 0; t < Ndat; ++t){
    Out << time[t] << " |rp ";
    for(int i = 0; i < 3; ++i){
      for(int j = 0; j < 3; ++j)
	Out << mshf[t][i][j] / norm[t] << " ";
      Out << "| ";
    }
    Out << "|P ";
    for(int i = 0; i < 3; ++i){
      for(int j = 0; j < 3; ++j)
	Out << pacf[t][i][j] / norm[t] << " ";
      Out << "| ";
    }
    Out << endl;
  }
  Out.close();
  return true;
}


int main(){
  int Ncomp;
  deque<int> np_tc;
  int mx_tcf;
  int Nfiles;
  int Ndat;
  double avP;
  double V;
  
  ifstream IN("tc_input.txt");
  if(! IN){
    cout<<"Could not read tc_input.txt\n";
    IN.close();
    return false;
  }
  IN >> Ncomp;
  for(int n = 0; n < Ncomp; ++n){
    int sp;
    IN >> sp;
    np_tc.push_back(sp);
  }
  IN >> mx_tcf;
  IN >> Nfiles;
  IN >> Ndat;
  IN >> avP;
  IN >> V;
  IN.close();
  
//   cout << "Enter number of species: ";
//   cin >> Ncomp;
//   
//   cout << "Enter number of particle per species.\n";
//   for(int n = 0; n < Ncomp; ++n){
//     int sp;
//     cout << "Species " << n << ": ";
//     cin >> sp;
//     np_tc.push_back(sp);
//   }
//   
//   cout << "Enter the maximum number of particles per species taken into account for the diffusion computation: ";
//   cin >> mx_tcf;
//   
//   cout << "Enter the number of files to read for computing diffusion: ";
//   cin >> Nfiles;
//   
//   cout << "Enter the number of rows in ptens.dat: ";
//   cin >> Ndat;
//   
//   cout << "Enter the average pressure in the system: ";
//   cin >> avP;
//   
//   cout << "Finally.. enter the volume of the system: ";
//   cin >> V;
  
  //Num. of particles taken into account for computing time correlations
  for(int i = 0; i < Ncomp; ++i){
    if(np_tc[i] > mx_tcf) np_tc[i] = mx_tcf;
  }
  
  //error
  ofstream Err("ErrorTC.dat");
  Err.close();
  
  if(! shear_visc(Ndat, V, avP)){
    Err.open("ErrorTC.dat", ios::app);
    Err << "problem opening viscodata.dat\n"; 
    Err.close();
    return 1;
  }
  
  Err.open("ErrorTC.dat", ios::app);
  Err << "Computing diffusion graphs...\n";
  Err.close();
  if(! diffusion(Nfiles, np_tc)){
    Err.open("ErrorTC.dat", ios::app);
    Err << "problem opening difudata.dat\n"; 
    Err.close();
    return 1;
  }
  Err.open("ErrorTC.dat", ios::app);
  Err << "Diffusion graphs ready!\n";
  Err.close();
  
  return 0;
}