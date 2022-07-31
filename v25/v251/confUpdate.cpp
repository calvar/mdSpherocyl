#include "functions.hpp"

bool print_conf_(const Dynvec<Dynvec<Spherocyl> >& part, Thermostat thermo[2], const deque<int>& NP, 
		const double L[3], bool fin){
  int NC = NP.size();
  char name[30];
  sprintf(name, "conf_.dat");
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
    int NN = NP[n];
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

bool chrg_conf_(Dynvec<Dynvec<Spherocyl> >& part, Thermostat thermo[2], const deque<int>& NP, 
		double L[3]){
  std::string line;
  int NC = NP.size();
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
    int NN = NP[n];
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
//       Con >> cl;
//       part[n][i].set_clust(cl);
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


int main(){
  int Ncomp = 1;
  deque<int> NP(Ncomp);
  NP[0] = 1500;
  double length[1] = {0.}; 
  double diameter[1] = {1.};
  
  //initialize particles
  Dynvec<Dynvec<Spherocyl> > part(Ncomp);
  for(int n = 0; n < Ncomp; ++n)
    part[n].set_size(NP[n]);
 
  int Nsum = 0;
  for(int n = 0; n < Ncomp; ++n){
    double sh[2] = {length[n], diameter[n]};
    int NN = NP[n];
    for(int i = 0; i < NN; ++i){
      part[n][i].set_id(Nsum+i);
      part[n][i].set_shape(sh);
      part[n][i].set_clust(0);
    }
    Nsum += NN;
  }
  
  //initialize thermostats
  Thermostat thermo[2]; //0. is for translation and 1. is for rotation
//   thermo[0].set_mass(Q);
//   thermo[0].set_shr(gamma);
//   thermo[0].set_seed(1/*time(NULL)*/);
//   thermo[1].set_mass(Q);
//   thermo[1].set_shr(gamma);
//   thermo[1].set_seed(2/*time(NULL)*/);
  
  double L[3] = {42.8239,42.8239,42.8239};
  
  //charge configuration
  if(! chrg_conf_(part, thermo, NP, L)){
    cout << "problem opening conf.dat for reading\n";
    return 1;
  }
  
  //write configuration//
  if(! print_conf_(part, thermo, NP, L, true)){
    cout << "problem opening conf_.dat\n"; 
    return 1;
  }
  
  return 0;
}