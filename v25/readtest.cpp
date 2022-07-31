#include "functions.hpp"

int main(){
  int Nsteps, pre_steps;
  int iprint;
  int signal;
  part_num NP;
  Dynvec<double> length; //length (L) of the spherocylinders
  Dynvec<double> diameter; //diameter (D) of the spherocylinders
  Dynvec<double> mass; //mass of the spherocylinders
  //new
  Dynvec<double> dip[5]; //particle dipoles
  double B[3]; //external field (B*=B*sqrt(sigma^3/(epsilon*relative_mag_permeability)))
  double dens; //density
  double temp; //initial temperature
  double sr_cut; //short range cutoff (long range uses minimum image convention)
  int K_max; //max K vectors
  double alpha; //convergence parameter
  double Dt; //time interval
  double bar_w; //histogram bar width
  double ArLx, ArLy; //aspect ratio of the box (Lx/Lz and Ly/Lz)
  int corrf; // correlation sum frequency
  int blk_size = 1000;
  //epsilon = 1 (energies normalized by epsilon)
  
  r_input(Nsteps, pre_steps, iprint, signal, NP, length, diameter, mass, dip, dens, temp, Dt, 
	  sr_cut, K_max, alpha, bar_w, ArLx, ArLy, corrf, B);
  
  int NC = NP.get_Ncomp();
  int NS = NP.get_Ncros();
  int Tpart = NP.get_Tpart();

  //error
  ofstream Err("Error.dat");
  Err.close();
  
  //check number of particles
  for(int n = 0; n < NC; ++n){
    if(NP.get_Npart(n) < 1){
      Err.open("Error.dat", ios::app);
      Err << "No. of particles per species must be >= 1" << endl;
      Err.close();
      return 1;
    }
  }
  
  //new
  //check taht dipolea are inside the particle
  for(int n = 0; n < NC; ++n){
    if(fabs(dip[0][n]) >= (length[n]+diameter[n])/2-0.01){
      Err.open("Error.dat", ios::app);
      Err << "dipole of species " << n << " is ouside the core!" << endl;
      Err.close();
      return 1;
    }
  }
  
  
   //find largest particle dimension
  double dmax = 0, dmin = 1.e+6;
  for(int n = 0; n < NC; ++n){
    dmax = max(dmax, length[n]+diameter[n]);
    dmin = min(dmin, diameter[n]);
  }
  
  //new
  //initialize particles
  Dynvec<Dynvec<Spherocyl> > part(NC);
  for(int n = 0; n < NC; ++n)
    part[n].set_size(NP.get_Npart(n));
 
  int Nsum = 0;
  for(int n = 0; n < NC; ++n){
    double sh[2] = {length[n], diameter[n]};
    double dp[5];
    for(int a = 0; a < 5; ++a) dp[a] = dip[a][n];
    int NN = NP.get_Npart(n);
    for(int i = 0; i < NN; ++i){
      part[n][i].set_id(Nsum+i);
      part[n][i].set_shape(sh);
      part[n][i].set_mass(mass[n]);
      part[n][i].set_inertia();
      part[n][i].set_dip(dp);  
    }
    Nsum += NN;
  }
  
   //box dimmensions
  double L[3], Lmin=0;
  
  for(int step = 0; step < pre_steps+Nsteps; ++step){
   
    //charge configuration
    if(! chrg_conf(part, NP, L)){
      Err.open("Error.dat", ios::app);
      Err << "problem opening conf.dat for reading\n";
      Err.close();
      return 1;
    }
    //new
    for(int n = 0; n < NC; ++n){
      double dp[5];
      for(int a = 0; a < 5; ++a) dp[a] = dip[a][n];
      int NN = NP.get_Npart(n);
      for(int i = 0; i < NN; ++i){
	part[n][i].set_mass(mass[n]);
	part[n][i].set_inertia();
	part[n][i].set_dip(dp);
      }
    }
    
    //write configuration//
    if(! print_conf(part, NP, L)){
      Err.open("Error.dat", ios::app);
      Err << "problem opening conf.dat\n"; 
      Err.close();
      return 1;
    }
    
  }
  
  return 0;
}