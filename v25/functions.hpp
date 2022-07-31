#ifndef __FUNCTIONS_HPP
#define __FUNCTIONS_HPP

#include <string>
#include <complex>
#include <deque>
#include "classes.hpp"

//get input parameters
bool r_input(int&, int&, int&, int& , part_num&, Dynvec<double>&, 
	     Dynvec<double>&, Dynvec<double>&, Dynvec<Matrix<double> >&, double&, double&, 
	     double&, double&, int&, double&, double&, double&, double&, int&, 
	     double&, double&, int&, int&, double[3], int&, double&);
// minumum dist. between spherocylinders
double spherocyl_dist(const Spherocyl&, const Spherocyl&, const double[3], const double[3], double[3], double[2]);
//initial configuration
int ini_pos(Dynvec<Dynvec<Spherocyl> >&, const double[3], const part_num&);
void ini_mom(Dynvec<Dynvec<Spherocyl> >&, const part_num&, double, const gsl_rng*);
//gaussian distribution
void gauss(double, double&, double&);
//find species and number from id
int calc_II(const part_num&, int, int);
void rev_II(const part_num&, int&, int&, int);
//print configuration
bool print_conf(const Dynvec<Dynvec<Spherocyl> >&, Thermostat[2], const part_num&, 
		const double[3], bool);
//charge configuration
bool chrg_conf(Dynvec<Dynvec<Spherocyl> >&, Thermostat[2], const part_num&, double[3]);
//neighbor list
void neigh_list(const Dynvec<Dynvec<Spherocyl> >&, const part_num&, 
		const double[3], deque<int>[], double);
bool savec_list(const Dynvec<Dynvec<Spherocyl> >&, const part_num&);
bool check_list(const Dynvec<Dynvec<Spherocyl> >&, const part_num&, double, 
		const double[3], double, bool&);
bool check_list(const Dynvec<Dynvec<Spherocyl> >&, const part_num&, const double[3], 
		double, double, bool&);
//compute forces
bool force_sc(Dynvec<Dynvec<Spherocyl> >&, const part_num&, const double[3], 
	   double, double&, double&, double[3][3], const deque<int>[]);
void force_dipR(Dynvec<Dynvec<Spherocyl> >&, const part_num&, const double[3], 
		const double[32768][3], double&, double[3][3]);
void force_dipK(Dynvec<Dynvec<Spherocyl> >&, const part_num&, const double[3], int, 
		double, const Dynvec<Matrix<double> >&, double&, double[3][3]);
void ext_field(Dynvec<Dynvec<Spherocyl> >&, const part_num&, const double[3], double&);
double self_ener(const Dynvec<Dynvec<Spherocyl> >&, const part_num&);
//tables
void tabul(double[32768][3], double);
void wave_v(Dynvec<Matrix<double> >&, const double[3], double, int);
//velocity Verlet
void moveX(Dynvec<Dynvec<Spherocyl> >&, const part_num&, double);
void moveP(Dynvec<Dynvec<Spherocyl> >&, Thermostat[2], const part_num&, double, 
	   double[2], double[3][3], bool);
//move thermo
void moveTher(Thermostat[2], double, double, const double[2], double, int);
//print/charge counters
bool save_cont(double, long, int, int, int, const accumulator&, const accumulator&);
bool chrg_cont(double&, long&, int&, int&, int&, accumulator&, accumulator&);
//static correlation function
void acc_correl(const Dynvec<Dynvec<Spherocyl> >&, Dynvec<Matrix<double> >&, 
		Dynvec<Matrix<double> >&, const double[3], double, const part_num&, 
		const cross[]);
bool nor_correl(const Dynvec<Matrix<double> >&, const Dynvec<Matrix<double> >&,
		double, double, const part_num&, double);
bool chrg_correl(double& cor_cont, Dynvec<Matrix<double> >&, Dynvec<Matrix<double> >&,
		 const part_num&);
//print values for computing time correlations
bool pos_vel(const Dynvec<Dynvec<Spherocyl> >&, const Dynvec<int>&, double, int);
bool print_ptens(const Dynvec<Dynvec<Spherocyl> >&, const part_num&, double, const double[3][3]);
bool print_magmom(double, const double[3]);
//orientational order
void oo_params(const Dynvec<Dynvec<Spherocyl> >&, const part_num&, double&, double&, double[3]);
// void oo_params(const Dynvec<Dynvec<Spherocyl> >&, const part_num&, double&, double&, double[3],
// 	       const double[3]);
//print values
bool print_out(double, const accumulator&, long, long, int, const double[6]);
bool print_final(double Time, double avK, double flK, double avU, double flU, double avP, 
		 double flP, double, double, double, double, 
		 double, double, double, double, 
		 double, double, double, double, double, double,
		 double[3][3], double[3][3],
		 double[3][3], double[3][3],
		 double[3][3], double[3][3],
		 double[3][3], double[3][3],
		 double[3][3], double[3][3],
		 double, double, double, double,
		 const Dynvec<double>&, const Dynvec<double>&, 
		 const Dynvec<double>&, const Dynvec<Matrix<double> >&, 
		 const double[3], const part_num&, double, double, 
		 double, const double B[3]);

//cluster***************
//choose a random seed particle
void rnd_part(const std::vector<std::vector<int> >&, int&, int&, int&, const gsl_rng*);
//test distance
bool dist_test(const Spherocyl&, const Spherocyl&, const double[3], double);
//create a cluster
int clust_gen(Dynvec<Dynvec<Spherocyl> >&, const part_num&, vector<vector<int> >&, 
	      Dynvec<bool>&, vector<vector<cl_element> >&, const double[3], double, 
	      const gsl_rng*, int);
//check for ghost bonds
void gb_check(const Dynvec<Dynvec<Spherocyl> >&, const part_num&, vector<cl_element>&, 
	      const double[3], double, int);
//check for percolation
bool percolation(const Dynvec<Dynvec<Spherocyl> >&, const part_num&, vector<cl_element>&,
		 const double[3], double, double);
//divide the configuration into clusters
bool clust_conf(Dynvec<Dynvec<Spherocyl> >&, const part_num&, const double[3], double, 
		double, const gsl_rng*, bool);

////////////////////////////
void GCAmix(Dynvec<Dynvec<Spherocyl> >&, const part_num&, const double[3], const gsl_rng*);

void dot(const double[3], const double[3], double&);

void cros(const double[3], const double[3], double[3]);

int delta(int, int);

void Qmat_3d(const Dynvec<Dynvec<Spherocyl> >&, const part_num&, float**);

/////////////////////////

//Matrix functions********************************************************
Matrix<double> transpose(const Matrix<double>& RESTRICT);
Matrix<double> LU_dec(const Matrix<double>& RESTRICT,
		      int* RESTRICT, int&, bool&);
Matrix<double> fwd_bkd(const Matrix<double>& RESTRICT, const double* RESTRICT);
Matrix<double> solve(const Matrix<double>& RESTRICT ,
		     const Matrix<double>& RESTRICT);
Matrix<double> inverse(const Matrix<double>& RESTRICT);

#endif