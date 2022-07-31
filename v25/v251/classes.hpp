#ifndef __CLASSES_HPP
#define __CLASSES_HPP

#include <cmath>
#include <ctime>
#include <sys/time.h>
#include <iomanip>
#include <vector>
#include <gsl/gsl_rng.h>
#include <omp.h>
#include "Matrix.hpp"

//part_num------
class part_num {
 protected:
  int m_Ncomp;
  int m_Ncros;
  int m_Tpart;
  Dynvec<int> m_Npart;
  
 public:
  part_num() : m_Ncomp(1), m_Ncros(0), m_Tpart(1) {}
  ~part_num() {}

  bool init(int, const Dynvec<int>&);
  int get_Ncomp() { return m_Ncomp; }
  int get_Ncomp() const{ return m_Ncomp; }
  int get_Ncros() { return m_Ncros; }
  int get_Ncros() const{ return m_Ncros; }
  int get_Tpart() { return m_Tpart; }
  int get_Tpart() const{ return m_Tpart; }
  int get_Npart(int comp) { return m_Npart[comp]; }
  int get_Npart(int comp) const{ return m_Npart[comp]; } 
  void get_Npart(Dynvec<int>& np);
  void operator=(const part_num&);
};

//Dipole
class Dipole {
 protected:
  double m_u; //dipole moment
  double m_r[3]; //position
  double m_o[3]; //orientation
  
 public:
  Dipole() : m_u(0.) { m_r[0] = m_r[1] = m_r[2] = 0.;
     m_o[0] = m_o[1] = 0.; m_o[2] = 1; }
  ~Dipole() {}
   
  void operator=(const Dipole&);
  double get_dip() { return m_u; }
  double get_dip() const{ return m_u; }
  void set_dip(double u) { m_u = u; }
  void get_pos(double r[3]) { r[0] = m_r[0]; r[1] = m_r[1]; r[2] = m_r[2]; }
  void get_pos(double r[3]) const{ r[0] = m_r[0]; r[1] = m_r[1]; r[2] = m_r[2]; }
  void set_pos(const double r[3]) { m_r[0] = r[0]; m_r[1] = r[1]; m_r[2] = r[2]; }
  void get_ori(double o[3]) { o[0] = m_o[0]; o[1] = m_o[1]; o[2] = m_o[2]; }
  void get_ori(double o[3]) const{ o[0] = m_o[0]; o[1] = m_o[1]; o[2] = m_o[2]; }
  void set_ori(const double o[3]) { m_o[0] = o[0]; m_o[1] = o[1]; m_o[2] = o[2]; }
};

//Spherocyl------
class Spherocyl {
 protected:
  int m_id;
  int m_clust; // cluster to which it belongs
  double m_mass; //mass
  double m_L; //length
  double m_D; //diameter
  double m_I[3]; //principal components of moment of inertia
  double m_r[3]; //position
  double m_p[3]; //momentum
  double m_F[3]; //force acting on this particle
  Matrix<double> m_ori; //orientation (director cosines matrix transposed)
  double m_l[3]; //angular momentum
  double m_T[3]; //torque acting on this particle
  //The dipole is constrained to be constant in the body reference frame, si we keep
  // this constraint in the following parameters
  int m_dips;
  Dynvec<double> m_pos_constr; //is the signed distance of the dip to the center of mass.
  Matrix<double> m_ori_constr; //gives the orientation in terms of the body axes (dip.u, dip.v, dip.w) 
  Dynvec<Dipole> m_u; //dipoles in the lab reference frame
  
 public:
  Spherocyl();
  Spherocyl(int id, double mass, double L, double D, int n_dip, 
	    const Matrix<double>& u);
  Spherocyl(const Spherocyl&);
  ~Spherocyl() {}
  
  void operator=(const Spherocyl&);
  int get_id() { return m_id; }
  int get_id() const{ return m_id; }
  void set_id(int id) { m_id = id; }
  int get_clust() { return m_clust; }
  int get_clust() const{ return m_clust; }
  void set_clust(int clust) { m_clust = clust; }
  double get_mass() { return m_mass; }
  double get_mass() const{ return m_mass; }
  void set_mass(double mass) { m_mass = mass; }  
  void get_shape(double sh[2]) { sh[0] = m_L; sh[1] = m_D; }
  void get_shape(double sh[2]) const{ sh[0] = m_L; sh[1] = m_D; }
  void set_shape(const double sh[2]) { m_L = sh[0]; m_D = sh[1]; }
  int get_ndips() { return m_dips; }
  int get_ndips() const{ return m_dips; }
  void set_ndips(int d) { m_dips = d; }
  Dipole get_dip(int i) { return m_u[i]; }
  Dipole get_dip(int i) const{ return m_u[i]; }
  void set_dip(const Matrix<double>& u);
  void set_inertia();
  void get_inertia(double I[3]) { I[0] = m_I[0]; I[1] = m_I[1]; I[2] = m_I[2]; }
  void get_inertia(double I[3]) const{ I[0] = m_I[0]; I[1] = m_I[1]; I[2] = m_I[2]; }
  void get_pos(double r[3]) { r[0] = m_r[0]; r[1] = m_r[1]; r[2] = m_r[2]; }
  void get_pos(double r[3]) const{ r[0] = m_r[0]; r[1] = m_r[1]; r[2] = m_r[2]; }
  void set_pos(const double r[3]) { m_r[0] = r[0]; m_r[1] = r[1]; m_r[2] = r[2]; }
  void get_mom(double p[3]) { p[0] = m_p[0]; p[1] = m_p[1]; p[2] = m_p[2]; }
  void get_mom(double p[3]) const{ p[0] = m_p[0]; p[1] = m_p[1]; p[2] = m_p[2]; }
  void set_mom(const double p[3]) { m_p[0] = p[0]; m_p[1] = p[1]; m_p[2] = p[2]; }
  void get_F(double f[3]) { f[0] = m_F[0]; f[1] = m_F[1]; f[2] = m_F[2]; }
  void get_F(double f[3]) const{ f[0] = m_F[0]; f[1] = m_F[1]; f[2] = m_F[2]; }
  void set_F(const double f[3]) { m_F[0] = f[0]; m_F[1] = f[1]; m_F[2] = f[2]; }
  void add_F(const double f[3]);
  void get_ori(double o[3]) { o[0] = m_ori[2][0]; o[1] = m_ori[2][1]; o[2] = m_ori[2][2]; }
  void get_ori(double o[3]) const{ o[0] = m_ori[2][0]; o[1] = m_ori[2][1]; o[2] = m_ori[2][2]; }
  void get_omat(Matrix<double>& o) { o = m_ori; }
  void get_omat(Matrix<double>& o) const{ o = m_ori; }
  void set_omat(const Matrix<double>& o) { m_ori = o; }
  void get_angm(double l[3]) { l[0] = m_l[0]; l[1] = m_l[1]; l[2] = m_l[2]; }
  void get_angm(double l[3]) const{ l[0] = m_l[0]; l[1] = m_l[1]; l[2] = m_l[2]; }
  void set_angm(const double l[3]) { m_l[0] = l[0]; m_l[1] = l[1]; m_l[2] = l[2]; }
  void get_T(double t[3]) { t[0] = m_T[0]; t[1] = m_T[1]; t[2] = m_T[2]; }
  void get_T(double t[3]) const{ t[0] = m_T[0]; t[1] = m_T[1]; t[2] = m_T[2]; }
  void set_T(const double t[3]) { m_T[0] = t[0]; m_T[1] = t[1]; m_T[2] = t[2]; }
  void add_T(const double t[3]);
  void get_lev(int i, double d[3]) { d[0] = m_pos_constr[i]*m_ori[2][0]; d[1] = m_pos_constr[i]*m_ori[2][1]; 
			      d[2] = m_pos_constr[i]*m_ori[2][2];}
  void get_lev(int i, double d[3]) const{ d[0] = m_pos_constr[i]*m_ori[2][0]; d[1] = m_pos_constr[i]*m_ori[2][1]; 
				   d[2] = m_pos_constr[i]*m_ori[2][2];}
  void translate(double dt);
  void rotate(double dt);
  void renorm();
};

//Thermostat-----
class Thermostat {
 protected:
  double m_Q; //"mass" (mass * area)
  double m_s; //generalized coord. (1. / time)
  double m_g; //"shear rate" on the thermostat (1. / time)
  gsl_rng* m_ran; //random number generator
  double m_strn; //stored normal random number 
  
 public:
  Thermostat() : m_Q(1.), m_s(1.) , m_g(0.) { 
    m_ran = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(m_ran, 1);
  }
  Thermostat(double Q, double g, long seed) : m_Q(Q), m_s(1.), m_g(g) { 
    m_ran = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(m_ran, seed);
  }
  Thermostat(const Thermostat&);
  ~Thermostat() { gsl_rng_free(m_ran); }
  
  void operator=(const Thermostat&);
  double get_mass() { return m_Q; }
  double get_mass() const{ return m_Q; }
  void set_mass(double Q) { m_Q = Q; }
  double get_coo() { return m_s; }
  double get_coo() const{ return m_s; }
  void set_coo(const double s) { m_s = s; }
  double get_shr() { return m_g; }
  double get_shr() const{ return m_g; }
  void set_shr(const double g) { m_g = g; }
  void set_seed(const long seed) { gsl_rng_set(m_ran, seed); }
  double get_rand() { return gsl_rng_uniform(m_ran); }
  void set_strn(const double r) { m_strn = r; }
  double get_strn() { return m_strn; }
  double get_strn() const{ return m_strn; }
  bool print_ran(unsigned n, bool fin);
  bool read_ran(unsigned n);
};

//accumulator----
class accumulator {
 protected:
  double m_K, m_U, m_P;
  double m_Uss, m_Pss;
  double m_UdipR, m_PdipR;
  double m_UdipK, m_PdipK;
  double m_Uext;
  double m_PT[3][3]; //total pressure tensor
  double m_PTss[3][3]; //pressure tensor short range
  double m_PTdipR[3][3]; //pressure tensor dipole R
  double m_PTdipK[3][3]; //pressure tensor dipole K
  double m_PTkin[3][3]; //pressure tensor kinetic
  double m_G1; //first order parameter
  double m_G2; //maier-saupe order parameter

 public:
  accumulator();
  ~accumulator() {}
  
  void set_K(double K) { m_K = K; }
  double get_K() { return m_K; }
  double get_K() const{ return m_K; }
  void set_U(double U) { m_U = U; }
  double get_U() { return m_U; }
  double get_U() const{ return m_U; }
  void set_P(double P) { m_P = P; }
  double get_P() { return m_P; }
  double get_P() const{ return m_P; }
  void set_Uss(double Uss) { m_Uss = Uss; }
  double get_Uss() { return m_Uss; }
  double get_Uss() const{ return m_Uss; }
  void set_Pss(double Pss) { m_Pss = Pss; }
  double get_Pss() { return m_Pss; }
  double get_Pss() const{ return m_Pss; }
  void set_UdipR(double UdipR) { m_UdipR = UdipR; }
  double get_UdipR() { return m_UdipR; }
  double get_UdipR() const{ return m_UdipR; }
  void set_PdipR(double PdipR) { m_PdipR = PdipR; }
  double get_PdipR() { return m_PdipR; }
  double get_PdipR() const{ return m_PdipR; }
  void set_UdipK(double UdipK) { m_UdipK = UdipK; }
  double get_UdipK() { return m_UdipK; }
  double get_UdipK() const{ return m_UdipK; }
  void set_PdipK(double PdipK) { m_PdipK = PdipK; }
  double get_PdipK() { return m_PdipK; }
  double get_PdipK() const{ return m_PdipK; }
  void set_Uext(double Uext) { m_Uext = Uext; }
  double get_Uext() { return m_Uext; }
  double get_Uext() const{ return m_Uext; }
  void set_PT(const double PT[3][3]) { 
    for(int i = 0; i < 3; ++i){ for(int j = 0; j < 3; ++j) m_PT[i][j] = PT[i][j];} }
  void get_PT(double PT[3][3]) {
    for(int i = 0; i < 3; ++i){ for(int j = 0; j < 3; ++j) PT[i][j] = m_PT[i][j];} }
  void get_PT(double PT[3][3]) const {
    for(int i = 0; i < 3; ++i){ for(int j = 0; j < 3; ++j) PT[i][j] = m_PT[i][j];} }
  void set_PTss(const double PT[3][3]) { 
    for(int i = 0; i < 3; ++i){ for(int j = 0; j < 3; ++j) m_PTss[i][j] = PT[i][j];} }
  void get_PTss(double PT[3][3]) {
    for(int i = 0; i < 3; ++i){ for(int j = 0; j < 3; ++j) PT[i][j] = m_PTss[i][j];} }
  void get_PTss(double PT[3][3]) const {
    for(int i = 0; i < 3; ++i){ for(int j = 0; j < 3; ++j) PT[i][j] = m_PTss[i][j];} }
  void set_PTdipR(const double PT[3][3]) { 
    for(int i = 0; i < 3; ++i){ for(int j = 0; j < 3; ++j) m_PTdipR[i][j] = PT[i][j];} }
  void get_PTdipR(double PT[3][3]) {
    for(int i = 0; i < 3; ++i){ for(int j = 0; j < 3; ++j) PT[i][j] = m_PTdipR[i][j];} }
  void get_PTdipR(double PT[3][3]) const {
    for(int i = 0; i < 3; ++i){ for(int j = 0; j < 3; ++j) PT[i][j] = m_PTdipR[i][j];} }
  void set_PTdipK(const double PT[3][3]) { 
    for(int i = 0; i < 3; ++i){ for(int j = 0; j < 3; ++j) m_PTdipK[i][j] = PT[i][j];} }
  void get_PTdipK(double PT[3][3]) {
    for(int i = 0; i < 3; ++i){ for(int j = 0; j < 3; ++j) PT[i][j] = m_PTdipK[i][j];} }
  void get_PTdipK(double PT[3][3]) const {
    for(int i = 0; i < 3; ++i){ for(int j = 0; j < 3; ++j) PT[i][j] = m_PTdipK[i][j];} }
  void set_PTkin(const double PT[3][3]) { 
    for(int i = 0; i < 3; ++i){ for(int j = 0; j < 3; ++j) m_PTkin[i][j] = PT[i][j];} }
  void get_PTkin(double PT[3][3]) {
    for(int i = 0; i < 3; ++i){ for(int j = 0; j < 3; ++j) PT[i][j] = m_PTkin[i][j];} }
  void get_PTkin(double PT[3][3]) const{
    for(int i = 0; i < 3; ++i){ for(int j = 0; j < 3; ++j) PT[i][j] = m_PTkin[i][j];} }
  void set_G1(double G1) { m_G1 = G1; }
  double get_G1() { return m_G1; }
  double get_G1() const{ return m_G1; }
  void set_G2(double G2) { m_G2 = G2; }
  double get_G2() { return m_G2; }
  double get_G2() const{ return m_G2; }
  void set_zero();
  void operator=(const accumulator&);
  void operator+=(const accumulator&);
  const accumulator operator*(const accumulator&);
  const accumulator operator/(long);
};

//for percolation checking-----
struct cl_element {
  int type; //species
  int num; //number within the species
  int gen; //round (generation) at which this particle was added 
//   bool end; //is it an end point of the cluster?
  bool state; //on/off depending if the particle has been visited or not
  double dis[3]; //net displacement from the initial node
//   int convty; //connectivity of the particle (node)
  std::vector<int> link; //to whom does it link within memb
};

class Walker {
 protected:
  double m_Dr[3]; //displacement relative to first node (0)
  int m_actual; //current node 
  int m_next; //next node to move to
   
 public:
  Walker() : m_actual(0), m_next(0) { m_Dr[0] = 0; m_Dr[1] = 0; m_Dr[2] = 0; }
  Walker(const Walker&);
  ~Walker() {}
  
  void operator=(const Walker&);
  void get_Dr(double D[3]) { D[0] = m_Dr[0]; D[1] = m_Dr[1]; D[2] = m_Dr[2]; }
  void get_Dr(double D[3]) const{ D[0] = m_Dr[0]; D[1] = m_Dr[1]; D[2] = m_Dr[2]; }
  void set_Dr(const double D[3]) { m_Dr[0] = D[0]; m_Dr[1] = D[1]; m_Dr[2] = D[2]; }
  int get_actual() { return m_actual; }
  int get_actual() const{ return m_actual; }
  void set_actual(int actual) { m_actual = actual; }
  int get_next() { return m_next; }
  int get_next() const{ return m_next; }
  void set_next(int next) { m_next = next; }
  double get_Drsq() { return m_Dr[0]*m_Dr[0] + m_Dr[1]*m_Dr[1] + m_Dr[2]*m_Dr[2]; }
  void move(const Dynvec<Dynvec<Spherocyl> >&, const double[3], const std::vector<cl_element>&);
};


//Timer--------
class Timer {
 private:
  timeval m_tim;
  long m_tcpu;
  
 public:
  Timer();
  ~Timer() {}
  
   void start();
   void stop(double&, long&);
};

//inter species
struct cross {
  Dynvec<int> C;
  Dynvec<int> P;
  Dynvec<bool> T;
};


#endif