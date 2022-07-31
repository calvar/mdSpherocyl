#include "classes.hpp"
#include "functions.hpp"

//part_num-----
bool part_num::init(int Ncomp, const Dynvec<int>& N) {
  m_Ncomp = Ncomp;
  if(N.get_size() != m_Ncomp ) return false;
  m_Ncros = m_Ncomp * (m_Ncomp - 1) / 2;
  m_Npart.set_size(m_Ncomp);
  m_Tpart = 0;
  for(int i = 0; i < m_Ncomp; ++i){
    m_Npart[i] = N[i];
    m_Tpart += N[i];
  }
  return true; 
}

void part_num::get_Npart(Dynvec<int>& np) {
  np.set_size(m_Ncomp);
  np = m_Npart;
}

void part_num::operator=(const part_num& NP) {
  m_Ncomp = NP.m_Ncomp;
  m_Ncros = NP.m_Ncros;
  m_Tpart = NP.m_Tpart;
  m_Npart.set_size(m_Ncomp);
  for(int i = 0; i < m_Ncomp; ++i){
    m_Npart[i] = NP.m_Npart[i];
  }
}

//Dipole----------
void Dipole::operator=(const Dipole& d) { 
  m_u = d.m_u;
  for(int i = 0; i < 3; ++i){
    m_r[i] = d.m_r[i]; 
    m_o[i] = d.m_o[i];
  } 
}

//Spherocyl-------
Spherocyl::Spherocyl() : m_id(0), m_clust(0), m_mass(1.), m_L(0.), m_D(0.5) {
  for(int i = 0; i < 3; ++i){
    m_r[i] = 0.;
    m_p[i] = 0.;
    m_F[i] = 0.;
    m_l[i] = 0.;
    m_T[i] = 0.;
  }
  double R = m_D / 2;
  double Rsq = R * R;
  double Rft = Rsq * Rsq;
  double Lsq = m_L * m_L;
  double vol = M_PI * Rsq * (m_L + 3*R/4); 
  double mdens = m_mass / vol;
  m_I[0] = mdens * M_PI * (Lsq*m_L*Rsq/12 + Lsq*Rsq*R/3 + m_L*Rft/4 + 8*Rft*R/15);
  m_I[1] = m_I[0];
  m_I[2] = mdens * M_PI * (m_L*Rft/2 + 8*Rft*R/15);
  m_ori.set_size(3, 3); //transpose of director cosines matrix
  m_ori[0][0] = 1.; m_ori[0][1] = 0.; m_ori[0][2] = 0.;
  m_ori[1][0] = 0.; m_ori[1][1] = 1.; m_ori[1][2] = 0.; 
  m_ori[2][0] = 0.; m_ori[2][1] = 0.; m_ori[2][2] = 1.;
  m_pos_constr.set_size(1);
  m_pos_constr[0] = 0.;
  m_ori_constr.set_size(1, 3);
  m_ori_constr[0][0] = 0.; m_ori_constr[0][1] = 0.; m_ori_constr[0][2] = 1.;
  m_u.set_size(1);
  double v[3] = {0., 0., m_pos_constr[0]};
  m_u[0].set_pos(v);
  v[0] = m_ori_constr[0][0]; v[1] = m_ori_constr[0][1]; v[2] = m_ori_constr[0][2];
  m_u[0].set_ori(v);
  m_u[0].set_dip(0.);
}

Spherocyl::Spherocyl(int id, double mass, double L, double D, int n_dip, 
		     const Matrix<double>& u) : m_id(id), m_clust(0), m_mass(mass), m_L(L), m_D(D),
		     m_dips(n_dip){
  for(int i = 0; i < 3; ++i){
    m_r[i] = 0.;
    m_p[i] = 0.;
    m_F[i] = 0.;
    m_l[i] = 0.;
    m_T[i] = 0.;
  }
  double R = m_D / 2;
  double Rsq = R * R;
  double Rft = Rsq * Rsq;
  double Lsq = m_L * m_L;
  double vol = M_PI * Rsq * (m_L + 3.*R/4); 
  double mdens = m_mass / vol;
  m_I[0] = mdens * M_PI * (Lsq*m_L*Rsq/12 + Lsq*Rsq*R/3 + m_L*Rft/4 + 8.*Rft*R/15);
  m_I[1] = m_I[0];
  m_I[2] = mdens * M_PI * (m_L*Rft/2 + 8.*Rft*R/15);
  m_ori.set_size(3, 3); //transpose of director cosines matrix
  m_ori[0][0] = 1.; m_ori[0][1] = 0.; m_ori[0][2] = 0.;
  m_ori[1][0] = 0.; m_ori[1][1] = 1.; m_ori[1][2] = 0.; 
  m_ori[2][0] = 0.; m_ori[2][1] = 0.; m_ori[2][2] = 1.;
  
  int a, b;
  u.get_size(a, b);
  if(a != m_dips || b != 5){
    ofstream Err("Error.dat");
    Err << "Trying to initialize dipoles with incorrect array size: "<<a<<"x"<<b<<endl;
    Err.close();
  }
  m_pos_constr.set_size(m_dips);
  m_ori_constr.set_size(m_dips, 3);
  double v[3];
  for(int i = 0; i < m_dips; ++i){
    m_pos_constr[i] = u[i][0];
    for(int j = 0; j < 3; ++j) 
      m_ori_constr[i][j] = u[i][j+1];
    
    v[0] = 0.; v[1] = 0.; v[2] = m_pos_constr[i];
    m_u[i].set_pos(v);
    v[0] = m_ori_constr[0][0]; v[1] = m_ori_constr[0][1]; v[2] = m_ori_constr[0][2];
    m_u[i].set_ori(v);
    m_u[i].set_dip(u[i][4]);
  }
}

Spherocyl::Spherocyl(const Spherocyl& P) {
  m_id = P.m_id;
  m_clust = P.m_clust;
  m_mass = P.m_mass;
  m_L = P.m_L;
  m_D = P.m_D;
  for(int i = 0; i < 3; ++i){
    m_I[i] = P.m_I[i];
    m_r[i] = P.m_r[i];
    m_p[i] = P.m_p[i];
    m_F[i] = P.m_F[i];
    m_l[i] = P.m_l[i];
    m_T[i] = P.m_T[i];
  }
  m_ori.set_size(3, 3);
  m_ori = P.m_ori;
  m_dips = P.m_dips;
  m_pos_constr.set_size(m_dips);
  m_pos_constr = P.m_pos_constr;
  m_ori_constr.set_size(m_dips, 3);
  m_ori_constr = P.m_ori_constr;
  m_u.set_size(m_dips);
  m_u = P.m_u;
}

void Spherocyl::set_dip(const Matrix<double>& u) { 
  int b;
  u.get_size(m_dips, b);
  
  m_pos_constr.set_size(m_dips);
  m_ori_constr.set_size(m_dips, 3);
  m_u.set_size(m_dips);
  double v[3];
  Matrix<double> l0(3,1), l1(3,1);
  for(int i = 0; i < m_dips; ++i){
    m_pos_constr[i] = u[i][0];
    for(int j = 0; j < 3; ++j) 
      m_ori_constr[i][j] = u[i][j+1];
    
    v[0] = m_r[0] + m_pos_constr[i]*m_ori[2][0]; 
    v[1] = m_r[1] + m_pos_constr[i]*m_ori[2][1]; 
    v[2] = m_r[2] + m_pos_constr[i]*m_ori[2][2];
    m_u[i].set_pos(v);
    for(int j = 0; j < 3; ++j) 
      l0[j][0] = m_ori_constr[i][j];
    l1 = transpose(m_ori) * l0;
    v[0] = l1[0][0]; v[1] = l1[1][0]; v[2] = l1[2][0];
    m_u[i].set_ori(v);
    m_u[i].set_dip(u[i][4]);
  }
}

void Spherocyl::set_inertia() {
  double R = m_D / 2;
  double Rsq = R * R;
  double Rft = Rsq * Rsq;
  double Lsq = m_L * m_L;
  double vol = M_PI * Rsq * (m_L + 3.*R/4); 
  double mdens = m_mass / vol;
  m_I[0] = mdens * M_PI * (Lsq*m_L*Rsq/12 + Lsq*Rsq*R/3 + m_L*Rft/4 + 8.*Rft*R/15);
  m_I[1] = m_I[0];
  m_I[2] = mdens * M_PI * (m_L*Rft/2 + 8.*Rft*R/15);
}

void Spherocyl::operator=(const Spherocyl& P) {
  m_id = P.m_id;
  m_clust = P.m_clust;
  m_mass = P.m_mass;
  m_L = P.m_L;
  m_D = P.m_D;
  for(int i = 0; i < 3; ++i){
    m_I[i] = P.m_I[i];
    m_r[i] = P.m_r[i];
    m_p[i] = P.m_p[i];
    m_F[i] = P.m_F[i];
    m_l[i] = P.m_l[i];
    m_T[i] = P.m_T[i];
  }
  m_ori = P.m_ori;
  m_pos_constr.set_size(m_dips);
  m_pos_constr = P.m_pos_constr;
  m_ori_constr.set_size(m_dips, 3);
  m_ori_constr = P.m_ori_constr;
  m_u.set_size(m_dips);
  m_u = P.m_u;
}

void Spherocyl::add_F(const double f[3]) {
  #pragma omp atomic
  m_F[0] += f[0]; 
  #pragma omp atomic
  m_F[1] += f[1]; 
  #pragma omp atomic
  m_F[2] += f[2];
}

void Spherocyl::add_T(const double t[3]) {
  #pragma omp atomic
  m_T[0] += t[0]; 
  #pragma omp atomic
  m_T[1] += t[1]; 
  #pragma omp atomic
  m_T[2] += t[2];
}

void Spherocyl::translate(double dt) { 
  for(int j = 0; j < 3; j++)
    m_r[j] += dt * m_p[j] / m_mass;
  //update dipole positions
  double v[3];
  for(int i = 0; i < m_dips; ++i){
    v[0] = m_r[0] + m_pos_constr[i]*m_ori[2][0];
    v[1] = m_r[1] + m_pos_constr[i]*m_ori[2][1];
    v[2] = m_r[2] + m_pos_constr[i]*m_ori[2][2];
    m_u[i].set_pos(v);
  }
}

void Spherocyl::rotate(double dt) { //(J.Chem.Phys 107, 5840)
  Matrix<double> R(3,3), Q(3,3), l0(3,1), l1(3,1);
  double ang, cosa, sina, dt2 = dt / 2;
  //take l to the body frame
  for(int i = 0; i < 3; ++i) l0[i][0] = m_l[i];
  l1 = m_ori * l0;
  //first rotation Rx
  ang = -l1[0][0] * dt2 / m_I[0];
  cosa = cos(ang); sina = sin(ang);
  R[0][0] = 1.;   R[0][1] = 0.;   R[0][2] = 0.;
  R[1][0] = 0.;   R[1][1] = cosa; R[1][2] = -sina;
  R[2][0] = 0.;   R[2][1] = sina; R[2][2] = cosa;
  l0 = R * l1;
  Q = R * m_ori;
  //second rotation Ry
  ang = -l0[1][0] * dt2 / m_I[1];
  cosa = cos(ang); sina = sin(ang);
  R[0][0] = cosa; R[0][1] = 0.;   R[0][2] = sina;
  R[1][0] = 0.;   R[1][1] = 1.;   R[1][2] = 0.;
  R[2][0] = -sina;R[2][1] = 0.;   R[2][2] = cosa;
  l1 = R * l0;
  m_ori = R * Q;
  //third rotation Rz
  ang = -l1[2][0] * dt / m_I[2];
  cosa = cos(ang); sina = sin(ang);
  R[0][0] = cosa; R[0][1] = -sina;R[0][2] = 0.;
  R[1][0] = sina; R[1][1] = cosa; R[1][2] = 0.;
  R[2][0] = 0.;   R[2][1] = 0.;   R[2][2] = 1.;
  l0 = R * l1;
  Q = R * m_ori;
  //forth rotation Ry
  ang = -l0[1][0] * dt2 / m_I[1];
  cosa = cos(ang); sina = sin(ang);
  R[0][0] = cosa; R[0][1] = 0.;   R[0][2] = sina;
  R[1][0] = 0.;   R[1][1] = 1.;   R[1][2] = 0.;
  R[2][0] = -sina;R[2][1] = 0.;   R[2][2] = cosa;
  l1 = R * l0;
  m_ori = R * Q;
  //fifth rotation Rx
  ang = -l1[0][0] * dt2 / m_I[0];
  cosa = cos(ang); sina = sin(ang);
  R[0][0] = 1.;   R[0][1] = 0.;   R[0][2] = 0.;
  R[1][0] = 0.;   R[1][1] = cosa; R[1][2] = -sina;
  R[2][0] = 0.;   R[2][1] = sina; R[2][2] = cosa;
  l0 = R * l1;
  Q = R * m_ori;
  
  m_ori = Q;
  
  //update dipole orientations and positions
  double v[3];
  for(int i = 0; i < m_dips; ++i){
    v[0] = m_r[0] + m_pos_constr[i]*m_ori[2][0];
    v[1] = m_r[1] + m_pos_constr[i]*m_ori[2][1];
    v[2] = m_r[2] + m_pos_constr[i]*m_ori[2][2];
    m_u[i].set_pos(v);
    
    for(int j = 0; j < 3; ++j) 
      l0[j][0] = m_ori_constr[i][j];
    l1 = transpose(m_ori) * l0;
    v[0] = l1[0][0]; v[1] = l1[1][0]; v[2] = l1[2][0];
    m_u[i].set_ori(v);
  }
}

//reortho-normalize the basis every n steps to correct numerical errors(effect on the dinamics?)
void Spherocyl::renorm() {
  //normalize u3
  double norm = sqrt(m_ori[2][0]*m_ori[2][0]+m_ori[2][1]*m_ori[2][1]+m_ori[2][2]*m_ori[2][2]);
  m_ori[2][0] /= norm; m_ori[2][1] /= norm; m_ori[2][2] /= norm;
  //find a vector in the u3xu2 plane between u3 and u2
  double Up[3] = {m_ori[2][0]+m_ori[1][0], m_ori[2][1]+m_ori[1][1], m_ori[2][2]+m_ori[1][2]};
  //find u2 (non-normalized)
  double dot = Up[0] * m_ori[2][0] + Up[1] * m_ori[2][1] + Up[2] * m_ori[2][2];
  m_ori[1][0] = Up[0] - dot * m_ori[2][0];
  m_ori[1][1] = Up[1] - dot * m_ori[2][1];
  m_ori[1][2] = Up[2] - dot * m_ori[2][2];
  //normalize u2
  norm = sqrt(m_ori[1][0]*m_ori[1][0]+m_ori[1][1]*m_ori[1][1]+m_ori[1][2]*m_ori[1][2]);
  m_ori[1][0] /= norm; m_ori[1][1] /= norm; m_ori[1][2] /= norm;
  //find u1 = u2 x u3
  m_ori[0][0] = m_ori[1][1]*m_ori[2][2] - m_ori[2][1]*m_ori[1][2];
  m_ori[0][1] = m_ori[1][2]*m_ori[2][0] - m_ori[2][2]*m_ori[1][0];
  m_ori[0][2] = m_ori[1][0]*m_ori[2][1] - m_ori[2][0]*m_ori[1][1];
}

//Thermostat-----
Thermostat::Thermostat(const Thermostat& T) {
  m_Q = T.m_Q;
  m_s = T.m_s;
  m_g = T.m_g;
  m_ran = T.m_ran;
  m_strn = T.m_strn;
}

void Thermostat::operator=(const Thermostat& T) {
  m_Q = T.m_Q;
  m_s = T.m_s;
  m_g = T.m_g;
  m_ran = T.m_ran;
  m_strn = T.m_strn;
}

bool Thermostat::print_ran(unsigned n, bool fin) {
  char name[30];
  char c;
  if(! fin) c = 'I';
  else c = 'F';
  sprintf(name, "ran_state%d%c.dat", n, c);
  FILE* stream = fopen(name, "w");
  if(stream == NULL)
    return false;
  gsl_rng_fwrite(stream, m_ran);
  fclose(stream);
  return true;
}

bool Thermostat::read_ran(unsigned n) {
  char name[30];
  sprintf(name, "ran_state%dF.dat", n);
  FILE* stream = fopen(name, "r");
  if(stream == NULL)
    return false;
  gsl_rng_fread(stream, m_ran);
  fclose(stream);
  return true;
}

//accumulator----
accumulator::accumulator() : m_K(0.), m_U(0.), m_P(0.), m_Uss(0.), m_Pss(0.),
		  m_UdipR(0.), m_PdipR(0.), m_UdipK(0.), m_PdipK(0.), m_Uext(0.){
  for(int i = 0; i < 3; ++i){
    for(int j = 0; j < 3; ++j){
      m_PT[i][j] = 0.;
      m_PTss[i][j] = 0.;
      m_PTdipR[i][j] = 0.;
      m_PTdipK[i][j] = 0.;
      m_PTkin[i][j] = 0.;
    }
  }
  m_G1 = 0.;
  m_G2 = 0.;
}

void accumulator::set_zero() {
  m_K = 0.; m_U = 0.; m_P = 0.;
  m_Uss = 0.; m_Pss = 0.; m_UdipR = 0.; m_PdipR = 0.; 
  m_UdipK = 0.; m_PdipK = 0.; m_Uext = 0.;
  for(int i = 0; i < 3; ++i){
    for(int j = 0; j < 3; ++j){
      m_PT[i][j] = 0.;
      m_PTss[i][j] = 0.;
      m_PTdipR[i][j] = 0.;
      m_PTdipK[i][j] = 0.;
      m_PTkin[i][j] = 0.;
    }
  }
  m_G1 = 0.;
  m_G2 = 0.;
}

void accumulator::operator=(const accumulator& acc) {
  m_K = acc.m_K;
  m_U = acc.m_U;
  m_P = acc.m_P;
  m_Uss = acc.m_Uss;
  m_Pss = acc.m_Pss;
  m_UdipR = acc.m_UdipR;
  m_PdipR = acc.m_PdipR;
  m_UdipK = acc.m_UdipK;
  m_PdipK = acc.m_PdipK;
  m_Uext = acc.m_Uext;
  for(int i = 0; i < 3; ++i){
    for(int j = 0; j < 3; ++j){
      m_PT[i][j] = acc.m_PT[i][j];
      m_PTss[i][j] = acc.m_PTss[i][j];
      m_PTdipR[i][j] = acc.m_PTdipR[i][j];
      m_PTdipK[i][j] = acc.m_PTdipK[i][j];
      m_PTkin[i][j] = acc.m_PTkin[i][j];
    }
  }
  m_G1 = acc.m_G1;
  m_G2 = acc.m_G2;
}

void accumulator::operator+=(const accumulator& acc) {
  m_K += acc.m_K;
  m_U += acc.m_U;
  m_P += acc.m_P;
  m_Uss += acc.m_Uss;
  m_Pss += acc.m_Pss;
  m_UdipR += acc.m_UdipR;
  m_PdipR += acc.m_PdipR;
  m_UdipK += acc.m_UdipK;
  m_PdipK += acc.m_PdipK;
  m_Uext = acc.m_Uext;
  for(int i = 0; i < 3; ++i){
    for(int j = 0; j < 3; ++j){
      m_PT[i][j] += acc.m_PT[i][j];
      m_PTss[i][j] += acc.m_PTss[i][j];
      m_PTdipR[i][j] += acc.m_PTdipR[i][j];
      m_PTdipK[i][j] += acc.m_PTdipK[i][j];
      m_PTkin[i][j] += acc.m_PTkin[i][j];
    }
  }
  m_G1 += acc.m_G1;
  m_G2 += acc.m_G2;
}
  
const accumulator accumulator::operator*(const accumulator& acc) {
  accumulator result;
  result.m_K = m_K * acc.m_K;
  result.m_U = m_U * acc.m_U;
  result.m_P = m_P * acc.m_P;
  result.m_Uss = m_Uss * acc.m_Uss;
  result.m_Pss = m_Pss * acc.m_Pss;
  result.m_UdipR = m_UdipR * acc.m_UdipR;
  result.m_PdipR = m_PdipR * acc.m_PdipR;
  result.m_UdipK = m_UdipK * acc.m_UdipK;
  result.m_PdipK = m_PdipK * acc.m_PdipK;
  result.m_Uext = m_Uext * acc.m_Uext;
  for(int i = 0; i < 3; ++i){
    for(int j = 0; j < 3; ++j){
      result.m_PT[i][j] = m_PT[i][j] * acc.m_PT[i][j];
      result.m_PTss[i][j] = m_PTss[i][j] * acc.m_PTss[i][j];
      result.m_PTdipR[i][j] = m_PTdipR[i][j] * acc.m_PTdipR[i][j];
      result.m_PTdipK[i][j] = m_PTdipK[i][j] * acc.m_PTdipK[i][j];
      result.m_PTkin[i][j] = m_PTkin[i][j] * acc.m_PTkin[i][j];
    }
  }
  result.m_G1 = m_G1 * acc.m_G1;
  result.m_G2 = m_G2 * acc.m_G2;
  return result;
}

const accumulator accumulator::operator/(long a) {
  accumulator result;
  result.m_K = m_K / a;
  result.m_U = m_U / a;
  result.m_P = m_P / a;
  result.m_Uss = m_Uss / a;
  result.m_Pss = m_Pss / a;
  result.m_UdipR = m_UdipR / a;
  result.m_PdipR = m_PdipR / a;
  result.m_UdipK = m_UdipK / a;
  result.m_PdipK = m_PdipK / a;
  result.m_Uext = m_Uext / a;
  for(int i = 0; i < 3; ++i){
    for(int j = 0; j < 3; ++j){
      result.m_PT[i][j] = m_PT[i][j] / a;
      result.m_PTss[i][j] = m_PTss[i][j] / a;
      result.m_PTdipR[i][j] = m_PTdipR[i][j] / a;
      result.m_PTdipK[i][j] = m_PTdipK[i][j] / a;
      result.m_PTkin[i][j] = m_PTkin[i][j] / a;
    }
  }
  result.m_G1 = m_G1 / a;
  result.m_G2 = m_G2 / a;
  return result;
}

// Walker----------------
Walker::Walker(const Walker& W) {
  for(int i = 0; i < 3; ++i)
    m_Dr[i] = W.m_Dr[i]; 
  m_actual = W.m_actual;
  m_next = W.m_next;
}

void Walker::operator=(const Walker& W) {
  for(int i = 0; i < 3; ++i)
    m_Dr[i] = W.m_Dr[i]; 
  m_actual = W.m_actual;
  m_next = W.m_next;
}

void Walker::move(const Dynvec<Dynvec<Spherocyl> >& part, const double L[3],
		  const vector<cl_element>& memb) {
  int m = memb[m_actual].type;
  int i = memb[m_actual].num;
//   int n = memb.at(m_next).type;
//   int j = memb.at(m_next).num;
  int n = memb[m_next].type;
  int j = memb[m_next].num;
  double rA[3], rB[3];
  
  part[m][i].get_pos(rA);
  part[n][j].get_pos(rB);
    
  //PBC applied so that the displacement doesnt change sign because of PBC
  double rij[3]/*, rsq=0*/;
  for(int a = 0; a < 3; ++a){
    rij[a] = rA[a] - rB[a];
    rij[a] -= L[a] * floor(rij[a] / L[a] + 0.5);
    m_Dr[a] += rij[a];
  }
  
  m_actual = m_next;
}

//Timer--------
Timer::Timer() {
  gettimeofday(&m_tim, NULL); 
  m_tcpu = clock();
}

void Timer::start() {
  gettimeofday(&m_tim, NULL); 
  m_tcpu = clock();
}

void Timer::stop(double& telap, long& twork) { 
  double t1 = m_tim.tv_sec * 1e+6 + m_tim.tv_usec;
  long start = m_tcpu;
  gettimeofday(&m_tim, NULL);
  double t2 = m_tim.tv_sec * 1e+6 + m_tim.tv_usec;
  m_tcpu = clock();
  long end = m_tcpu;
  telap = t2 - t1;
  twork = end - start;
}