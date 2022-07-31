#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>
#include <complex>
using namespace std;

#include "GL/gl.h"
#include "GL/glu.h"
#include "SDL/SDL.h"
#include "Matrix.hpp"
#include "gl2ps.h"

//Screen attributes
const int SCREEN_WIDTH = 480;
const int SCREEN_HEIGHT = 480;
const int SCREEN_BPP = 32;

//The frame rate
const int FRAMES_PER_SECOND = 20;

//Event handler
SDL_Event event;

double user_theta  = 0.1;
double user_phi = 0.;
double user_dist = 15.;
double usr_x = 0., usr_y = 0.01*user_dist, usr_z = 0.99*user_dist; //pos. from where im looking
double ref_x = 0, ref_y = 0, ref_z = 0; //position to look at
double usr_d = user_dist;
double d[3] = {0., 0., 0.}; 
double image[3] = {0., 0., 0.}; 
bool images = false;
bool h_res = false; //high resolution

double ma_x=1000, ma_y=1000, ma_z=1000; //section limits
double mi_x=-1000, mi_y=-1000, mi_z=-1000;
double cl;
bool ifcl = false;

int num = 0;

float tone_list[5] = {1.,0.4,0.6,0.2,0.8};

class Dipole_g {
 private:
  //position
  double x, y, z;
  //orientation
  double ox, oy, oz;
  
 public:
   Dipole_g() : x(0), y(0), z(0), ox(0), oy(0), oz(1) {}
   ~Dipole_g() {}
   
  void set_pos(const double p[3]) { x = p[0]; y = p[1]; z = p[2]; }
  void get_pos(double p[3]) { p[0] = x; p[1] = y; p[2] = z; }
  void set_ori(const double o[3]) { ox = o[0]; oy = o[1]; oz = o[2]; }
  void get_ori(double o[3]) { o[0] = ox; o[1] = oy; o[2] = oz; }
};

//The particles
class Spherocyl_g {
private:
  //id
  int id;
  //cluster
  int clust;
  //length
  double l;
  //radius
  double r;
  //position
  double x, y, z;
  //orientation
  Matrix<double> o;
  //rotation
  double nx, ny, nz, ang;
  //latitudes and longitudes
  int lats, longs;
  //vertices
  Matrix<double> vert_x;
  Matrix<double> vert_y;
  Matrix<double> vert_z;
  //alpha
  float A;
  //dipole
  int ndips;
  Dynvec<Dipole_g> dip;

public:
  Spherocyl_g();
  Spherocyl_g(int, int);
  ~Spherocyl_g() {}

  void set_id(int ID) { id = ID; }
  int get_id() { return id; }
  void set_clust(int cl) { clust = cl; }
  int get_clust() { return clust; }
  void set_l(double ll) { l = ll; }
  double get_l() { return l; }
  void set_r(double rr) { r = rr; }
  double get_r() { return r; }
  void set_pos(const double p[3]) { x = p[0]; y = p[1]; z = p[2]; }
  void get_pos(double p[3]) { p[0] = x; p[1] = y; p[2] = z; }
  void set_omat(const Matrix<double>& oo) { o = oo; }
  void get_omat(Matrix<double>& oo) { oo = o; }
  void set_lalo(const int l[2]) { lats = l[0]; longs = l[1]; }
  void mul_lalo() { lats *= 2; longs *= 3; }
  void div_lalo() { lats /= 2; longs /= 3; }
  void vert_calc();
  void set_A(float a) { A = a; }
  float get_A() { return A; }
  void set_ndips(int n) { ndips = n; dip.set_size(ndips); }
  int get_ndips() { return ndips; }
  void set_dip(int i, double p[3], double o[3]) { dip[i].set_pos(p); dip[i].set_ori(o); }
  Dipole_g get_dip(int i) { return dip[i]; }
  void draw(bool);
};

//The timer
class Timer {
    private:
    //The clock time when the timer started
    int startTicks;
    
    //The ticks stored when the timer was paused
    int pausedTicks;
    
    //The timer status
    bool paused;
    bool started;
    
    public:
    //Initializes variables
    Timer();
    
    //The various clock actions
    void start();
    void stop();
    void pause();
    void unpause();
    
    //Gets the timer's time
    int get_ticks();
    
    //Checks the status of the timer
    bool is_started();
    bool is_paused();    
};


Spherocyl_g::Spherocyl_g() : id(0), l(0.), r(0.5), x(0.), y(0.), z(0.), 
		       lats(10), longs(12), ndips(1) {
  o.set_size(3,3);
  o[0][0] = 1.; o[0][1] = 0.; o[0][2] = 0.;
  o[1][0] = 0.; o[1][1] = 1.; o[1][2] = 0.; 
  o[2][0] = 0.; o[2][1] = 0.; o[2][2] = 1.;
  A = 1.;
  vert_calc();
  dip.set_size(1);
  double v[3] = {0, 0, 0};
  dip[0].set_pos(v);
  v[2] = 1;
  dip[0].set_ori(v);
}

Spherocyl_g::Spherocyl_g(int lts, int lngs) : id(0), l(0.), r(0.5), x(0.), y(0.), z(0.), 
			  lats(lts), longs(lngs), ndips(1) {
  o.set_size(3,3);
  o[0][0] = 1.; o[0][1] = 0.; o[0][2] = 0.;
  o[1][0] = 0.; o[1][1] = 1.; o[1][2] = 0.; 
  o[2][0] = 0.; o[2][1] = 0.; o[2][2] = 1.;
  A = 1.;
  vert_calc();
  dip.set_size(1);
  double v[3] = {0, 0, 0};
  dip[0].set_pos(v);
  v[2] = 1;
  dip[0].set_ori(v);
}

void Spherocyl_g::vert_calc() {
  vert_x.set_size(longs+1, 1);
  vert_y.set_size(longs+1, 1);
  vert_z.set_size(lats+2, 2);

//   double rsq = r*r;
//   double p4 = 4. * M_PI;

  for(int i = 0; i <= lats+1; ++i) {
    double lat = M_PI * (-0.5 + static_cast<double>(i - 1) / lats);
    vert_z[i][0] = r * sin(lat);
    vert_z[i][1] = r * cos(lat);
  }
  for(int j = 0; j <= longs; ++j) {
    double lng = 2 * M_PI * static_cast<double>(j - 1) / longs;
    vert_x[j][0] = cos(lng);
    vert_y[j][0] = sin(lng);
  }
  //relation between euler parameters and director cosines matrix
  double trace = o[0][0] + o[1][1] + o[2][2];
  if(1+trace > 1e-12){
    double a0 = sqrt(1+trace) / 2;
    ang = acos(a0) * 360. / (M_PI);
    nx = (o[1][2]-o[2][1]) / (4*a0);
    ny = (o[2][0]-o[0][2]) / (4*a0);
    nz = (o[0][1]-o[1][0]) / (4*a0);
    //
    double norm = sqrt(nx*nx+ny*ny+nz*nz);
    nx /= norm; ny /= norm; nz /= norm;
    //
  }else if(1+2*o[0][0]-trace > 1e-12){
    nx = sqrt(1+2*o[0][0]-trace) / 2;
    ny = (o[1][0]+o[0][1]) / (4*nx);
    nz = (o[2][0]+o[0][2]) / (4*nx);
    double a0 = (o[1][2]-o[2][1]) / (4*nx);
    ang = acos(a0) * 360. / (M_PI);
    //
    double norm = sqrt(nx*nx+ny*ny+nz*nz);
    nx /= norm; ny /= norm; nz /= norm;
    //
  }else{
    ang = 180.;
    if(fabs(o[0][0]-1) < 1e-12){
      nx = 1.; ny = 0.; nz = 0.;
    }else if(fabs(o[1][1]-1) < 1e-12){
      nx = 0.; ny = 1.; nz = 0.;
    }else if(fabs(o[2][2]-1) < 1e-12){
      nx = 0.; ny = 0.; nz = 1.;
    }
  }
  
//   double trace = o[0][0] + o[1][1] + o[2][2];
//   if(trace != -1){
//     ang = acos((trace-1.)/2) * 360. / (2 * M_PI);
//     double e0 = sqrt(trace+1.) / 2;
//     nx = (o[1][2]-o[2][1])/(4*e0);
//     ny = (o[2][0]-o[0][2])/(4*e0);
//     nz = (o[0][1]-o[1][0])/(4*e0);         
//   }else{
//     ang = 180. / (2 * M_PI);
//     if(o[0][0] == 1){
//       nx = 1.; ny = 0.; nz = 0.;
//     }else if(o[1][1] == 1){
//       nx = 0.; ny = 1.; nz = 0.;
//     }else if(o[2][2] == 1){
//       nx = 0.; ny = 0.; nz = 1.;
//     }
//   }
//   cout<<ang<<" "<<nx<<","<<ny<<","<<nz<<endl;
}

void Spherocyl_g::draw(bool draw_dip) {
  double rsq = r * r;
  glPushMatrix();
  glTranslated(x+image[0]+d[0], y+image[1]+d[1], z+image[2]+d[2]);
  glRotated(ang, nx, ny, nz);
  
  //spherocylinder
  for(int i = 0; i <= lats; ++i) {
    glBegin(GL_QUAD_STRIP);
    for(int j = 0; j <= longs; ++j) {
      double x0 = vert_x[j][0] * vert_z[i][1];
      double x1 = vert_x[j][0] * vert_z[i+1][1];
      double y0 = vert_y[j][0] * vert_z[i][1];
      double y1 = vert_y[j][0] * vert_z[i+1][1];
      
      if(vert_z[i][0] >= 0){
	if(vert_z[i][0]/r > 0.8)
	  glColor4f(0, 0, 1, A);//glColor4f(fabs(o[2][0]), fabs(o[2][1]), fabs(o[2][2]), A);
	else
	  glColor4f(tone_list[id], tone_list[id], tone_list[id], A);
	glNormal3f(x0 / rsq, y0 / rsq, vert_z[i][0] / rsq);
	glVertex3f(x0, y0, vert_z[i][0] + l/2);
	glNormal3f(x1 / rsq, y1 / rsq, vert_z[i+1][0] / rsq);
	glVertex3f(x1, y1, vert_z[i+1][0] + l/2);
      }else{
	glColor4f(tone_list[id], tone_list[id], tone_list[id], A);
	glNormal3f(x0 / rsq, y0 / rsq, vert_z[i][0] / rsq);
	glVertex3f(x0, y0, vert_z[i][0] - l/2);
	glNormal3f(x1 / rsq, y1 / rsq, vert_z[i+1][0] / rsq);
	glVertex3f(x1, y1, vert_z[i+1][0] - l/2);
      }
    }
    glEnd();
  }
  glBegin(GL_QUAD_STRIP);
  for(int j = 0; j <= longs; ++j) {
    double x0 = r * vert_x[j][0];
    double y0 = r * vert_y[j][0];
    
    if(x0 == r)
      glColor4f(1, 0, 0, A);//glColor4f(fabs(o[0][0]), fabs(o[0][1]), fabs(o[0][2]), A);
    else if(y0 == r)
      glColor4f(0, 1, 0, A);//glColor4f(fabs(o[1][0]), fabs(o[1][1]), fabs(o[1][2]), A);
    else
      glColor4f(tone_list[id], tone_list[id], tone_list[id], A); 
    glNormal3f(x0 / rsq, y0 / rsq, 0);
    glVertex3f(x0, y0, l/2);
    glNormal3f(x0 / rsq, y0 / rsq, 0);
    glVertex3f(x0, y0, -l/2); 
  }
  glEnd();
  glPopMatrix();
  
  //dipole
  if(draw_dip){
    double p[3], o[3];
    double mlt = 1.5 * r;
    glLineWidth(2);
    for(int i = 0; i < ndips; ++i){
      dip[i].get_pos(p); dip[i].get_ori(o);
      glColor3f(0, 0, 1);
      glBegin(GL_LINES);
	glVertex3f(p[0]-mlt*o[0], p[1]-mlt*o[1], p[2]-mlt*o[2]);
	glVertex3f(p[0], p[1], p[2]);
      glEnd();
      glColor3f(1, 0, 0);
      glBegin(GL_LINES);
	glVertex3f(p[0], p[1], p[2]);
	glVertex3f(p[0]+mlt*o[0], p[1]+mlt*o[1], p[2]+mlt*o[2]);
      glEnd();
    }
    glLineWidth(1);
  }
}


Timer::Timer() {
  //Initialize the variables
  startTicks = 0;
  pausedTicks = 0;
  paused = false;
  started = false;    
}

void Timer::start() {
  //Start the timer
  started = true;
  
  //Unpause the timer
  paused = false;
  
  //Get the current clock time
  startTicks = SDL_GetTicks();    
}

void Timer::stop() {
  //Stop the timer
  started = false;
  
  //Unpause the timer
  paused = false;    
}

void Timer::pause() {
  //If the timer is running and isn't already paused
  if( ( started == true ) && ( paused == false ) ){
    //Pause the timer
    paused = true;
    
    //Calculate the paused ticks
    pausedTicks = SDL_GetTicks() - startTicks;
  }
}

void Timer::unpause() {
  //If the timer is paused
  if( paused == true ){
    //Unpause the timer
    paused = false;
    
    //Reset the starting ticks
    startTicks = SDL_GetTicks() - pausedTicks;
    
    //Reset the paused ticks
    pausedTicks = 0;
  }
}

int Timer::get_ticks() {
  //If the timer is running
  if( started == true ){
    //If the timer is paused
    if( paused == true ){
      //Return the number of ticks when the timer was paused
      return pausedTicks;
    }else{
      //Return the current time minus the start time
      return SDL_GetTicks() - startTicks;
    }    
  }
  
  //If the timer isn't running
  return 0;    
}

bool Timer::is_started() {
  return started;    
}

bool Timer::is_paused() {
  return paused;    
}




void computeLocation() {
  glMatrixMode(GL_PROJECTION); // Set projection parameters.
  //glPushMatrix();
  glLoadIdentity();
  gluPerspective(60., 1., 0.01, 1000.);
  glTranslatef(0., 0., - usr_d);
  gluLookAt(usr_x, usr_y, usr_z,  ref_x, ref_y, ref_z,  0, 0, 1);
  //glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
}

// Initializes information for drawing within OpenGL.
bool init_GL() {
  glClearColor(1.0, 1.0, 1.0, 0.0);   // Set window color to white.
  
  glEnable(GL_DEPTH_TEST);            // Draw only closest surfaces
  glEnable(GL_COLOR_MATERIAL);        // Configure glColor().
  glEnable(GL_LIGHTING);              // Set up ambient light.
  glEnable(GL_LIGHT0);                // Set up lamp.
//   glEnable(GL_LIGHT1);                // Set up sunlight.
  glEnable(GL_NORMALIZE);
  
  //glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
//   glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
//   glMaterialfv(GL_FRONT, GL_SHININESS, mat_shine);

  glEnable(GL_BLEND);
  glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  glEnable(GL_CLIP_PLANE0); glEnable(GL_CLIP_PLANE1);
  glEnable(GL_CLIP_PLANE2); glEnable(GL_CLIP_PLANE3);
  glEnable(GL_CLIP_PLANE4); glEnable(GL_CLIP_PLANE5);
  
  //If there was any errors
  if( glGetError() != GL_NO_ERROR ){
    return false;    
  }
  computeLocation();
  //If everything initialized
  return true;
}

//gets input
bool input(int& Ncomp, Dynvec<int>& Npart, Dynvec<double>& length, Dynvec<double>& diameter,
	   Dynvec<int>& ndip){
  std::string line;
  
  ifstream In("input.txt");
  if(! In){
    cout << "Couldn't open " << "input.txt" << endl;
    In.close();
    return false;
  }
  In >> line >> line >> line >> line;
  In >> line >> line >> line >> line;
  In >> line >> line;
  In >> line >> line >> line >> Ncomp;
  Npart.set_size(Ncomp);
  In >> line >> line >> line;
  for(int i = 0; i < Ncomp; ++i)
    In >> Npart[i];
  length.set_size(Ncomp);
  In >> line;
  for(int i = 0; i < Ncomp; ++i)
    In >> length[i];
  diameter.set_size(Ncomp);
  In >> line;
  for(int i = 0; i < Ncomp; ++i)
    In >> diameter[i];
  In >> line;
  for(int i = 0; i < Ncomp; ++i)
    In >> line;
  In >> line >> line;
  for(int i = 0; i < 3; ++i)
    In >> line;    
  In >> line >> line;
  In >> line >> line;
  In >> line >> line;
  In >> line >> line;
  In >> line >> line >> line >> line >> line;
  In >> line >> line; 
  In >> line >> line >> line >> line;
  In >> line >> line >> line;
  In >> line >> line >> line >> line;
  In >> line >> line >> line;
  In >> line >> line >> line >> line;
  In >> line >> line >> line >> line;
  In >> line >> line >> line;
  In >> line >> line >> line;
  In >> line >> line >> line;
  In >> line >> line >> line;
  In >> line >> line >> line;
  //read dipole number
  ndip.set_size(Ncomp);
  In >> line;
  In >> line >> line >> line >> line >> line;
  for(int i = 0; i < Ncomp; ++i){
    In >> line >> ndip[i];
    for(int j = 0; j < ndip[i]; ++j){
      for(int a = 0; a < 7; ++a)
	In >> line;
    }
  }
  In.close();
  return true;
}

//charges configuration
bool chrg_conf(Dynvec<Dynvec<Spherocyl_g> >& part, 
	       const Dynvec<int>& Npart, int Ncomp, 
	       double L[3]){
  std::string dummy;
  char name[30];
  sprintf(name, "conf.dat");
  ifstream InPos(name);
  if(! InPos){
    InPos.close();
    return false;
  }
  int cl;
  double pi[3], o[3];
  Matrix<double> oi(3, 3);
  InPos >> L[0] >> L[1] >> L[2];
  for(int n = 0; n < Ncomp; ++n){
    int NN = Npart[n];
    for(int i = 0; i < NN; ++i){
      InPos >> pi[0] >> pi[1] >> pi[2]
	    >> dummy >> dummy >> dummy;
      for(int a = 0; a < 3; ++a){
	for(int b = 0; b < 3; ++b){
	  InPos >> oi[a][b];
	}
      }
      for(int k = 0; k < 3; ++k)
	pi[k] -=  L[k] * floor(pi[k] / L[k] + 0.5);
      part[n][i].set_pos(pi);
      part[n][i].set_omat(oi);
      part[n][i].vert_calc();
      for(int a = 0; a < 5; ++a)
	InPos >> dummy;
      int nd = part[n][i].get_ndips();
      for(int l = 0; l < nd; ++l){
	InPos >> pi[0] >> pi[1] >> pi[2]
	      >> o[0] >> o[1] >> o[2];
// 	//
// 	    cout<<pi[0] <<" "<< pi[1] <<" "<< pi[2]
// 	      <<" "<< o[0] <<" "<< o[1] <<" "<< o[2]<<"    ";
// 	//
	for(int k = 0; k < 3; ++k)
	  pi[k] -=  L[k] * floor(pi[k] / L[k] + 0.5);
	part[n][i].set_dip(l, pi, o);
	InPos >> dummy >> dummy;
// 	//
// 	Dipole_g u = part[n][i].get_dip(l);
// 	u.get_pos(pi); u.get_ori(o);
// 	cout<<pi[0] <<" "<< pi[1] <<" "<< pi[2]
// 	      <<" "<< o[0] <<" "<< o[1] <<" "<< o[2]<<endl;
// 	//
      }
      InPos >> cl;
      part[n][i].set_clust(cl);
    }
  }
  InPos.close();
  return true;
}

void draw_box(const double L[3]){
  double L2[3] = {L[0]/2, L[1]/2, L[2]/2};
  glPushMatrix();
  glBegin(GL_LINE_LOOP);
  glVertex3f(-L2[0]+d[0], -L2[1]+d[1], -L2[2]+d[2]);
  glVertex3f(L2[0]+d[0], -L2[1]+d[1], -L2[2]+d[2]);
  glVertex3f(L2[0]+d[0], L2[1]+d[1], -L2[2]+d[2]);
  glVertex3f(-L2[0]+d[0], L2[1]+d[1], -L2[2]+d[2]);
  glEnd();
  glBegin(GL_LINE_LOOP);
  glVertex3f(-L2[0]+d[0], -L2[1]+d[1], L2[2]+d[2]);
  glVertex3f(L2[0]+d[0], -L2[1]+d[1], L2[2]+d[2]);
  glVertex3f(L2[0]+d[0], L2[1]+d[1], L2[2]+d[2]);
  glVertex3f(-L2[0]+d[0], L2[1]+d[1], L2[2]+d[2]);
  glEnd();
  glBegin(GL_LINES);
    glVertex3f(-L2[0]+d[0], -L2[1]+d[1], -L2[2]+d[2]);
    glVertex3f(-L2[0]+d[0], -L2[1]+d[1], L2[2]+d[2]);
  glEnd();
  glBegin(GL_LINES);
    glVertex3f(L2[0]+d[0], -L2[1]+d[1], -L2[2]+d[2]);
    glVertex3f(L2[0]+d[0], -L2[1]+d[1], L2[2]+d[2]);
  glEnd();
  glBegin(GL_LINES);
    glVertex3f(L2[0]+d[0], L2[1]+d[1], -L2[2]+d[2]);
    glVertex3f(L2[0]+d[0], L2[1]+d[1], L2[2]+d[2]);
  glEnd();
  glBegin(GL_LINES);
    glVertex3f(-L2[0]+d[0], L2[1]+d[1], -L2[2]+d[2]);
    glVertex3f(-L2[0]+d[0], L2[1]+d[1], L2[2]+d[2]);
  glEnd();
  glPopMatrix();
}

// Draws the current image.
void show(Dynvec<Dynvec<Spherocyl_g> >& part, 
	  const Dynvec<int>& Npart, int Ncomp, const double L[3], 
	  bool toggle[], bool draw_dip){
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // Clear window.
  glColor3f(1., 1., 1.);
  glShadeModel(/*GL_FLAT*/GL_SMOOTH);
  GLfloat ambientColor[] = {0.5f, 0.5f, 0.5f, 1.0f};
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambientColor);
  //directer light
  GLfloat lightColor0[] = {0.7f, 0.7f, 0.7f, 1.0f};
  GLfloat lightPos0[] = {usr_x, usr_y, usr_z, 0.0f};
  glLightfv(GL_LIGHT0, GL_DIFFUSE, lightColor0);
  glLightfv(GL_LIGHT0, GL_POSITION, lightPos0);
// //   //positioned ligth
// //   GLfloat lightColor1[] = {0.2f, 0.2f, 0.2f, 1.0f};
// //   GLfloat lightPos1[] = {0.0f, 0.0f, 0.0f, 1.0f};
// //   glLightfv(GL_LIGHT1, GL_DIFFUSE, lightColor1);
// //   glLightfv(GL_LIGHT1, GL_POSITION, lightPos1);
//   //clipping planes
//   double eqmi_x[] = { -1, 0, 0, mi_x};
//   double eqma_x[] = { 1, 0, 0, ma_x};
//   double eqmi_y[] = { 0, -1, 0, mi_y};
//   double eqma_y[] = { 0, 1, 0, ma_y};
//   double eqmi_z[] = { 0, 0, -1, mi_z};
//   double eqma_z[] = { 0, 0, 1, ma_z};
//   glClipPlane(GL_CLIP_PLANE0 , eqmi_x);
//   glClipPlane(GL_CLIP_PLANE1 , eqma_x);
//   glClipPlane(GL_CLIP_PLANE2 , eqmi_y);
//   glClipPlane(GL_CLIP_PLANE3 , eqma_y);
//   glClipPlane(GL_CLIP_PLANE4 , eqmi_z);
//   glClipPlane(GL_CLIP_PLANE5 , eqma_z);
  
  int mimag = 0;
  if(images) mimag = 2;
  for(int nz = 0; nz <= mimag; ++nz){
    for(int ny = 0; ny <= mimag; ++ny){
      for(int nx = 0; nx <= mimag; ++nx){
	if(images){
	  image[0] = L[0]*(nx - 1); 
	  image[1] = L[1]*(ny - 1);
	  image[2] = L[2]*(nz - 1);
	}
	for(int n = 0; n < Ncomp; ++n){
	  if(toggle[n]){
	    int NN = Npart[n];
	    double p[3];
	    for(int i = 0; i < NN; ++i){
	      if(!ifcl){
		part[n][i].get_pos(p);
		if(p[0] <= ma_x && p[1] <= ma_y && p[2] <= ma_z){
		  part[n][i].draw(draw_dip);
		}
	      }else if(part[n][i].get_clust() == cl){
		part[n][i].get_pos(p);
		if(p[0] <= ma_x && p[1] <= ma_y && p[2] <= ma_z){
		  part[n][i].draw(draw_dip);
		}
	      }
	    }
	  }
	}
      }
    }
  }
  image[0] = image[1] = image[2] = 0;
  glColor3f(0.6, 0.3, 0.3);
  draw_box(L);
}

void calc_usr(){
  double cost = cos(user_theta);
  double sint = sin(user_theta);
  double cosp = cos(user_phi);
  double sinp = sin(user_phi);
  usr_x = user_dist * sint * cosp;
  usr_y = user_dist * sint * sinp;
  usr_z = user_dist * cost;
  usr_d = sqrt(usr_x * usr_x + usr_y * usr_y + usr_z * usr_z); // distance to origin
  //Adjust the view
  computeLocation();
}

//writes an output ps
void writefile(int format, int sort, int options, int nbcol, char *filename, 
	       const char *extension, 
	       Dynvec<Dynvec<Spherocyl_g> >& part, 
	       const Dynvec<int>& Npart, int Ncomp, double L[3], 
	       bool toggle[], bool draw_dip){
  FILE *fp;
  char file[256];
  int state = GL2PS_OVERFLOW, buffsize = 0;
  GLint viewport[4];

  strcpy(file, filename);
  strcat(file, ".");
  strcat(file, extension);

  viewport[0] = 0;
  viewport[1] = 0;
  viewport[2] = SCREEN_WIDTH;
  viewport[3] = SCREEN_HEIGHT;

  fp = fopen(file, "wb");

  if(!fp){
    printf("Unable to open file %s for writing\n", file);
    exit(1);
  }

  printf("Saving image to file %s... ", file);
  fflush(stdout);

  while(state == GL2PS_OVERFLOW){
    buffsize += 1024*1024;
    gl2psBeginPage(file, "gl2psTest", viewport, format, sort, options,
                   GL_RGBA, 0, NULL, nbcol, nbcol, nbcol,
                   buffsize, fp, file);
    show(part, Npart, Ncomp, L, toggle, draw_dip);
    state = gl2psEndPage();
  }

  fclose(fp);

  printf("Done!\n");
  fflush(stdout);
}

void handle_input(Dynvec<Dynvec<Spherocyl_g> >& part, 
		  const Dynvec<int>& Npart, int Ncomp, double L[3], 
		  bool toggle[], bool& draw_dip){
  char ext[32];
  char name[32];
  float a, val;
  int NN, ll;
  //If a key was pressed
  if( event.type == SDL_KEYDOWN ){
    //Adjust the velocity
    switch( event.key.keysym.sym ){
    //translate view***************
    case SDLK_q:
      d[0] += 0.5;
      break;
    case SDLK_a:
      d[0] -= 0.5;
      break;
    case SDLK_w:
      d[1] += 0.5;
      break;
    case SDLK_s:
      d[1] -= 0.5;
      break;
    case SDLK_e:
      d[2] += 0.5;
      break;
    case SDLK_d:
      d[2] -= 0.5;
      break;
    //rotate view******************
    case SDLK_UP:
      user_theta += 0.1;
      if(user_theta > M_PI)
	user_theta = M_PI;
      calc_usr();
      break;
    case SDLK_DOWN:
      user_theta -= 0.1;
      if(user_theta < 0.1)
	user_theta = 0.1;
      calc_usr();
      break;
    case SDLK_LEFT: 
      user_phi  += 0.1;
      calc_usr();
      break;
    case SDLK_RIGHT:
      user_phi  -= 0.1;
      calc_usr();
      break;
    //zoom in/out********************
    case SDLK_PAGEUP:
      user_dist += 0.5;
      calc_usr();
      break;
    case SDLK_PAGEDOWN: 
      user_dist -= 0.5;
      calc_usr();
      break;
    //update view*******************
    case SDLK_SPACE:
      if(! chrg_conf(part, Npart, Ncomp, L))
	cout << "problem opening ptPosF.dat" << endl;
      break;
    //change mesh resolution********
    case SDLK_h:
      if(! h_res){
	for(int n = 0; n < Ncomp; ++n){
	  NN = Npart[n];
	  for(int i = 0; i < NN; ++i){
	    part[n][i].mul_lalo();
	    part[n][i].vert_calc();
	  }
	}
	h_res = true;
      }else{
	for(int n = 0; n < Ncomp; ++n){
	  NN = Npart[n];
	  for(int i = 0; i < NN; ++i){
	    part[n][i].div_lalo();
	    part[n][i].vert_calc();
	  }
	}
	h_res = false;
      }
      break;
    //draw dipoles
    case SDLK_TAB:
      if(draw_dip)
	draw_dip = false;
      else
	draw_dip = true;
      break;
    //show clusters********************
    case SDLK_c:
      if(ifcl)
	ifcl = false;
      else
	ifcl = true;
      if(ifcl){
	ifstream IN("cl_number.dat");
	if(!IN){
	  cout<<"proble opening cl_number.dat\n";
	  IN.close();
	  return;
	}
	IN >> cl;
	IN.close();
      }
      break;
    case SDLK_KP_MINUS:
      cl--;
      cout<<cl<<endl;
      break;
    case SDLK_KP_PLUS:
      cl++;
      cout<<cl<<endl;
      break;
    //show images
    case SDLK_t:
      if(images) images = false;
      else images = true;
      break;
    //toggle on/off components******
    case SDLK_KP0:
      if(Ncomp > 0){
	if(toggle[0])
	  toggle[0] = false; 
	else
	  toggle[0] = true;
      }
      break;
    case SDLK_KP1:
      if(Ncomp > 1){
	if(toggle[1])
	  toggle[1] = false;
	else
	  toggle[1] = true;
      }
      break;
    case SDLK_KP2:
      if(Ncomp > 2){
	if(toggle[2])
	  toggle[2] = false;
	else
	  toggle[2] = true;
      }
      break;
    case SDLK_KP3:
      if(Ncomp > 3){
	if(toggle[3])
	  toggle[3] = false;
	else
	  toggle[3] = true;
      }
      break;
    //print to a heavy eps file*******
    case SDLK_p:
      sprintf(name, "configuration%d", num);
      strcpy(ext, gl2psGetFileExtension(GL2PS_EPS));
      writefile(GL2PS_EPS, GL2PS_SIMPLE_SORT, GL2PS_NONE, 0, name, ext, 
		part, Npart, Ncomp, L, toggle, draw_dip);
      num++;
      break;
    //clip the system*****************
    case SDLK_i:
      if(event.key.keysym.mod & KMOD_SHIFT){
	if(ma_x > 0)
	  ma_x -= L[0]/20;
      }else{
	if(mi_x > 0)
	  mi_x -= L[0]/20;
      }
      break;
    case SDLK_j:
      if(event.key.keysym.mod & KMOD_SHIFT){
	if(ma_x < 3.*L[0]/2)
	  ma_x += L[0]/20;
      }else{
	if(mi_x < 3.*L[0]/2)
	  mi_x += L[0]/20;
      }
      break; 
    case SDLK_m:
      if(event.key.keysym.mod & KMOD_SHIFT){
	if(ma_y > 0)
	  ma_y -= L[1]/20;
      }else{
	if(mi_y > 0)
	  mi_y -= L[1]/20;
      }
      break;
    case SDLK_n:
      if(event.key.keysym.mod & KMOD_SHIFT){
	if(ma_y < 3.*L[1]/2)
	  ma_y += L[1]/20;
      }else{
	if(mi_y < 3.*L[1]/2)
	  mi_y += L[1]/20;
      }
      break;
    case SDLK_k:
      if(event.key.keysym.mod & KMOD_SHIFT){
	if(ma_z > 0)
	  ma_z -= L[2]/20;
      }else{
	if(mi_z > 0)
	  mi_z -= L[2]/20;
      }
      break;
    case SDLK_l:
      if(event.key.keysym.mod & KMOD_SHIFT){
	if(ma_z < 3.*L[2]/2)
	  ma_z += L[2]/20;
      }else{
	if(mi_z < 3.*L[2]/2)
	  mi_z += L[2]/20;
      }
      break;
    //make components transparent*****
    case SDLK_F1:
      if(Ncomp > 0){
	a = part[0][0].get_A();
	if(a < 1)
	  val = 1;
	else
	  val = 0.3;
	NN = Npart[0];
	for(int i = 0; i < NN ; ++i){
	  part[0][i].set_A(val);
	}
      }
      break;
    case SDLK_F2:
      if(Ncomp > 1){
	a = part[1][0].get_A();
	if(a < 1)
	  val = 1;
	else
	  val = 0.3;
	NN = Npart[1];
	for(int i = 0; i < NN ; ++i){
	  part[1][i].set_A(val);
	}
      }
      break;
    case SDLK_F3:
      if(Ncomp > 2){
	a = part[2][0].get_A();
	if(a < 1)
	  val = 1;
	else
	  val = 0.3;
	NN = Npart[2];
	for(int i = 0; i < NN ; ++i){
	  part[2][i].set_A(val);
	}
      }
      break;
    case SDLK_F4:
      if(Ncomp > 3){
	a = part[3][0].get_A();
	if(a < 1)
	  val = 1;
	else
	  val = 0.3;
	NN = Npart[3];
	for(int i = 0; i < NN ; ++i){
	  part[3][i].set_A(val);
	}
      }
      break;
    }
    //Show the scene////////////////////
    show(part, Npart, Ncomp, L, toggle, draw_dip);
  }
}

bool init(){
  //Initialize SDL
  if( SDL_Init( SDL_INIT_EVERYTHING ) < 0 ){
    return false;    
  }
  SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);//disable(0)-enable(1) double buffering
  //Create Window
  if( SDL_SetVideoMode( SCREEN_WIDTH, SCREEN_HEIGHT, SCREEN_BPP, 
			SDL_OPENGL ) == NULL ){
    return false;
  }
  //Initialize OpenGL
  if( init_GL() == false ){
    return false;    
  }
  //Set caption
  SDL_WM_SetCaption( "Visualization", NULL );
  return true;    
}

void clean_up(){
  //Quit SDL
  SDL_Quit();
}

int main( int argc, char *argv[] ){
  //Quit flag
  bool quit = false;
  //Initialize
  if( init() == false ){
    return 1;    
  }

  int Ncomp;
  Dynvec<int> Npart;
  double L[3] = {0.1, 0.1, 0.1};
  Dynvec<double> length;
  Dynvec<double> diameter;
  Dynvec<int> ndip;
  input(Ncomp, Npart, length, diameter, ndip);
//   //
//   cout << Ncomp << endl;
//   for(int i = 0; i < Ncomp; ++i)
//     cout << Npart[i] << " ";
//   cout << endl;
//   for(int i = 0; i < Ncomp; ++i)
//     cout << length[i] << " ";
//   cout << endl;
//   for(int i = 0; i < Ncomp; ++i)
//     cout << diameter[i] << " ";
//   cout << endl;
//   for(int i = 0; i < Ncomp; ++i)
//     cout << ndip[i] << " ";
//   cout << endl;
//   //

  //The frame rate regulator
  Timer fps;

  //the particles
  Dynvec<Dynvec<Spherocyl_g> > part(Ncomp);
  for(int n = 0; n < Ncomp; ++n)
    part[n].set_size(Npart[n]);
  double r;
  int l[2];
  bool toggle[Ncomp]; //for toggling components
  bool draw_dip = true; // for drawing dipoles 
  
  for(int n = 0; n < Ncomp; ++n){
    toggle[n] = true;
    int NN = Npart[n];
    l[0] = 10; l[1] = 12;
    for(int i = 0; i < NN; ++i){
      part[n][i].set_id(n);
      part[n][i].set_l(length[n]);
      part[n][i].set_r(diameter[n]/2);
      part[n][i].set_lalo(l);
      part[n][i].set_ndips(ndip[n]);
    }
  }
   
  //charge initiel configuration
  if(! chrg_conf(part, Npart, Ncomp, L))
    cout << "problem opening conf.dat" << endl;
//   //
//   Matrix<double> on(3,3);
//   part[0][1].get_omat(on);
//   on.print();
//   //
  //section limits
  ma_x = L[0]/2; ma_y = L[1]/2; ma_z = L[2]/2;   
  //Show the scene
  show(part, Npart, Ncomp, L, toggle, draw_dip);
  //Update screen
  SDL_GL_SwapBuffers();

  //Wait for user exit
  while( quit == false ){
    //Start the frame timer
    fps.start();
    //While there are events to handle
    while( SDL_PollEvent( &event ) ){
      //Handle key presses
      handle_input(part, Npart, Ncomp, L, toggle, draw_dip);
            
      if( event.type == SDL_QUIT ){
	quit = true;
      }
    }
    
    //Show the scene (animation) 
    if(! chrg_conf(part, Npart, Ncomp, L))
      cout << "problem opening ptPosF.dat" << endl;
    show(part, Npart, Ncomp, L, toggle, draw_dip);
    
    //Update screen
    SDL_GL_SwapBuffers();

    //Cap the frame rate
    if( fps.get_ticks() < 5000 / FRAMES_PER_SECOND ){
      SDL_Delay( ( 5000 / FRAMES_PER_SECOND ) - fps.get_ticks() );
    }
  }
  
  //Clean up
  clean_up();

  return 0;
}
