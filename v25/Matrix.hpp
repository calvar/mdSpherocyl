#ifndef __MATRIX_HPP
#define __MATRIX_HPP

#include <iostream>
#include <fstream>
#include <exception>
using namespace std;

#define RESTRICT __restrict__

class BadSize: public exception {
  virtual const char* what() const throw()
  {
    return "Error! Atempting a non allowed Matrix operation!\n";
  }
};

class BadSizeTa: public exception {
  virtual const char* what() const throw()
  {
    return "Error! Atempting a non allowed Tarray equality\n";
  }
};

class BadSizeV: public exception {
  virtual const char* what() const throw()
  {
    return "Error! Atempting a non allowed Dynvec equality\n";
  }
};

template <typename T>
class Matrix {
protected:
  int m_rows;
  int m_cols;
  T** m_data;
  
  friend Matrix<double> transpose(const Matrix<double>& RESTRICT);
  friend Matrix<double> LU_dec(const Matrix<double>& RESTRICT,
			       int* RESTRICT, int&, bool&);
  friend Matrix<double> fwd_bkd(const Matrix<double>& RESTRICT,
				const double* RESTRICT);
  friend Matrix<double> solve(const Matrix<double>& RESTRICT ,
			      const Matrix<double>& RESTRICT);
  friend Matrix<double> inverse(const Matrix<double>& RESTRICT);

public:
  Matrix();
  Matrix(int, int);
  Matrix(int, int, T);
  Matrix(const Matrix<T>& RESTRICT);
  ~Matrix();

  void get_size(int& rows, int& cols) {
    rows = m_rows;
    cols = m_cols; }
  void get_size(int& rows, int& cols) const {
    rows = m_rows;
    cols = m_cols; }
  void set_size(int, int);
  //operators
  T* operator[](int i) { return m_data[i]; }
  T const* operator[](int i) const{ return m_data[i]; }
  void operator=(const Matrix<T>& RESTRICT);
  void operator+=(const Matrix<T>& RESTRICT);
  void operator-=(const Matrix<T>& RESTRICT);
  const Matrix<T> operator+(const Matrix<T>& RESTRICT);
  const Matrix<T> operator+(const Matrix<T>& RESTRICT) const;
  const Matrix<T> operator-(const Matrix<T>& RESTRICT);
  const Matrix<T> operator-(const Matrix<T>& RESTRICT) const;
  bool operator==(const Matrix<T>& RESTRICT);
  bool operator==(const Matrix<T>& RESTRICT) const;
  bool operator!=(const Matrix<T>& RESTRICT);
  bool operator!=(const Matrix<T>& RESTRICT) const;
  const Matrix<T> operator*(const Matrix<T>& RESTRICT);
  const Matrix<T> operator*(const Matrix<T>& RESTRICT) const;
  const Matrix<T> operator*(T c);
  const Matrix<T> operator*(T c) const;
  const Matrix<T> transpose();
  void set_to(T);
  bool print();
  bool print() const;
};

template <typename T>
class Tarray {
protected:
  int m_size;
  T** m_data;
  
public:
  Tarray();
  Tarray(int);
  Tarray(int, T);
  Tarray(const Tarray<T>& RESTRICT);
  ~Tarray();

  int get_size() { return m_size; }
  int get_size() const { return m_size; }
  void set_size(int);
  //operators
  T* operator[](int i) { return m_data[i]; }
  T const* operator[](int i) const{ return m_data[i]; }
  void operator=(const Tarray<T>& RESTRICT);
  void set_to(T);
  bool print();
  bool print() const;
};

template <typename T>
class Dynvec {
 protected:
  int m_size;
  T* m_data;

 public:
  Dynvec();
  Dynvec(int);
  Dynvec(int, T);
  Dynvec(const Dynvec<T>& RESTRICT);
  ~Dynvec();

  int get_size() { return m_size; }
  int get_size() const { return m_size; }
  void set_size(int);
  T& operator[](int i) { return m_data[i]; }
  T const& operator[](int i) const{ return m_data[i]; }
  void operator()(int i, T val) { m_data[i] = val; }
  void operator=(const Dynvec<T>& RESTRICT);
  void set_to(T);
};

/* template <typename T> */
/* class Quat { */
/*  public: */
/*   Quat(); */
/*   Quat(T, T, T, T); */
/*   Quat(const Quat<T>&); */
/*   ~Quat() {} */

/*   Matrix<T> v; */
/*   //operators */
/*   void operator=(const Quat<T>&); */
/*   const Quat<T> operator+(const Quat<T>&); */
/*   const Quat<T> operator+(const Quat<T>&) const; */
/*   const Quat<T> operator-(const Quat<T>&); */
/*   const Quat<T> operator-(const Quat<T>&) const; */
/*   const Quat<T> operator*(const Quat<T>&); */
/*   const Quat<T> operator*(const Quat<T>&) const; */
/*   const Quat<T> conj(); */
/*   const Quat<T> conj() const; */
/* }; */

//MATRIX methods*********************************************************

//constructors---------------
template<typename T>
inline Matrix<T>::Matrix() try: m_rows(1), m_cols(1) {
  m_data = new T* [m_rows];
  for(int m = 0; m < m_rows; ++m)
    m_data[m] = new T[m_cols];
}catch(bad_alloc& ba){
  ofstream Err("Error.dat", ios::app);
  Err << "Error allocating memory in Matrix() constructor: "
      << ba.what() << "\n";
  Err.close();
}

template<typename T>
inline Matrix<T>::Matrix(int rows, int cols) try: m_rows(rows), m_cols(cols) {
  m_data = new T* [m_rows];
  for(int m = 0; m < m_rows; ++m)
    m_data[m] = new T[m_cols];
}catch(bad_alloc& ba){
  ofstream Err("Error.dat", ios::app);
  Err << "Error allocating memory in Matrix(int, int) constructor: "
      << ba.what() << "\n";
  Err.close();
}

template<typename T>
inline Matrix<T>::Matrix(int rows, int cols, T val) try: m_rows(rows), m_cols(cols) {
  m_data = new T* [m_rows];
  for(int m = 0; m < m_rows; ++m){
    m_data[m] = new T[m_cols];
    for(int n = 0; n < m_cols; ++n)
      m_data[m][n] = val;
  }
}catch(bad_alloc& ba){
  ofstream Err("Error.dat", ios::app);
  Err << "Error allocating memory in Matrix(int, int, T) constructor: "
      << ba.what() << "\n";
  Err.close();
}

template<typename T>
inline Matrix<T>::Matrix(const Matrix<T>& RESTRICT M) try{
  m_rows = M.m_rows;
  m_cols = M.m_cols;
  m_data = new T* [m_rows];
  for(int m = 0; m < m_rows; ++m){
    m_data[m] = new T[m_cols];
    for(int n = 0; n < m_cols; ++n){
      m_data[m][n] = M.m_data[m][n];
    }
  }
}catch(bad_alloc& ba){
  ofstream Err("Error.dat", ios::app);
  Err << "Error allocating memory in Matrix copy constructor: "
      << ba.what() << "\n";
  Err.close();
}

//destructor----------------
template<typename T>
inline Matrix<T>::~Matrix() {
  for(int m = 0; m < m_rows; ++m)
    delete[] m_data[m];
  delete[] m_data;
}

//set size------------------
template<typename T>
inline void Matrix<T>::set_size(int rows, int cols) try{
  for(int m = 0; m < m_rows; ++m)
    delete[] m_data[m];
  delete[] m_data;
  m_rows = rows;
  m_cols = cols;
  m_data = new T* [m_rows];
  for(int m = 0; m < m_rows; ++m)
    m_data[m] = new T[m_cols];
}catch(bad_alloc& ba){
  ofstream Err("Error.dat", ios::app);
  Err << "Error allocating memory in Matrix::set_size: "
      << ba.what() << "\n";
  Err.close();
}

//operators-----------------
template<typename T>
inline void Matrix<T>::operator=(const Matrix<T>& RESTRICT M) try{
  if(M.m_rows != m_rows || M.m_cols != m_cols){
    BadSize bads; throw bads; }
  for(int m = 0; m < m_rows; ++m){
    for(int n = 0; n < m_cols; ++n){
      m_data[m][n] = M.m_data[m][n];
    }
  }
}catch(exception& e){
  ofstream Err("Error.dat", ios::app);
  Err << e.what();
  Err.close();
}

template<typename T>
inline void Matrix<T>::operator+=(const Matrix<T>& RESTRICT M) try{
  if(M.m_rows != m_rows || M.m_cols != m_cols){
    BadSize bads; throw bads; }
  for(int m = 0; m < m_rows; ++m){
    for(int n = 0; n < m_cols; ++n){
      m_data[m][n] += M.m_data[m][n];
    }
  }
}catch(exception& e){
  ofstream Err("Error.dat", ios::app);
  Err << e.what();
  Err.close();
}

template<typename T>
inline void Matrix<T>::operator-=(const Matrix<T>& RESTRICT M) try{
  if(M.m_rows != m_rows || M.m_cols != m_cols){
    BadSize bads; throw bads; }
  for(int m = 0; m < m_rows; ++m){
    for(int n = 0; n < m_cols; ++n){
      m_data[m][n] -= M.m_data[m][n];
    }
  }
}catch(exception& e){
  ofstream Err("Error.dat", ios::app);
  Err << e.what();
  Err.close();
}

template<typename T>
inline const Matrix<T> Matrix<T>::operator+(const Matrix<T>& RESTRICT M) try{
  if(M.m_rows != m_rows || M.m_cols != m_cols){
    BadSize bads; throw bads; }
  Matrix<T> result(*this);
  result += M;
  return result;
}catch(exception& e){
  ofstream Err("Error.dat", ios::app);
  Err << e.what();
  Err.close();
}

template<typename T>
inline const Matrix<T> Matrix<T>::operator+(const Matrix<T>& RESTRICT M) const try{
  if(M.m_rows != m_rows || M.m_cols != m_cols){
    BadSize bads; throw bads; }
  Matrix<T> result(*this);
  result += M;
  return result;
}catch(exception& e){
  ofstream Err("Error.dat", ios::app);
  Err << e.what();
  Err.close();
}

template<typename T>
inline const Matrix<T> Matrix<T>::operator-(const Matrix<T>& RESTRICT M) try{
  if(M.m_rows != m_rows || M.m_cols != m_cols){
    BadSize bads; throw bads; }
  Matrix<T> result(*this);
  result -= M;
  return result;
}catch(exception& e){
  ofstream Err("Error.dat", ios::app);
  Err << e.what();
  Err.close();
}

template<typename T>
inline const Matrix<T> Matrix<T>::operator-(const Matrix<T>& RESTRICT M) const try{
  if(M.m_rows != m_rows || M.m_cols != m_cols){
    BadSize bads; throw bads; }
  Matrix<T> result(*this);
  result -= M;
  return result;
}catch(exception& e){
  ofstream Err("Error.dat", ios::app);
  Err << e.what();
  Err.close();
}


template<typename T>
inline bool Matrix<T>::operator==(const Matrix<T>& RESTRICT M) {
  if(M.m_rows != m_rows || M.m_cols != m_cols)
    return false;
  for(int m = 0; m < m_rows; ++m){
    for(int n = 0; n < m_cols; ++n){
      if(M.m_data[m][n] != m_data[m][n])
	return false;
    }
  }
  return true;
}

template<typename T>
inline bool Matrix<T>::operator==(const Matrix<T>& RESTRICT M) const{
  if(M.m_rows != m_rows || M.m_cols != m_cols)
    return false;
  for(int m = 0; m < m_rows; ++m){
    for(int n = 0; n < m_cols; ++n){
      if(M.m_data[m][n] != m_data[m][n])
	return false;
    }
  }
  return true;
}

template<typename T>
inline bool Matrix<T>::operator!=(const Matrix<T>& RESTRICT M) {
  if(M == this)
    return false;
  else
    return true;
}

template<typename T>
inline bool Matrix<T>::operator!=(const Matrix<T>& RESTRICT M) const{
  if(M == this)
    return false;
  else
    return true;
}

template<typename T>
inline const Matrix<T> Matrix<T>::operator*(const Matrix<T>& RESTRICT M) try{
  if(M.m_rows != m_cols){
    BadSize bads; throw bads; }
  Matrix<T> result(m_rows, M.m_cols, 0);
  for(int m = 0; m < m_rows; ++m){
    for(int n = 0; n < M.m_cols; ++n){
      for(int i = 0; i < m_cols; ++i)
	result.m_data[m][n] += m_data[m][i] * M.m_data[i][n];
    }
  }
  return result;
}catch(exception& e){
  ofstream Err("Error.dat", ios::app);
  Err << e.what();
  Err.close();
}

template<typename T>
inline const Matrix<T> Matrix<T>::operator*(const Matrix<T>& RESTRICT M) const try{
  if(M.m_rows != m_cols){
    BadSize bads; throw bads; }
  Matrix<T> result(m_rows, M.m_cols, 0);
  for(int m = 0; m < m_rows; ++m){
    for(int n = 0; n < M.m_cols; ++n){
      for(int i = 0; i < m_cols; ++i)
	result.m_data[m][n] += m_data[m][i] * M.m_data[i][n];
    }
  }
  return result;
}catch(exception& e){
  ofstream Err("Error.dat", ios::app);
  Err << e.what();
  Err.close();
}

template<typename T>
inline const Matrix<T> Matrix<T>::operator*(T c) {
  Matrix<T> result(m_rows, m_cols);
  for(int m = 0; m < m_rows; ++m){
    for(int n = 0; n < m_cols; ++n)
      result.m_data[m][n] = m_data[m][n] * c;
  }
  return result;
}

template<typename T>
inline const Matrix<T> Matrix<T>::operator*(T c) const{
  Matrix<T> result(m_rows, m_cols);
  for(int m = 0; m < m_rows; ++m){
    for(int n = 0; n < m_cols; ++n)
      result.m_data[m][n] = m_data[m][n] * c;
  }
  return result;
}

template<typename T>
inline void Matrix<T>::set_to(T val) {
  for(int m = 0; m < m_rows; ++m){
    for(int n = 0; n < m_cols; ++n){
      m_data[m][n] = val;
    }
  }
}

//print---------------------
template<typename T>
bool Matrix<T>::print() {
  std::ofstream Mat("mat.dat", ios::app);
  if(! Mat){
    Mat.close();
    return false;
  }
  for(int i = 0; i < m_rows; ++i){
    for(int j = 0; j < m_cols; ++j){
      Mat << m_data[i][j] << " ";
    }
    Mat << "\n";
  }
  Mat << endl;
  Mat.close();
  return true;
}

template<typename T>
bool Matrix<T>::print() const{
  std::ofstream Mat("mat.dat", ios::app);
  if(! Mat){
    Mat.close();
    return false;
  }
  for(int i = 0; i < m_rows; ++i){
    for(int j = 0; j < m_cols; ++j){
      Mat << m_data[i][j] << " ";
    }
    Mat << "\n";
  }
  Mat << endl;
  Mat.close();
  return true;
}


//TARRAY methods*********************************************************

//consructors---------------
template<typename T>
inline Tarray<T>::Tarray() try: m_size(1) {
  m_data = new T* [m_size];
  for(int n = 0; n < m_size; ++n)
    m_data[n] = new T[m_size-n];
}catch(bad_alloc& ba){
  ofstream Err("Error.dat", ios::app);
  Err << "Error allocating memory in Tarray() constructor: "
      << ba.what() << "\n";
  Err.close();
}

template<typename T>
inline Tarray<T>::Tarray(int size) try: m_size(size) {
  m_data = new T* [m_size];
  for(int n = 0; n < m_size; ++n)
    m_data[n] = new T[m_size-n];
}catch(bad_alloc& ba){
  ofstream Err("Error.dat", ios::app);
  Err << "Error allocating memory in Tarray(int) constructor: "
      << ba.what() << "\n";
  Err.close();
}

template<typename T>
inline Tarray<T>::Tarray(int size, T val) try: m_size(size) {
  m_data = new T* [m_size];
  for(int n = 0; n < m_size; ++n){
    m_data[n] = new T[m_size-n];
    for(int i = 0; i < m_size-n; ++i)
      m_data[n][i] = val;
  }
}catch(bad_alloc& ba){
  ofstream Err("Error.dat", ios::app);
  Err << "Error allocating memory in Tarray(int, T) constructor: "
      << ba.what() << "\n";
  Err.close();
}
  
template<typename T>
inline Tarray<T>::Tarray(const Tarray<T>& RESTRICT Ta) try{
  m_size = Ta.m_size;
  m_data = new T* [m_size];
  for(int n = 0; n < m_size; ++n){
    m_data[n] = new T[m_size-n];
    for(int i = 0; i < m_size-n; ++i){
      m_data[n][i] = Ta.m_data[n][i];
    }
  }
}catch(bad_alloc& ba){
  ofstream Err("Error.dat", ios::app);
  Err << "Error allocating memory in Tarray copy constructor: "
      << ba.what() << "\n";
  Err.close();
}

//destructor-----------------
template<typename T>
inline Tarray<T>::~Tarray() {
  for(int n = 0; n < m_size; ++n)
    delete[] m_data[n];
  delete[] m_data;
}

//set size------------------
template<typename T>
inline void Tarray<T>::set_size(int size) try{
  for(int n = 0; n < m_size; ++n)
    delete[] m_data[n];
  delete[] m_data;
  m_size = size;
  m_data = new T* [m_size];
  for(int n = 0; n < m_size; ++n)
    m_data[n] = new T[m_size-n];
}catch(bad_alloc& ba){
  ofstream Err("Error.dat", ios::app);
  Err << "Error allocating memory in Tarray::set_size: "
      << ba.what() << "\n";
  Err.close();
}

//operators-----------------
template<typename T>
inline void Tarray<T>::operator=(const Tarray<T>& RESTRICT Ta) try{
  if(Ta.m_size != m_size){
    BadSizeTa badt; throw badt; }
  for(int n = 0; n < m_size; ++n){
    for(int i = 0; i < m_size-n; ++i){
      m_data[n][i] = Ta.m_data[n][i];
    }
  }
}catch(exception& e){
  ofstream Err("Error.dat", ios::app);
  Err << e.what();
  Err.close();
}

template<typename T>
inline void Tarray<T>::set_to(T val) {
  for(int n = 0; n < m_size; ++n){
    for(int i = 0; i < m_size-n; ++i){
      m_data[n][i] = val;
    }
  }
}

//print---------------------
template<typename T>
bool Tarray<T>::print() {
  std::ofstream Mat("mat.dat", ios::app);
  if(! Mat){
    Mat.close();
    return false;
  }
  for(int n = 0; n < m_size; ++n){
    for(int i = 0; i < m_size-n; ++i){
      Mat << m_data[n][i] << " ";
    }
    Mat << "\n";
  }
  Mat << endl;
  Mat.close();
  return true;
}

template<typename T>
bool Tarray<T>::print() const{
  std::ofstream Mat("mat.dat", ios::app);
  if(! Mat){
    Mat.close();
    return false;
  }
  for(int n = 0; n < m_size; ++n){
    for(int i = 0; i < m_size-n; ++i){
      Mat << m_data[n][i] << " ";
    }
    Mat << "\n";
  }
  Mat << endl;
  Mat.close();
  return true;
}

//DYNVEC methods*********************************************************

//consructors---------------
template<typename T>
inline Dynvec<T>::Dynvec() try: m_size(1) {
  m_data = new T[m_size];
}catch(bad_alloc& ba){
  ofstream Err("Error.dat", ios::app);
  Err << "Error allocating memory in Dynvec() constructor: "
      << ba.what() << "\n";
  Err.close();
}

template<typename T>
inline Dynvec<T>::Dynvec(int size) try: m_size(size) {
  m_data = new T[m_size];
}catch(bad_alloc& ba){
  ofstream Err("Error.dat", ios::app);
  Err << "Error allocating memory in Dynvec(int) constructor: "
      << ba.what() << "\n";
  Err.close();
}

template<typename T>
inline Dynvec<T>::Dynvec(int size, T val) try: m_size(size) {
  m_data = new T[m_size];
  for(int i = 0; i < m_size; ++i)
    m_data[i] = val;
}catch(bad_alloc& ba){
  ofstream Err("Error.dat", ios::app);
  Err << "Error allocating memory in Dynvec(int, T) constructor: "
      << ba.what() << "\n";
  Err.close();
}
  
template<typename T>
inline Dynvec<T>::Dynvec(const Dynvec<T>& RESTRICT V) try{
  m_size = V.m_size;
  m_data = new T[m_size];
  for(int i = 0; i < m_size; ++i){
    m_data[i] = V.m_data[i];
  }
}catch(bad_alloc& ba){
  ofstream Err("Error.dat", ios::app);
  Err << "Error allocating memory in Dynvec copy constructor: "
      << ba.what() << "\n";
  Err.close();
}

//destructor-----------------
template<typename T>
inline Dynvec<T>::~Dynvec() {
  delete[] m_data;
}

//set size------------------
template<typename T>
inline void Dynvec<T>::set_size(int size) try{
  delete[] m_data;
  m_size = size;
  m_data = new T[m_size];
}catch(bad_alloc& ba){
  ofstream Err("Error.dat", ios::app);
  Err << "Error allocating memory in Dynvec::set_size: "
      << ba.what() << "\n";
  Err.close();
}

//operators-----------------
template<typename T>
inline void Dynvec<T>::operator=(const Dynvec<T>& RESTRICT V) try{
  if(V.m_size != m_size){
    BadSizeV badv; throw badv; }
  for(int i = 0; i < m_size; ++i){
    m_data[i] = V.m_data[i];
  }
}catch(exception& e){
  ofstream Err("Error.dat", ios::app);
  Err << e.what();
  Err.close();
}

template<typename T>
inline void Dynvec<T>::set_to(T val) {
  for(int i = 0; i < m_size; ++i){
    m_data[i] = val;
  }
}


//QUATERNION methods******************************************************

/* //constructors */
/* template<typename T> */
/* inline Quat<T>::Quat() { */
/*   v.set_size(4, 0); */
/*   for(int i = 0; i < 4; ++i) */
/*     v.m_data[i][0] = 0; */
/* } */

/* template<typename T> */
/* inline Quat<T>::Quat(T q1, T q2, T q3, T q4) { */
/*   v.set_size(4, 0); */
/*   v.m_data[0][0] = q1; */
/*   v.m_data[1][0] = q2; */
/*   v.m_data[2][0] = q3; */
/*   v.m_data[3][0] = q4; */
/* } */

/* template<typename T> */
/* inline Quat<T>::Quat(const Quat<T>& Q) { */
/*   v.set_size(4, 0); */
/*   for(int i = 0; i < 4; ++i) */
/*     v.m_data[i][0] = Q.m_data[i][0]; */
/* } */

/* //operators */
/* template<typename T> */
/* inline void Quat<T>::operator=(const Quat<T>& Q) { */
/*   for(int i = 0; i < 4; ++i) */
/*     v.m_data[i][0] = Q.m_data[i][0]; */
/* } */

/* template<typename T> */
/* inline const Quat<T> Quat<T>::operator+(const Quat<T>& Q) { */
/*   Quat<T> result(*this); */
/*   result.v += Q.v; */
/*   return result; */
/* } */

/* template<typename T> */
/* inline const Quat<T> Quat<T>::operator+(const Quat<T>& Q) const{ */
/*   Quat<T> result(*this); */
/*   result.v += Q.v; */
/*   return result; */
/* } */

/* template<typename T> */
/* inline const Quat<T> Quat<T>::operator-(const Quat<T>& Q) { */
/*   Quat<T> result(*this); */
/*   result.v -= Q.v; */
/*   return result; */
/* } */

/* template<typename T> */
/* inline const Quat<T> Quat<T>::operator-(const Quat<T>& Q) const{ */
/*   Quat<T> result(*this); */
/*   result.v -= Q.v; */
/*   return result; */
/* } */

/* template<typename T> */
/* inline const Quat<T> Quat<T>::operator*(const Quat<T>& Q) { */
/*   Quat<T> result; */
/*   //-- */
/*   Matrix<T> M(4, 4); */

/*   M.m_data[0][0] = v.[3][0]; */
/*   M.m_data[0][1] = -v.[2][0]; */
/*   M.m_data[0][2] = v.[1][0]; */
/*   M.m_data[0][3] = v.[0][0]; */

/*   M.m_data[1][0] = v.[2][0]; */
/*   M.m_data[1][1] = v.[3][0]; */
/*   M.m_data[1][2] = -v.[0][0]; */
/*   M.m_data[1][3] = v.[1][0]; */

/*   M.m_data[2][0] = -v.[1][0]; */
/*   M.m_data[2][1] = v.[0][0]; */
/*   M.m_data[2][2] = v.[3][0]; */
/*   M.m_data[2][3] = v.[2][0]; */

/*   M.m_data[3][0] = -v.[0][0]; */
/*   M.m_data[3][1] = -v.[1][0]; */
/*   M.m_data[3][2] = -v.[2][0]; */
/*   M.m_data[3][3] = v.[3][0]; */
/*   //-- */
/*   result.v = M * Q.v; */
/*   return result; */
/* } */

/* template<typename T> */
/* inline const Quat<T> Quat<T>::operator*(const Quat<T>& Q) const{ */
/*   Quat<T> result; */
/*   //-- */
/*   Matrix<T> M(4, 4); */

/*   M.m_data[0][0] = v.[3][0]; */
/*   M.m_data[0][1] = -v.[2][0]; */
/*   M.m_data[0][2] = v.[1][0]; */
/*   M.m_data[0][3] = v.[0][0]; */

/*   M.m_data[1][0] = v.[2][0]; */
/*   M.m_data[1][1] = v.[3][0]; */
/*   M.m_data[1][2] = -v.[0][0]; */
/*   M.m_data[1][3] = v.[1][0]; */

/*   M.m_data[2][0] = -v.[1][0]; */
/*   M.m_data[2][1] = v.[0][0]; */
/*   M.m_data[2][2] = v.[3][0]; */
/*   M.m_data[2][3] = v.[2][0]; */

/*   M.m_data[3][0] = -v.[0][0]; */
/*   M.m_data[3][1] = -v.[1][0]; */
/*   M.m_data[3][2] = -v.[2][0]; */
/*   M.m_data[3][3] = v.[3][0]; */
/*   //-- */
/*   result.v = M * Q.v; */
/*   return result; */
/* } */

/* template<typename T> */
/* inline const Quat<T> Quat<T>::conj() { */
/*   Quat<T> result; */
/*   for(int i = 0; i < 3; ++i) */
/*     result.v.m_data[i][0] = -v.m_data[i][0]; */
/*   result.v.m_data[3][0] = v.m_data[3][0]; */
/*   return result; */
/* } */

/* template<typename T> */
/* inline const Quat<T> Quat<T>::conj() const{ */
/*   Quat<T> result; */
/*   for(int i = 0; i < 3; ++i) */
/*     result.v.m_data[i][0] = -v.m_data[i][0]; */
/*   result.v.m_data[3][0] = v.m_data[3][0]; */
/*   return result; */
/* } */

#endif
