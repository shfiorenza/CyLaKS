#ifndef _CYLAKS_SYSTEM_DEFINITIONS_HPP_
#define _CYLAKS_SYSTEM_DEFINITIONS_HPP_
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstring>
#include <functional>
#include <iostream>
#include <map>
#include <unordered_map>
#include <vector>

enum Ligand { NONE, ATP, ADPP, ADP };
inline static const double _max_weight{1e3};

/* Physical constants */
inline static const size_t _n_dims_max{2};
inline static const size_t _n_neighbs_max{2};
/* Lambda values used in Boltzmann factors for various energy-based FX */
inline static const double _lambda_lattice{0.5};
inline static const double _lambda_spring{0.5};
inline static const double _lambda_neighb{1.0};
/* Protein species IDs */
inline static const size_t _id_site{0};
inline static const size_t _id_xlink{1};
inline static const size_t _id_motor{2};
/* Protein size constants; in nm */
inline static const double _r_site{8.0};
inline static const double _r_xlink_head{4.0};
inline static const double _r_motor_head{4.0};
/* Lab frame coordinate vectors */
inline static const std::vector<double> _x_hat{1.0, 0.0};
inline static const std::vector<double> _y_hat{0.0, 1.0};

/* Stylistic stuff */
using Str = std::string;
using SysClock = std::chrono::steady_clock;
using SysTimepoint = SysClock::time_point;
template <typename DATA_T> using Vec = std::vector<DATA_T>;
template <typename DATA_T> using Vec2D = Vec<Vec<DATA_T>>;
template <typename DATA_T> using Vec3D = Vec<Vec<Vec<DATA_T>>>;
template <typename DATA_T> using Vec4D = Vec<Vec<Vec<Vec<DATA_T>>>>;
template <typename DATA_T> using Vec5D = Vec<Vec<Vec<Vec<Vec<DATA_T>>>>>;
template <typename DATA_T> using Fn = std::function<DATA_T>;
template <typename T1, typename T2> using Map = std::map<T1, T2>;
template <typename T1, typename T2> using UMap = std::unordered_map<T1, T2>;
template <typename T1, typename T2> using Pair = std::pair<T1, T2>;

/* Common macros */
inline double Square(double x) { return x * x; }
inline double Cube(double x) { return x * x * x; }
template <typename DATA_T> inline DATA_T Pow(DATA_T x, size_t n) {
  DATA_T result{1.0};
  while (n > 0) {
    result *= x;
    n--;
  }
  return result;
}
template <typename... ARGS> inline double Avg(ARGS... vals) {
  double sum{0.0};
  for (auto const &val : {vals...}) {
    sum += val;
  }
  return sum / sizeof...(vals);
}
inline double Dot(Vec<double> a, Vec<double> b) {
  assert(a.size() == b.size());
  assert(a.size() <= _n_dims_max);
  double dotprod{0.0};
  for (int i_dim{0}; i_dim < a.size(); i_dim++) {
    dotprod += a[i_dim] * b[i_dim];
  }
  return dotprod;
}
inline double Dot(Vec<double> a, int i_dim) {
  switch (i_dim) {
  case 0:
    return Dot(a, _x_hat);
  case 1:
    return Dot(a, _y_hat);
  }
  return NAN;
}
// Pseudo cross-product in 2-D (torque is a scalar; in/out of page)
inline double Cross(Vec<double> a, Vec<double> b) {
  return a[0] * b[1] - a[1] * b[0];
}
inline Vec<double> Cross(double a, Vec<double> b) {
  return {-a * b[1], a * b[0]};
}
inline Vec<double> Cross(Vec<double> a, double b) {
  return {b * a[1], -b * a[0]};
}
inline Vec2D<double> Outer(Vec<double> a, Vec<double> b) {
  Vec2D<double> matrix(a.size(), Vec<double>(b.size(), 0.0));
  for (int i{0}; i < a.size(); i++) {
    for (int j{0}; j < b.size(); j++) {
      matrix[i][j] = a[i] * b[j];
    }
  }
  return matrix;
}
inline Vec2D<double> GetProjectionMatrix(Vec<double> a) { return Outer(a, a); }
inline Vec2D<double> GetOrthonormalBasis(Vec<double> a) {
  return {{a[0], a[1]}, {a[1], -a[0]}};
}
#endif