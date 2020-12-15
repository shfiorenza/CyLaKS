#ifndef _CYLAKS_DEFINITIONS_HPP_
#define _CYLAKS_DEFINITIONS_HPP_
#include <cassert>
#include <chrono>
#include <cmath>
#include <functional>
#include <map>
#include <unordered_map>
#include <vector>

#define GetVarName(Variable) (#Variable)

/* Lab coordinate vectors */
inline static const std::vector<double> _x_hat{1.0, 0.0};
inline static const std::vector<double> _y_hat{0.0, 1.0};
/* Physical constants */
inline static const size_t _n_dims_max{2};
inline static const size_t _n_neighbs_max{2};
/* Protein species IDs */
inline static const size_t _id_site{0};
inline static const size_t _id_motor{1};
inline static const size_t _id_xlink{2};
/* Protein size constants; in nm */
inline static const double _r_site{8.0};
inline static const double _r_motor_head{4.0};
inline static const double _r_xlink_head{4.0};

/* Stylistic stuff */
using SysClock = std::chrono::steady_clock;
using SysTimepoint = SysClock::time_point;
using Str = std::string;
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
    break;
  case 1:
    return Dot(a, _y_hat);
    break;
  }
}
// Pseudo cross-product in 2-D (torque is a scalar; in/out of page)
inline Vec<double> Cross(double tq, Vec<double> u) {
  return {-tq * u[1], tq * u[0]};
}
#endif