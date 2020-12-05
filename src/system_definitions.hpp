#ifndef _CYLAKS_SYSTEM_DEFINITIONS_HPP_
#define _CYLAKS_SYSTEM_DEFINITIONS_HPP_
#include <chrono>
#include <cmath>
#include <functional>
#include <map>
#include <unordered_map>
#include <vector>

#define GetVarName(Variable) (#Variable)

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
template <typename T1, typename T2> using Pair = std::pair<T1, T2>;
template <typename T1, typename T2> using UMap = std::unordered_map<T1, T2>;
template <typename T1, typename T2> using Map = std::map<T1, T2>;

struct Sys {

  /* Physical constants */
  static const size_t _n_dims_max{2};
  static const size_t _n_neighbs_max{2};
  /* Protein species IDs */
  static const size_t _id_site{0};
  static const size_t _id_motor{1};
  static const size_t _id_xlink{2};
  /* Protein size constants; in nm */
  static constexpr double _r_site{8.0};
  static constexpr double _r_motor_head{4.0};
  static constexpr double _r_xlink_head{4.0};
};

#endif