#ifndef _CYLAKS_DEFINITIONS_HPP_
#define _CYLAKS_DEFINITIONS_HPP_
#include <chrono>
#include <cmath>
#include <functional>
#include <map>
#include <unordered_map>
#include <vector>

#define GetVarName(Variable) (#Variable)

/* Stylistic stuff */
using Str = std::string;
using SysClock = std::chrono::steady_clock;
using SysTimepoint = SysClock::time_point;
template <typename DATA_T> using Vec = std::vector<DATA_T>;
template <typename DATA_T> using Vec2D = Vec<Vec<DATA_T>>;
template <typename DATA_T> using Vec3D = Vec<Vec<Vec<DATA_T>>>;
template <typename DATA_T> using Vec4D = Vec<Vec<Vec<Vec<DATA_T>>>>;
template <typename DATA_T> using Vec5D = Vec<Vec<Vec<Vec<Vec<DATA_T>>>>>;
template <typename T1, typename T2> using Map = std::map<T1, T2>;
template <typename T1, typename T2> using UMap = std::unordered_map<T1, T2>;
template <typename T1, typename T2> using Pair = std::pair<T1, T2>;
template <typename T1, typename... ARGS> using Fn = std::function<T1(ARGS...)>;

/* Physical constants */
extern const size_t _n_dims_max{2};
extern const size_t _n_neighbs_max{2};

/* Protein species IDs */
extern const size_t _id_site{0};
extern const size_t _id_motor{1};
extern const size_t _id_xlink{2};

/* Protein size constants; in nm */
extern const double _r_generic{1.0};
extern const double _r_site{8.0};
extern const double _r_motor{4.0};
extern const double _r_xlink{4.0};

#endif