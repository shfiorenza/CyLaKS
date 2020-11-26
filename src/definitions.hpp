#ifndef _CYLAKS_DEFINITIONS_HPP_
#define _CYLAKS_DEFINITIONS_HPP_
#include <chrono>
#include <cmath>
#include <functional>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <map>
#include <unordered_map>
#include <vector>
#include <yaml-cpp/yaml.h>

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
extern size_t _n_dims_max{2};
extern size_t _n_neighbs_max{2};

/* Protein species IDs */
extern size_t _id_tubulin{0};
extern size_t _id_kinesin{1};
extern size_t _id_crosslinker{2};

/* Protein size constants; in nm */
extern double _r_generic{1.0};
extern double _r_tubulin{8.0};
extern double _r_kinesin{4.0};
extern double _r_crosslinker{4.0};

#endif
