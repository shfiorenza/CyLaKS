#ifndef _CYLAKS_SYSTEM_PARAMETERS_HPP_
#define _CYLAKS_SYSTEM_PARAMETERS_HPP_
#include <vector>

namespace Params {
inline size_t seed; // Random number seed
inline double kbT;  // boltzmann constant * temp, in pN * nm
inline double eta;  // Viscosity of liquid; in (pN/um^2)*s
inline double dt;
inline double t_run;
inline double t_equil;
inline double t_snapshot;
inline double dynamic_equil_window; // Set to 0 or negative value to disable
inline size_t verbosity;
namespace Filaments {
inline size_t count;                  // Number of filaments to simulate
inline double radius;                 // Radius of rod (or barrel for MTs); nm
inline double site_size;              // Length of filament binding site; nm
inline size_t n_bd_per_kmc;           // BD iterations done per KMC iteration
inline std::vector<size_t> n_sites;   // Length of each filament; n_sites
inline std::vector<size_t> polarity;  // 0 (1): plus-end at i = 0 (n_sites - 1)
inline std::vector<double> x_initial; // Starting x-coord of filament center; nm
inline std::vector<double> y_initial; // Starting y-coord of filament center; nm
inline std::vector<double> immobile_until; // Time at which filament can move; s
inline std::vector<bool> translation_enabled; // translational movement in x,y
inline bool rotation_enabled; // rotational movement within the x-y plane

}; // namespace Filaments
namespace Motors {
inline size_t n_runs_to_exit;
inline size_t gaussian_range;
inline double gaussian_amp_solo;
inline double gaussian_ceiling_bulk;
inline double neighb_neighb_energy; // absolute value of interaction energy
inline double t_active;             // Time at which motors/ATP is flowed in
inline double k_on;                 // Binding rate of ADP heads to MT; 1/(nM*s)
inline double c_bulk;               // Bulk concentration of motors; nM
inline double c_eff_bind;           // For 2nd ADP head when 1st is bound; nM
inline double k_on_ATP;             // ATP binding rate; 1/(micromolar*s)
inline double c_ATP;                // Bulk concentration of ATP; micromolar
inline double k_hydrolyze;          // Rate that ATP->ADPP occurs in motors; 1/s
inline double k_off_i;              // Unbinding rate for ADPP-bound heads; 1/s
inline double k_off_ii;             // Unbinding rate for ADPP-bound heads; 1/s
inline double applied_force;        // Perpetually applied force on motors; pN
inline double internal_force;       // internal necklinker tensionpN
inline double sigma_ATP;            // nm
inline double sigma_off_i;          // nm
inline double sigma_off_ii;         // nm
inline bool endpausing_active;
inline bool tethers_active;
inline double k_tether;     // Tethering rate; 1/(nM*s)
inline double c_eff_tether; // Effective concentration of free_teth motors
inline double k_untether;   // Untethering rate when bound; 1/s
inline double r_0;          // Rest length of stalk (or tether); nm
inline double k_spring;     // Spring constant of tether; pN/nm
inline double k_slack;      // '' but when shorter than rest length

}; // namespace Motors
namespace Xlinks {
inline double neighb_neighb_energy; // Energy between two PRC1 neighbors
inline double t_active;
inline double k_on;       // Bulk binding rate;  1/(nM*s)
inline double c_bulk;     // Bulk concentration; nM
inline double c_eff_bind; // Effective conc. of 2nd head binding; nM
inline double k_off_i;    // Unbinding rate while singly-bound;  1/s
inline double k_off_ii;   // Unbinding rate while doubly-bound;  1/s
inline double d_i;        // Diffusion coefficient; um^2/s
inline double d_ii;       // Diffusion coefficient; um^2/s
inline double r_0;        // Rest length of coiled-coil domain; nm
inline double k_spring;   // Spring constant of CC-domain; pN/nm
inline double theta_0;
inline double k_rot;
}; // namespace Xlinks
}; // namespace Params

#endif