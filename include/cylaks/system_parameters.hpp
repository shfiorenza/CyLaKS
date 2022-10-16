#ifndef _CYLAKS_SYSTEM_PARAMETERS_HPP_
#define _CYLAKS_SYSTEM_PARAMETERS_HPP_
#include <vector>

// Params namespace: used to easily store parameters across CyLaKS classes
namespace Params {
inline size_t seed;       // Seed of random number generator
inline double kbT;        // Boltzmann constant multiplied by temperature; pN*nm
inline double eta;        // Viscosity of liquid; (pN/um^2)*s
inline double dt;         // Time elapsed each kmc-bd step; s
inline double t_run;      // Time to run post-equilibration; s
inline double t_equil;    // Minimum time to equilibrate; s
inline double t_snapshot; // Time between each data output; s
inline double dynamic_equil_window; // Set to <=0 to disable dynamic equil; s
inline size_t verbosity;            // Verbosity level; 0 (quiet) to 3 (max)
namespace Filaments {
inline size_t count;          // Number of filaments in simulation
inline size_t n_subfilaments; // For multiple protofilaments in a MT
inline bool periodic_barrel;  // whether or not barrel wraps around completelty
inline double radius;         // Radius of rod (or barrel for MTs); nm
inline double site_size;      // Length of each binding site; nm
inline size_t n_bd_per_kmc;   // BD iterations done per KMC iteration
inline std::vector<size_t> n_sites;  // Length of each filament; n_sites
inline std::vector<size_t> polarity; // 0 (1) sets plus-end to i=0 (n_sites - 1)
inline std::vector<double> x_initial; // Starting x-coord of filament COM; nm
inline std::vector<double> y_initial; // Starting y-coord of filament COM; nm
inline std::vector<double> immobile_until; // Time at which filament can move; s
inline std::vector<double> f_applied;      // Force applied to filament COM; pN
inline std::vector<bool> translation_enabled; // Toggle translation in x,y
inline bool rotation_enabled;      // Toggle rotation within  x-y plane
inline bool wca_potential_enabled; // Toggle WCA potential between filaments

}; // namespace Filaments
namespace Motors {
inline size_t n_runs_to_exit; // Number of post-equil. runs to trigger an exit
inline size_t gaussian_range; // Range of Gaussian interaction; n_sites
inline double gaussian_amp_solo; // Amplitude of interaction for one motor; kBT
inline double gaussian_ceiling_bulk; // Ceiling of amp. for many motors; kBT
inline bool gaussian_stepping_coop;  // Do long-rane FX apply to stepping?
inline double neighb_neighb_energy;  // Short-range interaction magnitude; kBT
inline double t_active;              // Time at which motors/ATP flows in; s
inline double k_on;            // Binding rate of ADP heads to MT; 1/(nM*s)
inline double c_bulk;          // Bulk concentration of motors; nM
inline double c_eff_bind;      // For 2nd ADP head when 1st is bound; nM
inline double k_on_ATP;        // ATP binding rate; 1/(micromolar*s)
inline double c_ATP;           // Bulk concentration of ATP; micromolar
inline double k_hydrolyze;     // Rate that ATP->ADPP occurs in motors; 1/s
inline double k_off_i;         // Unbinding rate for ADPP-bound heads; 1/s
inline double k_off_ii;        // Unbinding rate for ADPP-bound heads; 1/s
inline double applied_force;   // Perpetually applied force on motors; pN
inline double internal_force;  // internal necklinker tension; pN
inline double sigma_ATP;       // ATP-binding distance parameter; nm
inline double sigma_off_i;     // Singly-bound off-rate distance parameter; nm
inline double sigma_off_ii;    // Doubly-bound off-rate distance parameter; nm
inline bool endpausing_active; // Toggles if motor step off plus-ends
inline bool tethers_active;    // Toggles tethers (Kif4A stalks or 'tails')
inline double k_tether;        // Tethering rate; 1/(nM*s)
inline double c_eff_tether;    // Effective concentration of free_teth motors
inline double k_untether;      // Untethering rate when bound; 1/s
inline double r_0;             // Rest length of stalk (or tether); nm
inline double k_spring;        // Spring constant of tether; pN/nm
inline double k_slack;         // '' but when shorter than rest length

}; // namespace Motors
namespace Xlinks {
inline double neighb_neighb_energy; // Energy between two PRC1 neighbors
inline double t_active;             // Time at which crosslinkers flow in; s
inline double k_on;                 // Bulk binding rate;  1/(nM*s)
inline double c_bulk;               // Bulk concentration; nM
inline double c_eff_bind;           // Effective conc. of 2nd head binding; nM
inline double k_off_i;              // Unbinding rate while singly-bound;  1/s
inline double k_off_ii;             // Unbinding rate while doubly-bound;  1/s
inline double d_i;                  // Diffusion coefficient; um^2/s
inline double d_ii;                 // Diffusion coefficient; um^2/s
inline double d_side;   // Diffusion coefficient; um^2/s (side-stepping)
inline double r_0;      // Rest length of coiled-coil domain; nm
inline double k_spring; // Spring constant of CC-domain; pN/nm
inline double theta_0;  // Rotational rest angle
inline double k_rot;    // Rotational spring constant;
};                      // namespace Xlinks
};                      // namespace Params

#endif