#ifndef _CYLAKS_PARAMETERS_HPP_
#define _CYLAKS_PARAMETERS_HPP_

struct SysParameters {

  // Simulation stuff
  size_t seed; /* Random number seed */
  double kbT;  /* boltzmann constant * temp, in pN * nm */
  double eta;  /* Viscosity of liquid; in (pN/um^2)*s */
  double dt;
  double t_tot;
  double t_equil;
  double t_snapshot;
  bool dynamic_equil;

  struct FilamentParameters {
    int count;                  // Number of MTs in simulation
    Vec<int> length;            // Length of each MT in sites (tubulin dimers)
    double y_dist;              // Vertical distance between each MT
    double site_size;           // Length of tubulin dimer; nm
    double radius;              // Outer radius of MT barrel; nm
    double elevation;           // Distance above glass slide; nm
    Vec<double> start_coord;    // Coords of each MT at sim start; nm
    Vec<double> immobile_until; // MTs immobilzed until this time; s
    double applied_force;       // Perpetually applied force on each MT (pN)
    bool printout_on;           // Whether or not ASCII printout is enabled
    bool diffusion_on;          // Whether or not MT diffusion is enabled
  } filaments;
  struct MotorParameters {
    size_t n_runs_desired;
    int lattice_coop_range;
    double lattice_coop_Emax_solo;
    double lattice_coop_Emax_bulk;
    double interaction_energy; // absolute value of interaction energy
    double t_active;           // Time at which motors/ATP is flowed in
    double k_on;               // Binding rate of ADP heads to MT; 1/(nM*s)
    double c_bulk;             // Bulk concentration of motors; nM
    double c_eff_bind;         // For 2nd ADP head when 1st is bound; nM
    double k_on_ATP;           // ATP binding rate to mots; 1/(micromolar*s)
    double c_ATP;              // Bulk concentration of ATP; micromolar
    double k_hydrolyze;        // Rate that ATP->ADPP occurs in motors; 1/s
    double k_off_i;            // Unbinding rate for ADPP-bound heads; 1/s
    double k_off_ii;           // Unbinding rate for ADPP-bound heads; 1/s
    double applied_force;      // Perpetually applied force on motors; pN
    double internal_force;     // internal necklinker tensionpN
    double sigma_ATP;          // nm
    double sigma_off_i;        // nm
    double sigma_off_ii;       // nm
    double k_tether;           // Tethering rate; 1/(nM*s)
    double c_eff_tether;       // Effective concentration of free_teth motors
    double k_untether;         // Untethering rate when bound; 1/s
    double r_0;                // Rest length of stalk (or tether); nm
    double k_spring;           // Spring constant of tether; pN/nm
    double k_slack;            // '' but when shorter than rest length
    bool tethers_active;
    bool endpausing_active;
  } motors;
  struct XlinkParameters {
    double k_on;               // Bulk binding rate;  1/(nM*s)
    double c_bulk;             // Bulk concentration; nM
    double c_eff_bind;         // Effective conc. of 2nd head binding; nM
    double k_off_i;            // Unbinding rate while singly-bound;  1/s
    double k_off_ii;           // Unbinding rate while doubly-bound;  1/s
    double r_0;                // Rest length of coiled-coil domain; nm
    double k_spring;           // Spring constant of CC-domain; pN/nm
    double diffu_coeff_i;      // Diffusion coefficient; um^2/s
    double diffu_coeff_ii;     // Diffusion coefficient; um^2/s
    double interaction_energy; // Energy between two PRC1 neighbors
  } xlinks;
};
#endif
