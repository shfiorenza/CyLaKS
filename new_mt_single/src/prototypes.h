// Functions for kMC.cpp
void motors_bind(system_parameters *parameters, microtubule *mt_array, std::vector <motor> &motor_list, std::vector<motor*> &bound_list, std::vector<tubulin*> &unbound_list, gsl_rng *rng);
void motors_unbind(system_parameters *parameters, microtubule *mt_array, std::vector <motor> &motor_list, std::vector<motor*> &bound_list, std::vector<tubulin*> &unbound_list, gsl_rng *rng);
void motors_switch(system_parameters *parameters, microtubule *mt_array, std::vector <motor> &motor_list, std::vector<motor*> &bound_list, std::vector<tubulin*> &unbound_list, gsl_rng *rng);
void motors_move(system_parameters *parameters, microtubule *mt_array, std::vector <motor> &motor_list, std::vector<motor*> &bound_list, std::vector<tubulin*> &unbound_list, gsl_rng *rng);
void motors_boundaries(system_parameters *parameters, microtubule *mt_array, std::vector <motor> &motor_list, std::vector<motor*> &bound_list, gsl_rng *rng, int n_events);

// Functions for output.cpp
void output_data(system_parameters *parameters, microtubule *mt_array, FILE *output_file, FILE *ID_file);
void print_microtubules(system_parameters *parameters, microtubule *mt_array);

// Functions for parse_parameters.cpp
void parse_parameters(char *param_file, system_parameters *parameters);
int parse_tokens(char *line, char ***token);

// "Graceful" file open function for gfopen.cpp
FILE *gfopen(const char *file_name, const char *type);

