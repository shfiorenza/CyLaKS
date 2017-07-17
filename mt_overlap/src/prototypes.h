// Functions for kMC file
void motors_bind(system_parameters *parameters, microtubule *mt_array, std::vector< std::vector<int> > &bound_list, std::vector< std::vector<int> > &unbound_list, int mt_index, gsl_rng *rng);
void motors_unbind(system_parameters *parameters, microtubule *mt_array, std::vector< std::vector<int> > &bound_list, std::vector< std::vector<int> > &unbound_list, int mt_index, gsl_rng *rng);
void motors_switch(system_parameters *parameters, microtubule *mt_array, std::vector< std::vector<int> > &bound_list, std::vector< std::vector<int> > &unbound_list, int mt_index, gsl_rng *rng);
void motors_move(system_parameters *parameters, microtubule *mt_array, std::vector< std::vector<int> > &bound_list, std::vector< std::vector<int> > &unbound_list, int mt_index, gsl_rng *rng);
void motors_boundaries(system_parameters *parameters, microtubule *mt_array, std::vector< std::vector<int> > &bound_list, std::vector< std::vector<int> > &unbound_list, gsl_rng *rng);

// Functions for allocate
void allocate_memory(system_parameters *parameters, microtubule *mt_array);

// Functions for output file
void output_data(system_parameters *parameters, microtubule *mt_array, FILE *output_file);
void print_microtubules(system_parameters *parameters, microtubule *mt_array);

// Functions for parameter file
void parse_parameters(char *param_file, system_parameters *parameters);
int parse_tokens(char *line, char ***token);

// "Graceful" file open function
FILE *gfopen(const char *file_name, const char *type);

// "Graceful" allocation/free functions -- currently not used
void *allocate_1d_array(size_t n, size_t size);
void **allocate_2d_array(size_t n1, size_t n2, size_t size);
void ***allocate_3d_array(size_t n1, size_t n2, size_t n3, size_t size);
void ****allocate_4d_array(size_t n1, size_t n2, size_t n3, size_t n4, size_t size);
void free_1d_array(void *ptr);
void free_2d_array(void **ptr, size_t n1);
void free_3d_array(void ***ptr, size_t n1, size_t n2);
void free_4d_array(void ****ptr, size_t n1, size_t n2, size_t n3);
void *gmalloc(size_t size);
void *gcalloc(size_t n, size_t size);
void *grealloc(void *ptr, size_t size);

//Some weird in-house RNG
double ran3(long *idum);
