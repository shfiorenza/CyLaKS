/* Function prototypes for motors. */

void overlap_growth(system_parameters * parameters, int **microtubule_right, int **microtubule_left);
void motors_move(system_parameters * parameters, int **microtubule_right, int **microtubule_left);
void output(int ****microtubule_final, int time, int n_protofilaments, FILE *output_file, FILE *detail_file);
void parse_parameters(char *param_file, system_parameters *parameters);

int parse_tokens(char *line, char ***token);

FILE *gfopen(const char *file_name, const char *type);

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

double ran3(long *idum);
