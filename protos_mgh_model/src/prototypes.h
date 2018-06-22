// Functions for parse_parameters.cpp
void parse_parameters(char *param_file, system_parameters *parameters);
int parse_tokens(char *line, char ***token);

// "Graceful" file open function for gfopen.cpp
FILE *gfopen(const char *file_name, const char *type);

