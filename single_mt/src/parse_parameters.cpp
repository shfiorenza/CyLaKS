/* This routine assigns values to parameters referenced in a parameter file.
   The functional body of this routine is generated automatically from a
   parameter configuration file using the utility routine configure_parameters.

   Input: name of the parameter file (param_file)
          pointer to parameters structure (parameters)

   Output: the values of parameters referenced in the parameter file are
           modified on output */
#include <yaml-cpp/yaml.h>
#include "master_header.h"

void parse_parameters(char *param_file, system_parameters *parameters) {
    /* Print message to standard output. */
    fprintf(stdout, "\nReading parameters from %s:\n\n", param_file);

    /* Break parameter file into YAML parser node */
    YAML::Node node = YAML::LoadFile(param_file);

    for(YAML::iterator it=node.begin(); it!=node.end(); ++it) {
        std::string key = it->first.as<std::string>();
        std::string value = it->second.as<std::string>();

        const char* const param_name = key.c_str();
        const char* const param_value = value.c_str();

#include "parse_parameters_body.h"

    }
    
    fprintf(stdout, "\n");
    fflush(stdout);

    return;
}

#undef __LINE_MAX
