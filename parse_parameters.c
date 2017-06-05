/* This routine assigns values to parameters referenced in a parameter file.
   The functional body of this routine is generated automatically from a
   parameter configuration file using the utility routine configure_parameters.

   Input: name of the parameter file (param_file)
          pointer to parameters structure (parameters)

   Output: the values of parameters referenced in the parameter file are
           modified on output */

#include "hskuan.h"

#define LINE_MAX 2048

void parse_parameters(char *param_file, system_parameters *parameters)
{
   FILE *f_param;
   char line[LINE_MAX], line_tmp[LINE_MAX], *param_name, *param_value, **token = NULL;
   int n_lines, n_tokens;

   /* Print message to standard output. */
   fprintf(stdout, "reading parameters from %s:\n\n", param_file);

   /* Open parameter file. */
   f_param = gfopen(param_file, "r");

   /* Read through file. */
   n_lines = 0;
   while (fgets(line, LINE_MAX, f_param) != NULL) {

      /* Increment line counter. */
      ++n_lines;

      /* Copy line into scratch space. */
      strcpy(line_tmp, line);

      /* Parse line into tokens. */
      n_tokens = parse_tokens(line_tmp, &token);

      /* Skip blank lines and comment lines. */
      if (n_tokens > 0 && strncmp(line_tmp, "#", 1) != 0) {

         /* Set value of parameter. */
         if (n_tokens == 3 && strcmp(token[1], "=") == 0) {
            param_name = token[0];
            param_value = token[2];

#include "parse_parameters_body.c"

         }
         else {
            fprintf(stderr, "error parsing parameter file %s on line %d:\n", param_file, n_lines);
            fprintf(stderr, "%s", line);
            exit(1);
         }
      }
   }
   fprintf(stdout, "\n");
   fflush(stdout);

   /* Close parameter file. */
   fclose(f_param);

   return;
}

#undef LINE_MAX
