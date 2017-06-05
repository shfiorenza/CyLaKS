/* Error handling routine. */

#include "hskuan.h"

void error_exit(char *error_msg)
{

   /* Print error message to standard output and exit. */
   fprintf(stderr, "%s\n", error_msg);

   exit(1);
}
