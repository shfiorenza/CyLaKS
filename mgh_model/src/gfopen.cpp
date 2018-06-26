/* Graceful fopen function.

Input: name of file to open (file_name)
type of file access (type)

Output: pointer to file (return value) */

#include <iostream>

/* Graceful fopen routine. */
FILE *gfopen(const char *file_name, const char *type)
{
    FILE *ptr;

    if ((ptr = fopen(file_name, type)) == NULL) {
		fprintf(stderr, "Cannot open %s in gfopen\n", file_name);
		exit(1);
    }

    return ptr;
}
