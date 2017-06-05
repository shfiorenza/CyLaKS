/* This routine parses a string containing tokens separated by white space (spaces and tabs),
   and returns the number of tokens and an array of strings containing the tokens.

   input: the string to be parsed into tokens (line)

   output: number of tokens (return value)
           pointer to array of tokens (token) */

#include "hskuan.h"

int parse_tokens(char *line, char ***token)
{
   int n_tokens = 0, length;
   char *token_tmp;

   /* Parse string into array of tokens. */
   token_tmp = strtok(line, " \t\n");
   while (token_tmp != NULL) {
      (*token) = grealloc((*token), (n_tokens+1) * sizeof(char*));
      (*token)[n_tokens] = gmalloc((strlen(token_tmp)+1) * sizeof(char));
      strcpy((*token)[n_tokens], token_tmp);
      ++n_tokens;
      token_tmp = strtok(NULL, " \t\n");
   }

   return n_tokens;
}
