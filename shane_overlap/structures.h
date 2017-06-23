typedef struct node {

  int coord;

  struct node *next;

} node;

typedef struct site {

  int *occupancy; /*0 for nothing, 1 for static cross-linker, 2 for motor*/

  int *coord; /*coordinate relative to origin*/

} site; 

typedef struct microtubule {

  int *polarity;   /*0 for plus-end on right; 1 for plus-end on left */

  int *coord;      /*coordinate of MT edge to 'origin' */ 

  int *n_bound;  /*total number of motors bound; 1 species for now */

  site *track;

} microtubule;

typedef struct environment {

  int *time;

  microtubule *mt_array;


} environment; 
