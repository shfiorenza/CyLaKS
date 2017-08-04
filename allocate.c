#include "hskuan.h"

/* Allocate 1-dimensional array. */
void *allocate_1d_array(size_t n, size_t size)
{
   void *ptr;

   ptr = gcalloc(n, size);

   return ptr;
}

/* Allocate 2-dimensional array. */
void **allocate_2d_array(size_t n1, size_t n2, size_t size)
{
   int i;
   void **ptr;

   ptr = gcalloc(n1, sizeof(void*));
   for (i = 0; i < n1; ++i)
      ptr[i] = gcalloc(n2, size);

   return ptr;
}

/* Allocate 3-dimensional array. */
void ***allocate_3d_array(size_t n1, size_t n2, size_t n3, size_t size)
{
   int i, j;
   void ***ptr;

   ptr = gcalloc(n1, sizeof(void**));
   for (i = 0; i < n1; ++i) {
      ptr[i] = gcalloc(n2, sizeof(void*));
      for (j = 0; j < n2; ++j)
         ptr[i][j] = gcalloc(n3, size);
   }

   return ptr;
}

/* Allocate 4-dimensional array. */
void ****allocate_4d_array(size_t n1, size_t n2, size_t n3, size_t n4, size_t size)
{
   int i, j, k;
   void ****ptr;

   ptr = gcalloc(n1, sizeof(void***));
   for (i = 0; i < n1; ++i) {
      ptr[i] = gcalloc(n2, sizeof(void**));
      for (j = 0; j < n2; ++j) {
         ptr[i][j] = gcalloc(n3, sizeof(void*));
         for (k = 0; k < n3; ++k)
            ptr[i][j][k] = gcalloc(n4, size);
      }
   }

   return ptr;
}

/* Free 1-dimensional array. */
void free_1d_array(void *ptr)
{
   free(ptr);

   return;
}

/* Free 2-dimensional array. */
void free_2d_array(void **ptr, size_t n1)
{
   int i;

   for (i = 0; i < n1; ++i)
      free(ptr[i]);;
   free(ptr);

   return;
}

/* Free 3-dimensional array. */
void free_3d_array(void ***ptr, size_t n1, size_t n2)
{
   int i, j;

   for (i = 0; i < n1; ++i) {
      for (j = 0; j < n2; ++j)
         free(ptr[i][j]);
      free(ptr[i]);
   }
   free(ptr);

   return;
}

/* Free 4-dimensional array. */
void free_4d_array(void ****ptr, size_t n1, size_t n2, size_t n3)
{
   int i, j, k;

   for (i = 0; i < n1; ++i) {
      for (j = 0; j < n2; ++j) {
         for (k = 0; k < n3; ++k)
            free(ptr[i][j][k]);
         free(ptr[i][j]);
      }
      free(ptr[i]);
   }
   free(ptr);

   return;
}

/* Graceful malloc routine. */
void *gmalloc(size_t size)
{
   void *ptr;

   if ((ptr = malloc(size)) == NULL)
      error_exit("Cannot allocate memory in gmalloc");

   return ptr;
}

/* Graceful calloc routine. */
void *gcalloc(size_t n, size_t size)
{
   void *ptr;

   if ((ptr = calloc(n, size)) == NULL)
      error_exit("Cannot allocate memory in gcalloc");

   return ptr;
}

/* Graceful realloc routine. */
void *grealloc(void *ptr, size_t size)
{
   if (ptr == NULL)
      ptr = malloc(size);
   else
      ptr = realloc(ptr, size);
   if (ptr == NULL)
      error_exit("Cannot allocate memory in grealloc");

   return ptr;
}
