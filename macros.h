
#define SQR(x)             ((x) * (x))
#define CUBE(x)            ((x) * (x) * (x))
#define NINT(x)            ((x) < 0.0 ? (int) ((x) - 0.5) : (int) ((x) + 0.5))
#define ABS(x)             ((x) < 0 ? -(x) : (x))
#define MAX(x,y)           ((x) > (y) ? (x) : (y))
#define MIN(x,y)           ((x) < (y) ? (x) : (y))
#define SIGN(x)		   ( (x/fabs(x)) )
