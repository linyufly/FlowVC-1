#ifndef INC_MACROS_H
#define INC_MACROS_H

/* Put all macros here */

#define LONGSTRING 200 

//#define TINY 1.0e-10
#define UNLAUNCHED 0
#define LAUNCHED 1
#define COMPLETE 2
#define CARTESIAN 0
#define UNSTRUCTURED 1
#define ANALYTIC 2

#define INITIALIZED 0
#define ASCII_LIST 1
#define VTK_POLYDATA 2
#define VTK_UNSTRUCTURED 3
#define BINARY_LIST 4

#define STAGGERED 1
#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define fmin(a, b) ((a) < (b) ? (a) : (b))
#define fmax(a, b) ((a) > (b) ? (a) : (b))
#define SHORTSTRING 100
#define LONGSTRING  200 

#define MAX_SOURCE_SIZE (0x100000)

#define CL_CHECK(_expr)                                                         \
   do {                                                                         \
     cl_int _err = _expr;                                                       \
     if (_err == CL_SUCCESS)                                                    \
       break;                                                                   \
     fprintf(stderr, "OpenCL Error: '%s' returned %d!\n", #_expr, (int)_err);   \
     abort();                                                                   \
   } while (0)
   

/*
#define CL_CHECK_ERR(_expr)                                                     \
   ({                                                                           \
     cl_int _err = CL_INVALID_VALUE;                                            \
     typeof(_expr) _ret = _expr;                                                \
     if (_err != CL_SUCCESS) {                                                  \
       fprintf(stderr, "OpenCL Error: '%s' returned %d!\n", #_expr, (int)_err); \
       abort();                                                                 \
     }                                                                          \
     _ret;                                                                      \
   })
*/

#endif
