#ifndef INC_CARTESIAN_H
#define INC_CARTESIAN_H



// 2D Cartesian kernels 

extern __global__ void initialize_Cart2D(Point Tracer_dev, double *posx, double *posy, double *xn0, double *xn1, int ss, int offset);


extern __global__ void compute_points_Cartesian2D_1(double *posx, double *posy, double *tloc, double *xn0, double *xn1, double *x, double *y, int *ld, VelData_double vel, double *h, int *integrate, int stage, int ss);


extern __global__ void compute_points_Cartesian2D_2(double *posx, double *posy, double *tloc, double *xn0, double *xn1, double *x, double *y, int *ld, VelData_double vel, double *h, int *integrate, int ss);


// 3D Cartesian  kernels

extern __global__ void initialize_Cart3D(Point Tracer_dev, double *posx, double *posy, double *posz, double *xn0, double *xn1, double *xn2, int ss, int offset);


extern __global__ void compute_points_Cartesian3D_1(double *posx, double *posy, double *posz, double *tloc, double *xn0, double *xn1, double *xn2, double *x, double *y, double *z, int *ld, VelData_double vel, double *h, int *integrate, int stage, int ss);


extern __global__ void compute_points_Cartesian3D_2(double *posx, double *posy, double *posz, double *tloc, double *xn0, double *xn1, double *xn2, double *x, double *y, double *z, int *ld, VelData_double vel, double *h, int *integrate, int ss);


#endif
