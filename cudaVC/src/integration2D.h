#ifndef INC_INTEGRATION2D_H
#define INC_INTEGRATION2D_H


extern __global__ void compute_points_Cartesian2D(double *posx, double *posy, double tloc, double *xn0, double *xn1, double *x, double *y, VelData_double vel, double stage, double stage2, int output, int ss, int offset, int *integrate, int *ld);

extern void ComputePoints_Unstructured2D(double tmin, double tmax, int N1, double Data_TMin, double Data_TMax, int Num_launch);

extern __global__ void compute_points_Unstructure2D(double *posx, double *posy, double tloc, double *xn0, double *xn1, double *x, double *y, int *eid, double *r, double *s, VelData_double vel, Element MeshElementArray_device, double stage, double stage2, int output, int ss, int offset, int *integrate);

extern void ComputePoints_Cartesian2D(double tmin, double tmax, int N1, double Data_TMin, double Data_TMax, int Num_launch);

#endif
