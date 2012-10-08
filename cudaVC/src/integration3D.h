#ifndef INC_INTEGRATION3D_H
#define INC_INTEGRATION3D_H


extern __global__ void compute_points_Cartesian3D(double *posx, double *posy, double *posz, double tloc, double *xn0, double *xn1, double *xn2, double *x, double *y, double *z, VelData_double vel, double stage, double stage2, int output, int ss, int offset, int *integrate, int *ld);


extern void ComputePoints_Cartesian3D(double tmin, double tmax, int N1, double Data_TMin, double Data_TMax, int Num_launch);


extern void ComputePoints_Unstructured3D(double tmin, double tmax, int N1, double Data_TMin, double Data_TMax, int Num_launch);


extern __global__ void compute_points_Unstructure3D(double *posx, double *posy, double *posz, double tloc, double *xn0, double *xn1, double *xn2, double *x, double *y, double *z, int *eid, double *r, double *s, double *t, VelData_double vel, Element MeshElementArray_device, int output, int ss, int *integrate, int offset);


extern __global__ void check_int(Point Tracer_dev, int *integrate, double t1, double t2, int ss, int offset);


#endif
