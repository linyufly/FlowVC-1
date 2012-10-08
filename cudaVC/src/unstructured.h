#ifndef INC_UNSTRUCTURED_H
#define INC_UNSTRUCTURED_H

// 2D Unstructured Kernels

extern __global__ void initialize_timestep2D(Point Tracer_dev, double *posx, double *posy, double *xn0, double *xn1, int ss, int offset);


extern __global__ void compute_points_Unstructure2D_2(double *posx, double *posy, double *tloc, double *xn0, double *xn1, double *x, double *y, int *eid, double *r, double *s, VelData_double vel, Element MeshElementArray_device, double *h, int *integrate, int ss);


extern __global__ void compute_points_Unstructure2D_1(double *posx, double *posy, double *tloc, double *xn0, double *xn1, double *x, double *y, int *eid, double *r, double *s, VelData_double vel, Element MeshElementArray_device, double *h, int *integrate, int stage, int ss);


// 3D Unstructured Kernels
extern __global__ void initialize_timestep3D(Point Tracer_dev, double *posx, double *posy, double *posz, double *xn0, double *xn1, double *xn2, int ss, int offset);


extern __global__ void compute_points_Unstructure3D_1(double *posx, double *posy, double *posz, double *tloc, double *xn0, double *xn1, double *xn2, double *x, double *y, double *z, int *eid, double *r, double *s, double *t, VelData_double vel, Element MeshElementArray_device, double *h, int *integrate, int stage, int ss);


extern __global__ void compute_points_Unstructure3D_2(double *posx, double *posy, double *posz, double *tloc, double *xn0, double *xn1, double *xn2, double *x, double *y, double *z, int *eid, double *r, double *s, double *t, VelData_double vel, Element MeshElementArray_device, double *h, int *integrate, int ss);


// Kernels Used by both 2D and 3D cases

extern __global__ void check_int_new(Point Tracer_dev, int *integrate, double *t1, double *t2, int ss, int offset);


#endif
