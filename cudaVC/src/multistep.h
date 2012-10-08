#ifndef INC_MULTISTEP_H
#define INC_MULTISTEP_H

extern void computePoints_AdamsBashford2(double tmin, double tmax, int N1, double Data_TMin, double Data_TMax, int Num_launch);


extern __device__ int TestDomain(double x, double y);

extern __device__ int TestDomain3D(double x, double y, double z);



extern __global__ void AdamsBashford_2Order_2D_Cartesian( double *x, double *y, VelData_double vel, int ss,  int *integrate, int *ld , double h, int offset ); 

extern __global__ void AdamsBashford_2Order_3D_Cartesian( double *x, double *y, double *z, VelData_double vel, int ss,  int *integrate, int *ld , double h, int offset );



extern __global__ void AdamsBashford_2D_Unstructured(double *x, double *y, VelData_double vel, Element MeshElementArray_device, int *eid, double *r, double *s, Node_double MeshNodeArray_double_device, int ss,  int *integrate, double h, int offset );

extern __global__ void AdamsBashford_3D_Unstructured(double *x, double *y, double *z, VelData_double vel, Element MeshElementArray_device, int *eid, double *r, double *s, double *t, Node_double MeshNodeArray_double_device, int ss,  int *integrate, double h, int offset );


//extern void AdamsBashford_3DUnstructured_optimized(int N1, double h);


//extern __global__ void Adams_Integrate3D(double *x, double *y, double *z, int *eid, double *xn0, double *xn1, double *xn2, VelData_double vel, Element MeshElementArray_device, double *r, double *s, double *t, int ss,  int *integrate, double h, int option, int offset );


#endif
