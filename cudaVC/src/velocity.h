#ifndef INC_VELOCITY_H
#define INC_VELOCITY_H

# include "structs.h"


// Loading velocity parameters
extern void readveldata(int ss);

// Compute Cartesian velocity
extern __device__ void GetVel_cartesian2D(double posx, double posy, double tloc, double *dXdt, VelData_double  vel, int ld );

extern __device__ void GetVel_cartesian3D(double posx, double posy, double posz, double tloc, double *dXdt, VelData_double  vel, int ld );

extern void GetVelocity_Cartesian(double tq, LagrangianPoint *pt, double *dXdt);


// Test if outside domain (for cartesian grid )

extern __device__ int TestOutDomain_dev(double x, double y);

extern __device__ int TestOutDomain3D_dev(double x, double y, double z);

extern int TestOutsideDomain_host(double point[3]);

extern int TestOutsideCartVelDomain(double X[3]);


// Local Search

extern __device__ inline int get_local_search_2D(double x, double y, Element MeshElementArray_device, Node_double MeshNodeArray_double_device, int eid, double r, double s);

extern __device__ inline int get_local_search_3D(double x, double y, double z, Element MeshElementArray_device, Node_double MeshNodeArray_double_device, int eid, double r, double s, double t);

extern __global__ void LocalSearch2D(double *posx, double *posy, Element MeshElementArray_device, Node_double MeshNodeArray_double_device, int *eid, int ss, int offset, double *r, double *s, int *integrate);

extern __global__ void LocalSearch3D(double *posx, double *posy, double *posz, Element MeshElementArray_device, Node_double MeshNodeArray_double_device, int *eid, int ss, int offset, double *r, double *s, double *t, int *integrate);


// Compute velocity for unstructure grid

extern __device__ inline void getvel_unstruct_2D(double *output, double t, double *dXdt, VelData_double vel, double TMin, double TMax, Element MeshElementArray_device, Node_double MeshNodeArray_double_device, int eid);

extern __device__  inline void GetVel_unstruct3D(double tloc, double *dXdt, VelData_double vel, Element MeshElementArray_device, int eid, double r, double s, double t );

extern __device__  inline void GetVel_unstruct2D(double tloc, double *dXdt, VelData_double vel, Element MeshElementArray_device, int eid, double r, double s );

extern void GetVelocity_Unstructured(const double tq, LagrangianPoint *pt, double *dXdt);




#endif
