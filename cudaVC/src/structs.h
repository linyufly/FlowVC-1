#ifndef INC_STRUCTS_H
#define INC_STRUCTS_H

#include <stdio.h>
//#include "macros.h"

typedef struct _VelData {

	float *u0; 	// U velocity related to first data file
	float *u1;
	float *v0; 	//  V velocity values related to first data file
	float *v1;
	float *w0; 	// W velocity values related to first data files
	float *w1;

	//float X[3];	// Position of the correspoding velocity data

			// Time interval of a data file
	float *time0;
	float *time1;

} VelData ;


typedef struct _VelData_double {

	double *u0; 	// U velocity related to first data file
	double *u1;
	double *v0; 	//  V velocity values related to first data file
	double *v1;
	double *w0; 	// W velocity values related to first data files
	double *w1;

	//float X[3];	// Position of the correspoding velocity data

			// Time interval of a data file
	double *time0;
	double *time1;

} VelData_double ;



typedef struct _Launch {

	double Start_time; // Start time for a frame
	double Stop_time; // Stop time of a frame
	int Status;
} Launch ;


typedef struct _Point {

	// Coordinates of a point
	double *x;
	double *y;
	double *z;

	// For unstructured: position in a mesh
	int *ElementIndex;
	
	// Integration time
	double *Start_time;
	double *Stop_time;
	
	int *Status;

	// Flag for cartesian mesh to check if the point is in the domain
	int *LeftDomain;


} Point ;
// Total memory for one element of Point for now: 42 bytes

// Structure defining mesh parameters ....

typedef struct _CartMesh {

	double XMin;
	double XMax;
	double YMin;
	double YMax;
	double ZMin;
	double ZMax;

	int XRes;
	int YRes;
	int ZRes;
	
	double XDelta;
	double YDelta;
	double ZDelta;


} CartMesh;


// For Optimization purpose: structure of array is used ...
typedef struct _Element {

	int *Node1;
	int *Node2;
	int *Node3;
	int *Node4;

	int *Neighborindex1;
	int *Neighborindex2;
	int *Neighborindex3;
	int *Neighborindex4;

} Element;

// Memory stats: Element1 takes 16 bytes of memory for one element....
// Node1: It depends on whether used as float or double. If used as float, memory used is 12 bytes... for double: 24 bytes
// Keep in mind that double uses twice the number of registers .....


typedef struct _Node {

	// Position of Node
	float *x;
	float *y;
	float *z;

} Node;

typedef struct _Node_double {

	// Position of Node
	double *x;
	double *y;
	double *z;

} Node_double;

typedef struct _LagrangianPoint { 
	double X[3];
	double V[3];
	int    ElementIndex;
	int    AuxElementIndex;
	int    LeftDomain;
	double LeftDomainTime;
	double Scalar;
} LagrangianPoint;


typedef struct _ReleaseLocation { 
	
	Launch slide;
	LagrangianPoint pt;
	struct _ReleaseLocation *next;

} ReleaseLocation;



#endif
