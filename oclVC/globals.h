#ifndef INC_GLOBALS_H
#define INC_GLOBALS_H

#include <CL/cl.h>
#include "structs.h"
#include "settings.h"

extern cl_mem					XXYYZZ; 	// 9
extern cl_mem					RES;		// 3

extern cl_mem					tlocation;	// 4
extern cl_mem					XXYYZZ_DATA;	// 9
extern cl_mem					RES_DATA;	// 3
extern cl_mem					stageconst;	// 2
extern cl_mem					stageconst1; // this will be 1 int

extern cl_platform_id 			platform_id;
extern cl_uint 					num_platforms;
extern cl_device_id 			device_id;
extern cl_context 				context;
extern cl_command_queue 		command_queue;

extern cl_program 				program;
extern cl_program 				program2;

extern cl_kernel 				kernel1;
extern cl_kernel 				kernel2;
extern cl_kernel 				kernel3;
extern cl_kernel 				kernel4;
extern cl_kernel 				kernel5;

extern cl_uint 					num_devices;
extern cl_int 					ret;

extern cl_mem					x_dev;
extern cl_mem					y_dev;
extern cl_mem 					z_dev;
extern cl_mem					Start_time_dev;
extern cl_mem 					Stop_time_dev;
extern cl_mem 					ElementIndex_dev;
extern cl_mem 					LeftDomain_dev;

extern cl_mem					posx;
extern cl_mem					posy;
extern cl_mem					posz;
extern cl_mem					xn0;
extern cl_mem					xn1;
extern cl_mem					xn2;
extern cl_mem					r;
extern cl_mem					s;
extern cl_mem					t;
extern cl_mem					eid;
extern cl_mem					integrate;

// Read Only memory
extern cl_mem					Mesh_Node_x;
extern cl_mem					Mesh_Node_y;
extern cl_mem					Mesh_Node_z;

extern cl_mem					Mesh_Element_Node1;
extern cl_mem					Mesh_Element_Node2;
extern cl_mem					Mesh_Element_Node3;
extern cl_mem					Mesh_Element_Node4;
extern cl_mem					Mesh_Element_Neighborindex1;
extern cl_mem 					Mesh_Element_Neighborindex2;
extern cl_mem					Mesh_Element_Neighborindex3;
extern cl_mem 					Mesh_Element_Neighborindex4;

extern cl_mem					Vel_U0;
extern cl_mem					Vel_U1;
extern cl_mem					Vel_V0;
extern cl_mem					Vel_V1;
extern cl_mem					Vel_W0;
extern cl_mem					Vel_W1;

extern int 						frames_data;

extern int						*index1;



extern int						Searching;
extern int						LocalSearchChecking;

extern int						Generate_Frame;
extern int 						Keep_Tempfile;


extern int						Dimensions;
extern int 						Num_Constant_Memory;
extern int 						Output_Release;
// This is just performance related variable which refers to the constant memory used during the global search process.
extern int 						Memory_Usage_GS;
extern int						Memory_Usage_Tracer;
// This variables stores the number of points loaded onto kernel function for global search purpose.
extern int 						Number_Points_loaded;
extern int 						MAIN_LOOP;
extern int						Device_ID;


extern int						N_Present;

// Globally defined file names 

extern char						Path_Data[LONGSTRING];
extern char 					Data_InFilePrefix[LONGSTRING];
extern char 					Trace_OutFilePrefix[LONGSTRING];
extern char						Path_Output[LONGSTRING];
extern char						Temp_OutFilePrefix[LONGSTRING];
extern char						Trace_InFile[LONGSTRING];


// Globals related to velocity data.

extern int						Data_MeshType;
extern int						Data_SuffixTMin;
extern int						Data_SuffixTDelta;
extern int						Data_TRes;
extern double					Data_TMin;
extern double					Data_TDelta;
extern double 					Data_TMax;
extern int 						Data_TPeriodic;
extern CartMesh					Data_MeshBounds;
//CartMesh						Data_MeshBounds_device;
extern CartMesh					Vel_CartMesh;
extern Launch					*DataTime1;


extern CartMesh					Trace_CartMesh;
extern int						Data_FirstFrame;
extern int						Data_LastFrame;
extern VelData_double			velocity;
//VelData_double				velocity_dev;
extern double					Data_Loaded_TMin;
extern double 					Data_Loaded_TMax;

// Used for staggered release
extern double 					Data_LoadedTMin;
extern double					Data_LoadedTMax;


// Globals related to how much output to generate 

extern double					Output_TStart;
extern double					Output_TEnd; 
extern double					Output_TDelta;
extern int						Output_TRes;


// Globals related to integrator settings

extern int 						Int_Type;
extern double					Int_TimeStep;
extern int						Int_TimeDirection;
extern int						Int_Extrapolate;

// Globals related to tracer computation

extern int 						Trace_Compute;
extern double					Trace_LaunchTimeSpacing;
extern int             			Trace_NumLaunchTimes;
extern double	        		Trace_IntTLength;
extern int             			Trace_AlwaysOutput;
extern int						Trace_GenerateMesh;
extern int 						Release_FrameDelta;
extern int						frames;
extern int 						N_Launches;
extern int 						Trace_ReleaseStrategy;
extern double					Trace_ReleaseTMax;
extern int						Trace_MultipleInFiles;
extern int						Trace_InFileFormat;

extern int 						Trace_NumTracers;

// Temp variables for staggered release.
extern int             			Trace_NumReleaseLocations;
extern ReleaseLocation 			**Trace_ReleaseList;
extern LagrangianPoint 			*Trace_ReleasePoints;

// Variables regarding total number of points

extern int						N_Total; // Total number of points including all release frames.
extern int						N_Frame; // Total number of points in a single frame.
extern int						N_Vel; // Total number of points in a velocity mesh.


// Cuda related variables

//cudaError_t				err;
extern int						nblocks;
extern int						nthreads;
extern int 						nblocks1;


// Tracer Launch time buffer ...

extern double					*Launch_time;
extern double					*Output_time;



// Unstructured Mesh Parameters ......

extern int						Vel_MeshNumElements;
extern double					**Vel_MeshNodeArray;
extern int 						Vel_MeshNumNodes;


// Optimized Velocity Mesh variables .....

extern Element					MeshElementArray;
//Element					MeshElementArray_device;

//extern Node					MeshNodeArray;
//Node					MeshNodeArray_device;

extern Node_double				MeshNodeArray_double;
//Node_double				MeshNodeArray_double_device;


// Optimized Tracer related variables

extern Point					Tracer;
extern Point					Tracer1;



#endif



