

#include "structs.h"
#include "settings.h"

// MISC globals


//double 					*tempx;
//double					*tempy;
double 					*posx;
double					*posy;
double					*posz;
//int						*eid;
double					*r;
double					*s;
double					*t;
double					*xn0;
double 					*xn1;
double 					*xn2;
int 					*integrate;


int 					frames_data;

int						*index1;



int						Searching;
int						LocalSearchChecking;

int						Generate_Frame;
int 					Keep_Tempfile;


int						Dimensions;
int 					Num_Constant_Memory;
int 					Output_Release;
// This is just performance related variable which refers to the constant memory used during the global search process.
int 					Memory_Usage_GS;
int						Memory_Usage_Tracer;
// This variables stores the number of points loaded onto kernel function for global search purpose.
int 					Number_Points_loaded;
int 					MAIN_LOOP = 0;
int						Device_ID;


int						N_Present = 0;

// Globally defined file names 

char					Path_Data[LONGSTRING];
char 					Data_InFilePrefix[LONGSTRING];
char 					Trace_OutFilePrefix[LONGSTRING];
char					Path_Output[LONGSTRING];
char					Temp_OutFilePrefix[LONGSTRING];
char					Trace_InFile[LONGSTRING];


// Globals related to velocity data.

int						Data_MeshType;
int						Data_SuffixTMin;
int						Data_SuffixTDelta;
int						Data_TRes;
double					Data_TMin;
double					Data_TDelta;
double 					Data_TMax;
int 					Data_TPeriodic;
CartMesh				Data_MeshBounds;
CartMesh				Data_MeshBounds_device;
CartMesh				Vel_CartMesh;
Launch					*DataTime1;


CartMesh				Trace_CartMesh;
int						Data_FirstFrame;
int						Data_LastFrame;
VelData_double			velocity;
VelData_double			velocity_dev;
double					Data_Loaded_TMin;
double 					Data_Loaded_TMax;

// Used for staggered release
double 					Data_LoadedTMin;
double					Data_LoadedTMax;


// Globals related to how much output to generate 

double					Output_TStart;
double					Output_TEnd; 
double					Output_TDelta;
int						Output_TRes;


// Globals related to integrator settings

int 					Int_Type;
double					Int_TimeStep;
int						Int_TimeDirection;
int						Int_Extrapolate;

// Globals related to tracer computation

int 					Trace_Compute;
double					Trace_LaunchTimeSpacing;
int             		Trace_NumLaunchTimes;
double	        		Trace_IntTLength;
int             		Trace_AlwaysOutput;
int						Trace_GenerateMesh;
int 					Release_FrameDelta;
int						frames;
int 					N_Launches;
int 					Trace_ReleaseStrategy;
double					Trace_ReleaseTMax;
int						Trace_MultipleInFiles;
int						Trace_InFileFormat;

int 					Trace_NumTracers;

// Temp variables for staggered release.
int             		Trace_NumReleaseLocations = 0;
ReleaseLocation 		**Trace_ReleaseList = NULL;
LagrangianPoint 		*Trace_ReleasePoints = NULL;

// Variables regarding total number of points

int						N_Total; // Total number of points including all release frames.
int						N_Frame; // Total number of points in a single frame.
int						N_Vel; // Total number of points in a velocity mesh.


// Cuda related variables

cudaError_t				err;
int						nblocks;
int						nthreads;
int 					nblocks1;


// Tracer Launch time buffer ...

double					*Launch_time;
double					*Output_time;



// Unstructured Mesh Parameters ......

int						Vel_MeshNumElements;
double					**Vel_MeshNodeArray;
int 					Vel_MeshNumNodes;


// Optimized Velocity Mesh variables .....

Element					MeshElementArray;
Element					MeshElementArray_device;

Node					MeshNodeArray;
Node					MeshNodeArray_device;

Node_double				MeshNodeArray_double;
Node_double				MeshNodeArray_double_device;


// Optimized Tracer related variables

Point					Tracer;
Point					Tracer1;
Point					Tracer_dev;



// temp is the temporary buffer used to store the launch related information.
//Point					temp;



//outputid is just an array which will be used to store elementid's for global search.....
int 					*outputid;
int 					*outputid_dev;


// These three variables are used to store x,y and z position of our tracer which will be copied to constant memory for quick use.
float 					*x_host; 
float					*y_host; 
float 					*z_host;


// Constant memory on GPU 
	__constant__ float	x_dev[CONSTANT_MEMORY];
	__constant__ float 	y_dev[CONSTANT_MEMORY];
	__constant__ float 	z_dev[CONSTANT_MEMORY];

	__constant__ double	stageconst[2];
	__constant__ int stageconst1[1];
	
////////////////////////////////////////////////////////////////////

// 			XX[0] -> XMIN
//			XX[1] -> XMAX
//			XX[2] -> XDELTA
// Same applies to YY and ZZ
//
//
//
//		RES[0] -> XRES
//		RES[1] -> YRES
//		RES[2] -> ZRES
//
/////////////////////////////////////////////////////////////////////

	__constant__ double XX[3];
	__constant__ double YY[3];
	__constant__ double ZZ[3];
	
	__constant__ int RES[3];
	
	__constant__ double XX_Data[3];
	__constant__ double YY_Data[3];
	__constant__ double ZZ_Data[3];
	
	__constant__ int RES_Data[3];	
	
	__constant__ double DataLoaded[3];
	////////////////////////////////////////////
	//
	//	DataLoaded[0] = Data_TMin
	//	DataLoaded[1] = Data_TMax
	//	DataLoaded[2] = 1/(Data_TMax - Data_TMin)
	//
	/////////////////////////////////////////////


	__constant__ double tlocation[4];



