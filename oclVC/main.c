#include <stdio.h>
#include <stdlib.h>

#include "settings.h"

#include <CL/cl.h>
#pragma OPENCL EXTENSION cl_amd_fp64 : enable
#include <stdarg.h>

#include <time.h>
#include "error.h"

#include "structs.h"

#include "globals.h"
#include "openclfuncs.h"
#include "parameters.h"
#include "mesh.h"
#include "allocate.h"
#include "clean.h"
#include "memcopy.h"
#include "fileoutput.h"
#include "mymath.h"
#include "index.h"
#include "GlobalSearch.h"
#include "tracers.h"
#include "velocity.h"
#include "integration3D.h"




int main(int argc, const char *argv[]) {	
	
	// Read in the data file, check the settings and calculated some of the derived parameters....	
	ReadInParameters(argc,argv);
	CheckParameters(); 
	SetDerivedParameters();
	
	// Initialize OpenCL environment.
	initopencl();
	
	// Load Mesh data....
	LoadMeshData();
	
	// Allocate Host variables (Variables on CPU)
	allocate_host();
	
	// Staggered release ...
	// Generate Staggered release
	if (Trace_ReleaseStrategy == 1) {
			
			GenerateStaggeredRelease();
			for (int kk = 0; kk < N_Total; kk++) {
			
				if (Tracer.Start_time[kk] == Output_TStart)
					Tracer.Status[kk] = LAUNCHED;
			
			}
			copyresult2vtk_staggered(0, N_Total, Output_TStart, 0.00, 0.00);
	}
		
	// Allocate Device
	allocate_device();
	
	// host to device some constant variables
	host2device_const();
	
	// Initializing our Tracer
	initialize_new();
	
	// Start Integrating
	printf("Integrating ..... \n\n");

	int i,j;
	int num_integrated = 0;
	
	
	if (Trace_Compute) {
		
				// Main loop
				for(i = 0; i < frames; i++) {
				
					
					
					MAIN_LOOP = 1;   //This flag will memcopy some parameters needed during integrating function .
					
					int N1;
					
					// Time interval for Output array
					Data_Loaded_TMin = Output_TStart + i*Output_TDelta;
					Data_Loaded_TMax = Data_Loaded_TMin + Output_TDelta;
					
					// Flags for counting intervals into which data will be divided
					int sum = 0;
					int counter;
					
					double tmin, tmax;
		
					for(j = 0; j < frames_data; j++) {
					
				
						// If our output start time is 
						if ( ( (Data_Loaded_TMin >= DataTime1[j].Start_time) && (Data_Loaded_TMin < (DataTime1[j].Stop_time - TINY)) ) || ( (Data_Loaded_TMax > (DataTime1[j].Start_time + TINY)) && (Data_Loaded_TMax <= (DataTime1[j].Stop_time + TINY)) || ( (Data_Loaded_TMin < (DataTime1[j].Start_time)) && (Data_Loaded_TMax > (DataTime1[j].Stop_time + TINY ) ) ) ) ) {
							if (sum == 0) // first data time range .....
								counter = j; // this will be first data file set which will be loaded for integration
							sum++;
						}
						
					
					
					}
					
					
					for (j = 0; j < sum; j++) { // Intermediate Main Loop ....
					
						// assigning tmin first ... tmin represents the starting time of integration for current loaded data set. 
						if (Data_Loaded_TMin > DataTime1[counter + j].Start_time) { 					
							tmin = 	Data_Loaded_TMin;
						} else {				

							tmin = DataTime1[counter + j].Start_time;
						}
						//Assigning value of tmax .... tmax represents the ending time of integration for current loaded data set.
						if (Data_Loaded_TMax > (DataTime1[counter + j].Stop_time + TINY)) {
							tmax = DataTime1[counter + j].Stop_time;
						} else {
							tmax = Data_Loaded_TMax;
						}

						
						if (Trace_ReleaseStrategy == 0) // Normal release
							// Launch Tracers if they were to be launched between time interval spanned by tmin and tmax ....
							N_Launches = Launch_Tracer(tmin,tmax, N_Launches );
					
						
						
						// Read Velocity data from velocity data files....
						readveldata(counter + j);
						
						
						
						
						printf("Data loaded in memory spans from t = %f to %f \n", tmin, tmax );
						
						
	
						if (Trace_ReleaseStrategy == 0) {// Normal release
							// Number of points to be integrated
							N1 = (N_Launches ) * N_Frame; // N1 represents the number of points to be integrated .........
						} else
							N1 = N_Total;
						
						
							printf("N = %d points loaded ....\n\n", N1);
						
						// Number of Threads per block.
						nthreads = POINTS_BLOCKSIZE_MAIN;
					
					
						host2device(N1);
					
						// Calculation of Block size for single chunk of integration kernel ...
						if ( N1 % nthreads == 0 ) {
							nblocks = (int)(N1/nthreads);
						}
						else {
							nblocks = (int)((N1/nthreads) + 1);
						}
						//printf("tmin = %f and tmax = %f \n",tmin,tmax);
						//printf("Number of blocks = %d and number of threads per block = %d \n\n",nblocks,nthreads);
						
						if (Int_Type == 1) {
						
							if (Data_MeshType == UNSTRUCTURED) {
								// Our main integration kernel ......
								//compute_points(tmin, tmax, N1);
								if (Dimensions == 2) {
								
										//ComputePoints(tmin, tmax, N1, DataTime1[counter + j].Start_time, DataTime1[counter + j].Stop_time);
										
								
								} else { // 3D integration
							
									ComputePoints3D_new(tmin, tmax, N1, DataTime1[counter + j].Start_time, DataTime1[counter + j].Stop_time);
								}
							} else { // Cartesian mesh 
						
								if (Dimensions == 2) {
								
									//ComputePoints_Cartesian(tmin, tmax, N1, DataTime1[counter + j].Start_time, DataTime1[counter + j].Stop_time);
								
								} else { // 3D integration
									
									//ComputePoints_Cartesian3D(tmin, tmax, N1, DataTime1[counter + j].Start_time, DataTime1[counter + j].Stop_time);
								}
						
						
						
							}
						
						} else if (Int_Type == 3) {
						
							//computePoints_AdamsBashford2(tmin, tmax, N1, DataTime1[counter + j].Start_time, DataTime1[counter + j].Stop_time);
						
						}
						//printf("Loop %d completed .... \n",b);
						// Its time for copying results back into our host arrays ......
						device2host(N1);
						//printf("N1 = %d \n", N1);
						
						int ii;
						num_integrated = 0;
						if (Trace_ReleaseStrategy == 0) {
							// Crude manner to check for number of points integrated. (This Logic does not give the exact number of points integrated. Also this number will deviate from the number from flowVC file.)
							for(ii = 0; ii < N1; ii++) {
						
							
								if (Data_MeshType == UNSTRUCTURED) {
									// Check element index for unstructured
									if (Tracer.ElementIndex[ii] != -1) {
										num_integrated++;
										//printf("x = %f ,  y = %f and z = %f \n", Tracer.x[ii], Tracer.y[ii], Tracer.z[ii]);
									}
								} else {
									// For catresian check if a point is in the domain.
									if (Tracer.LeftDomain[ii] != 1)
											num_integrated++;
							
								}
							
							}
							printf("N = %d points integrated ....\n\n", num_integrated);
						
						} else {
					
							for(ii = 0; ii < N1; ii++) {
						
							
								if (Data_MeshType == UNSTRUCTURED) {
									// Check element index for unstructured
									if (Tracer.ElementIndex[ii] != -1) {
										if (Tracer.Stop_time[ii] < tmin || Tracer.Start_time[ii] > tmax) {
											continue;
										} else {
											num_integrated++;
										}
										//printf("x = %f ,  y = %f and z = %f \n", Tracer.x[ii], Tracer.y[ii], Tracer.z[ii]);
									}
								} else {
									// For catresian check if a point is in the domain.
									if (Tracer.LeftDomain[ii] != 1)
											num_integrated++;
							
								}
							
							}
							printf("N = %d points integrated ....\n\n", num_integrated);
						
						
						
						}
						
					}
					
					
					
					if (Trace_ReleaseStrategy == 0) {
						copyresult2vtk(i + 1, (N_Launches ), Output_time[i + 1]);
						copyresult2bin(i + 1, (N_Launches ), Output_time[i + 1]);
						printf("Output %d released ....\n\n\n", i+1);
					} else {
					
						copyresult2vtk_staggered(i + 1, N_Total, Output_time[i + 1], Data_Loaded_TMin, Data_Loaded_TMax);
						//copyresult2bin(i + 1, num_integrated, Output_time[i + 1]);
						printf("Output %d released ....\n\n\n", i+1);
					
					
					}

				} // End of main loop


		}
	
	clean_all();
	
	return 0;

}
