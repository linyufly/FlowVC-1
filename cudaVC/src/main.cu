# include <stdio.h>
# include <math.h>
# include <stdlib.h>
# include <cuda.h>
# include <time.h>

# include "cuPrintf.cu"
# include "settings.h"
# include "error.cu"
# include "structs.h"
# include "globals.cu"
# include "allocate.cu"
# include "fileoutput.cu"



# include "parameters.cu"
# include "memcopy.cu"
# include "mesh.cu"
# include "mymath.cu"
# include "GlobalSearch.cu"
# include "tracers.cu"
# include "velocity.cu"

# include "unstructured.cu"
# include "cartesian.cu"
# include "integration.cu"

# include "integration3D.cu"
# include "integration2D.cu"
# include "multistep.cu"
# include "clean.cu"

// Main function ////



int main(int argc, const char *argv[]) {
		
		float elapsedTime;
		
		// Read in the data file, check the settings and calculated some of the derived parameters....	
		ReadInParameters(argc,argv);
		CheckParameters(); 
		SetDerivedParameters();
		
		
		// Time counters which use GPU to time our entire program....
		cudaEvent_t start, stop, stop_initial;
		cudaEventCreate( &start );
		cudaEventCreate( &stop );
		cudaEventCreate( &stop_initial );
		cudaEventRecord( start, 0 );
		



		// Load Mesh data....
		LoadMeshData();

		// Allocate Host variables (Variables on CPU)
		allocate_host();
		
		// Generate Staggered release
		if (Trace_ReleaseStrategy == 1) {
			
			GenerateStaggeredRelease();
			for (int kk = 0; kk < N_Total; kk++) {
			
				if (Tracer.Start_time[kk] == Output_TStart)
					Tracer.Status[kk] = LAUNCHED;
			
			}
			copyresult2vtk_staggered(0, N_Total, Output_TStart, 0.00, 0.00);
		}

		// Allocating Device variables (These variables are allocated on GPU)
		allocate_device();
		
		// Allocate MISC variables
		alloc_misc();
		
		// Memcpy some constant parameters . These are mesh parameters for unstructured mesh.
		host2device_const();

		
			
		// Initializing our Tracer
		initialize_new();
			
		

		// Calculate Performance
		cudaEventRecord( stop_initial, 0 );
		cudaEventSynchronize( stop_initial );
		cudaEventElapsedTime( &elapsedTime, start, stop_initial );
		printf( "Time for initialization: %3.2f ms\n", elapsedTime );
		cudaEventDestroy( stop_initial ) ;	

		// Start Integrating
		printf("Integrating ..... \n\n");

		int i,j;
		int num_integrated = 0;
		
		
		if (Trace_Compute) {
				
				// Main loop for unstructured mesh
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
					
					// Check how many velocity data file lies between each output
					// Check for all Data times ....
					for(j = 0; j < frames_data; j++) {
					
						
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
						
						
							printf("N = %d points loaded ....\n", N1);
						
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
						
						if (Int_Type == 1) { // RK4 Routine
							
							// Integrate all points from tmin to tmax on GPU using RK4 routine.
							if (Trace_ReleaseStrategy == 1) { // Staggered release where tloc is different for particles.
								
								ComputePoints(tmin, tmax, N1, DataTime1[counter + j].Start_time, DataTime1[counter + j].Stop_time);
							
							} else if (Trace_ReleaseStrategy == 0) { // For Normal Releases, kernels in which value of tloc remains the same.
								
								if (Data_MeshType == UNSTRUCTURED) {
									
									if (Dimensions == 2) {
								
										ComputePoints_Unstructured2D(tmin, tmax, N_Frame, DataTime1[counter + j].Start_time, DataTime1[counter + j].Stop_time, N_Launches);
								
									} else { // 3D integration
							
										ComputePoints_Unstructured3D(tmin, tmax, N_Frame, DataTime1[counter + j].Start_time, DataTime1[counter + j].Stop_time, N_Launches);
										
									}
								} else { // Cartesian mesh 
						
									if (Dimensions == 2) { // 2D cartesian 
								
										ComputePoints_Cartesian2D(tmin, tmax, N_Frame, DataTime1[counter + j].Start_time, DataTime1[counter + j].Stop_time, N_Launches);
								
									} else { // 3D cartesian integration
									
										ComputePoints_Cartesian3D(tmin, tmax, N_Frame, DataTime1[counter + j].Start_time, DataTime1[counter + j].Stop_time, N_Launches);
									}
										
								}
								
							}							
						
						} else if (Int_Type == 3) { // Second order Adams-Bashford method.
						
							if (Trace_ReleaseStrategy == 1) { // Staggered release where tloc is different for particles.
							
								FatalError("Multistep Methods are currently not supported for staggered release. \n\n");
							
							} else if (Trace_ReleaseStrategy == 0) { // Normal release
							
								// Compute point function for normal release in multistep function
								computePoints_AdamsBashford2(tmin, tmax, N_Frame, DataTime1[counter + j].Start_time, DataTime1[counter + j].Stop_time, N_Launches);
							
							}
						
						}
					
						// copy results back into host arrays ......
						device2host(N1);
						
						
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
					
						// Copy results to Output file.
						copyresult2vtk(i + 1, (N_Launches ), Output_time[i + 1]);
						copyresult2bin(i + 1, (N_Launches ), Output_time[i + 1]);
						printf("Output %d released ....\n\n\n", i+1);
					} else {
						
						// Copy results to Output file.
						copyresult2vtk_staggered(i + 1, N_Total, Output_time[i + 1], Data_Loaded_TMin, Data_Loaded_TMax);
						copyresult2bin_staggered(i + 1, N_Total, Output_time[i + 1], Data_Loaded_TMin, Data_Loaded_TMax);
						
						printf("Output %d released ....\n\n\n", i+1);
					
					}

				} // End of main loop for unstructured mesh 


		}
		// Calculate Performance
		cudaEventRecord( stop, 0 );
		cudaEventSynchronize( stop );
		
		cudaEventElapsedTime( &elapsedTime, start, stop );
		printf( "Time to complete entire program: %3.2f ms\n", elapsedTime );
		cudaEventDestroy( start ) ;
		cudaEventDestroy( stop ) ;


		// Cleaning memory space.
		clean_all();

		return 0;

}


