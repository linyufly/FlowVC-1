

# include "integration3D.h"

// File with kernel functions of 3D cartesian and unstructured integration in normal release.


__global__ void compute_points_Cartesian3D(double *posx, double *posy, double *posz, double tloc, double *xn0, double *xn1, double *xn2, double *x, double *y, double *z, VelData_double vel, double stage, double stage2, int output, int ss, int offset, int *integrate, int *ld) {

	int tid;
	
	// Get thread index
	tid=(blockIdx.x*blockDim.x)+threadIdx.x; 
	if( tid >= ss) { // Redundant thread ...
		return; 
	} else {
		
		tid = tid + offset;
		
		// Only integrate this point if it needs to be integrated
		if (integrate[tid] == 1) {
		
			double k[3];
		
			// Get velocity at current position given by posx, and posy.
			GetVel_cartesian3D(posx[tid], posy[tid], posz[tid], tloc, k, vel, ld[tid] );
	
			// Add k to xn buffer
			xn0[tid] = xn0[tid] + stage*k[0];
			xn1[tid] = xn1[tid] + stage*k[1];
			xn2[tid] = xn2[tid] + stage*k[2];
		
			if (output == 0) {
				// Update intermediate position 
				posx[tid] = x[tid] + stage2*k[0];
				posy[tid] = y[tid] + stage2*k[1];
				posz[tid] = z[tid] + stage2*k[2];
				
				// Check if Point Left Domain
				ld[tid] = TestOutDomain3D_dev(posx[tid], posy[tid], posz[tid]);
		
			} else { // // For updating position . Done only after calculating k4 
			
				x[tid] = x[tid] + xn0[tid];
				y[tid] = y[tid] + xn1[tid];
				z[tid] = z[tid] + xn2[tid];
				
				// Check if Point Left Domain
				ld[tid] = TestOutDomain3D_dev(x[tid], y[tid], z[tid]);
		
			}
		
		}
	
	}

	return;
}




void ComputePoints_Cartesian3D(double tmin, double tmax, int N1, double Data_TMin, double Data_TMax, int Num_launch) {
	
	
	double t_int, h; // Integration time from tmin to tmax.
	t_int = tmax - tmin;
	int n_loops;
	
	float elapsedTime;
	
	if ( fmod(t_int , Int_TimeStep) == 0 ) {
		n_loops = fabs(t_int/Int_TimeStep); // Total number of loops for integrating particles over a single time step.
	} else {
	
		n_loops = fabs(t_int/Int_TimeStep) + 1;
	}
	// Calculation of Block size for single chunk of integration kernel ...
	if ( N1 % nthreads == 0 ) {
		nblocks = (int)(N1/nthreads);
	}
	else {
		nblocks = (int)((N1/nthreads) + 1);
	}
	
	int a; // counter variable.
	double t1,t2, tloc, tcurrent, stage, stage2;
	int releases, offset;
	
	// Time counters which use GPU to time our entire program....
	cudaEvent_t start_int,stop_int;
	cudaEventCreate( &start_int );
	cudaEventCreate( &stop_int );
	cudaEventRecord( start_int, 0 );
	
	for (releases = 0; releases < Num_launch; releases++) { // For all releases

		printf("Integrating Release %d from %f to %f \n", releases, tmin, tmax);
		
		offset = releases * N1;
	
		for (a = 0; a < n_loops; a++) {
	
			t1 = fmax(Tracer.Start_time[offset], (tmin + a*Int_TimeStep));
			t2 = t1 + Int_TimeStep;
		
			if (t2 > tmax) {
				t2 = tmax;
				h = tmax - t1;
			} else {
				h = Int_TimeStep;
			}
	
			// Check which point needs to be integrated ...
			cudaPrintfInit(); 
			check_int<<<nblocks, nthreads>>>(Tracer_dev, integrate, t1, t2, N1, offset);
			cudaPrintfDisplay(stdout, true);
			err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
	
			if(err != cudaSuccess) {
				fprintf(stderr, "Error in check_int.  \n");
				fprintf(stderr,"CUDA Error: %s \n\n", cudaGetErrorString(err));
				exit(1);
			}
		
			// Initialize  Misc variables used for integration
			cudaPrintfInit(); 
			initialize_Cart3D<<<nblocks, nthreads>>>(Tracer_dev, posx, posy, posz, xn0, xn1, xn2, N1, offset);
			cudaPrintfDisplay(stdout, true);
			err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
			if(err != cudaSuccess) {
				fprintf(stderr, " Error in initialize_Cart3D.  \n");
				fprintf(stderr,"CUDA Error: %s \n\n", cudaGetErrorString(err));
				exit(1);
			}
		
			tcurrent = t1;
			tloc = (tcurrent - Data_TMin)/(Data_TMax - Data_TMin);
			stage = (h/6.0);
			stage2 = (0.5 * h);
		
			// Compute k1//
			cudaPrintfInit(); 
			compute_points_Cartesian3D<<<nblocks, nthreads>>>(posx, posy, posz, tloc, xn0, xn1,xn2, Tracer_dev.x, Tracer_dev.y, Tracer_dev.z, velocity_dev, stage, stage2, 0, N1, offset, integrate, Tracer_dev.LeftDomain);
			cudaPrintfDisplay(stdout, true);
			err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
			if(err != cudaSuccess) {
				fprintf(stderr, "Error in compute_points_Cartesian3D-1.  \n");
				fprintf(stderr,"CUDA Error: %s \n\n" , cudaGetErrorString(err));
				exit(1);
			}
		
			tcurrent = t1 + stage2;
			tloc = (tcurrent - Data_TMin)/(Data_TMax - Data_TMin);
			stage = (h/3.0);
			stage2 = (0.5 * h);
		
			// Compute k2//
			cudaPrintfInit(); 
			compute_points_Cartesian3D<<<nblocks, nthreads>>>(posx, posy, posz, tloc, xn0, xn1,xn2, Tracer_dev.x, Tracer_dev.y, Tracer_dev.z, velocity_dev, stage, stage2, 0, N1, offset, integrate, Tracer_dev.LeftDomain);
		
			cudaPrintfDisplay(stdout, true);
			err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
			if(err != cudaSuccess) {
				fprintf(stderr, "Error in compute_points_Cartesian3D-2.  \n");
				fprintf(stderr,"CUDA Error: %s \n\n" , cudaGetErrorString(err));
				exit(1);
			}
		
			tcurrent = t1 + stage2;
			tloc = (tcurrent - Data_TMin)/(Data_TMax - Data_TMin);
			stage = (h/3.0);
			stage2 = (h);
		
			// Compute k3//
			cudaPrintfInit(); 
			compute_points_Cartesian3D<<<nblocks, nthreads>>>(posx, posy, posz, tloc, xn0, xn1,xn2, Tracer_dev.x, Tracer_dev.y, Tracer_dev.z, velocity_dev, stage, stage2, 0, N1, offset, integrate, Tracer_dev.LeftDomain);
		
			cudaPrintfDisplay(stdout, true);
			err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
			if(err != cudaSuccess) {
				fprintf(stderr, "Error in compute_points_Cartesian3D-3.  \n");
				fprintf(stderr,"CUDA Error: %s \n\n" , cudaGetErrorString(err));
				exit(1);
			}
		
			tcurrent = t1 + stage2;
			tloc = (tcurrent - Data_TMin)/(Data_TMax - Data_TMin);
			stage = (h/6.0);
			stage2 = (h);
		
			// Compute k4 and Update Tracer_dev//
			cudaPrintfInit(); 
			compute_points_Cartesian3D<<<nblocks, nthreads>>>(posx, posy, posz, tloc, xn0, xn1,xn2, Tracer_dev.x, Tracer_dev.y, Tracer_dev.z, velocity_dev, stage, stage2, 0, N1, offset, integrate, Tracer_dev.LeftDomain);
		
			cudaPrintfDisplay(stdout, true);
			err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
			if(err != cudaSuccess) {
				fprintf(stderr, "Error in compute_points_Cartesian3D-4.  \n");
				fprintf(stderr,"CUDA Error: %s \n\n" , cudaGetErrorString(err));
				exit(1);
			}
		}
		printf("Release %d integrated \n",releases);
	}
	
	cudaEventRecord( stop_int, 0 );
	cudaEventSynchronize( stop_int );
	cudaEventElapsedTime( &elapsedTime, start_int, stop_int );
	printf( "Time for integration: %3.2f ms\n", elapsedTime );
	cudaEventDestroy( stop_int ) ;	
	cudaEventDestroy( start_int ) ;	

}


void ComputePoints_Unstructured3D(double tmin, double tmax, int N1, double Data_TMin, double Data_TMax, int Num_launch) {

	double t_int, h; // Integration time from tmin to tmax.
	t_int = tmax - tmin;
	int n_loops;
	
	float elapsedTime;
	
	// Total number of loops for integrating particles from tmin to tmax.
	n_loops = (fmod(t_int , Int_TimeStep) == 0) ? (fabs(t_int/Int_TimeStep)) : (fabs(t_int/Int_TimeStep) + 1);
	
	// Calculation of Block size for single chunk of integration kernel ...
	nblocks = ( N1 % nthreads == 0 ) ? (int)(N1/nthreads) : (int)((N1/nthreads) + 1);
	
	int a; // counter variable.
	double t1,t2, tloc, tcurrent, stage[2];
	int releases, offset;

	// Time counters which use GPU to time our entire program....
	cudaEvent_t start_int,stop_int;
	cudaEventCreate( &start_int );
	cudaEventCreate( &stop_int );
	cudaEventRecord( start_int, 0 );
	
	for (releases = 0; releases < Num_launch; releases++) {

		printf("Integrating Release %d from %f to %f \n", releases, tmin, tmax);
		
		offset = releases * N1;
	
		for (a = 0; a < n_loops; a++) {

			t1 = fmax(Tracer.Start_time[offset], (tmin + a*Int_TimeStep));
			t2 = t1 + Int_TimeStep;
		
			if (t2 > tmax) {
				t2 = tmax;
				h = tmax - t1;
			} else {
				h = Int_TimeStep;
			}

			cudaPrintfInit(); 
			check_int<<<nblocks, nthreads>>>(Tracer_dev, integrate, t1, t2, N1, offset);
			cudaPrintfDisplay(stdout, true);
			err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
	
			if(err != cudaSuccess) {
				fprintf(stderr, "Error in check_int.  \n");
				fprintf(stderr,"CUDA Error: %s \n\n", cudaGetErrorString(err));
				exit(1);
			}
		
			// Initialize MISC variables used during integration
			cudaPrintfInit(); 
			initialize_timestep3D<<<nblocks, nthreads>>>(Tracer_dev, posx, posy, posz, xn0, xn1, xn2, N1, offset);
			cudaPrintfDisplay(stdout, true);
			err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
	
			if(err != cudaSuccess) {
				fprintf(stderr, "Error in initialize_timestep.  \n");
				fprintf(stderr,"CUDA Error: %s \n\n", cudaGetErrorString(err));
				exit(1);
			}
		
		
			// Do local search .....
			cudaPrintfInit(); 
			LocalSearch3D<<<nblocks, nthreads>>>(posx, posy, posz, MeshElementArray_device, MeshNodeArray_double_device, Tracer_dev.ElementIndex, N1, offset, r, s, t, integrate);
			cudaPrintfDisplay(stdout, true);

			err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
		
			if(err != cudaSuccess) {
				fprintf(stderr, " Error in LocalSearch3D-1.  \n");
				fprintf(stderr,"CUDA Error: %s \n\n" , cudaGetErrorString(err));
				exit(1);
			}
		
			tcurrent = t1;
		
			tloc = (tcurrent - Data_TMin)/ (Data_TMax - Data_TMin);
			stage[0] = (h /6.0);
			stage[1] = (0.5 * h);
		
			// Memcpy of constant memory
			err = cudaMemcpyToSymbol( "stageconst", stage, 2*sizeof(double), 0, cudaMemcpyHostToDevice);
			if(err != cudaSuccess) {
					fprintf(stderr, "Memory copy from Host to Device constant stage failed\n");
					printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
					fprintf(stderr, " \n\n");
					exit(1);
			}
		
			cudaPrintfInit(); 
			// Compute k1//
			compute_points_Unstructure3D<<<nblocks, nthreads>>>(posx, posy, posz, tloc, xn0, xn1, xn2, Tracer_dev.x, Tracer_dev.y, Tracer_dev.z, Tracer_dev.ElementIndex, r, s, t, velocity_dev, MeshElementArray_device, 0, N1, integrate, offset);
			cudaPrintfDisplay(stdout, true);
			err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
			if(err != cudaSuccess) {
				fprintf(stderr, "Error in compute_points_Unstructure3D.  \n");
				fprintf(stderr,"CUDA Error: %s \n\n" , cudaGetErrorString(err));
				exit(1);
			}
		
			// Do local search .....
			LocalSearch3D<<<nblocks, nthreads>>>(posx, posy, posz, MeshElementArray_device, MeshNodeArray_double_device, Tracer_dev.ElementIndex, N1, offset, r, s, t, integrate);
		
			err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
			if(err != cudaSuccess) {
				fprintf(stderr, "Error in LocalSearch3D-2.  %i\n", nblocks);
				fprintf(stderr,"CUDA Error: %s \n\n" , cudaGetErrorString(err));
				exit(1);
			}
		
			tcurrent = t1 + stage[1];
			tloc = (tcurrent - Data_TMin)/ (Data_TMax - Data_TMin);
			stage[0] = (h/3.0);
			stage[1] = (0.5 * h);
	
			// Memcpy of constant memory
			err = cudaMemcpyToSymbol( "stageconst", stage, 2*sizeof(double), 0, cudaMemcpyHostToDevice);
			if(err != cudaSuccess) {
					fprintf(stderr, "Memory copy from Host to Device constant stage failed\n");
					printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
					fprintf(stderr, " \n\n");
					exit(1);
			}
		
		
			cudaPrintfInit(); 
			// Compute k2//
			compute_points_Unstructure3D<<<nblocks, nthreads>>>(posx, posy, posz, tloc, xn0, xn1, xn2, Tracer_dev.x, Tracer_dev.y, Tracer_dev.z, Tracer_dev.ElementIndex, r, s, t, velocity_dev, MeshElementArray_device, 0, N1, integrate, offset);

			cudaPrintfDisplay(stdout, true);
			err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
			if(err != cudaSuccess) {
				fprintf(stderr, "Error in compute_points_Unstructure3D.  \n");
				fprintf(stderr,"CUDA Error: %s \n\n" , cudaGetErrorString(err));
				exit(1);
			}
	
			// Do local search .....
			LocalSearch3D<<<nblocks, nthreads>>>(posx, posy, posz, MeshElementArray_device, MeshNodeArray_double_device, Tracer_dev.ElementIndex, N1, offset, r, s, t, integrate);
		

			err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
			if(err != cudaSuccess) {
				fprintf(stderr, "Error in LocalSearch3D-3.  \n");
				fprintf(stderr,"CUDA Error: %s \n\n" , cudaGetErrorString(err));
				exit(1);
			}
		
			tcurrent = t1 + stage[1];
			tloc = (tcurrent - Data_TMin)/ (Data_TMax - Data_TMin);
			stage[0] = (h/3.0);
			stage[1] = (h);
		
			// Memcpy of constant memory
			err = cudaMemcpyToSymbol( "stageconst", stage, 2*sizeof(double), 0, cudaMemcpyHostToDevice);
			if(err != cudaSuccess) {
					fprintf(stderr, "Memory copy from Host to Device constant stage failed\n");
					printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
					fprintf(stderr, " \n\n");
					exit(1);
			}
		
		
			cudaPrintfInit(); 
			// Compute k3//
			compute_points_Unstructure3D<<<nblocks, nthreads>>>(posx, posy, posz, tloc, xn0, xn1, xn2, Tracer_dev.x, Tracer_dev.y, Tracer_dev.z, Tracer_dev.ElementIndex, r, s, t, velocity_dev, MeshElementArray_device, 0, N1, integrate, offset);
		
			cudaPrintfDisplay(stdout, true);
			err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
			if(err != cudaSuccess) {
				fprintf(stderr, "Error in compute_points_Unstructure3D.  \n");
				fprintf(stderr,"CUDA Error: %s \n\n" , cudaGetErrorString(err));
				exit(1);
			}
		
			// Do local search .....
			LocalSearch3D<<<nblocks, nthreads>>>(posx, posy, posz, MeshElementArray_device, MeshNodeArray_double_device, Tracer_dev.ElementIndex, N1, offset, r, s, t, integrate);
	
		
			err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
			if(err != cudaSuccess) {
				fprintf(stderr, "Error in LocalSearch3D-4.  \n");
				fprintf(stderr,"CUDA Error: %s \n\n" , cudaGetErrorString(err));
				exit(1);
			}
		
			tcurrent = t1 + stage[1];
			tloc = (tcurrent - Data_TMin)/ (Data_TMax - Data_TMin);
			stage[0] = (h/6.0);
			stage[1] = (h);
		
			// Memcpy of constant memory
			err = cudaMemcpyToSymbol( "stageconst", stage, 2*sizeof(double), 0, cudaMemcpyHostToDevice);
			if(err != cudaSuccess) {
					fprintf(stderr, "Memory copy from Host to Device constant stage failed\n");
					printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
					fprintf(stderr, " \n\n");
					exit(1);
			}
		
			cudaPrintfInit(); 
			// Compute k4 and Update the position of Tracer_dev//
			compute_points_Unstructure3D<<<nblocks, nthreads>>>(posx, posy, posz, tloc, xn0, xn1, xn2, Tracer_dev.x, Tracer_dev.y, Tracer_dev.z, Tracer_dev.ElementIndex, r, s, t, velocity_dev, MeshElementArray_device, 1, N1, integrate, offset);
		
			cudaPrintfDisplay(stdout, true);
			err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
			if(err != cudaSuccess) {
				fprintf(stderr, "Error in compute_points_Unstructure3D.  \n");
				fprintf(stderr,"CUDA Error: %s \n\n" , cudaGetErrorString(err));
				exit(1);
			}
		
			// Do local search .....
			LocalSearch3D<<<nblocks, nthreads>>>(Tracer_dev.x, Tracer_dev.y, Tracer_dev.z, MeshElementArray_device, MeshNodeArray_double_device, Tracer_dev.ElementIndex, N1, offset, r, s, t, integrate);
		
		
			err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
			if(err != cudaSuccess) {
				fprintf(stderr, "Error in LocalSearch3D-5.  \n");
				fprintf(stderr,"CUDA Error: %s \n\n" , cudaGetErrorString(err));
				exit(1);
			}
		
		}
		printf("Release %d integrated \n",releases);
	}
	
	cudaEventRecord( stop_int, 0 );
	cudaEventSynchronize( stop_int );
	cudaEventElapsedTime( &elapsedTime, start_int, stop_int );
	printf( "Time for integration: %3.2f ms\n", elapsedTime );
	cudaEventDestroy( stop_int ) ;	
	cudaEventDestroy( start_int ) ;

}

__global__ void check_int(Point Tracer_dev, int *integrate, double t1, double t2, int ss, int offset) {

	int tid;
	int index;
	// Get thread index
	tid=(2*blockIdx.x*blockDim.x)+threadIdx.x; 
	
	if( tid >= ss) { // Redundant thread ...
		return; 
	} else {
	
		index = tid + offset;
		
		// Load global memory into register
		double tstart,tstop;
		double tstart1,tstop1;
		tstart = Tracer_dev.Start_time[index];
		tstop = Tracer_dev.Stop_time[index];
		if ((tid + blockDim.x) < ss) {
			tstart1 = Tracer_dev.Start_time[index + blockDim.x];
			tstop1 = Tracer_dev.Stop_time[index + blockDim.x];
		}
		
		// Check if this point needs to be integrated
		if ( (tstop < t1) || (tstart > t2) ) {
			integrate[index] = 0;
		} else {
			integrate[index] = 1;
		}
		
		if ((tid + blockDim.x) < ss) {
			
			if ( (tstop1 < t1) || (tstart1 > t2) ) {
				integrate[index+ blockDim.x] = 0;
			} else {
				integrate[index+ blockDim.x] = 1;
			}
		}
		
		return;
	}

}



__global__ void compute_points_Unstructure3D(double *posx, double *posy, double *posz, double tloc, double *xn0, double *xn1, double *xn2, double *x, double *y, double *z, int *eid, double *r, double *s, double *t, VelData_double vel, Element MeshElementArray_device, int output, int ss, int *integrate, int offset) {

	int tid;
	
	// Get thread index
	tid=(blockIdx.x*blockDim.x)+threadIdx.x; 
	if( tid >= ss) { // Redundant thread ...
		return; 
	} else {
		tid = tid + offset;
		if (integrate[tid] == 1) {
			double k[3];
		
			// Get velocity at current position given by posx, and posy.
			GetVel_unstruct3D(tloc, k, vel, MeshElementArray_device, eid[tid], r[tid], s[tid], t[tid] );
	
			// Add k to xn buffer
			xn0[tid] = xn0[tid] + stageconst[0]*k[0];
			xn1[tid] = xn1[tid] + stageconst[0]*k[1];
			xn2[tid] = xn2[tid] + stageconst[0]*k[2];
		
			if (output == 0) {
				
				// Update intermediate position 
				posx[tid] = x[tid] + stageconst[1]*k[0];
				posy[tid] = y[tid] + stageconst[1]*k[1];
				posz[tid] = z[tid] + stageconst[1]*k[2];
		
			} else { // // For updating position . Done only after calculating k4 
			
				// Update position
				x[tid] = x[tid] + xn0[tid];
				y[tid] = y[tid] + xn1[tid];
				z[tid] = z[tid] + xn2[tid];
		
			} 
		
		}
	}

	return;
}
