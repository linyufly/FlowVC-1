

# include "integration2D.h"

// File with kernel functions for 2D cartesian and unstructured integration in normal release.


__global__ void compute_points_Cartesian2D(double *posx, double *posy, double tloc, double *xn0, double *xn1, double *x, double *y, VelData_double vel, double stage, double stage2, int output, int ss, int offset, int *integrate, int *ld) {

	int tid;
	
	// Get thread index
	tid=(blockIdx.x*blockDim.x)+threadIdx.x; 
	if( tid >= ss) { // Redundant thread ...
		return; 
	} else {
		
		tid = tid + offset;
		
		if (integrate[tid] == 1) {
		
			double k[2];
		
			// Get velocity at current position given by posx, and posy.
			GetVel_cartesian2D(posx[tid], posy[tid], tloc, k, vel, ld[tid] );
			
			// Add k to the xn buffer
			xn0[tid] = xn0[tid] + stage*k[0];
			xn1[tid] = xn1[tid] + stage*k[1];
		
			if (output == 0) {
		
				posx[tid] = x[tid] + stage2*k[0];
				posy[tid] = y[tid] + stage2*k[1];
				
				// Check if Point Left Domain
				ld[tid] = TestOutDomain_dev(posx[tid], posy[tid]);
			
			} else { // For updating position. Done only after calculating k4 
			
				x[tid] = x[tid] + xn0[tid];
				y[tid] = y[tid] + xn1[tid];
				
				// Check if Point Left Domain
				ld[tid] = TestOutDomain_dev(x[tid], y[tid]);
		
			}
		
		}
	
	}

	return;
}




void ComputePoints_Cartesian2D(double tmin, double tmax, int N1, double Data_TMin, double Data_TMax, int Num_launch) {
	
	
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
	int releases;
	int offset;
	// Time counters which use GPU to time our entire program....
	cudaEvent_t start_int,stop_int;
	cudaEventCreate( &start_int );
	cudaEventCreate( &stop_int );
	cudaEventRecord( start_int, 0 );
	
	
	for (releases = 0; releases < Num_launch; releases++) {
	
		printf("Integrating Release %d from %f to %f \n", releases, tmin, tmax);
		
		offset = releases * N1;
		// Integrate all points for this release
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
		
			// Initialize the misc variables used for integration
			cudaPrintfInit(); 
			initialize_Cart2D<<<nblocks, nthreads>>>(Tracer_dev, posx, posy, xn0, xn1, N1, offset);
			cudaPrintfDisplay(stdout, true);
			err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
			if(err != cudaSuccess) {
				fprintf(stderr, " Error in initialize_Cart2D kernel.  \n");
				fprintf(stderr,"CUDA Error: %s \n\n", cudaGetErrorString(err));
				exit(1);
			}
		
			tcurrent = t1;
			tloc = (tcurrent - Data_TMin)/(Data_TMax - Data_TMin);
			stage = (h/6.0);
			stage2 = (0.5 * h);
		
			// compute k1 //
			cudaPrintfInit(); 
			compute_points_Cartesian2D<<<nblocks, nthreads>>>(posx, posy, tloc, xn0, xn1, Tracer_dev.x, Tracer_dev.y, velocity_dev, stage, stage2, 0, N1, offset, integrate, Tracer_dev.LeftDomain);
			cudaPrintfDisplay(stdout, true);
			err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
			if(err != cudaSuccess) {
				fprintf(stderr, "Error in compute_points_Cartesian2D-1.  \n");
				fprintf(stderr,"CUDA Error: %s \n\n" , cudaGetErrorString(err));
				exit(1);
			}
		
			tcurrent = t1 + stage2;
			tloc = (tcurrent - Data_TMin)/(Data_TMax - Data_TMin);
			stage = (h/3.0);
			stage2 = (0.5 * h);
		
			// compute k2//
			cudaPrintfInit(); 
			compute_points_Cartesian2D<<<nblocks, nthreads>>>(posx, posy, tloc, xn0, xn1, Tracer_dev.x, Tracer_dev.y, velocity_dev, stage, stage2, 0, N1, offset, integrate, Tracer_dev.LeftDomain);
			cudaPrintfDisplay(stdout, true);
			err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
			if(err != cudaSuccess) {
				fprintf(stderr, "Error in compute_points_Cartesian2D-2.  \n");
				fprintf(stderr,"CUDA Error: %s \n\n" , cudaGetErrorString(err));
				exit(1);
			}
	
			tcurrent = t1 + stage2;
			tloc = (tcurrent - Data_TMin)/(Data_TMax - Data_TMin);
			stage = (h/3.0);
			stage2 = (h);
		
			// compute k3//
			cudaPrintfInit(); 
			compute_points_Cartesian2D<<<nblocks, nthreads>>>(posx, posy, tloc, xn0, xn1, Tracer_dev.x, Tracer_dev.y, velocity_dev, stage, stage2, 0, N1, offset, integrate, Tracer_dev.LeftDomain);
			cudaPrintfDisplay(stdout, true);
			err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
			if(err != cudaSuccess) {
				fprintf(stderr, "Error in compute_points_Cartesian2D-3.  \n");
				fprintf(stderr,"CUDA Error: %s \n\n" , cudaGetErrorString(err));
				exit(1);
			}
		
			tcurrent = t1 + stage2;
			tloc = (tcurrent - Data_TMin)/(Data_TMax - Data_TMin);
			stage = (h/6.0);
			stage2 = (h);
		
			// Compute k4 and update the position//
			cudaPrintfInit(); 
			compute_points_Cartesian2D<<<nblocks, nthreads>>>(posx, posy, tloc, xn0, xn1, Tracer_dev.x, Tracer_dev.y, velocity_dev, stage, stage2, 1, N1, offset, integrate, Tracer_dev.LeftDomain);
			cudaPrintfDisplay(stdout, true);
			err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
			if(err != cudaSuccess) {
				fprintf(stderr, "Error in compute_points_Cartesian2D-4.  \n");
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



void ComputePoints_Unstructured2D(double tmin, double tmax, int N1	, double Data_TMin, double Data_TMax, int Num_launch) {

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
	int releases;
	int offset;
	double t1,t2, tloc, tcurrent, stage, stage2;
	
	// Time counters which use GPU to time our entire program....
	cudaEvent_t start_int,stop_int;
	cudaEventCreate( &start_int );
	cudaEventCreate( &stop_int );
	cudaEventRecord( start_int, 0 );
	
	for (releases = 0; releases < Num_launch; releases++) { // For all current Releases

		printf("Integrating Release %d from %f to %f \n", releases, tmin, tmax);
		
		// offset for index.
		offset = releases * N1;
	
		// Loop over all points in one releases
		for (a = 0; a < n_loops; a++) {
	
			t1 = fmax(Tracer.Start_time[offset], (tmin + a*Int_TimeStep));
			t2 = t1 + Int_TimeStep;
		
			if (t2 > tmax) {
				t2 = tmax;
				h = tmax - t1;
			} else {
				h = t2 - t1;
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
			
			// Initialize Misc variables used for integration
			cudaPrintfInit(); 
			initialize_timestep2D<<<nblocks, nthreads>>>(Tracer_dev, posx, posy, xn0, xn1, N1, offset);
			cudaPrintfDisplay(stdout, true);
			err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
			if(err != cudaSuccess) {
				fprintf(stderr, "Something went terribly wrong in initialize_timestep2D kernel.  %i\n", nblocks);
				fprintf(stderr,"CUDA Error: %s \n\n", cudaGetErrorString(err));
				exit(1);
			}

			cudaPrintfInit(); 
			// Do local search .....
			LocalSearch2D<<<nblocks, nthreads>>>(posx, posy, MeshElementArray_device, MeshNodeArray_double_device, Tracer_dev.ElementIndex, N1, offset, r, s, integrate);
			cudaPrintfDisplay(stdout, true);
			err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
			if(err != cudaSuccess) {
				fprintf(stderr, "Something went terribly wrong in LocalSearch2D.  %i\n", nblocks);
				fprintf(stderr,"CUDA Error: %s \n\n" , cudaGetErrorString(err));
				exit(1);
			}
			
			tcurrent = t1;
			tloc = (tcurrent - Data_TMin)/ (Data_TMax - Data_TMin);
			stage = (h/6.0);
			stage2 = (0.5 * h);
		
			cudaPrintfInit(); 
			// Compute k1//
			compute_points_Unstructure2D<<<nblocks, nthreads>>>(posx, posy, tloc, xn0, xn1, Tracer_dev.x, Tracer_dev.y, Tracer_dev.ElementIndex, r, s, velocity_dev, MeshElementArray_device, stage, stage2, 0, N1, offset, integrate);
			cudaPrintfDisplay(stdout, true);
			err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
			if(err != cudaSuccess) {
				fprintf(stderr, "Something went terribly wrong in compute_points_Unstructure2D.  %i\n", nblocks);
				fprintf(stderr,"CUDA Error: %s \n\n" , cudaGetErrorString(err));
				exit(1);
			}

			// Do local search .....
			LocalSearch2D<<<nblocks, nthreads>>>(posx, posy, MeshElementArray_device, MeshNodeArray_double_device, Tracer_dev.ElementIndex, N1, offset, r, s, integrate);
		
			err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
			if(err != cudaSuccess) {
				fprintf(stderr, "Something went terribly wrong in LocalSearch2D.  %i\n", nblocks);
				fprintf(stderr,"CUDA Error: %s \n\n" , cudaGetErrorString(err));
				exit(1);
			}
		
			tcurrent = t1 + stage2;
			tloc = (tcurrent - Data_TMin)/ (Data_TMax - Data_TMin);
			stage = (h/3.0);
			stage2 = (0.5 * h);
	
	
			cudaPrintfInit(); 
			// Compute k2//
			compute_points_Unstructure2D<<<nblocks, nthreads>>>(posx, posy, tloc, xn0, xn1, Tracer_dev.x, Tracer_dev.y, Tracer_dev.ElementIndex, r, s, velocity_dev, MeshElementArray_device, stage, stage2, 0, N1, offset, integrate);
			cudaPrintfDisplay(stdout, true);
			err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
			if(err != cudaSuccess) {
				fprintf(stderr, "Something went terribly wrong in compute_points_Unstructure2D.  %i\n", nblocks);
				fprintf(stderr,"CUDA Error: %s \n\n" , cudaGetErrorString(err));
				exit(1);
			}
	
			// Do local search .....
			LocalSearch2D<<<nblocks, nthreads>>>(posx, posy, MeshElementArray_device, MeshNodeArray_double_device, Tracer_dev.ElementIndex, N1, offset, r, s, integrate);

			err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
			if(err != cudaSuccess) {
				fprintf(stderr, "Something went terribly wrong in LocalSearch2D.  %i\n", nblocks);
				fprintf(stderr,"CUDA Error: %s \n\n" , cudaGetErrorString(err));
				exit(1);
			}
		
			tcurrent = t1 + stage2;
			tloc = (tcurrent - Data_TMin)/ (Data_TMax - Data_TMin);
			stage = (h/3.0);
			stage2 = (h);

			cudaPrintfInit(); 
			// Compute k3//
			compute_points_Unstructure2D<<<nblocks, nthreads>>>(posx, posy, tloc, xn0, xn1, Tracer_dev.x, Tracer_dev.y, Tracer_dev.ElementIndex, r, s, velocity_dev, MeshElementArray_device, stage, stage2, 0, N1, offset, integrate);
			cudaPrintfDisplay(stdout, true);
			err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
			if(err != cudaSuccess) {
				fprintf(stderr, "Something went terribly wrong in compute_points_Unstructure2D.  %i\n", nblocks);
				fprintf(stderr,"CUDA Error: %s \n\n" , cudaGetErrorString(err));
				exit(1);
			}
		
			// Do local search .....
			LocalSearch2D<<<nblocks, nthreads>>>(posx, posy, MeshElementArray_device, MeshNodeArray_double_device, Tracer_dev.ElementIndex, N1, offset, r, s, integrate);
		
			err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
			if(err != cudaSuccess) {
				fprintf(stderr, "Something went terribly wrong in LocalSearch2D.  %i\n", nblocks);
				fprintf(stderr,"CUDA Error: %s \n\n" , cudaGetErrorString(err));
				exit(1);
			}
		
			tcurrent = t1 + stage2;
			tloc = (tcurrent - Data_TMin)/ (Data_TMax - Data_TMin);
			stage = (h/6.0);
			stage2 = (h);

			cudaPrintfInit(); 
			// Compute k4 and update the position(Tracer_dev)//
			compute_points_Unstructure2D<<<nblocks, nthreads>>>(posx, posy, tloc, xn0, xn1, Tracer_dev.x, Tracer_dev.y, Tracer_dev.ElementIndex, r, s, velocity_dev, MeshElementArray_device, stage, stage2, 1, N1, offset, integrate);
			cudaPrintfDisplay(stdout, true);
			err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
			if(err != cudaSuccess) {
				fprintf(stderr, "Something went terribly wrong in compute_points_Unstructure2D.  %i\n", nblocks);
				fprintf(stderr,"CUDA Error: %s \n\n" , cudaGetErrorString(err));
				exit(1);
			}
		
			// Do local search .....
			LocalSearch2D<<<nblocks, nthreads>>>(Tracer_dev.x, Tracer_dev.y, MeshElementArray_device, MeshNodeArray_double_device, Tracer_dev.ElementIndex, N1, offset, r, s, integrate);
		
			err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
			if(err != cudaSuccess) {
				fprintf(stderr, "Something went terribly wrong in LocalSearch2D.  %i\n", nblocks);
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


__global__ void compute_points_Unstructure2D(double *posx, double *posy, double tloc, double *xn0, double *xn1, double *x, double *y, int *eid, double *r, double *s, VelData_double vel, Element MeshElementArray_device, double stage, double stage2, int output, int ss, int offset, int *integrate) {

	int tid;
	
	// Get thread index
	tid=(blockIdx.x*blockDim.x)+threadIdx.x; 
	if( tid >= ss) { // Redundant thread ...
		return; 
	} else {
		
		tid = tid + offset;
		
		if (integrate[tid] == 1) {
		
			double k[2];
		
			// Get velocity at current position given by posx, and posy.
			GetVel_unstruct2D(tloc, k, vel, MeshElementArray_device, eid[tid], r[tid], s[tid] );
			
			// Add k to the xn buffer
			xn0[tid] = xn0[tid] + stage*k[0];
			xn1[tid] = xn1[tid] + stage*k[1];
		
			if (output == 0) {
		
				posx[tid] = x[tid] + stage2*k[0];
				posy[tid] = y[tid] + stage2*k[1];
		
			} else { // For updating position . Done only after calculating k4 
			
				x[tid] = x[tid] + xn0[tid];
				y[tid] = y[tid] + xn1[tid];
		
			}
		
		}
	
	}

	return;
}


