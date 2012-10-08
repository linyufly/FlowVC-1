
# include "multistep.h"


void computePoints_AdamsBashford2(double tmin, double tmax, int N1, double Data_TMin, double Data_TMax, int Num_launch) {

	//////////////////////////////////////////////////////////////////////////////////////
	// 	2nd order Adams-Bashford Method
	// 	For Given Initial Value Problem:
	//	(d/dt)y = f(y,t)
	//	y(n) = constant (For our application)
	//	equation is y(n+2) = y(n+1) + ( 3/2 * h * f( t(n+1), y(n+1) ) ) - ( 1/2 * h * f( t(n), y(n) ) )
	//	y(n) is given.
	//	First y(n+1) is obtained from the Euler method.
	//	Simple Euler Method used is y(n+1) = y(n) + ( h * f( t(n), y(n) ) )
	//////////////////////////////////////////////////////////////////////////////////////
	

	double h; // Integration time from tmin to tmax.

	int n_loops;
	
	float elapsedTime;
	
	// Total number of loops for integrating particles over a single time step.
	n_loops = ( fmod((tmax - tmin) , Int_TimeStep) == 0 ) ? fabs((tmax - tmin)/Int_TimeStep) : (fabs((tmax - tmin)/Int_TimeStep) + 1);
	
	// Calculation of Block size for single chunk of integration kernel ...
	if ( N1 % nthreads == 0 ) {
		nblocks = (int)(N1/nthreads);
	}
	else {
		nblocks = (int)((N1/nthreads) + 1);
	}
	
	double t1, t2, tloc[2], tcurrent;

	int releases, offset;
	
	
	// Time counters which use GPU to time our entire program....
	cudaEvent_t start_int,stop_int;
	cudaEventCreate( &start_int );
	cudaEventCreate( &stop_int );
	cudaEventRecord( start_int, 0 );
	
	
	for (releases = 0; releases < Num_launch; releases++) {
	
		printf("Integrating Release %d from %f to %f \n", releases, tmin, tmax);
		
		offset = releases * N1;
		// Integrate all points for this release
		for (int a = 0; a < n_loops; a++) {

			t1 = fmax (( tmin + a*Int_TimeStep ), Tracer.Start_time[offset]);
			t2 = t1 + Int_TimeStep;
		
			if (t2 > tmax) {
				t2 = tmax;
				h = t2 - t1;
			
			} else {
				h = t2 - t1;
			
			}
		
			tcurrent = t1;
	
			tloc[0] = (tcurrent - Data_TMin) / (Data_TMax - Data_TMin);
		
			tcurrent = t1 + (h/2);
	
			tloc[1] = (tcurrent - Data_TMin) / (Data_TMax - Data_TMin);

		
			// Memcpy of constant memory
			err = cudaMemcpyToSymbol( "tlocation", tloc, 2*sizeof(double), 0, cudaMemcpyHostToDevice);
			if(err != cudaSuccess) {
					fprintf(stderr, "Memory copy of constant variable (tloc) from Host to Device failed\n");
					printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
					fprintf(stderr, " \n\n");
					exit(1);
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
		
			if (Dimensions == 3) { // 3-D case.
		
				if (Data_MeshType == UNSTRUCTURED) {
			
					//AdamsBashford_3DUnstructured_optimized(N1, h);
				
					// This function is faster than the optimized version
					AdamsBashford_3D_Unstructured<<<nblocks, nthreads>>>(Tracer_dev.x, Tracer_dev.y, Tracer_dev.z, velocity_dev, MeshElementArray_device, Tracer_dev.ElementIndex, r, s, t, MeshNodeArray_double_device, N1, integrate, h, offset );
			
				} else { // Cartesian Mesh
			
					AdamsBashford_2Order_3D_Cartesian<<<nblocks, nthreads>>>( Tracer_dev.x, Tracer_dev.y, Tracer_dev.z, velocity_dev, N1,  integrate, Tracer_dev.LeftDomain , h, offset );
		
				}
		
			} else { // 2-D case.
		
				if (Data_MeshType == UNSTRUCTURED) {
			
					AdamsBashford_2D_Unstructured<<<nblocks, nthreads>>>(Tracer_dev.x, Tracer_dev.y, velocity_dev, MeshElementArray_device, Tracer_dev.ElementIndex, r, s, MeshNodeArray_double_device, N1, integrate, h, offset );
			
				} else { // Cartesian Mesh
		
	 				AdamsBashford_2Order_2D_Cartesian<<<nblocks, nthreads>>>( Tracer_dev.x, Tracer_dev.y, velocity_dev, N1,  integrate, Tracer_dev.LeftDomain , h, offset );
		
				}
			}
			err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
			if(err != cudaSuccess) {
				fprintf(stderr, "Error in AdamsBashford kernel. \n");
				fprintf(stderr,"CUDA Error: %s \n\n", cudaGetErrorString(err));
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

/*

void AdamsBashford_3DUnstructured_optimized(int N1, double h) {

	// Do local search .....
	LocalSearch3D<<<nblocks, nthreads>>>(Tracer_dev.x, Tracer_dev.y, Tracer_dev.z, MeshElementArray_device, MeshNodeArray_double_device, Tracer_dev.ElementIndex, N1, 0, r, s, t, integrate);
	err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
	if(err != cudaSuccess) {
			fprintf(stderr, "Error in LocalSearch3D.  \n");
			fprintf(stderr,"CUDA Error: %s \n\n" , cudaGetErrorString(err));
			exit(1);
	}

	// Get Velocity ..
	Adams_Integrate3D<<<nblocks, nthreads>>>(Tracer_dev.x, Tracer_dev.y, Tracer_dev.z, Tracer_dev.ElementIndex, xn0, xn1, xn2, velocity_dev, MeshElementArray_device, r, s, t, N1, integrate, h, 0, 0);
	err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
	if(err != cudaSuccess) {
			fprintf(stderr, "Error in Euler3D.  \n");
			fprintf(stderr,"CUDA Error: %s \n\n" , cudaGetErrorString(err));
			exit(1);
	}
	
	// Do local search .....
	LocalSearch3D<<<nblocks, nthreads>>>(Tracer_dev.x, Tracer_dev.y, Tracer_dev.z, MeshElementArray_device, MeshNodeArray_double_device, Tracer_dev.ElementIndex, N1, 0, r, s, t, integrate);
	err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
	if(err != cudaSuccess) {
			fprintf(stderr, "Error in LocalSearch3D.  \n");
			fprintf(stderr,"CUDA Error: %s \n\n" , cudaGetErrorString(err));
			exit(1);
	}
	
	// Get Velocity ..
	Adams_Integrate3D<<<nblocks, nthreads>>>(Tracer_dev.x, Tracer_dev.y, Tracer_dev.z, Tracer_dev.ElementIndex, xn0, xn1, xn2, velocity_dev, MeshElementArray_device, r, s, t, N1, integrate, h, 1, 0);
	
	

}

__global__ void Adams_Integrate3D(double *x, double *y, double *z, int *eid, double *xn0, double *xn1, double *xn2, VelData_double vel, Element MeshElementArray_device, double *r, double *s, double *t, int ss,  int *integrate, double h, int option, int offset ) {


	// Option = 0 => Euler method
	
	int tid;
	
	// Get thread index
	tid=(blockIdx.x*blockDim.x)+threadIdx.x; 
	if( tid >= ss) { // Redundant thread ...
		return; 
	} else {
	
		tid = tid + offset;
		
		if (integrate[tid] == 1) {
		
			double k[3];
			
			if (option == 0 ) {
			
				// Get velocity
				GetVel_unstruct3D( tlocation[0], k, vel, MeshElementArray_device, eid[tid], r[tid], s[tid], t[tid] );
			
				// This is Euler method with step of h = h/2
				//	x1 = x1 + (h/2)*k[0]; 
				//	y1 = y1 + (h/2)*k[1];
				//	z1 = z1 + (h/2)*k[2];
			
				x[tid] = x[tid] + (0.5*h)*k[0]; 
				y[tid] = y[tid] + (0.5*h)*k[1];
				z[tid] = z[tid] + (0.5*h)*k[2];
				
				xn0[tid] = k[0];
				xn1[tid] = k[1];
				xn2[tid] = k[2]; 
				
			} else if (option == 1) {
			
				// Get velocity
				GetVel_unstruct3D( tlocation[1], k, vel, MeshElementArray_device, eid[tid], r[tid], s[tid], t[tid] );
			
				// Corrector Adambashford 2nd order equation.
			//	x1 = x1 + ( (3/2) * (h/2) * k[0] ) - ( (1/2) * (h/2) * xn0[tid] );
			//	y1 = y1 + ( (3/2) * (h/2) * k[1] ) - ( (1/2) * (h/2) * xn1[tid] );
			//	z1 = z1 + ( (3/2) * (h/2) * k[2] ) - ( (1/2) * (h/2) * xn2[tid] );
				
				x[tid] = x[tid] + ( (1.5) * (0.5*h) * k[0] ) - ( (0.5) * (0.5*h) * xn0[tid] );
				y[tid] = y[tid] + ( (1.5) * (0.5*h) * k[1] ) - ( (0.5) * (0.5*h) * xn1[tid] );
				z[tid] = z[tid] + ( (1.5) * (0.5*h) * k[2] ) - ( (0.5) * (0.5*h) * xn2[tid] );
				
			
			}
			
		
			
			return;
		}
	
	}

} */

__global__ void AdamsBashford_3D_Unstructured(double *x, double *y, double *z, VelData_double vel, Element MeshElementArray_device, int *eid, double *r, double *s, double *t, Node_double MeshNodeArray_double_device, int ss,  int *integrate, double h, int offset ) {

	int tid;
	
	// Get thread index
	tid=(blockIdx.x*blockDim.x)+threadIdx.x; 
	if( tid >= ss) { // Redundant thread ...
		return; 
	} else {
		
		tid = tid + offset;
		
		if (integrate[tid] == 1) {
		
			double k1[3], k2[3];
			double x1,y1,z1;
			int ID = eid[tid];
		
			x1 = x[tid];
			y1 = y[tid];
			z1 = z[tid];
			
			// Do Local Search 
			ID = get_local_search_3D(x1, y1, z1, MeshElementArray_device, MeshNodeArray_double_device, ID, r[tid], s[tid], t[tid]);
			
			// Get velocity
			GetVel_unstruct3D( tlocation[0], k1, vel, MeshElementArray_device, ID, r[tid], s[tid], t[tid] );
		
			// This is Euler method with step of h = h/2
			x1 = x1 + (h*0.5)*k1[0]; 
			y1 = y1 + (h*0.5)*k1[1];
			z1 = z1 + (h*0.5)*k1[2];
		
			// Do Local Search 
			ID = get_local_search_3D(x1, y1, z1, MeshElementArray_device, MeshNodeArray_double_device, ID, r[tid], s[tid], t[tid]);
		
			// Get velocity
			GetVel_unstruct3D( tlocation[1], k2, vel, MeshElementArray_device, ID, r[tid], s[tid], t[tid] );
			
			// Corrector Adambashford 2nd order equation.
			x1 = x1 + ( (1.5) * (h*0.5) * k2[0] ) - ( (0.5) * (h*0.5) * k1[0] );
			y1 = y1 + ( (1.5) * (h*0.5) * k2[1] ) - ( (0.5) * (h*0.5) * k1[1] );
			z1 = z1 + ( (1.5) * (h*0.5) * k2[2] ) - ( (0.5) * (h*0.5) * k1[2] );
			
			x[tid] = x1;
			y[tid] = y1;
			z[tid] = z1;
			eid[tid] = ID;
		
		}
		
		return;

	}

}



__global__ void AdamsBashford_2D_Unstructured(double *x, double *y, VelData_double vel, Element MeshElementArray_device, int *eid, double *r, double *s, Node_double MeshNodeArray_double_device, int ss,  int *integrate, double h, int offset ) {

	int tid;
	
	// Get thread index
	tid=(blockIdx.x*blockDim.x)+threadIdx.x; 
	if( tid >= ss) { // Redundant thread ...
		return; 
	} else {
		
		tid = tid + offset;
		
		if (integrate[tid] == 1) {
		
			double k1[2], k2[2];
			double x1,y1;
			int ID = eid[tid];
		
			x1 = x[tid];
			y1 = y[tid];
			
			// Do Local Search 
			ID = get_local_search_2D(x1, y1,MeshElementArray_device, MeshNodeArray_double_device, ID, r[tid], s[tid]);
			
			// Get velocity
			GetVel_unstruct2D( tlocation[0], k1, vel, MeshElementArray_device, ID, r[tid], s[tid] );
		
			// This is Euler method with step of h = h/2
			x1 = x1 + (h/2)*k1[0]; 
			y1 = y1 + (h/2)*k1[1];
		
			// Do Local Search 
			ID = get_local_search_2D(x1, y1, MeshElementArray_device, MeshNodeArray_double_device, ID, r[tid], s[tid]);
		
			// Get velocity
			GetVel_unstruct2D( tlocation[1], k2, vel, MeshElementArray_device, ID, r[tid], s[tid] );
			
			// Corrector Adambashford 2nd order equation.
			x1 = x1 + ( (3/2) * (h/2) * k2[0] ) - ( (1/2) * (h/2) * k1[0] );
			y1 = y1 + ( (3/2) * (h/2) * k2[1] ) - ( (1/2) * (h/2) * k1[1] );
			
			x[tid] = x1;
			y[tid] = y1;
			eid[tid] = ID;
		
		}
		
		return;

	}

}


__global__ void AdamsBashford_2Order_2D_Cartesian( double *x, double *y, VelData_double vel, int ss,  int *integrate, int *ld , double h, int offset  ) {
	
	//extern __shared__ cache[2*POINTS_BLOCKSIZE_MAIN];
	// Each thread will use their own cache element twice. Order is thread 0 -> cache[0] and cache[1]
	
	// tloc will be in constant memory (Add this feature if doing multiple segments)
	
	// intervals represents how many timesteps will be evaluated by a single kernel.(Add this feature if doing multiple segments)
	
	// Since the first step involved is Euler, we will carryout integration over half a time step for first case. Then we will discard it completely. This is done for the first step run only.
	
	int tid;
	
	// Get thread index
	tid=(blockIdx.x*blockDim.x)+threadIdx.x; 
	if( tid >= ss) { // Redundant thread ...
		return; 
	} else {
	
		tid = tid + offset;
		
		if (integrate[tid] == 1) {
		
			double k1[2], k2[2];
			double x1,y1;
		
			x1 = x[tid];
			y1 = y[tid];
		
			// Get velocity for x and y position at time t1
		
			GetVel_cartesian2D(x1, y1, tlocation[0], k1, vel, ld[tid] );
		
			//cache[2*threadIdx.x] = k[0]; 
			//cache[2*threadIdx.x + 1] = k[1];
		
			// This is Euler method with step of h = h/2
			x1 = x1 + (h/2)*k1[0]; 
			y1 = y1 + (h/2)*k1[1];
		
			// Test if point is outside domain.
			ld[tid] = TestDomain(x1, y1);
		
		
			// Get velocity at new x and y positions
			GetVel_cartesian2D( x1, y1, tlocation[1], k2, vel, ld[tid] );
		
			// Corrector Adambashford 2nd order equation.
			x1 = x1 + ( (3/2) * (h/2) * k2[0] ) - ( (1/2) * (h/2) * k1[0] );
			y1 = y1 + ( (3/2) * (h/2) * k2[1] ) - ( (1/2) * (h/2) * k1[1] );
		
			// Test if point is outside domain.
		
			ld[tid] = TestDomain(x1, y1);
		
			x[tid] = x1;
			y[tid] = y1;
	
		}
	}


}

__global__ void AdamsBashford_2Order_3D_Cartesian( double *x, double *y, double *z, VelData_double vel, int ss,  int *integrate, int *ld , double h, int offset) {
	
	
	// Each thread will use their own cache element twice. Order is thread 0 -> cache[0] and cache[1]
	
	// tloc will be in constant memory (Add this feature if doing multiple segments)
	
	// intervals represents how many timesteps will be evaluated by a single kernel.(Add this feature if doing multiple segments)
	
	// Since the first step involved is Euler, we will carryout integration over half a time step for first case. Then we will discard it completely. This is done for the first step run only.
	
	int tid;
	
	// Get thread index
	tid=(blockIdx.x*blockDim.x)+threadIdx.x; 
	if( tid >= ss) { // Redundant thread ...
		return; 
	} else {
		
		tid = tid + offset;
		
		if (integrate[tid] == 1) {
		
			double k1[3], k2[3];
			double x1,y1,z1;
		
			x1 = x[tid];
			y1 = y[tid];
			z1 = z[tid];
		
			// Get velocity for x, y and z position at time t1
		
			GetVel_cartesian3D(x1, y1, z1, tlocation[0], k1, vel, ld[tid] );
		
			// This is Euler method with step of h = h/2
			x1 = x1 + (h/2)*k1[0]; 
			y1 = y1 + (h/2)*k1[1];
			z1 = z1 + (h/2)*k1[2];
		
			// Test if point is outside domain.
			ld[tid] = TestDomain3D(x1, y1, z1);
		
		
			// Get velocity at new x, y and z positions
			GetVel_cartesian3D( x1, y1, z1, tlocation[1], k2, vel, ld[tid] );
		
			// Corrector Adambashford 2nd order equation.
			x1 = x1 + ( (3/2) * (h/2) * k2[0] ) - ( (1/2) * (h/2) * k1[0] );
			y1 = y1 + ( (3/2) * (h/2) * k2[1] ) - ( (1/2) * (h/2) * k1[1] );
			z1 = z1 + ( (3/2) * (h/2) * k2[2] ) - ( (1/2) * (h/2) * k1[2] );
		
			// Test if point is outside domain.
		
			ld[tid] = TestDomain3D(x1, y1, z1);
		
			x[tid] = x1;
			y[tid] = y1;
			z[tid] = z1;
	
		}
		
		return;
	}


}



__device__ int TestDomain(double x, double y) {

	// This function uses constant memory.
	// Test if point is outside domain.
		if(x < (XX_Data[0] - TINY) || x > (XX_Data[1] + TINY)) {

			return(1);
						
		} else if(y < (YY_Data[0] - TINY) || y > (YY_Data[1] + TINY)) {
			
			return(1);
		
		} else {
		
			return(0);
		
		}


}

__device__ int TestDomain3D(double x, double y, double z) {

	// This function uses constant memory.
	// Test if point is outside domain.
		if(x < (XX_Data[0] - TINY) || x > (XX_Data[1] + TINY)) {

			return(1);
						
		} else if(y < (YY_Data[0] - TINY) || y > (YY_Data[1] + TINY)) {
			
			return(1);
			
		} else if (z < (ZZ_Data[0] - TINY) || z > (ZZ_Data[1] + TINY)) {
		
			return(1);
		
		} else {
		
			return(0);
		
		}


}


