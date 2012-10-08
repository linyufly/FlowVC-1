# include "integration.h"

// This Function integrates all points on GPU between two velocity files.
// ComputePoints is a host function that calls different kernels at execution time.
// All kernels integrate on GPU for one time step. This processes is repeated until all points are integrated from time tmin to time tmax. 


void ComputePoints(double tmin, double tmax, int N1, double Data_TMin, double Data_TMax) {

	double t_int; // Integration time from tmin to tmax.
	t_int = tmax - tmin;
	int n_loops;
	
	float elapsedTime;
	
	// Total number of loops for integrating particles over a single time step.
	n_loops = ( fmod(t_int , Int_TimeStep) == 0 ) ? fabs(t_int/Int_TimeStep) : ( fabs(t_int/Int_TimeStep) + 1 ); 

	// Calculation of Block size for single chunk of integration kernel ...	
	nblocks = ( N1 % nthreads == 0 ) ? (int)(N1/nthreads) : (int)((N1/nthreads) + 1);
	
	double *tloc_new, *tloc_device;
	double *t11, *t22;
	double *h_new, *h_device;
	
	// Allocating memory for host tloc_new
	err = cudaHostAlloc( (void**)&tloc_new, N1 * sizeof(double ), cudaHostAllocDefault );
	if(err != cudaSuccess) {
		fprintf(stderr, "Something went terribly wrong in allocating Host Memory for tloc_new .\n");
		printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
		exit(1);
	}
			
	// Allocating memory for host tloc_new
	err = cudaHostAlloc( (void**)&h_new, N1 * sizeof(double ), cudaHostAllocDefault );
	if(err != cudaSuccess) {
		fprintf(stderr, "Something went terribly wrong in allocating Host Memory for tloc_new .\n");
		printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
		exit(1);
	}
	
	// Allocating memory for host t11
	err = cudaHostAlloc( (void**)&t11, N1 * sizeof(double ), cudaHostAllocDefault );
	if(err != cudaSuccess) {
		fprintf(stderr, "Something went terribly wrong in allocating Host Memory for tloc_new .\n");
		printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
		exit(1);
	}
	
	// Allocating memory for host t22
	err = cudaHostAlloc( (void**)&t22, N1 * sizeof(double ), cudaHostAllocDefault );
	if(err != cudaSuccess) {
		fprintf(stderr, "Something went terribly wrong in allocating Host Memory for tloc_new .\n");
		printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
		exit(1);
	}
			
	// Allocating memory for device tloc_device
	err = cudaMalloc((double **)&tloc_device, N1*sizeof(double));
	if(err != cudaSuccess) {
		fprintf(stderr, "Something went terribly wrong in allocating tloc_device.\n");
		printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
		exit(1);
	}
	
	err = cudaMalloc((double **)&h_device, N1*sizeof(double));
	if(err != cudaSuccess) {
		fprintf(stderr, "Something went terribly wrong in allocating h_device.\n");
		printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
		exit(1);
	}
	
	
	int a; // counter variable.

	
	for(int k = 0; k < N1; k++) {
	
		t11[k] = tmin;
	
	}
	
	double *t11_device, *t22_device;
	
	// Allocating memory for device t11_device and t22_device
	err = cudaMalloc((double **)&t11_device, N1*sizeof(double));
	if(err != cudaSuccess) {
		fprintf(stderr, "Something went terribly wrong in allocating t11_device.\n");
		printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
		exit(1);
	}
	
	err = cudaMalloc((double **)&t22_device, N1*sizeof(double));
	if(err != cudaSuccess) {
		fprintf(stderr, "Something went terribly wrong in allocating t22_device.\n");
		printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
		exit(1);
	}

	// Time counters which use GPU to time our entire program....
	cudaEvent_t start_int,stop_int;
	cudaEventCreate( &start_int );
	cudaEventCreate( &stop_int );
	cudaEventRecord( start_int, 0 );
	
	for (a = 0; a < n_loops; a++) {

		for (int k = 0; k < N1; k++) {
		
			t11[k] = fmax(t11[k], Tracer.Start_time[k]);
			
			t22[k] = t11[k] + Int_TimeStep;
			
			// set Timestep of integration for this particular point
			if (t22[k] > tmax) {
				t22[k] = tmax;
				h_new[k] = tmax - t11[k];
			} else 
				h_new[k] = Int_TimeStep;
		
		}
		
		// Transfer Host Memory to GPU Memory. (h_new, t11, t22 variables to device ) 
		err = cudaMemcpy(h_device, h_new , N1 * sizeof(double) , cudaMemcpyHostToDevice);
		if(err != cudaSuccess) {
			fprintf(stderr, "Memory copy from Host structure array h_new to Device structure array h_device failed\n");
			printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
			exit(1);
		}	
		
		err = cudaMemcpy(t11_device, t11 , N1 * sizeof(double) , cudaMemcpyHostToDevice);
		if(err != cudaSuccess) {
			fprintf(stderr, "Memory copy from Host structure array t11 to Device structure array t11_device failed\n");
			printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
			exit(1);
		}	
		err = cudaMemcpy(t22_device, t22 , N1 * sizeof(double) , cudaMemcpyHostToDevice);
		if(err != cudaSuccess) {
			fprintf(stderr, "Memory copy from Host structure array t22 to Device structure array t22_device failed\n");
			printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
			exit(1);
		}	
		
		// Check which thread needs to be integrated between t11 and t22.
		cudaPrintfInit(); 
		check_int_new<<<nblocks, nthreads>>>(Tracer_dev, integrate, t11_device, t22_device, N1, 0);
		cudaPrintfDisplay(stdout, true);
		err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
	
		if(err != cudaSuccess) {
			fprintf(stderr, "Error in check_int.  \n");
			fprintf(stderr,"CUDA Error: %s \n\n", cudaGetErrorString(err));
			exit(1);
		}
		
		// Initialize global variables on GPU
		if (Data_MeshType == UNSTRUCTURED)  { // Unstructured Mesh
		
			if (Dimensions == 3) { // 3D case 
				cudaPrintfInit();
				initialize_timestep3D<<<nblocks, nthreads>>>(Tracer_dev, posx, posy, posz, xn0, xn1, xn2, N1, 0);
				cudaPrintfDisplay(stdout, true);
				err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
	
				if(err != cudaSuccess) {
					fprintf(stderr, "Error in initialize_timestep3D Kernel.  \n");
					fprintf(stderr,"CUDA Error: %s \n\n", cudaGetErrorString(err));
					exit(1);
				}
			
				// Do 3D Local Search
				cudaPrintfInit(); 
				LocalSearch3D<<<nblocks, nthreads>>>(posx, posy, posz, MeshElementArray_device, MeshNodeArray_double_device, Tracer_dev.ElementIndex, N1, 0, r, s, t, integrate);
				cudaPrintfDisplay(stdout, true);

				err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
		
				if(err != cudaSuccess) {
					fprintf(stderr, " Error in LocalSearch3D Kernel. ( Stage - 1 )  \n");
					fprintf(stderr,"CUDA Error: %s \n\n" , cudaGetErrorString(err));
					exit(1);
				}
		
			} else { // 2D case 
		
				cudaPrintfInit();
				initialize_timestep2D<<<nblocks, nthreads>>>(Tracer_dev, posx, posy, xn0, xn1, N1, 0);
				cudaPrintfDisplay(stdout, true);
				err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
	
				if(err != cudaSuccess) {
					fprintf(stderr, "Error in initialize_timestep2D Kernel.  \n");
					fprintf(stderr,"CUDA Error: %s \n\n", cudaGetErrorString(err));
					exit(1);
				}
			
			
				// Do 2D Local Search
				cudaPrintfInit(); 
				LocalSearch2D<<<nblocks, nthreads>>>(posx, posy, MeshElementArray_device, MeshNodeArray_double_device, Tracer_dev.ElementIndex, N1, 0, r, s, integrate);
				cudaPrintfDisplay(stdout, true);

				err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
		
				if(err != cudaSuccess) {
					fprintf(stderr, " Error in LocalSearch2D Kernel. ( Stage - 1 )  \n");
					fprintf(stderr,"CUDA Error: %s \n\n" , cudaGetErrorString(err));
					exit(1);
				}	
		
			}
		
		} else if (Data_MeshType == CARTESIAN) { // Cartesian Mesh 
		
			if (Dimensions == 3) { // 3D case 
				cudaPrintfInit();
				initialize_Cart3D<<<nblocks, nthreads>>>(Tracer_dev, posx, posy, posz, xn0, xn1, xn2, N1, 0);
				cudaPrintfDisplay(stdout, true);
				err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
	
				if(err != cudaSuccess) {
					fprintf(stderr, "Error in initialize_timestep3D Kernel.  \n");
					fprintf(stderr,"CUDA Error: %s \n\n", cudaGetErrorString(err));
					exit(1);
				}
		
			} else { // 2D case 
		
				cudaPrintfInit();
				initialize_Cart2D<<<nblocks, nthreads>>>(Tracer_dev, posx, posy, xn0, xn1, N1, 0);				
				cudaPrintfDisplay(stdout, true);
				
				err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
	
				if(err != cudaSuccess) {
					fprintf(stderr, "Error in initialize_timestep2D kernel.  \n");
					fprintf(stderr,"CUDA Error: %s \n\n", cudaGetErrorString(err));
					exit(1);
				}
					
			}
			
		}

		// Set location of all particles in time.
		for (int k =0; k < N1; k++) {
			t11[k] = fmax(Tracer.Start_time[k], t11[k]);
			tloc_new[k] = (t11[k] - Data_TMin)/ (Data_TMax - Data_TMin);
			
		}
		
		// Transfer tloc_new variable to device.
		err = cudaMemcpy(tloc_device, tloc_new , N1 * sizeof(double) , cudaMemcpyHostToDevice);
		if(err != cudaSuccess) {
			fprintf(stderr, "Memory copy from Host structure array tloc_new to Device structure array tloc_device failed\n");
			printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
			exit(1);
		}	
		

		// Compute k1 //
		if (Data_MeshType == UNSTRUCTURED)  { // Unstructured Mesh
		
			if (Dimensions == 3) { // 3D case 
		
				cudaPrintfInit(); 
			
				compute_points_Unstructure3D_1<<<nblocks, nthreads>>>(posx, posy, posz, tloc_device, xn0, xn1, xn2, Tracer_dev.x, Tracer_dev.y, Tracer_dev.z, Tracer_dev.ElementIndex, r, s, t, velocity_dev, MeshElementArray_device, h_device, integrate, 0, N1);
			
				cudaPrintfDisplay(stdout, true);
			
				err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
				if(err != cudaSuccess) {
					fprintf(stderr, "Error in compute_points_Unstructure3D_1 Kernel ( Stage - 1 )  \n");
					fprintf(stderr,"CUDA Error: %s \n\n" , cudaGetErrorString(err));
					exit(1);
				}
			
				// Do 3D local search .....
				cudaPrintfInit(); 
				LocalSearch3D<<<nblocks, nthreads>>>(posx, posy, posz, MeshElementArray_device, MeshNodeArray_double_device, Tracer_dev.ElementIndex, N1, 0, r, s, t, integrate);
				cudaPrintfDisplay(stdout, true);

				err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
		
				if(err != cudaSuccess) {
					fprintf(stderr, " Error in LocalSearch3D Kernel ( Stage - 2 )  \n");
					fprintf(stderr,"CUDA Error: %s \n\n" , cudaGetErrorString(err));
					exit(1);
				}
			
			
			} else { // 2D case
		
				cudaPrintfInit(); 
			
				compute_points_Unstructure2D_1<<<nblocks, nthreads>>>(posx, posy, tloc_device, xn0, xn1, Tracer_dev.x, Tracer_dev.y, Tracer_dev.ElementIndex, r, s, velocity_dev, MeshElementArray_device, h_device, integrate, 0, N1);
			
				cudaPrintfDisplay(stdout, true);
			
				err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
				if(err != cudaSuccess) {
					fprintf(stderr, "Error in compute_points_Unstructure2D_1 Kernel ( Stage - 1 )  \n");
					fprintf(stderr,"CUDA Error: %s \n\n" , cudaGetErrorString(err));
					exit(1);
				}
		
				// Do 2D local search .....
				cudaPrintfInit(); 
				LocalSearch2D<<<nblocks, nthreads>>>(posx, posy, MeshElementArray_device, MeshNodeArray_double_device, Tracer_dev.ElementIndex, N1, 0, r, s, integrate);
				cudaPrintfDisplay(stdout, true);

				err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
		
				if(err != cudaSuccess) {
					fprintf(stderr, " Error in LocalSearch2D Kernel ( Stage - 2 )  \n");
					fprintf(stderr,"CUDA Error: %s \n\n" , cudaGetErrorString(err));
					exit(1);
				}		
			}

		} else if (Data_MeshType == CARTESIAN) {
		
			if (Dimensions == 3) { // 3D case 
		
				cudaPrintfInit(); 
				
				compute_points_Cartesian3D_1<<<nblocks, nthreads>>>(posx, posy, posz, tloc_device, xn0, xn1, xn2, Tracer_dev.x, Tracer_dev.y, Tracer_dev.z, Tracer_dev.LeftDomain, velocity_dev, h_device, integrate, 0, N1);
			
				cudaPrintfDisplay(stdout, true);
			
				err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
				if(err != cudaSuccess) {
					fprintf(stderr, "Error in compute_points_Cartesian3D_1 Kernel ( Stage - 1 )  \n");
					fprintf(stderr,"CUDA Error: %s \n\n" , cudaGetErrorString(err));
					exit(1);
				}
						
			} else { // 2D case
		
				cudaPrintfInit(); 
				
				compute_points_Cartesian2D_1<<<nblocks, nthreads>>>(posx, posy, tloc_device, xn0, xn1, Tracer_dev.x, Tracer_dev.y, Tracer_dev.LeftDomain, velocity_dev, h_device, integrate, 0, N1);
			
				cudaPrintfDisplay(stdout, true);
			
				err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
				if(err != cudaSuccess) {
					fprintf(stderr, "Error in compute_points_Cartesian2D_1 Kernel ( Stage - 1 )  \n");
					fprintf(stderr,"CUDA Error: %s \n\n" , cudaGetErrorString(err));
					exit(1);
				}	
			}
		
		}
		// Set location of all particles in time.
		for (int k =0; k < N1; k++) {
			t11[k] = t11[k]+ (0.5 * h_new[k]);
			tloc_new[k] = (t11[k] - Data_TMin)/ (Data_TMax - Data_TMin);
					
		}
		
		// Transfer tloc_new variable to device.
		err = cudaMemcpy(tloc_device, tloc_new , N1 * sizeof(double) , cudaMemcpyHostToDevice);
		if(err != cudaSuccess) {
			fprintf(stderr, "Memory copy from Host structure array tloc_new to Device structure array tloc_device failed\n");
			printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
			exit(1);
		}	
		

		
		
		// Compute k2 //
		
		if (Data_MeshType == UNSTRUCTURED)  { // Unstructured Mesh

			if (Dimensions == 3) { // 3D case 
			
				cudaPrintfInit(); 
		
				compute_points_Unstructure3D_1<<<nblocks, nthreads>>>(posx, posy, posz, tloc_device, xn0, xn1, xn2, Tracer_dev.x, Tracer_dev.y, Tracer_dev.z, Tracer_dev.ElementIndex, r, s, t, velocity_dev, MeshElementArray_device, h_device, integrate, 1, N1);

				cudaPrintfDisplay(stdout, true);
				err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
				if(err != cudaSuccess) {
					fprintf(stderr, "Error in compute_points_Unstructure3D_1 Kernel ( Stage - 2 )  \n");
					fprintf(stderr,"CUDA Error: %s \n\n" , cudaGetErrorString(err));
					exit(1);
				}
			
			
				// Do local search .....
				cudaPrintfInit(); 
				LocalSearch3D<<<nblocks, nthreads>>>(posx, posy, posz, MeshElementArray_device, MeshNodeArray_double_device, Tracer_dev.ElementIndex, N1, 0, r, s, t, integrate);
				cudaPrintfDisplay(stdout, true);

				err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
		
				if(err != cudaSuccess) {
					fprintf(stderr, " Error in LocalSearch3D Kernel ( Stage - 3 )  \n");
					fprintf(stderr,"CUDA Error: %s \n\n" , cudaGetErrorString(err));
					exit(1);
				}
			
			} else { // 2D case
		
				cudaPrintfInit(); 
		
				compute_points_Unstructure2D_1<<<nblocks, nthreads>>>(posx, posy, tloc_device, xn0, xn1, Tracer_dev.x, Tracer_dev.y, Tracer_dev.ElementIndex, r, s, velocity_dev, MeshElementArray_device, h_device, integrate, 1, N1);

				cudaPrintfDisplay(stdout, true);
				err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
				if(err != cudaSuccess) {
					fprintf(stderr, "Error in compute_points_Unstructure2D_1 Kernel ( Stage - 2 )  \n");
					fprintf(stderr,"CUDA Error: %s \n\n" , cudaGetErrorString(err));
					exit(1);
				}
			
				// Do local search .....
				cudaPrintfInit(); 
				LocalSearch2D<<<nblocks, nthreads>>>(posx, posy, MeshElementArray_device, MeshNodeArray_double_device, Tracer_dev.ElementIndex, N1, 0, r, s, integrate);
				cudaPrintfDisplay(stdout, true);

				err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
		
				if(err != cudaSuccess) {
					fprintf(stderr, " Error in LocalSearch2D Kernel ( Stage - 3 )  \n");
					fprintf(stderr,"CUDA Error: %s \n\n" , cudaGetErrorString(err));
					exit(1);
				}		
			}
		
		} else if (Data_MeshType == CARTESIAN) {
		
			if (Dimensions == 3) { // 3D case 
		
				cudaPrintfInit(); 
				
				compute_points_Cartesian3D_1<<<nblocks, nthreads>>>(posx, posy, posz, tloc_device, xn0, xn1, xn2, Tracer_dev.x, Tracer_dev.y, Tracer_dev.z, Tracer_dev.LeftDomain, velocity_dev, h_device, integrate, 1, N1);
			
				cudaPrintfDisplay(stdout, true);
			
				err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
				if(err != cudaSuccess) {
					fprintf(stderr, "Error in compute_points_Cartesian3D_1 Kernel ( Stage - 2 )  \n");
					fprintf(stderr,"CUDA Error: %s \n\n" , cudaGetErrorString(err));
					exit(1);
				}			
			
			} else { // 2D case
		
				cudaPrintfInit(); 
				
				compute_points_Cartesian2D_1<<<nblocks, nthreads>>>(posx, posy, tloc_device, xn0, xn1, Tracer_dev.x, Tracer_dev.y, Tracer_dev.LeftDomain, velocity_dev, h_device, integrate, 1, N1);
			
				cudaPrintfDisplay(stdout, true);
			
				err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
				if(err != cudaSuccess) {
					fprintf(stderr, "Error in compute_points_Cartesian2D_1 Kernel ( Stage - 2 )  \n");
					fprintf(stderr,"CUDA Error: %s \n\n" , cudaGetErrorString(err));
					exit(1);
				}
			}
		
		}
				
		for (int k =0; k < N1; k++) {
		
			t11[k] = t11[k];
		
			tloc_new[k] = (t11[k] - Data_TMin)/ (Data_TMax - Data_TMin);	
		
		}
		
		err = cudaMemcpy(tloc_device, tloc_new , N1 * sizeof(double) , cudaMemcpyHostToDevice);
		if(err != cudaSuccess) {
			fprintf(stderr, "Memory copy from Host structure array tloc_new to Device structure array tloc_device failed\n");
			printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
			exit(1);
		}	
		

		// Compute k3 //
		if (Data_MeshType == UNSTRUCTURED)  { // Unstructured Mesh
			
			if (Dimensions == 3) { // 3D case
		
				cudaPrintfInit(); 
				compute_points_Unstructure3D_1<<<nblocks, nthreads>>>(posx, posy, posz, tloc_device, xn0, xn1, xn2, Tracer_dev.x, Tracer_dev.y, Tracer_dev.z, Tracer_dev.ElementIndex, r, s, t, velocity_dev, MeshElementArray_device, h_device, integrate, 2, N1);
		
				cudaPrintfDisplay(stdout, true);
				err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
				if(err != cudaSuccess) {
					fprintf(stderr, "Error in compute_points_Unstructure3D_1 Kernel. ( Stage - 3 )  \n");
					fprintf(stderr,"CUDA Error: %s \n\n" , cudaGetErrorString(err));
					exit(1);
				} 
		
		
				// Do 3D local search .....
				cudaPrintfInit(); 
				LocalSearch3D<<<nblocks, nthreads>>>(posx, posy, posz, MeshElementArray_device, MeshNodeArray_double_device, Tracer_dev.ElementIndex, N1, 0, r, s, t, integrate);
				cudaPrintfDisplay(stdout, true);

				err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
		
				if(err != cudaSuccess) {
					fprintf(stderr, " Error in LocalSearch3D Kernel. (Stage - 4)  \n");
					fprintf(stderr,"CUDA Error: %s \n\n" , cudaGetErrorString(err));
					exit(1);
				}
			
			} else { // 2D case
		
				cudaPrintfInit(); 
				compute_points_Unstructure2D_1<<<nblocks, nthreads>>>(posx, posy, tloc_device, xn0, xn1, Tracer_dev.x, Tracer_dev.y, Tracer_dev.ElementIndex, r, s, velocity_dev, MeshElementArray_device, h_device, integrate, 2, N1);
		
				cudaPrintfDisplay(stdout, true);
				err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
				if(err != cudaSuccess) {
					fprintf(stderr, "Error in compute_points_Unstructure2D_1 Kernel. (Stage - 3)  \n");
					fprintf(stderr,"CUDA Error: %s \n\n" , cudaGetErrorString(err));
					exit(1);
				} 
		
			
				// Do 2D local search .....
				cudaPrintfInit(); 
				LocalSearch2D<<<nblocks, nthreads>>>(posx, posy, MeshElementArray_device, MeshNodeArray_double_device, Tracer_dev.ElementIndex, N1, 0, r, s, integrate);
				cudaPrintfDisplay(stdout, true);

				err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
		
				if(err != cudaSuccess) {
					fprintf(stderr, " Error in LocalSearch2D Kernel (Stage - 4)  \n");
					fprintf(stderr,"CUDA Error: %s \n\n" , cudaGetErrorString(err));
					exit(1);
				}
					
			}
		
		} else if (Data_MeshType == CARTESIAN) {
		
			if (Dimensions == 3) { // 3D case 
		
				cudaPrintfInit(); 
				
				compute_points_Cartesian3D_1<<<nblocks, nthreads>>>(posx, posy, posz, tloc_device, xn0, xn1, xn2, Tracer_dev.x, Tracer_dev.y, Tracer_dev.z, Tracer_dev.LeftDomain, velocity_dev, h_device, integrate, 2, N1);
			
				cudaPrintfDisplay(stdout, true);
			
				err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
				if(err != cudaSuccess) {
					fprintf(stderr, "Error in compute_points_Cartesian3D_1 Kernel ( Stage - 3 )  \n");
					fprintf(stderr,"CUDA Error: %s \n\n" , cudaGetErrorString(err));
					exit(1);
				}
			
			} else { // 2D case
		
				cudaPrintfInit(); 
				
				compute_points_Cartesian2D_1<<<nblocks, nthreads>>>(posx, posy, tloc_device, xn0, xn1, Tracer_dev.x, Tracer_dev.y, Tracer_dev.LeftDomain, velocity_dev, h_device, integrate, 2, N1);
			
				cudaPrintfDisplay(stdout, true);
			
				err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
				if(err != cudaSuccess) {
					fprintf(stderr, "Error in compute_points_Cartesian2D_1 Kernel ( Stage - 3 )  \n");
					fprintf(stderr,"CUDA Error: %s \n\n" , cudaGetErrorString(err));
					exit(1);
				}
			}
		}
	
		for (int k =0; k < N1; k++) {
			
			t11[k] = t11[k] + 0.5*h_new[k];
			tloc_new[k] = (t11[k] - Data_TMin)/ (Data_TMax - Data_TMin);	
		
		}
		
		err = cudaMemcpy(tloc_device, tloc_new , N1 * sizeof(double) , cudaMemcpyHostToDevice);
		if(err != cudaSuccess) {
			fprintf(stderr, "Memory copy from Host structure array tloc_new to Device structure array tloc_device failed\n");
			printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
			exit(1);
		}	
		
		//  Compute k4  //
		if (Data_MeshType == UNSTRUCTURED)  { // Unstructured Mesh
			
			if (Dimensions == 3) { // 3D case 
		
				cudaPrintfInit(); 
		
				compute_points_Unstructure3D_2<<<nblocks, nthreads>>>(posx, posy, posz, tloc_device, xn0, xn1, xn2, Tracer_dev.x, Tracer_dev.y, Tracer_dev.z, Tracer_dev.ElementIndex, r, s, t, velocity_dev, MeshElementArray_device, h_device, integrate, N1);
		
				cudaPrintfDisplay(stdout, true);
				err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
				if(err != cudaSuccess) {
					fprintf(stderr, "Error in compute_points_Unstructure3D_2 Kernel.  \n");
					fprintf(stderr,"CUDA Error: %s \n\n" , cudaGetErrorString(err));
					exit(1);
				}
			
				// Do 3D local search .....
				cudaPrintfInit(); 
				LocalSearch3D<<<nblocks, nthreads>>>(Tracer_dev.x, Tracer_dev.y, Tracer_dev.z, MeshElementArray_device, MeshNodeArray_double_device, Tracer_dev.ElementIndex, N1, 0, r, s, t, integrate);
				cudaPrintfDisplay(stdout, true);

				err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
		
				if(err != cudaSuccess) {
					fprintf(stderr, " Error in LocalSearch3D Kernel ( Stage - 5 )  \n");
					fprintf(stderr,"CUDA Error: %s \n\n" , cudaGetErrorString(err));
					exit(1);
				}
			
			} else { // 2D case
		
				cudaPrintfInit(); 
		
				compute_points_Unstructure2D_2<<<nblocks, nthreads>>>(posx, posy, tloc_device, xn0, xn1, Tracer_dev.x, Tracer_dev.y, Tracer_dev.ElementIndex, r, s, velocity_dev, MeshElementArray_device, h_device, integrate, N1);
		
				cudaPrintfDisplay(stdout, true);
				err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
				if(err != cudaSuccess) {
					fprintf(stderr, "Error in compute_points_Unstructure2D_2 Kernel.  \n");
					fprintf(stderr,"CUDA Error: %s \n\n" , cudaGetErrorString(err));
					exit(1);
				}
			
		
				// Do 2D local search .....
				cudaPrintfInit(); 
				LocalSearch2D<<<nblocks, nthreads>>>(Tracer_dev.x, Tracer_dev.y, MeshElementArray_device, MeshNodeArray_double_device, Tracer_dev.ElementIndex, N1, 0, r, s, integrate);
				cudaPrintfDisplay(stdout, true);

				err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
		
				if(err != cudaSuccess) {
					fprintf(stderr, " Error in LocalSearch2D Kernel ( Stage - 5 )  \n");
					fprintf(stderr,"CUDA Error: %s \n\n" , cudaGetErrorString(err));
					exit(1);
				}		
			}
			
		} else if (Data_MeshType == CARTESIAN) {
		
			if (Dimensions == 3) { // 3D case 
		
				cudaPrintfInit(); 
				
				compute_points_Cartesian3D_2<<<nblocks, nthreads>>>(posx, posy, posz, tloc_device, xn0, xn1, xn2, Tracer_dev.x, Tracer_dev.y, Tracer_dev.z, Tracer_dev.LeftDomain, velocity_dev, h_device, integrate, N1);
			
				cudaPrintfDisplay(stdout, true);
			
				err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
				if(err != cudaSuccess) {
					fprintf(stderr, "Error in compute_points_Cartesian3D_2 Kernel ( Stage - 4 )  \n");
					fprintf(stderr,"CUDA Error: %s \n\n" , cudaGetErrorString(err));
					exit(1);
				}
			
			} else { // 2D case
		
				cudaPrintfInit(); 
				
				compute_points_Cartesian2D_2<<<nblocks, nthreads>>>(posx, posy, tloc_device, xn0, xn1, Tracer_dev.x, Tracer_dev.y, Tracer_dev.LeftDomain, velocity_dev, h_device, integrate, N1);
			
				cudaPrintfDisplay(stdout, true);
			
				err = cudaThreadSynchronize(); // Synchronize function that will wait till all threads are done computing ....
				if(err != cudaSuccess) {
					fprintf(stderr, "Error in compute_points_Cartesian2D_2 Kernel ( Stage - 4 )  \n");
					fprintf(stderr,"CUDA Error: %s \n\n" , cudaGetErrorString(err));
					exit(1);
				}
		
			}
		
		}
		

		
	}
	
	cudaFree(tloc_device);
	cudaFree(h_device);
	cudaFreeHost(tloc_new);
	cudaFreeHost(h_new);
	
	cudaFree(t11_device);
	cudaFree(t22_device);
	cudaFreeHost(t11);
	cudaFreeHost(t22);
	
	
	cudaEventRecord( stop_int, 0 );
	cudaEventSynchronize( stop_int );
	cudaEventElapsedTime( &elapsedTime, start_int, stop_int );
	printf( "Time for integration: %3.2f ms\n", elapsedTime );
	cudaEventDestroy( stop_int ) ;	
	cudaEventDestroy( start_int ) ;

}


