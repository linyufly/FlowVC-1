
# include "allocate.h"

// This File defines all the functions related to memory allocation of variables on Host (CPU) and Device (GPU) respectively.



// This function will allocate parameters used by integration algorithm.
void alloc_misc(void) {

	printf("Allocating misc varibles .. \n\n");
	
	// Allocating misc variables onto GPU
	err = cudaMalloc((double **)&posx, N_Total*sizeof(double));
	if(err != cudaSuccess) {
		fprintf(stderr, "Something went terribly wrong in allocating posx.\n");
		printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
		exit(1);
	}
	
	err = cudaMalloc((double **)&posy, N_Total*sizeof(double));
	if(err != cudaSuccess) {
		fprintf(stderr, "Something went terribly wrong in allocating posy.\n");
		printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
		exit(1);
	}


	if (Dimensions == 3) {
	
		err = cudaMalloc((double **)&posz, N_Total*sizeof(double));
		if(err != cudaSuccess) {
			fprintf(stderr, "Something went terribly wrong in allocating posz.\n");
			printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
			exit(1);
		}
	}				
	if (Data_MeshType == UNSTRUCTURED) {	
	
		err = cudaMalloc((double **)&r, N_Total*sizeof(double));
		if(err != cudaSuccess) {
			fprintf(stderr, "Something went terribly wrong in allocating r.\n");
			printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
			exit(1);
		}
	
		err = cudaMalloc((double **)&s, N_Total*sizeof(double));
		if(err != cudaSuccess) {
			fprintf(stderr, "Something went terribly wrong in allocating s.\n");
			printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
			exit(1);
		}
	
		err = cudaMalloc((double **)&t, N_Total*sizeof(double));
		if(err != cudaSuccess) {
			fprintf(stderr, "Something went terribly wrong in allocating t.\n");
			printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
			exit(1);
		}	
		
	} 
	
	err = cudaMalloc((double **)&xn0, N_Total*sizeof(double));
	if(err != cudaSuccess) {
		fprintf(stderr, "Something went terribly wrong in allocating xn0.\n");
		printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
		exit(1);
	}
	
	err = cudaMalloc((double **)&xn1, N_Total*sizeof(double));
	if(err != cudaSuccess) {
		fprintf(stderr, "Something went terribly wrong in allocating xn1.\n");
		printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
		exit(1);
	}
	
	if (Dimensions == 3) {
		err = cudaMalloc((double **)&xn2, N_Total*sizeof(double));
		if(err != cudaSuccess) {
			fprintf(stderr, "Something went terribly wrong in allocating xn2.\n");
			printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
			exit(1);
		}
	}
	
	err = cudaMalloc((int **)&integrate, N_Total*sizeof(int));
	if(err != cudaSuccess) {
		fprintf(stderr, "Memory Allocation of integrate failed .\n");
		printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
		exit(1);
	}
	
	printf("Misc Variables Allocated .... \n\n");
}

// This function allocates all host variables. Page-locked memory Usage is enabled by default. It can be disabled from the input file.
void allocate_host(void)
{

	printf("Allocating Host Variables .... \n\n");

	int num;

	if (Data_MeshType == CARTESIAN)
		num = N_Vel;
	else
		num = Vel_MeshNumNodes;

	if (Memory_Usage_Tracer == Use_Pinned_Memory) { // Use Pinned-Memory
	
		// Launch_time array keeps track of different launch times of tracer grid .....
		err = cudaHostAlloc( (void**)&Launch_time, Trace_NumLaunchTimes * sizeof(double ), cudaHostAllocDefault );
		if(err != cudaSuccess) {
			fprintf(stderr, "Something went terribly wrong in allocating Host Memory for Launch_time.\n");
			printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
			exit(1);
		}
		
		// Output_time keeps track of output time stamps for output of the tracer grid .....
		err = cudaHostAlloc( (void**)&Output_time, Output_TRes * sizeof(double ), cudaHostAllocDefault );
		if(err != cudaSuccess) {
			fprintf(stderr, "Something went terribly wrong in allocating Host Memory for Output_time.\n");
			printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
			exit(1);
		}
		
		// Data_Time keeps track of the time stamps of our data.
		err = cudaHostAlloc( (void**)&DataTime1, (Data_TRes - 1) * sizeof(Launch ), cudaHostAllocDefault );
		if(err != cudaSuccess) {
			fprintf(stderr, "Something went terribly wrong in allocating Host Memory for Data_Time.\n");
			printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
			exit(1);
		}
		
		
		if (Trace_ReleaseStrategy == 0) { // Traditional Release.
		
			// Tracer is an structure of array used to store position and time related information.
			// x position
			err = cudaHostAlloc( (void**)&Tracer.x, N_Total * sizeof(double ), cudaHostAllocDefault );
			if(err != cudaSuccess) {
				fprintf(stderr, "Something went terribly wrong in allocating Host Memory for Tracer.x .\n");
				printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
				exit(1);
			}
		
			// y-position
			err = cudaHostAlloc( (void**)&Tracer.y, N_Total * sizeof(double ), cudaHostAllocDefault );
			if(err != cudaSuccess) {
				fprintf(stderr, "Something went terribly wrong in allocating Host Memory for Tracer.y .\n");
				printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
				exit(1);
			}
		
	
				// z-position
				err = cudaHostAlloc( (void**)&Tracer.z, N_Total * sizeof(double ), cudaHostAllocDefault );
				if(err != cudaSuccess) {
					fprintf(stderr, "Something went terribly wrong in allocating Host Memory for Tracer.z .\n");
					printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
					exit(1);
				}
		
			// ElementIndex
			err = cudaHostAlloc( (void**)&Tracer.ElementIndex, N_Total * sizeof(int ), cudaHostAllocDefault );
			if(err != cudaSuccess) {
				fprintf(stderr, "Something went terribly wrong in allocating Host Memory for Tracer.ElementIndex .\n");
				printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
				exit(1);
			}
		
			// LeftDomain
			err = cudaHostAlloc( (void**)&Tracer.LeftDomain, N_Total * sizeof(int ), cudaHostAllocDefault );
			if(err != cudaSuccess) {
				fprintf(stderr, "Something went terribly wrong in allocating Host Memory for Tracer.LeftDomain .\n");
				printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
				exit(1);
			}
		
			// Start Time
			err = cudaHostAlloc( (void**)&Tracer.Start_time, N_Total * sizeof(double ), cudaHostAllocDefault );
			if(err != cudaSuccess) {
				fprintf(stderr, "Something went terribly wrong in allocating Host Memory for Tracer.Start_time .\n");
				printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
				exit(1);
			}
		
			// Stop time
			err = cudaHostAlloc( (void**)&Tracer.Stop_time, N_Total * sizeof(double ), cudaHostAllocDefault );
			if(err != cudaSuccess) {
				fprintf(stderr, "Something went terribly wrong in allocating Host Memory for Tracer.Stop_time .\n");
				printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
				exit(1);
			}
		
		
		}
		// Allocating velocity related data variable
		err = cudaHostAlloc( (void**)&velocity.u0, num * sizeof(double ), cudaHostAllocDefault );
		if(err != cudaSuccess) {
			fprintf(stderr, "Something went terribly wrong in allocating Host Memory for velocity.u0 .\n");
			printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
			exit(1);
		}
		
		err = cudaHostAlloc( (void**)&velocity.u1, num * sizeof(double ), cudaHostAllocDefault );
		if(err != cudaSuccess) {
			fprintf(stderr, "Something went terribly wrong in allocating Host Memory for velocity.u1 .\n");
			printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
			exit(1);
		}
		
		err = cudaHostAlloc( (void**)&velocity.v0, num * sizeof(double ), cudaHostAllocDefault );
		if(err != cudaSuccess) {
			fprintf(stderr, "Something went terribly wrong in allocating Host Memory for velocity.v0 .\n");
			printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
			exit(1);
		}
		
		err = cudaHostAlloc( (void**)&velocity.v1, num * sizeof(double ), cudaHostAllocDefault );
		if(err != cudaSuccess) {
			fprintf(stderr, "Something went terribly wrong in allocating Host Memory for velocity.v1 .\n");
			printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
			exit(1);
		}
		
		
		err = cudaHostAlloc( (void**)&velocity.w0, num * sizeof(double ), cudaHostAllocDefault );
		if(err != cudaSuccess) {
			fprintf(stderr, "Something went terribly wrong in allocating Host Memory for velocity.w0 .\n");
			printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
			exit(1);
		}
		
		err = cudaHostAlloc( (void**)&velocity.w1, num * sizeof(double ), cudaHostAllocDefault );
		if(err != cudaSuccess) {
			fprintf(stderr, "Something went terribly wrong in allocating Host Memory for velocity.w1 .\n");
			printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
			exit(1);
		}
		
		err = cudaHostAlloc( (void**)&velocity.time0, num * sizeof(double ), cudaHostAllocDefault );
		if(err != cudaSuccess) {
			fprintf(stderr, "Something went terribly wrong in allocating Host Memory for velocity.time0 .\n");
			printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
			exit(1);
		}
		
		err = cudaHostAlloc( (void**)&velocity.time1, num * sizeof(double ), cudaHostAllocDefault );
		if(err != cudaSuccess) {
			fprintf(stderr, "Something went terribly wrong in allocating Host Memory for velocity.time1 .\n");
			printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
			exit(1);
		}
		
		
		
		
	}
	else {
		// Allocating some time related arrays ....
		// Launch_time array keeps track of different launch times of tracer grid .....
		if((Launch_time = (double *)malloc(Trace_NumLaunchTimes*sizeof(double))) == NULL) {
			fprintf(stderr, "Malloc failed for Launch_time buffer \n");
			exit(1);
		}
	
		// Output_time keeps track of output time stamps for output of the tracer grid .....
		if((Output_time = (double *)malloc(Output_TRes*sizeof(double))) == NULL) {
			fprintf(stderr, "Malloc failed for Output_time buffer \n");
			exit(1);
		}

		// Data_Time keeps track of the time stamps of our data.
		if((DataTime1 = (Launch *)malloc((Data_TRes - 1)*sizeof(Launch))) == NULL) {
			fprintf(stderr, "Malloc failed for data_time buffer \n");
			exit(1);
		}

		if (Trace_ReleaseStrategy == 0) {
			// Allocation of our varibles used for optimization .....

			// Tracer is an structure of array used to store position and time related information.
			// Intuitively x,y,z represents the x, the y and th z position of the tracer. This approach is SOA(Structure of Arrays) thus we have allocated every variables of the structure. ElementIndex represents the elementindex of a particular tracer point with respect to unstructured grid.

			if((Tracer.x = (double *)malloc(N_Total*sizeof(double))) == NULL) {
				fprintf(stderr, "Malloc failed for Tracer.x \n");
				exit(1);
			}

			if((Tracer.y = (double *)malloc(N_Total*sizeof(double))) == NULL) {
				fprintf(stderr, "Malloc failed for Tracer.y \n");
				exit(1);
			}

		
				if((Tracer.z = (double *)malloc(N_Total*sizeof(double))) == NULL) {
					fprintf(stderr, "Malloc failed for Tracer.z \n");
					exit(1);
				}
		
			if((Tracer.ElementIndex = (int *)malloc(N_Total*sizeof(int))) == NULL) {
				fprintf(stderr, "Malloc failed for Tracer.ElementIndex \n");
				exit(1);
			}
	
			if((Tracer.LeftDomain = (int *)malloc(N_Total*sizeof(int))) == NULL) {
				fprintf(stderr, "Malloc failed for Tracer.leftDomain \n");
				exit(1);
			}

			if((Tracer.Start_time = (double *)malloc(N_Total*sizeof(double))) == NULL) {
				fprintf(stderr, "Malloc failed for Tracer.Start_time \n");
				exit(1);
			}

			if((Tracer.Stop_time = (double *)malloc(N_Total*sizeof(double))) == NULL) {
				fprintf(stderr, "Malloc failed for Tracer.Stop_time \n");
				exit(1);
			}
		
		}
		// Allocating velocity related data variable

		if((velocity.u0 = (double *)malloc(num*sizeof(double))) == NULL) {
			fprintf(stderr, "Malloc failed for velocity.u0 \n");
			exit(1);
		}
	
		if((velocity.u1 = (double *)malloc(num*sizeof(double))) == NULL) {
			fprintf(stderr, "Malloc failed for velocity.u1 \n");
			exit(1);
		}
	
		if((velocity.v0 = (double *)malloc(num*sizeof(double))) == NULL) {
			fprintf(stderr, "Malloc failed for velocity.v0 \n");
			exit(1);
		}
	
		if((velocity.v1 = (double *)malloc(num*sizeof(double))) == NULL) {
			fprintf(stderr, "Malloc failed for velocity.v1 \n");
			exit(1);
		}

		if((velocity.w0 = (double *)malloc(num*sizeof(double))) == NULL) {
			fprintf(stderr, "Malloc failed for velocity.w0 \n");
			exit(1);
		}
	
		if((velocity.w1 = (double *)malloc(num*sizeof(double))) == NULL) {
			fprintf(stderr, "Malloc failed for velocity.w1 \n");
			exit(1);
		}

		if((velocity.time0 = (double *)malloc(num*sizeof(double))) == NULL) {
			fprintf(stderr, "Malloc failed for velocity.time0 \n");
			exit(1);
		}
	
		if((velocity.time1 = (double *)malloc(num*sizeof(double))) == NULL) {
			fprintf(stderr, "Malloc failed for velocity.time1 \n");
			exit(1);
		}


	}

	if(Data_MeshType == UNSTRUCTURED) { 
		// Some variables for global search (Required only for unstructured mesh) ........

		Number_Points_loaded = CONSTANT_MEMORY;

		if((outputid = (int *)malloc(Number_Points_loaded*sizeof(int))) == NULL) {
			fprintf(stderr, "Malloc failed for outputid1\n");
			exit(1);
		}

		// allocate memory for constant memory on host
		// x_host, y_host and z_host are constant memories that will be copied to GPU ....
		if((x_host = (float *)malloc(Number_Points_loaded*sizeof(float))) == NULL) {
			fprintf(stderr, "Malloc failed for x_host \n");
			exit(1);
		}

		if((y_host = (float *)malloc(Number_Points_loaded*sizeof(float))) == NULL) {
			fprintf(stderr, "Malloc failed for y_host \n");
			exit(1);
		}

		if(Dimensions == 3) {

			if((z_host = (float *)malloc(Number_Points_loaded*sizeof(float))) == NULL) {
				fprintf(stderr, "Malloc failed for z_host \n");
				exit(1);
			}
			
		}

	}


	printf("Host Variables Allocated .... \n\n");
	
}


// Allocating variables on Device (GPU)
void allocate_device(void)
{

	printf("Allocating Device Variables .... \n\n");

	int num;

	if (Data_MeshType == CARTESIAN)
		num = N_Vel;
	else
		num = Vel_MeshNumNodes;


	// Allocating velocity data on GPU
	err = cudaMalloc((double **)&velocity_dev.u0, num*sizeof(double));
	if(err != cudaSuccess) {
		fprintf(stderr, "Something went terribly wrong in allocating mesh velocity data (u0).\n");
		printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
		exit(1);
	}

	err = cudaMalloc((double **)&velocity_dev.u1, num*sizeof(double));
	if(err != cudaSuccess) {
		fprintf(stderr, "Something went terribly wrong in allocating Mesh velocity data (u1).\n");
		printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
		exit(1);
	}

	err = cudaMalloc((double **)&velocity_dev.v0, num*sizeof(double));
	if(err != cudaSuccess) {
		fprintf(stderr, "Something went terribly wrong in allocating mesh velocity data (v0).\n");
		printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
		exit(1);
	}

	err = cudaMalloc((double **)&velocity_dev.v1, num*sizeof(double));
	if(err != cudaSuccess) {
		fprintf(stderr, "Something went terribly wrong in allocating Mesh velocity data (v1).\n");
		printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
		exit(1);
	}


	err = cudaMalloc((double **)&velocity_dev.w0, num*sizeof(double));
	if(err != cudaSuccess) {
		fprintf(stderr, "Something went terribly wrong in allocating mesh velocity data (w0).\n");
		printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
		exit(1);
	}

	err = cudaMalloc((double **)&velocity_dev.w1, num*sizeof(double));
	if(err != cudaSuccess) {
		fprintf(stderr, "Something went terribly wrong in allocating Mesh velocity data (w1).\n");
		printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
		exit(1);
	}

	err = cudaMalloc((double **)&velocity_dev.time0, num*sizeof(double));
	if(err != cudaSuccess) {
		fprintf(stderr, "Something went terribly wrong in allocating mesh velocity data (time0).\n");
		printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
		exit(1);
	}

	err = cudaMalloc((double **)&velocity_dev.time1, num*sizeof(double));
	if(err != cudaSuccess) {
		fprintf(stderr, "Something went terribly wrong in allocating Mesh velocity data (time1).\n");
		printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
		exit(1);
	}


	// Allocating Tracer variable onto GPU
	err = cudaMalloc((double **)&Tracer_dev.x, N_Total*sizeof(double));
	if(err != cudaSuccess) {
		fprintf(stderr, "Something went terribly wrong in allocating tracer (x).\n");
		printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
		exit(1);
	}

	err = cudaMalloc((double **)&Tracer_dev.y, N_Total*sizeof(double));
	if(err != cudaSuccess) {
		fprintf(stderr, "Something went terribly wrong in allocating tracer (y).\n");
		printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
		exit(1);
	}

	if (Dimensions == 3) {
		err = cudaMalloc((double **)&Tracer_dev.z, N_Total*sizeof(double));
		if(err != cudaSuccess) {
			fprintf(stderr, "Something went terribly wrong in allocating tracer (z).\n");
			printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
			exit(1);
		}
	}
	err = cudaMalloc((int **)&Tracer_dev.ElementIndex, N_Total*sizeof(int));
	if(err != cudaSuccess) {
		fprintf(stderr, "Something went terribly wrong in allocating tracer (Elementid).\n");
		printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
		exit(1);
	}
	
	err = cudaMalloc((int **)&Tracer_dev.LeftDomain, N_Total*sizeof(int));
	if(err != cudaSuccess) {
		fprintf(stderr, "Something went terribly wrong in allocating tracer (LeftDomain).\n");
		printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
		exit(1);
	}
	
	err = cudaMalloc((double **)&Tracer_dev.Start_time, N_Total*sizeof(double));
	if(err != cudaSuccess) {
		fprintf(stderr, "Something went terribly wrong in allocating tracer (start time).\n");
		printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
		exit(1);
	}

	err = cudaMalloc((double **)&Tracer_dev.Stop_time, N_Total*sizeof(double));
	if(err != cudaSuccess) {
		fprintf(stderr, "Something went terribly wrong in allocating tracer (stop time).\n");
		printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
		exit(1);
	}




	if (Data_MeshType == UNSTRUCTURED) {

		Number_Points_loaded = CONSTANT_MEMORY;


		// Allocating Element data array on GPU
		err = cudaMalloc((int **)&MeshElementArray_device.Node1, Vel_MeshNumElements*sizeof(int));
		if(err != cudaSuccess) {
			fprintf(stderr, "Something went terribly wrong in allocating Mesh Element Array (Node 1).\n");
			printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
			exit(1);
		  }

		err = cudaMalloc((int **)&MeshElementArray_device.Node2, Vel_MeshNumElements*sizeof(int));
		if(err != cudaSuccess) {
			fprintf(stderr, "Something went terribly wrong in allocating Mesh Element Array(Node 2).\n");
			printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
			exit(1);
		  }

		err = cudaMalloc((int **)&MeshElementArray_device.Node3, Vel_MeshNumElements*sizeof(int));
		if(err != cudaSuccess) {
			fprintf(stderr, "Something went terribly wrong in allocating Mesh Element Array(Node 3).\n");
			printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
			exit(1);
		  }

		err = cudaMalloc((int **)&MeshElementArray_device.Node4, Vel_MeshNumElements*sizeof(int));
		if(err != cudaSuccess) {
			fprintf(stderr, "Something went terribly wrong in allocating Mesh Element Array(Node 4).\n");
			printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
			exit(1);
		  }


		err = cudaMalloc((int **)&MeshElementArray_device.Neighborindex1, Vel_MeshNumElements*sizeof(int));
		if(err != cudaSuccess) {
			fprintf(stderr, "Something went terribly wrong in allocating Mesh Element Array(Neighborindex1).\n");
			//fprintf("CUDA Error: %s \n\n", cudaGetErrorString(err));
			exit(1);
		  }
		
		err = cudaMalloc((int **)&MeshElementArray_device.Neighborindex2, Vel_MeshNumElements*sizeof(int));
		if(err != cudaSuccess) {
			fprintf(stderr, "Something went terribly wrong in allocating Mesh Element Array(Neighborindex2).\n");
			//fprintf("CUDA Error: %s \n\n", cudaGetErrorString(err));
			exit(1);
		  }

		err = cudaMalloc((int **)&MeshElementArray_device.Neighborindex3, Vel_MeshNumElements*sizeof(int));
		if(err != cudaSuccess) {
			fprintf(stderr, "Something went terribly wrong in allocating Mesh Element Array(Neighborindex3).\n");
			//fprintf("CUDA Error: %s \n\n", cudaGetErrorString(err));
			exit(1);
		  }

		err = cudaMalloc((int **)&MeshElementArray_device.Neighborindex4, Vel_MeshNumElements*sizeof(int));
		if(err != cudaSuccess) {
			fprintf(stderr, "Something went terribly wrong in allocating Mesh Element Array(Neighborindex4).\n");
			//fprintf("CUDA Error: %s \n\n", cudaGetErrorString(err));
			exit(1);
		  }

		// Allocating Node array data on GPU
		err = cudaMalloc((float **)&MeshNodeArray_device.x, Vel_MeshNumNodes*sizeof(float));
		if(err != cudaSuccess) {
			fprintf(stderr, "Something went terribly wrong in allocating Mesh Node Array(x position).\n");
			printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
			exit(1);
		  }

		err = cudaMalloc((float **)&MeshNodeArray_device.y, Vel_MeshNumNodes*sizeof(float));
		if(err != cudaSuccess) {
			fprintf(stderr, "Something went terribly wrong in allocating Mesh Node Array(y position).\n");
			printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
			exit(1);
		  }

		err = cudaMalloc((float **)&MeshNodeArray_device.z, Vel_MeshNumNodes*sizeof(float));
		if(err != cudaSuccess) {
			fprintf(stderr, "Something went terribly wrong in allocating Mesh Node Array(z position).\n");
			printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
			exit(1);
		  }

		// Allocating memory for MeshNodeArray_double_device
		err = cudaMalloc((double **)&MeshNodeArray_double_device.x, Vel_MeshNumNodes*sizeof(double));
		if(err != cudaSuccess) {
			fprintf(stderr, "Something went terribly wrong in allocating Mesh Node Array(x position).\n");
			printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
			exit(1);
		  }

		err = cudaMalloc((double **)&MeshNodeArray_double_device.y, Vel_MeshNumNodes*sizeof(double));
		if(err != cudaSuccess) {
			fprintf(stderr, "Something went terribly wrong in allocating Mesh Node Array(y position).\n");
			printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
			exit(1);
		  }

		err = cudaMalloc((double **)&MeshNodeArray_double_device.z, Vel_MeshNumNodes*sizeof(double));
		if(err != cudaSuccess) {
			fprintf(stderr, "Something went terribly wrong in allocating Mesh Node Array(z position).\n");
			printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
			exit(1);
		  }
	
		err = cudaMalloc((int **)&outputid_dev, Number_Points_loaded *sizeof(int));
		if(err != cudaSuccess) {
			fprintf(stderr, "Something went terribly wrong in allocating outputid data structure.\n");
			printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
			exit(1);
		}

	


	}

	printf("Device variables allocated ... \n\n");

}


