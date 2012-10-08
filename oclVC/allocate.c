
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <CL/cl.h>

# include "allocate.h"
# include "globals.h"
# include "structs.h"
# include "settings.h"

void allocate_device(void) {

	printf("Allocating Device Variables .... \n\n");
	int num;

	if (Data_MeshType == CARTESIAN)
		num = N_Vel;
	else
		num = Vel_MeshNumNodes;
		
	// Allocating Misc Variables.
		
	posx = clCreateBuffer(context, CL_MEM_READ_WRITE, N_Total * sizeof(double), NULL, &ret);
	if (ret != CL_SUCCESS) {	
		fprintf(stderr, "Error Allocating memory for posx. \n\n");
		exit(1);
	}
	
	posy = clCreateBuffer(context, CL_MEM_READ_WRITE, N_Total * sizeof(double), NULL, &ret);
	if (ret != CL_SUCCESS) {	
		fprintf(stderr, "Error Allocating memory for posy. \n\n");
		exit(1);
	}
	
	xn0 = clCreateBuffer(context, CL_MEM_READ_WRITE, N_Total * sizeof(double), NULL, &ret);
	if (ret != CL_SUCCESS) {	
		fprintf(stderr, "Error Allocating memory for xn0. \n\n");
		exit(1);
	}
	
	xn1 = clCreateBuffer(context, CL_MEM_READ_WRITE, N_Total * sizeof(double), NULL, &ret);
	if (ret != CL_SUCCESS) {	
		fprintf(stderr, "Error Allocating memory for xn1. \n\n");
		exit(1);
	}
	
	integrate = clCreateBuffer(context, CL_MEM_READ_WRITE, N_Total * sizeof(int), NULL, &ret);
	if (ret != CL_SUCCESS) {	
		fprintf(stderr, "Error Allocating memory for integrate. \n\n");
		exit(1);
	}
	
	if (Dimensions == 3) {
	
		posz = clCreateBuffer(context, CL_MEM_READ_WRITE, N_Total * sizeof(double), NULL, &ret);
		if (ret != CL_SUCCESS) {	
			fprintf(stderr, "Error Allocating memory for posz. \n\n");
			exit(1);
		}
		
		xn2 = clCreateBuffer(context, CL_MEM_READ_WRITE, N_Total * sizeof(double), NULL, &ret);
		if (ret != CL_SUCCESS) {	
			fprintf(stderr, "Error Allocating memory for xn2. \n\n");
			exit(1);
		}
		
	}
	
	if (Data_MeshType == UNSTRUCTURED) {
	
		r = clCreateBuffer(context, CL_MEM_READ_WRITE, N_Total * sizeof(double), NULL, &ret);
		if (ret != CL_SUCCESS) {	
			fprintf(stderr, "Error Allocating memory for r. \n\n");
			exit(1);
		}
	
		s = clCreateBuffer(context, CL_MEM_READ_WRITE, N_Total * sizeof(double), NULL, &ret);
		if (ret != CL_SUCCESS) {	
			fprintf(stderr, "Error Allocating memory for s. \n\n");
			exit(1);
		}
	
		t = clCreateBuffer(context, CL_MEM_READ_WRITE, N_Total * sizeof(double), NULL, &ret);
		if (ret != CL_SUCCESS) {	
			fprintf(stderr, "Error Allocating memory for t. \n\n");
			exit(1);
		}
	
		eid = clCreateBuffer(context, CL_MEM_READ_WRITE, N_Total * sizeof(int), NULL, &ret);
		if (ret != CL_SUCCESS) {	
			fprintf(stderr, "Error Allocating memory for eid. \n\n");
			exit(1);
		}
	
	}
	
	// Allocating Tracer variables
	
	Vel_U0 = clCreateBuffer(context, CL_MEM_READ_ONLY , num * sizeof(double), NULL, &ret);
	if (ret != CL_SUCCESS) {	
		fprintf(stderr, "Error Allocating memory for Vel_U0. \n\n");
		exit(1);
	}
	
	Vel_U1 = clCreateBuffer(context, CL_MEM_READ_ONLY , num * sizeof(double), NULL, &ret);
	if (ret != CL_SUCCESS) {	
		fprintf(stderr, "Error Allocating memory for Vel_U1. \n\n");
		exit(1);
	}
	
	Vel_V0 = clCreateBuffer(context, CL_MEM_READ_ONLY , num * sizeof(double), NULL, &ret);
	if (ret != CL_SUCCESS) {	
		fprintf(stderr, "Error Allocating memory for Vel_V0. \n\n");
		exit(1);
	}
	
	Vel_V1 = clCreateBuffer(context, CL_MEM_READ_ONLY, num * sizeof(double), NULL, &ret);
	if (ret != CL_SUCCESS) {	
		fprintf(stderr, "Error Allocating memory for Vel_V1. \n\n");
		exit(1);
	}
	
	Vel_W0 = clCreateBuffer(context, CL_MEM_READ_ONLY , num * sizeof(double), NULL, &ret);
	if (ret != CL_SUCCESS) {	
		fprintf(stderr, "Error Allocating memory for Vel_W0. \n\n");
		exit(1);
	}
	
	Vel_W1 = clCreateBuffer(context, CL_MEM_READ_ONLY , num * sizeof(double), NULL, &ret);
	if (ret != CL_SUCCESS) {	
		fprintf(stderr, "Error Allocating memory for Vel_W1. \n\n");
		exit(1);
	}
	
	x_dev = clCreateBuffer(context, CL_MEM_READ_WRITE, N_Total * sizeof(double), NULL, &ret);
	if (ret != CL_SUCCESS) {	
		fprintf(stderr, "Error Allocating memory for x_dev. \n\n");
		exit(1);
	}
	
	y_dev = clCreateBuffer(context, CL_MEM_READ_WRITE, N_Total * sizeof(double), NULL, &ret);
	if (ret != CL_SUCCESS) {	
		fprintf(stderr, "Error Allocating memory for y_dev. \n\n");
		exit(1);
	}
	
	if (Dimensions == 3) {
		z_dev = clCreateBuffer(context, CL_MEM_READ_WRITE, N_Total * sizeof(double), NULL, &ret);
		if (ret != CL_SUCCESS) {	
			fprintf(stderr, "Error Allocating memory for z_dev. \n\n");
			exit(1);
		}
	}
	
	Start_time_dev = clCreateBuffer(context, CL_MEM_READ_WRITE, N_Total * sizeof(double), NULL, &ret);
	if (ret != CL_SUCCESS) {	
		fprintf(stderr, "Error Allocating memory for Start_time_dev. \n\n");
		exit(1);
	}
	
	Stop_time_dev = clCreateBuffer(context, CL_MEM_READ_WRITE, N_Total * sizeof(double), NULL, &ret);
	if (ret != CL_SUCCESS) {	
		fprintf(stderr, "Error Allocating memory for Stop_time_dev. \n\n");
		exit(1);
	}
	
	ElementIndex_dev = clCreateBuffer(context, CL_MEM_READ_WRITE, N_Total * sizeof(int), NULL, &ret);
	if (ret != CL_SUCCESS) {	
		fprintf(stderr, "Error Allocating memory for ElementIndex_dev. \n\n");
		exit(1);
	}
	
	LeftDomain_dev = clCreateBuffer(context, CL_MEM_READ_WRITE, N_Total * sizeof(int), NULL, &ret);
	if (ret != CL_SUCCESS) {	
		fprintf(stderr, "Error Allocating memory for LeftDomain_dev. \n\n");
		exit(1);
	}
	//MeshElementArray.Node1
	if (Data_MeshType == UNSTRUCTURED) {
	
		Mesh_Element_Node1 = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, Vel_MeshNumElements * sizeof(int), MeshElementArray.Node1, &ret);
		if (ret != CL_SUCCESS) {	
			fprintf(stderr, "Error Allocating memory for Mesh_Element_Node1. \n\n");
			exit(1);
		}
		Mesh_Element_Node2 = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, Vel_MeshNumElements * sizeof(int), MeshElementArray.Node2, &ret);
		if (ret != CL_SUCCESS) {	
			fprintf(stderr, "Error Allocating memory for Mesh_Element_Node2. \n\n");
			exit(1);
		}
		Mesh_Element_Node3 = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, Vel_MeshNumElements * sizeof(int), MeshElementArray.Node3, &ret);
		if (ret != CL_SUCCESS) {	
			fprintf(stderr, "Error Allocating memory for Mesh_Element_Node3. \n\n");
			exit(1);
		}
		Mesh_Element_Node4 = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, Vel_MeshNumElements * sizeof(int), MeshElementArray.Node4, &ret);
		if (ret != CL_SUCCESS) {	
			fprintf(stderr, "Error Allocating memory for Mesh_Element_Node4. \n\n");
			exit(1);
		}
		
		Mesh_Element_Neighborindex1 = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, Vel_MeshNumElements * sizeof(int), MeshElementArray.Neighborindex1, &ret);
		if (ret != CL_SUCCESS) {	
			fprintf(stderr, "Error Allocating memory for Mesh_Element_Neighborindex1. \n\n");
			exit(1);
		}
		
		Mesh_Element_Neighborindex2 = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, Vel_MeshNumElements * sizeof(int), MeshElementArray.Neighborindex2, &ret);
		if (ret != CL_SUCCESS) {	
			fprintf(stderr, "Error Allocating memory for Mesh_Element_Neighborindex2. \n\n");
			exit(1);
		}
		
		Mesh_Element_Neighborindex3 = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, Vel_MeshNumElements * sizeof(int), MeshElementArray.Neighborindex3, &ret);
		if (ret != CL_SUCCESS) {	
			fprintf(stderr, "Error Allocating memory for Mesh_Element_Neighborindex3. \n\n");
			exit(1);
		}
		
		Mesh_Element_Neighborindex4 = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, Vel_MeshNumElements * sizeof(int), MeshElementArray.Neighborindex4, &ret);
		if (ret != CL_SUCCESS) {	
			fprintf(stderr, "Error Allocating memory for Mesh_Element_Neighborindex4. \n\n");
			exit(1);
		}
	
		Mesh_Node_x = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, Vel_MeshNumNodes * sizeof(double), MeshNodeArray_double.x, &ret);
		if (ret != CL_SUCCESS) {	
			fprintf(stderr, "Error Allocating memory for Mesh_Node_x. \n\n");
			exit(1);
		}
		
		Mesh_Node_y = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, Vel_MeshNumNodes * sizeof(double), MeshNodeArray_double.y, &ret);
		if (ret != CL_SUCCESS) {	
			fprintf(stderr, "Error Allocating memory for Mesh_Node_y. \n\n");
			exit(1);
		}
	
		Mesh_Node_z = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, Vel_MeshNumNodes * sizeof(double), MeshNodeArray_double.z, &ret);
		if (ret != CL_SUCCESS) {	
			fprintf(stderr, "Error Allocating memory for Mesh_Node_z. \n\n");
			exit(1);
		}
	
	
	}
		
	printf("Device memory allocated successfully .. \n\n");
}





void allocate_host(void) {

	printf("Allocating Host Variables .... \n\n");

	int num;

	if (Data_MeshType == CARTESIAN)
		num = N_Vel;
	else
		num = Vel_MeshNumNodes;


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


	printf("Host Variables Allocated .... \n\n");
	
}

