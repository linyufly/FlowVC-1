#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <CL/cl.h>
#include "structs.h"
#include "settings.h"
#pragma OPENCL EXTENSION cl_amd_fp64 : enable
#include "globals.h"
//#include <shrQATest.h>
//#include <oclUtils.h>
#pragma OPENCL EXTENSION cl_amd_fp64 : enable

# include "integration3D.h"


void ComputePoints3D_new(double tmin, double tmax, int N1, double Data_TMin, double Data_TMax) {

	double t_int; // Integration time from tmin to tmax.
	t_int = tmax - tmin;
	int n_loops;
	
	float elapsedTime;
	
	if ( fmod(t_int , Int_TimeStep) == 0 ) {
		n_loops = fabs(t_int/Int_TimeStep); // Total number of loops for integrating particles over a single time step.
	} else {
	
		n_loops = fabs(t_int/Int_TimeStep) + 1;
	}
	printf("nloops = %d \n", n_loops);

	// Calculation of Block size for single chunk of integration kernel ...
	if ( N1 % nthreads == 0 ) {
		nblocks = (int)(N1/nthreads);
	}
	else {
		nblocks = (int)((N1/nthreads) + 1);
	}
	
	double *tloc_new;
	double *t11, *t22;
	double *h_new;
	
	cl_mem 	h_device;
	cl_mem	tloc_device;
	
	// Allocating memory for host tloc_new, t11 and h_new
	if((tloc_new = (double *)malloc(N1 * sizeof(double))) == NULL) {
			fprintf(stderr, "Malloc failed for tloc_new buffer \n");
			exit(1);
	}
	
	if((h_new = (double *)malloc(N1 * sizeof(double))) == NULL) {
			fprintf(stderr, "Malloc failed for h_new buffer \n");
			exit(1);
	}
	
	if((t11 = (double *)malloc(N1 * sizeof(double))) == NULL) {
			fprintf(stderr, "Malloc failed for t11 buffer \n");
			exit(1);
	}
	
	if((t22 = (double *)malloc(N1 * sizeof(double))) == NULL) {
			fprintf(stderr, "Malloc failed for t22 buffer \n");
			exit(1);
	}
		
	int a; // counter variable.
	
	for(int k = 0; k < N1; k++) {
	
		t11[k] = tmin;
	
	}
	
	cl_mem t11_device, t22_device;

	// Load the kernel source code into the array source_str
/*	FILE *fp;
	char *source_str;
	size_t source_size;
	
	fp = fopen("integrate.cl", "r");
	if (!fp) {
	    fprintf(stderr, "Failed to load kernel.\n");
	    exit(1);
	}
	
	source_str = (char*)malloc(MAX_SOURCE_SIZE);
	source_size = fread( source_str, 1, MAX_SOURCE_SIZE, fp);
	fclose( fp );	
	cl_program program1;
	// Create a program from the kernel source
    program1 = clCreateProgramWithSource(context, 1, 
    	(const char **)&source_str, (const size_t *)&source_size, &ret);	
    if (ret != CL_SUCCESS) {	
			fprintf(stderr, "Error creating a program for integration3D. \n\n");
			exit(1);
	}
    // Build the program
    
    ret = clBuildProgram(program1, 1, &device_id, NULL, NULL, NULL); 
    if (ret != CL_SUCCESS)
    {
    
    	size_t length;
    	char buffer[1024000];
    	//clGetProgramBuildInfo(program, 1, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &length);
    	fprintf(stderr, "Error reurned %d. \n\n", (int)ret);
    	printf("Error Log: \n\n %s \n\n", buffer);
    	exit(0);
    }
	
   // Create the OpenCL kernel (compute_points_Unstructure3D_1)
    cl_kernel kernel1 = clCreateKernel(program, "compute_points_Unstructure3D_1", &ret);
	if (ret != CL_SUCCESS) {	
			fprintf(stderr, "Error creating a kernel for compute_points_Unstructure3D_1. \n\n");
			exit(1);
	}
		
	// Create the OpenCL kernel (check_int)
    cl_kernel kernel2 = clCreateKernel(program1, "check_int", &ret);
	if (ret != CL_SUCCESS) {	
			fprintf(stderr, "Error creating a kernel for check_int. %d \n\n", (int)ret);
			exit(1);
	}
	
	// Create the OpenCL kernel (initialize_timestep3D)
    cl_kernel kernel3 = clCreateKernel(program, "initialize_timestep3D", &ret);
	if (ret != CL_SUCCESS) {	
			fprintf(stderr, "Error creating a kernel for initialize_timestep3D. \n\n");
			exit(1);
	}
	
*/	int stage;
	for (a = 0; a < n_loops; a++) {

		for (int k = 0; k < N1; k++) {
		
			t11[k] = fmax(t11[k], Tracer.Start_time[k]);
			
			t22[k] = t11[k] + Int_TimeStep;
			
			// Timestep of our integration for this particular point
			if (t22[k] > tmax) {
				t22[k] = tmax;
				h_new[k] = tmax - t11[k];
			} else 
				h_new[k] = Int_TimeStep;
		
		}
		
		// Allocate and memcopy constant device objects (t11_device, t22_device, h_device)
		t11_device = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, N1 * sizeof(double), t11, &ret);
		if (ret != CL_SUCCESS) {	
			fprintf(stderr, "Error Allocating memory for t11_device. \n\n");
			exit(1);
		}
		
		t22_device = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, N1 * sizeof(double), t22, &ret);
		if (ret != CL_SUCCESS) {	
			fprintf(stderr, "Error Allocating memory for t22_device. \n\n");
			exit(1);
		}
		
		h_device = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, N1 * sizeof(double), h_new, &ret);
		if (ret != CL_SUCCESS) {	
			fprintf(stderr, "Error Allocating memory for h_device. \n\n");
			exit(1);
		}
		
		// Check if a point needs to be integrated
		
		// Set the arguments of the kernel (check_int)
		CL_CHECK(clSetKernelArg(kernel2, 0, sizeof(cl_mem), (void *)&Start_time_dev));
		CL_CHECK(clSetKernelArg(kernel2, 1, sizeof(cl_mem), (void *)&Stop_time_dev));
		CL_CHECK(clSetKernelArg(kernel2, 2, sizeof(cl_mem), (void *)&integrate));
		CL_CHECK(clSetKernelArg(kernel2, 3, sizeof(cl_mem), (void *)&t11_device));
		CL_CHECK(clSetKernelArg(kernel2, 4, sizeof(cl_mem), (void *)&t22_device));
		CL_CHECK(clSetKernelArg(kernel2, 5, sizeof(int), (void *)&N1));
		
		// Execute the OpenCL kernel on the list
		int threads = (N1 % 512 == 0 ? N1 : ((int)(N1/512)+1)*512 );
		const size_t global_item_size = threads; // Process the entire lists
		const size_t local_item_size = 512; // Process one item at a time
		CL_CHECK(clEnqueueNDRangeKernel(command_queue, kernel2, 1, NULL, &global_item_size, &local_item_size, 0, NULL, NULL));
		
		// Initialize
		// Set the arguments of the kernel (Initialize_timestep3D)
		CL_CHECK(clSetKernelArg(kernel3, 0, sizeof(cl_mem), (void *)&x_dev));
		CL_CHECK(clSetKernelArg(kernel3, 1, sizeof(cl_mem), (void *)&y_dev));
		CL_CHECK(clSetKernelArg(kernel3, 2, sizeof(cl_mem), (void *)&z_dev));
		CL_CHECK(clSetKernelArg(kernel3, 3, sizeof(cl_mem), (void *)&posx));
		CL_CHECK(clSetKernelArg(kernel3, 4, sizeof(cl_mem), (void *)&posy));
		CL_CHECK(clSetKernelArg(kernel3, 5, sizeof(cl_mem), (void *)&posz));
		CL_CHECK(clSetKernelArg(kernel3, 6, sizeof(cl_mem), (void *)&xn0));
		CL_CHECK(clSetKernelArg(kernel3, 7, sizeof(cl_mem), (void *)&xn1));
		CL_CHECK(clSetKernelArg(kernel3, 8, sizeof(cl_mem), (void *)&xn2));
		CL_CHECK(clSetKernelArg(kernel3, 9, sizeof(cl_mem), (void *)&eid));
		CL_CHECK(clSetKernelArg(kernel3, 10, sizeof(cl_mem), (void *)&ElementIndex_dev));
		CL_CHECK(clSetKernelArg(kernel3, 11, sizeof(int), (void *)&N1));
		
		// Execute the OpenCL kernel on the list
		//global_item_size = N1; // Process the entire lists
		//local_item_size = 512; // Process one item at a time
		CL_CHECK(clEnqueueNDRangeKernel(command_queue, kernel3, 1, NULL, &global_item_size, &local_item_size, 0, NULL, NULL));
	

		
		
		// Do local search .....
		// Set the arguments of the kernel (Initialize_timestep3D)
		CL_CHECK(clSetKernelArg(kernel4, 0, sizeof(cl_mem), (void *)&posx));
		CL_CHECK(clSetKernelArg(kernel4, 1, sizeof(cl_mem), (void *)&posy));
		CL_CHECK(clSetKernelArg(kernel4, 2, sizeof(cl_mem), (void *)&posz));
		CL_CHECK(clSetKernelArg(kernel4, 3, sizeof(cl_mem), (void *)&Mesh_Element_Node1));
		CL_CHECK(clSetKernelArg(kernel4, 4, sizeof(cl_mem), (void *)&Mesh_Element_Node2));
		CL_CHECK(clSetKernelArg(kernel4, 5, sizeof(cl_mem), (void *)&Mesh_Element_Node3));
		CL_CHECK(clSetKernelArg(kernel4, 6, sizeof(cl_mem), (void *)&Mesh_Element_Node4));
		CL_CHECK(clSetKernelArg(kernel4, 7, sizeof(cl_mem), (void *)&Mesh_Element_Neighborindex1));
		CL_CHECK(clSetKernelArg(kernel4, 8, sizeof(cl_mem), (void *)&Mesh_Element_Neighborindex2));
		CL_CHECK(clSetKernelArg(kernel4, 9, sizeof(cl_mem), (void *)&Mesh_Element_Neighborindex3));
		CL_CHECK(clSetKernelArg(kernel4, 10, sizeof(cl_mem), (void *)&Mesh_Element_Neighborindex4));
		CL_CHECK(clSetKernelArg(kernel4, 11, sizeof(cl_mem), (void *)&Mesh_Node_x));
		CL_CHECK(clSetKernelArg(kernel4, 12, sizeof(cl_mem), (void *)&Mesh_Node_y));
		CL_CHECK(clSetKernelArg(kernel4, 13, sizeof(cl_mem), (void *)&Mesh_Node_z));
		CL_CHECK(clSetKernelArg(kernel4, 14, sizeof(cl_mem), (void *)&eid));
		CL_CHECK(clSetKernelArg(kernel4, 15, sizeof(int), (void *)&N1));
		CL_CHECK(clSetKernelArg(kernel4, 16, sizeof(cl_mem), (void *)&r));
		CL_CHECK(clSetKernelArg(kernel4, 17, sizeof(cl_mem), (void *)&s));
		CL_CHECK(clSetKernelArg(kernel4, 18, sizeof(cl_mem), (void *)&t));
		CL_CHECK(clSetKernelArg(kernel4, 19, sizeof(cl_mem), (void *)&integrate));
		
		// Execute the OpenCL kernel on the list
		CL_CHECK(clEnqueueNDRangeKernel(command_queue, kernel4, 1, NULL, &global_item_size, &local_item_size, 0, NULL, NULL));
		
		
		

		for (int k =0; k < N1; k++) {
			t11[k] = fmax(Tracer.Start_time[k], t11[k]);
			
			tloc_new[k] = (t11[k] - Data_TMin)/ (Data_TMax - Data_TMin);
			
		}
		// Allocate and memcopy tloc_device
		tloc_device = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, N1 * sizeof(double), tloc_new, &ret);
		if (ret != CL_SUCCESS) {	
			fprintf(stderr, "Error Allocating memory for tloc_device. \n\n");
			exit(1);
		}
		
		//k1//
		stage = 0;
		// Set the arguments of the kernel1 (compute_points_Unstructure3D_1)
		CL_CHECK(clSetKernelArg(kernel1, 0, sizeof(cl_mem), (void *)&posx));
		CL_CHECK(clSetKernelArg(kernel1, 1, sizeof(cl_mem), (void *)&posy));
		CL_CHECK(clSetKernelArg(kernel1, 2, sizeof(cl_mem), (void *)&posz));
		CL_CHECK(clSetKernelArg(kernel1, 3, sizeof(cl_mem), (void *)&tloc_device));
		CL_CHECK(clSetKernelArg(kernel1, 4, sizeof(cl_mem), (void *)&xn0));
		CL_CHECK(clSetKernelArg(kernel1, 5, sizeof(cl_mem), (void *)&xn1));
		CL_CHECK(clSetKernelArg(kernel1, 6, sizeof(cl_mem), (void *)&xn2));
		CL_CHECK(clSetKernelArg(kernel1, 7, sizeof(cl_mem), (void *)&x_dev));
		CL_CHECK(clSetKernelArg(kernel1, 8, sizeof(cl_mem), (void *)&y_dev));
		CL_CHECK(clSetKernelArg(kernel1, 9, sizeof(cl_mem), (void *)&z_dev));
		CL_CHECK(clSetKernelArg(kernel1, 10, sizeof(cl_mem), (void *)&eid));
		CL_CHECK(clSetKernelArg(kernel1, 11, sizeof(cl_mem), (void *)&r));
		CL_CHECK(clSetKernelArg(kernel1, 12, sizeof(cl_mem), (void *)&s));
		CL_CHECK(clSetKernelArg(kernel1, 13, sizeof(cl_mem), (void *)&t));
		CL_CHECK(clSetKernelArg(kernel1, 14, sizeof(cl_mem), (void *)&Vel_U0));
		CL_CHECK(clSetKernelArg(kernel1, 15, sizeof(cl_mem), (void *)&Vel_U1));
		CL_CHECK(clSetKernelArg(kernel1, 16, sizeof(cl_mem), (void *)&Vel_V0));
		CL_CHECK(clSetKernelArg(kernel1, 17, sizeof(cl_mem), (void *)&Vel_V1));
		CL_CHECK(clSetKernelArg(kernel1, 18, sizeof(cl_mem), (void *)&Vel_W0));
		CL_CHECK(clSetKernelArg(kernel1, 19, sizeof(cl_mem), (void *)&Vel_W1));
		CL_CHECK(clSetKernelArg(kernel1, 20, sizeof(cl_mem), (void *)&Mesh_Element_Node1));
		CL_CHECK(clSetKernelArg(kernel1, 21, sizeof(cl_mem), (void *)&Mesh_Element_Node2));
		CL_CHECK(clSetKernelArg(kernel1, 22, sizeof(cl_mem), (void *)&Mesh_Element_Node3));
		CL_CHECK(clSetKernelArg(kernel1, 23, sizeof(cl_mem), (void *)&Mesh_Element_Node4));
		CL_CHECK(clSetKernelArg(kernel1, 24, sizeof(cl_mem), (void *)&Mesh_Element_Neighborindex1));
		CL_CHECK(clSetKernelArg(kernel1, 25, sizeof(cl_mem), (void *)&Mesh_Element_Neighborindex2));
		CL_CHECK(clSetKernelArg(kernel1, 26, sizeof(cl_mem), (void *)&Mesh_Element_Neighborindex3));
		CL_CHECK(clSetKernelArg(kernel1, 27, sizeof(cl_mem), (void *)&Mesh_Element_Neighborindex4));
		CL_CHECK(clSetKernelArg(kernel1, 28, sizeof(cl_mem), (void *)&h_device));
		CL_CHECK(clSetKernelArg(kernel1, 29, sizeof(cl_mem), (void *)&integrate));
		CL_CHECK(clSetKernelArg(kernel1, 30, sizeof(int), (void *)&stage));
		CL_CHECK(clSetKernelArg(kernel1, 31, sizeof(int), (void *)&N1));
		
		// Execute the OpenCL kernel on the list
		//threads = (N1 % 512 == 0 ? N1 : ((int)(N1/512)+1)*512 );
		//global_item_size = threads; // Process the entire lists
		//const size_t local_item_size = 512; // Process one item at a time
		CL_CHECK(clEnqueueNDRangeKernel(command_queue, kernel1, 1, NULL, &global_item_size, &local_item_size, 0, NULL, NULL));
		// Do local search .....
		// Do local search .....
		// Set the arguments of the kernel (Initialize_timestep3D)
		CL_CHECK(clSetKernelArg(kernel4, 0, sizeof(cl_mem), (void *)&posx));
		CL_CHECK(clSetKernelArg(kernel4, 1, sizeof(cl_mem), (void *)&posy));
		CL_CHECK(clSetKernelArg(kernel4, 2, sizeof(cl_mem), (void *)&posz));
		CL_CHECK(clSetKernelArg(kernel4, 3, sizeof(cl_mem), (void *)&Mesh_Element_Node1));
		CL_CHECK(clSetKernelArg(kernel4, 4, sizeof(cl_mem), (void *)&Mesh_Element_Node2));
		CL_CHECK(clSetKernelArg(kernel4, 5, sizeof(cl_mem), (void *)&Mesh_Element_Node3));
		CL_CHECK(clSetKernelArg(kernel4, 6, sizeof(cl_mem), (void *)&Mesh_Element_Node4));
		CL_CHECK(clSetKernelArg(kernel4, 7, sizeof(cl_mem), (void *)&Mesh_Element_Neighborindex1));
		CL_CHECK(clSetKernelArg(kernel4, 8, sizeof(cl_mem), (void *)&Mesh_Element_Neighborindex2));
		CL_CHECK(clSetKernelArg(kernel4, 9, sizeof(cl_mem), (void *)&Mesh_Element_Neighborindex3));
		CL_CHECK(clSetKernelArg(kernel4, 10, sizeof(cl_mem), (void *)&Mesh_Element_Neighborindex4));
		CL_CHECK(clSetKernelArg(kernel4, 11, sizeof(cl_mem), (void *)&Mesh_Node_x));
		CL_CHECK(clSetKernelArg(kernel4, 12, sizeof(cl_mem), (void *)&Mesh_Node_y));
		CL_CHECK(clSetKernelArg(kernel4, 13, sizeof(cl_mem), (void *)&Mesh_Node_z));
		CL_CHECK(clSetKernelArg(kernel4, 14, sizeof(cl_mem), (void *)&eid));
		CL_CHECK(clSetKernelArg(kernel4, 15, sizeof(int), (void *)&N1));
		CL_CHECK(clSetKernelArg(kernel4, 16, sizeof(cl_mem), (void *)&r));
		CL_CHECK(clSetKernelArg(kernel4, 17, sizeof(cl_mem), (void *)&s));
		CL_CHECK(clSetKernelArg(kernel4, 18, sizeof(cl_mem), (void *)&t));
		CL_CHECK(clSetKernelArg(kernel4, 19, sizeof(cl_mem), (void *)&integrate));
		
		// Execute the OpenCL kernel on the list
		CL_CHECK(clEnqueueNDRangeKernel(command_queue, kernel4, 1, NULL, &global_item_size, &local_item_size, 0, NULL, NULL));
		
		
		// release tloc_device
		clReleaseMemObject(tloc_device);
		
		for (int k =0; k < N1; k++) {
			t11[k] = t11[k]+ (0.5 * h_new[k]);
			//tcurrent = fmax(Tracer.Start_time[k], t1) + (0.5 * h_new[k]);
			tloc_new[k] = (t11[k] - Data_TMin)/ (Data_TMax - Data_TMin);
			//if (k == 2)
			//	printf("tloc[%d] = %0.9f \t tq = %0.9f \n", k, tloc_new[k], t11[k]);		
		}
		// Allocate and memcopy tloc_device
		tloc_device = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, N1 * sizeof(double), tloc_new, &ret);
		if (ret != CL_SUCCESS) {	
			fprintf(stderr, "Error Allocating memory for tloc_device. \n\n");
			exit(1);
		}
		
		
		
		
		
		//k2//
		stage = 1;
		// Set the arguments of the kernel1 (compute_points_Unstructure3D_1)
		CL_CHECK(clSetKernelArg(kernel1, 0, sizeof(cl_mem), (void *)&posx));
		CL_CHECK(clSetKernelArg(kernel1, 1, sizeof(cl_mem), (void *)&posy));
		CL_CHECK(clSetKernelArg(kernel1, 2, sizeof(cl_mem), (void *)&posz));
		CL_CHECK(clSetKernelArg(kernel1, 3, sizeof(cl_mem), (void *)&tloc_device));
		CL_CHECK(clSetKernelArg(kernel1, 4, sizeof(cl_mem), (void *)&xn0));
		CL_CHECK(clSetKernelArg(kernel1, 5, sizeof(cl_mem), (void *)&xn1));
		CL_CHECK(clSetKernelArg(kernel1, 6, sizeof(cl_mem), (void *)&xn2));
		CL_CHECK(clSetKernelArg(kernel1, 7, sizeof(cl_mem), (void *)&x_dev));
		CL_CHECK(clSetKernelArg(kernel1, 8, sizeof(cl_mem), (void *)&y_dev));
		CL_CHECK(clSetKernelArg(kernel1, 9, sizeof(cl_mem), (void *)&z_dev));
		CL_CHECK(clSetKernelArg(kernel1, 10, sizeof(cl_mem), (void *)&eid));
		CL_CHECK(clSetKernelArg(kernel1, 11, sizeof(cl_mem), (void *)&r));
		CL_CHECK(clSetKernelArg(kernel1, 12, sizeof(cl_mem), (void *)&s));
		CL_CHECK(clSetKernelArg(kernel1, 13, sizeof(cl_mem), (void *)&t));
		CL_CHECK(clSetKernelArg(kernel1, 14, sizeof(cl_mem), (void *)&Vel_U0));
		CL_CHECK(clSetKernelArg(kernel1, 15, sizeof(cl_mem), (void *)&Vel_U1));
		CL_CHECK(clSetKernelArg(kernel1, 16, sizeof(cl_mem), (void *)&Vel_V0));
		CL_CHECK(clSetKernelArg(kernel1, 17, sizeof(cl_mem), (void *)&Vel_V1));
		CL_CHECK(clSetKernelArg(kernel1, 18, sizeof(cl_mem), (void *)&Vel_W0));
		CL_CHECK(clSetKernelArg(kernel1, 19, sizeof(cl_mem), (void *)&Vel_W1));
		CL_CHECK(clSetKernelArg(kernel1, 20, sizeof(cl_mem), (void *)&Mesh_Element_Node1));
		CL_CHECK(clSetKernelArg(kernel1, 21, sizeof(cl_mem), (void *)&Mesh_Element_Node2));
		CL_CHECK(clSetKernelArg(kernel1, 22, sizeof(cl_mem), (void *)&Mesh_Element_Node3));
		CL_CHECK(clSetKernelArg(kernel1, 23, sizeof(cl_mem), (void *)&Mesh_Element_Node4));
		CL_CHECK(clSetKernelArg(kernel1, 24, sizeof(cl_mem), (void *)&Mesh_Element_Neighborindex1));
		CL_CHECK(clSetKernelArg(kernel1, 25, sizeof(cl_mem), (void *)&Mesh_Element_Neighborindex2));
		CL_CHECK(clSetKernelArg(kernel1, 26, sizeof(cl_mem), (void *)&Mesh_Element_Neighborindex3));
		CL_CHECK(clSetKernelArg(kernel1, 27, sizeof(cl_mem), (void *)&Mesh_Element_Neighborindex4));
		CL_CHECK(clSetKernelArg(kernel1, 28, sizeof(cl_mem), (void *)&h_device));
		CL_CHECK(clSetKernelArg(kernel1, 29, sizeof(cl_mem), (void *)&integrate));
		CL_CHECK(clSetKernelArg(kernel1, 30, sizeof(int), (void *)&stage));
		CL_CHECK(clSetKernelArg(kernel1, 31, sizeof(int), (void *)&N1));
		
		// Execute the OpenCL kernel on the list
		CL_CHECK(clEnqueueNDRangeKernel(command_queue, kernel1, 1, NULL, &global_item_size, &local_item_size, 0, NULL, NULL));
		
		
		
		// Do local search .....
		// Do local search .....
		// Set the arguments of the kernel (Initialize_timestep3D)
		CL_CHECK(clSetKernelArg(kernel4, 0, sizeof(cl_mem), (void *)&posx));
		CL_CHECK(clSetKernelArg(kernel4, 1, sizeof(cl_mem), (void *)&posy));
		CL_CHECK(clSetKernelArg(kernel4, 2, sizeof(cl_mem), (void *)&posz));
		CL_CHECK(clSetKernelArg(kernel4, 3, sizeof(cl_mem), (void *)&Mesh_Element_Node1));
		CL_CHECK(clSetKernelArg(kernel4, 4, sizeof(cl_mem), (void *)&Mesh_Element_Node2));
		CL_CHECK(clSetKernelArg(kernel4, 5, sizeof(cl_mem), (void *)&Mesh_Element_Node3));
		CL_CHECK(clSetKernelArg(kernel4, 6, sizeof(cl_mem), (void *)&Mesh_Element_Node4));
		CL_CHECK(clSetKernelArg(kernel4, 7, sizeof(cl_mem), (void *)&Mesh_Element_Neighborindex1));
		CL_CHECK(clSetKernelArg(kernel4, 8, sizeof(cl_mem), (void *)&Mesh_Element_Neighborindex2));
		CL_CHECK(clSetKernelArg(kernel4, 9, sizeof(cl_mem), (void *)&Mesh_Element_Neighborindex3));
		CL_CHECK(clSetKernelArg(kernel4, 10, sizeof(cl_mem), (void *)&Mesh_Element_Neighborindex4));
		CL_CHECK(clSetKernelArg(kernel4, 11, sizeof(cl_mem), (void *)&Mesh_Node_x));
		CL_CHECK(clSetKernelArg(kernel4, 12, sizeof(cl_mem), (void *)&Mesh_Node_y));
		CL_CHECK(clSetKernelArg(kernel4, 13, sizeof(cl_mem), (void *)&Mesh_Node_z));
		CL_CHECK(clSetKernelArg(kernel4, 14, sizeof(cl_mem), (void *)&eid));
		CL_CHECK(clSetKernelArg(kernel4, 15, sizeof(int), (void *)&N1));
		CL_CHECK(clSetKernelArg(kernel4, 16, sizeof(cl_mem), (void *)&r));
		CL_CHECK(clSetKernelArg(kernel4, 17, sizeof(cl_mem), (void *)&s));
		CL_CHECK(clSetKernelArg(kernel4, 18, sizeof(cl_mem), (void *)&t));
		CL_CHECK(clSetKernelArg(kernel4, 19, sizeof(cl_mem), (void *)&integrate));
		
		// Execute the OpenCL kernel on the list
		CL_CHECK(clEnqueueNDRangeKernel(command_queue, kernel4, 1, NULL, &global_item_size, &local_item_size, 0, NULL, NULL));
		
		// release tloc_device
		clReleaseMemObject(tloc_device);
		
		for (int k =0; k < N1; k++) {
		
			t11[k] = t11[k];
			
			tloc_new[k] = (t11[k] - Data_TMin)/ (Data_TMax - Data_TMin);	
		
		}
		// Allocate and memcopy tloc_device
		tloc_device = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, N1 * sizeof(double), tloc_new, &ret);
		if (ret != CL_SUCCESS) {	
			fprintf(stderr, "Error Allocating memory for tloc_device. \n\n");
			exit(1);
		}
		
		
		
		
		//k3//
		stage = 2;
		// Set the arguments of the kernel1 (compute_points_Unstructure3D_1)
		CL_CHECK(clSetKernelArg(kernel1, 0, sizeof(cl_mem), (void *)&posx));
		CL_CHECK(clSetKernelArg(kernel1, 1, sizeof(cl_mem), (void *)&posy));
		CL_CHECK(clSetKernelArg(kernel1, 2, sizeof(cl_mem), (void *)&posz));
		CL_CHECK(clSetKernelArg(kernel1, 3, sizeof(cl_mem), (void *)&tloc_device));
		CL_CHECK(clSetKernelArg(kernel1, 4, sizeof(cl_mem), (void *)&xn0));
		CL_CHECK(clSetKernelArg(kernel1, 5, sizeof(cl_mem), (void *)&xn1));
		CL_CHECK(clSetKernelArg(kernel1, 6, sizeof(cl_mem), (void *)&xn2));
		CL_CHECK(clSetKernelArg(kernel1, 7, sizeof(cl_mem), (void *)&x_dev));
		CL_CHECK(clSetKernelArg(kernel1, 8, sizeof(cl_mem), (void *)&y_dev));
		CL_CHECK(clSetKernelArg(kernel1, 9, sizeof(cl_mem), (void *)&z_dev));
		CL_CHECK(clSetKernelArg(kernel1, 10, sizeof(cl_mem), (void *)&eid));
		CL_CHECK(clSetKernelArg(kernel1, 11, sizeof(cl_mem), (void *)&r));
		CL_CHECK(clSetKernelArg(kernel1, 12, sizeof(cl_mem), (void *)&s));
		CL_CHECK(clSetKernelArg(kernel1, 13, sizeof(cl_mem), (void *)&t));
		CL_CHECK(clSetKernelArg(kernel1, 14, sizeof(cl_mem), (void *)&Vel_U0));
		CL_CHECK(clSetKernelArg(kernel1, 15, sizeof(cl_mem), (void *)&Vel_U1));
		CL_CHECK(clSetKernelArg(kernel1, 16, sizeof(cl_mem), (void *)&Vel_V0));
		CL_CHECK(clSetKernelArg(kernel1, 17, sizeof(cl_mem), (void *)&Vel_V1));
		CL_CHECK(clSetKernelArg(kernel1, 18, sizeof(cl_mem), (void *)&Vel_W0));
		CL_CHECK(clSetKernelArg(kernel1, 19, sizeof(cl_mem), (void *)&Vel_W1));
		CL_CHECK(clSetKernelArg(kernel1, 20, sizeof(cl_mem), (void *)&Mesh_Element_Node1));
		CL_CHECK(clSetKernelArg(kernel1, 21, sizeof(cl_mem), (void *)&Mesh_Element_Node2));
		CL_CHECK(clSetKernelArg(kernel1, 22, sizeof(cl_mem), (void *)&Mesh_Element_Node3));
		CL_CHECK(clSetKernelArg(kernel1, 23, sizeof(cl_mem), (void *)&Mesh_Element_Node4));
		CL_CHECK(clSetKernelArg(kernel1, 24, sizeof(cl_mem), (void *)&Mesh_Element_Neighborindex1));
		CL_CHECK(clSetKernelArg(kernel1, 25, sizeof(cl_mem), (void *)&Mesh_Element_Neighborindex2));
		CL_CHECK(clSetKernelArg(kernel1, 26, sizeof(cl_mem), (void *)&Mesh_Element_Neighborindex3));
		CL_CHECK(clSetKernelArg(kernel1, 27, sizeof(cl_mem), (void *)&Mesh_Element_Neighborindex4));
		CL_CHECK(clSetKernelArg(kernel1, 28, sizeof(cl_mem), (void *)&h_device));
		CL_CHECK(clSetKernelArg(kernel1, 29, sizeof(cl_mem), (void *)&integrate));
		CL_CHECK(clSetKernelArg(kernel1, 30, sizeof(int), (void *)&stage));
		CL_CHECK(clSetKernelArg(kernel1, 31, sizeof(int), (void *)&N1));
		
		// Execute the OpenCL kernel on the list
		
		CL_CHECK(clEnqueueNDRangeKernel(command_queue, kernel1, 1, NULL, &global_item_size, &local_item_size, 0, NULL, NULL));
		
		
		// Do local search .....
		// Do local search .....
		// Set the arguments of the kernel (Initialize_timestep3D)
		CL_CHECK(clSetKernelArg(kernel4, 0, sizeof(cl_mem), (void *)&posx));
		CL_CHECK(clSetKernelArg(kernel4, 1, sizeof(cl_mem), (void *)&posy));
		CL_CHECK(clSetKernelArg(kernel4, 2, sizeof(cl_mem), (void *)&posz));
		CL_CHECK(clSetKernelArg(kernel4, 3, sizeof(cl_mem), (void *)&Mesh_Element_Node1));
		CL_CHECK(clSetKernelArg(kernel4, 4, sizeof(cl_mem), (void *)&Mesh_Element_Node2));
		CL_CHECK(clSetKernelArg(kernel4, 5, sizeof(cl_mem), (void *)&Mesh_Element_Node3));
		CL_CHECK(clSetKernelArg(kernel4, 6, sizeof(cl_mem), (void *)&Mesh_Element_Node4));
		CL_CHECK(clSetKernelArg(kernel4, 7, sizeof(cl_mem), (void *)&Mesh_Element_Neighborindex1));
		CL_CHECK(clSetKernelArg(kernel4, 8, sizeof(cl_mem), (void *)&Mesh_Element_Neighborindex2));
		CL_CHECK(clSetKernelArg(kernel4, 9, sizeof(cl_mem), (void *)&Mesh_Element_Neighborindex3));
		CL_CHECK(clSetKernelArg(kernel4, 10, sizeof(cl_mem), (void *)&Mesh_Element_Neighborindex4));
		CL_CHECK(clSetKernelArg(kernel4, 11, sizeof(cl_mem), (void *)&Mesh_Node_x));
		CL_CHECK(clSetKernelArg(kernel4, 12, sizeof(cl_mem), (void *)&Mesh_Node_y));
		CL_CHECK(clSetKernelArg(kernel4, 13, sizeof(cl_mem), (void *)&Mesh_Node_z));
		CL_CHECK(clSetKernelArg(kernel4, 14, sizeof(cl_mem), (void *)&eid));
		CL_CHECK(clSetKernelArg(kernel4, 15, sizeof(int), (void *)&N1));
		CL_CHECK(clSetKernelArg(kernel4, 16, sizeof(cl_mem), (void *)&r));
		CL_CHECK(clSetKernelArg(kernel4, 17, sizeof(cl_mem), (void *)&s));
		CL_CHECK(clSetKernelArg(kernel4, 18, sizeof(cl_mem), (void *)&t));
		CL_CHECK(clSetKernelArg(kernel4, 19, sizeof(cl_mem), (void *)&integrate));
		
		// Execute the OpenCL kernel on the list
		CL_CHECK(clEnqueueNDRangeKernel(command_queue, kernel4, 1, NULL, &global_item_size, &local_item_size, 0, NULL, NULL));
		// release tloc_device
		clReleaseMemObject(tloc_device);
		for (int k =0; k < N1; k++) {
			
			t11[k] = t11[k] + 0.5*h_new[k];
			
			tloc_new[k] = (t11[k] - Data_TMin)/ (Data_TMax - Data_TMin);	
				
		}
		
		// Allocate and memcopy tloc_device
		tloc_device = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, N1 * sizeof(double), tloc_new, &ret);
		if (ret != CL_SUCCESS) {	
			fprintf(stderr, "Error Allocating memory for tloc_device. \n\n");
			exit(1);
		}
		
		
		
		//k4//
		// Set the arguments of the kernel5 (compute_points_Unstructure3D_2)
		CL_CHECK(clSetKernelArg(kernel5, 0, sizeof(cl_mem), (void *)&posx));
		CL_CHECK(clSetKernelArg(kernel5, 1, sizeof(cl_mem), (void *)&posy));
		CL_CHECK(clSetKernelArg(kernel5, 2, sizeof(cl_mem), (void *)&posz));
		CL_CHECK(clSetKernelArg(kernel5, 3, sizeof(cl_mem), (void *)&tloc_device));
		CL_CHECK(clSetKernelArg(kernel5, 4, sizeof(cl_mem), (void *)&xn0));
		CL_CHECK(clSetKernelArg(kernel5, 5, sizeof(cl_mem), (void *)&xn1));
		CL_CHECK(clSetKernelArg(kernel5, 6, sizeof(cl_mem), (void *)&xn2));
		CL_CHECK(clSetKernelArg(kernel5, 7, sizeof(cl_mem), (void *)&x_dev));
		CL_CHECK(clSetKernelArg(kernel5, 8, sizeof(cl_mem), (void *)&y_dev));
		CL_CHECK(clSetKernelArg(kernel5, 9, sizeof(cl_mem), (void *)&z_dev));
		CL_CHECK(clSetKernelArg(kernel5, 10, sizeof(cl_mem), (void *)&eid));
		CL_CHECK(clSetKernelArg(kernel5, 11, sizeof(cl_mem), (void *)&r));
		CL_CHECK(clSetKernelArg(kernel5, 12, sizeof(cl_mem), (void *)&s));
		CL_CHECK(clSetKernelArg(kernel5, 13, sizeof(cl_mem), (void *)&t));
		CL_CHECK(clSetKernelArg(kernel5, 14, sizeof(cl_mem), (void *)&Vel_U0));
		CL_CHECK(clSetKernelArg(kernel5, 15, sizeof(cl_mem), (void *)&Vel_U1));
		CL_CHECK(clSetKernelArg(kernel5, 16, sizeof(cl_mem), (void *)&Vel_V0));
		CL_CHECK(clSetKernelArg(kernel5, 17, sizeof(cl_mem), (void *)&Vel_V1));
		CL_CHECK(clSetKernelArg(kernel5, 18, sizeof(cl_mem), (void *)&Vel_W0));
		CL_CHECK(clSetKernelArg(kernel5, 19, sizeof(cl_mem), (void *)&Vel_W1));
		CL_CHECK(clSetKernelArg(kernel5, 20, sizeof(cl_mem), (void *)&Mesh_Element_Node1));
		CL_CHECK(clSetKernelArg(kernel5, 21, sizeof(cl_mem), (void *)&Mesh_Element_Node2));
		CL_CHECK(clSetKernelArg(kernel5, 22, sizeof(cl_mem), (void *)&Mesh_Element_Node3));
		CL_CHECK(clSetKernelArg(kernel5, 23, sizeof(cl_mem), (void *)&Mesh_Element_Node4));
		CL_CHECK(clSetKernelArg(kernel5, 24, sizeof(cl_mem), (void *)&Mesh_Element_Neighborindex1));
		CL_CHECK(clSetKernelArg(kernel5, 25, sizeof(cl_mem), (void *)&Mesh_Element_Neighborindex2));
		CL_CHECK(clSetKernelArg(kernel5, 26, sizeof(cl_mem), (void *)&Mesh_Element_Neighborindex3));
		CL_CHECK(clSetKernelArg(kernel5, 27, sizeof(cl_mem), (void *)&Mesh_Element_Neighborindex4));
		CL_CHECK(clSetKernelArg(kernel5, 28, sizeof(cl_mem), (void *)&h_device));
		CL_CHECK(clSetKernelArg(kernel5, 29, sizeof(cl_mem), (void *)&integrate));
		CL_CHECK(clSetKernelArg(kernel5, 30, sizeof(int), (void *)&N1));
		
		// Execute the OpenCL kernel on the list
		CL_CHECK(clEnqueueNDRangeKernel(command_queue, kernel5, 1, NULL, &global_item_size, &local_item_size, 0, NULL, NULL));
		
		
		// Do local search .....
		// Do local search .....
		// Set the arguments of the kernel (Initialize_timestep3D)
		CL_CHECK(clSetKernelArg(kernel4, 0, sizeof(cl_mem), (void *)&posx));
		CL_CHECK(clSetKernelArg(kernel4, 1, sizeof(cl_mem), (void *)&posy));
		CL_CHECK(clSetKernelArg(kernel4, 2, sizeof(cl_mem), (void *)&posz));
		CL_CHECK(clSetKernelArg(kernel4, 3, sizeof(cl_mem), (void *)&Mesh_Element_Node1));
		CL_CHECK(clSetKernelArg(kernel4, 4, sizeof(cl_mem), (void *)&Mesh_Element_Node2));
		CL_CHECK(clSetKernelArg(kernel4, 5, sizeof(cl_mem), (void *)&Mesh_Element_Node3));
		CL_CHECK(clSetKernelArg(kernel4, 6, sizeof(cl_mem), (void *)&Mesh_Element_Node4));
		CL_CHECK(clSetKernelArg(kernel4, 7, sizeof(cl_mem), (void *)&Mesh_Element_Neighborindex1));
		CL_CHECK(clSetKernelArg(kernel4, 8, sizeof(cl_mem), (void *)&Mesh_Element_Neighborindex2));
		CL_CHECK(clSetKernelArg(kernel4, 9, sizeof(cl_mem), (void *)&Mesh_Element_Neighborindex3));
		CL_CHECK(clSetKernelArg(kernel4, 10, sizeof(cl_mem), (void *)&Mesh_Element_Neighborindex4));
		CL_CHECK(clSetKernelArg(kernel4, 11, sizeof(cl_mem), (void *)&Mesh_Node_x));
		CL_CHECK(clSetKernelArg(kernel4, 12, sizeof(cl_mem), (void *)&Mesh_Node_y));
		CL_CHECK(clSetKernelArg(kernel4, 13, sizeof(cl_mem), (void *)&Mesh_Node_z));
		CL_CHECK(clSetKernelArg(kernel4, 14, sizeof(cl_mem), (void *)&eid));
		CL_CHECK(clSetKernelArg(kernel4, 15, sizeof(int), (void *)&N1));
		CL_CHECK(clSetKernelArg(kernel4, 16, sizeof(cl_mem), (void *)&r));
		CL_CHECK(clSetKernelArg(kernel4, 17, sizeof(cl_mem), (void *)&s));
		CL_CHECK(clSetKernelArg(kernel4, 18, sizeof(cl_mem), (void *)&t));
		CL_CHECK(clSetKernelArg(kernel4, 19, sizeof(cl_mem), (void *)&integrate));
		
		// Execute the OpenCL kernel on the list
		CL_CHECK(clEnqueueNDRangeKernel(command_queue, kernel4, 1, NULL, &global_item_size, &local_item_size, 0, NULL, NULL));
		
		// release tloc_device
		clReleaseMemObject(tloc_device);
		
		// Copy value of eid into Tracer_dev.eid
		
		
		
		
		// Release device constant memory objects 
		clReleaseMemObject(t11_device);
		clReleaseMemObject(t22_device);
		clReleaseMemObject(h_device);
	} // End of For loop
	
}

