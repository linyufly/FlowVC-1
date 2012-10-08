
#include "openclfuncs.h"

#include <stdlib.h>
#include <stdio.h>
#include <CL/cl.h>

#include "settings.h"
#include "globals.h"





void initopencl(void) {
	
	int i;


	// Get Platform and Device Info
	CL_CHECK(clGetPlatformIDs(1, &platform_id, &num_platforms));
	
	// Currently this program only runs on a SINGLE GPU.
	CL_CHECK(clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_GPU, 1, &device_id, &num_devices));
	printf("=== %d OpenCL platform(s) found: ===\n", num_platforms);
	printf("=== %d OpenCL device(s) found on platform:\n", num_devices);
	
	
	char buffer[10240];
	cl_uint buf_uint;
	cl_ulong buf_ulong;
	printf("  -- %d --\n", i);
	CL_CHECK(clGetDeviceInfo(device_id, CL_DEVICE_NAME, sizeof(buffer), buffer, NULL));
	printf("  DEVICE_NAME = %s\n", buffer);
	CL_CHECK(clGetDeviceInfo(device_id, CL_DEVICE_VENDOR, sizeof(buffer), buffer, NULL));
	printf("  DEVICE_VENDOR = %s\n", buffer);
	CL_CHECK(clGetDeviceInfo(device_id, CL_DEVICE_VERSION, sizeof(buffer), buffer, NULL));
	printf("  DEVICE_VERSION = %s\n", buffer);
	CL_CHECK(clGetDeviceInfo(device_id, CL_DRIVER_VERSION, sizeof(buffer), buffer, NULL));
	printf("  DRIVER_VERSION = %s\n", buffer);
	CL_CHECK(clGetDeviceInfo(device_id, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(buf_uint), &buf_uint, NULL));
	printf("  DEVICE_MAX_COMPUTE_UNITS = %u\n", (unsigned int)buf_uint);
	CL_CHECK(clGetDeviceInfo(device_id, CL_DEVICE_MAX_CLOCK_FREQUENCY, sizeof(buf_uint), &buf_uint, NULL));
	printf("  DEVICE_MAX_CLOCK_FREQUENCY = %u\n", (unsigned int)buf_uint);
	CL_CHECK(clGetDeviceInfo(device_id, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(buf_ulong), &buf_ulong, NULL));
	printf("  DEVICE_GLOBAL_MEM_SIZE = %llu\n", (unsigned long long)buf_ulong);

	if (num_devices == 0)
	{	
		fprintf(stderr, "No Devices found that can run OpenCL.");
		exit(0);	
	}
	// Create OpenCL context
	context = clCreateContext(NULL, 1, &device_id, NULL, NULL, &ret);
	if (ret != CL_SUCCESS) {
		
		fprintf(stderr, "Error creating context: Function returned %d \n\n", ret);
		exit(1);
	
	}
	// Create Command Queue
	command_queue = clCreateCommandQueue(context, device_id, 0, &ret);
	if (ret != CL_SUCCESS) {
		
		fprintf(stderr, "Error creating command Queue: Function returned %d \n\n", ret);
		exit(1);
	
	}
	
	// Load the kernel source code into the array source_str
	FILE *fp;
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
	
	
	// Create a program from the kernel source
    program = clCreateProgramWithSource(context, 1, (const char **)&source_str, (const size_t *)&source_size, &ret);	
    if (ret != CL_SUCCESS) {	
			fprintf(stderr, "Error creating a program for integration3D. %d \n\n", (int)ret);
			exit(1);
	}
    // Build the program
    
    ret = clBuildProgram(program, 1, &device_id, "-DUSE_DOUBLE=1", NULL, NULL); 
    if (ret != CL_SUCCESS)
    {
    
    	size_t length;
    	char buffer[10240];
    	clGetProgramBuildInfo(program, 1, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &length);
    	fprintf(stderr, "Error returned %d. \n\n", (int)ret);
    	printf("Error Log: \n\n %s \n\n", buffer);
    	exit(0);
    }
	
/*    // Create the OpenCL kernel (compute_points_Unstructure3D_1)
    kernel1 = clCreateKernel(program, "compute_points_Unstructure3D_1", &ret);
	if (ret != CL_SUCCESS) {	
			fprintf(stderr, "Error creating a kernel for compute_points_Unstructure3D_1. \n\n");
			exit(1);
	}
*/	
	// Create the OpenCL kernel (check_int)
    kernel2 = clCreateKernel(program, "check_int", &ret);
	if (ret != CL_SUCCESS) {	
			fprintf(stderr, "Error creating a kernel for check_int. %d \n\n", (int)ret);
			exit(1);
	}
	
    
    // Create the OpenCL kernel (compute_points_Unstructure3D_1)
    kernel1 = clCreateKernel(program, "compute_points_Unstructure3D_1", &ret);
	if (ret != CL_SUCCESS) {	
			fprintf(stderr, "Error creating a kernel for compute_points_Unstructure3D_1. \n\n");
			exit(1);
	}
	
	// Create the OpenCL kernel (initialize_timestep3D)
	kernel3 = clCreateKernel(program, "initialize_timestep3D", &ret);
	if (ret != CL_SUCCESS) {	
			fprintf(stderr, "Error creating a kernel for initialize_timestep3D. \n\n");
			exit(1);
	}
	
	// Create the OpenCL kernel (initialize_timestep3D)
    kernel4 = clCreateKernel(program, "LocalSearch3D", &ret);
	if (ret != CL_SUCCESS) {	
			fprintf(stderr, "Error creating a kernel for LocalSearch3D. \n\n");
			exit(1);
	}
	
	// Create the OpenCL kernel (initialize_timestep3D)
    kernel5 = clCreateKernel(program, "compute_points_Unstructure3D_2", &ret);
	if (ret != CL_SUCCESS) {	
			fprintf(stderr, "Error creating a kernel for LocalSearch3D. \n\n");
			exit(1);
	}
	
	
	
	printf("\n\n");
}




