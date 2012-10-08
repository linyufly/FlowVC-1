# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <CL/cl.h>



# include "memcopy.h"
# include "globals.h"

# include "settings.h"


void host2device_const(void) {

	if (Data_MeshType == UNSTRUCTURED) {
	
	} else {
		
		double XYZ_host[9];
		int Res_host[3];
		
		// Memcopy Velocity mesh data to constant memory on device..
			XYZ_host[0] = Vel_CartMesh.XMin;
			XYZ_host[1] = Vel_CartMesh.XMax;
			XYZ_host[2] = (1/Vel_CartMesh.XDelta);
			
			XYZ_host[3] = Vel_CartMesh.YMin;
			XYZ_host[4] = Vel_CartMesh.YMax;
			XYZ_host[5] = (1/Vel_CartMesh.YDelta);
	
			XYZ_host[6] = Vel_CartMesh.ZMin;
			XYZ_host[7] = Vel_CartMesh.ZMax;
			XYZ_host[8] = (1/Vel_CartMesh.ZDelta);
			
			Res_host[0] = Vel_CartMesh.XRes;
			Res_host[1] = Vel_CartMesh.YRes;
			Res_host[2] = Vel_CartMesh.ZRes;
			
			XXYYZZ = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, 9 * sizeof(double), XYZ_host, &ret);
			if (ret != CL_SUCCESS) {	
				fprintf(stderr, "Error Allocating memory for XXYYZZ. \n\n");
				exit(1);
			}
			
			RES = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, 3 * sizeof(int), Res_host, &ret);
			if (ret != CL_SUCCESS) {	
				fprintf(stderr, "Error Allocating memory for RES. \n\n");
				exit(1);
			}
			
			XYZ_host[0] = Data_MeshBounds.XMin ;
			XYZ_host[1] = Data_MeshBounds.XMax;
			XYZ_host[2] = Data_MeshBounds.XDelta;
			XYZ_host[3] = Data_MeshBounds.YMin ;
			XYZ_host[4] = Data_MeshBounds.YMax;
			XYZ_host[5] = Data_MeshBounds.YDelta;
			
			XYZ_host[6] = Data_MeshBounds.ZMin ;
			XYZ_host[7] = Data_MeshBounds.ZMax;
			XYZ_host[8] = Data_MeshBounds.ZDelta;
		
			Res_host[0] = Data_MeshBounds.XRes;
			Res_host[1] = Data_MeshBounds.YRes;
			Res_host[2] = Data_MeshBounds.ZRes;
	
			XXYYZZ_DATA = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, 9 * sizeof(double), XYZ_host, &ret);
			if (ret != CL_SUCCESS) {	
				fprintf(stderr, "Error Allocating memory for XXYYZZ_DATA. \n\n");
				exit(1);
			}
			
			RES_DATA = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, 3 * sizeof(int), Res_host, &ret);
			if (ret != CL_SUCCESS) {	
				fprintf(stderr, "Error Allocating memory for RES_DATA. \n\n");
				exit(1);
			}
	
	
	
	}
}



void host2device(int ss) {

	CL_CHECK(clEnqueueWriteBuffer(command_queue, x_dev, CL_TRUE, 0, ss*sizeof(double), Tracer.x , 0, NULL, NULL));
	CL_CHECK(clEnqueueWriteBuffer(command_queue, y_dev, CL_TRUE, 0, ss*sizeof(double), Tracer.y , 0, NULL, NULL));
	if (Dimensions == 3)
		CL_CHECK(clEnqueueWriteBuffer(command_queue, z_dev, CL_TRUE, 0, ss*sizeof(double), Tracer.z , 0, NULL, NULL));
	
	if (Data_MeshType == UNSTRUCTURED) 
		CL_CHECK(clEnqueueWriteBuffer(command_queue, ElementIndex_dev, CL_TRUE, 0, ss*sizeof(int), Tracer.ElementIndex , 0, NULL, NULL));
	
	CL_CHECK(clEnqueueWriteBuffer(command_queue, LeftDomain_dev, CL_TRUE, 0, ss*sizeof(int), Tracer.LeftDomain , 0, NULL, NULL));
	
	CL_CHECK(clEnqueueWriteBuffer(command_queue, Start_time_dev, CL_TRUE, 0, ss*sizeof(double), Tracer.Start_time , 0, NULL, NULL));
	CL_CHECK(clEnqueueWriteBuffer(command_queue, Stop_time_dev, CL_TRUE, 0, ss*sizeof(double), Tracer.Stop_time , 0, NULL, NULL));
	
	int num;

		if (Data_MeshType == CARTESIAN)
			num = N_Vel;
		else
			num = Vel_MeshNumNodes;
			
	
	CL_CHECK(clEnqueueWriteBuffer(command_queue, Vel_U0, CL_TRUE, 0, num*sizeof(double), velocity.u0 , 0, NULL, NULL));
	CL_CHECK(clEnqueueWriteBuffer(command_queue, Vel_U1, CL_TRUE, 0, num*sizeof(double), velocity.u1 , 0, NULL, NULL));
	CL_CHECK(clEnqueueWriteBuffer(command_queue, Vel_V0, CL_TRUE, 0, num*sizeof(double), velocity.v0 , 0, NULL, NULL));
	CL_CHECK(clEnqueueWriteBuffer(command_queue, Vel_V1, CL_TRUE, 0, num*sizeof(double), velocity.v1 , 0, NULL, NULL));
	CL_CHECK(clEnqueueWriteBuffer(command_queue, Vel_W0, CL_TRUE, 0, num*sizeof(double), velocity.w0 , 0, NULL, NULL));
	CL_CHECK(clEnqueueWriteBuffer(command_queue, Vel_W1, CL_TRUE, 0, num*sizeof(double), velocity.w1 , 0, NULL, NULL));
	
}



void device2host(int ss) {

	CL_CHECK(clEnqueueReadBuffer(command_queue, x_dev, CL_TRUE, 0, ss*sizeof(double), Tracer.x, 0, NULL, NULL)); 
	CL_CHECK(clEnqueueReadBuffer(command_queue, y_dev, CL_TRUE, 0, ss*sizeof(double), Tracer.y, 0, NULL, NULL)); 
	
	if (Dimensions == 3)
		CL_CHECK(clEnqueueReadBuffer(command_queue, z_dev, CL_TRUE, 0, ss*sizeof(double), Tracer.z, 0, NULL, NULL)); 

	CL_CHECK(clEnqueueReadBuffer(command_queue, ElementIndex_dev, CL_TRUE, 0, ss*sizeof(int), Tracer.ElementIndex, 0, NULL, NULL)); 

	CL_CHECK(clEnqueueReadBuffer(command_queue, LeftDomain_dev, CL_TRUE, 0, ss*sizeof(int), Tracer.LeftDomain, 0, NULL, NULL)); 
}


















