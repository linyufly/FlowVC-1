#include <stdlib.h>
#include <stdio.h>
#include <CL/cl.h>


# include "clean.h"
#include "structs.h"
#include "globals.h"

void clean_all(void) {

		printf("Cleaning Variables ... \n\n");
		
		// Opencl environment variables
		clReleaseCommandQueue(command_queue);
		clReleaseContext(context);
		
		
		// Release all memory allocated

		
		if (Data_MeshType == UNSTRUCTURED) {
	
			// Mesh Variables
			free(MeshElementArray.Node1);
			free(MeshElementArray.Node2);
			free(MeshElementArray.Node3);
			free(MeshElementArray.Node4);
			
			free(MeshNodeArray_double.x);
			free(MeshNodeArray_double.y);
			free(MeshNodeArray_double.z);


			free(MeshElementArray.Neighborindex1);
			free(MeshElementArray.Neighborindex2);
			free(MeshElementArray.Neighborindex3);
			free(MeshElementArray.Neighborindex4);
		
			clReleaseMemObject(Mesh_Node_x);
			clReleaseMemObject(Mesh_Node_y);
			clReleaseMemObject(Mesh_Node_z);
			
			clReleaseMemObject(Mesh_Element_Node1);
			clReleaseMemObject(Mesh_Element_Node2);
			clReleaseMemObject(Mesh_Element_Node3);
			clReleaseMemObject(Mesh_Element_Node4);
			
			clReleaseMemObject(Mesh_Element_Neighborindex1);
			clReleaseMemObject(Mesh_Element_Neighborindex2);
			clReleaseMemObject(Mesh_Element_Neighborindex3);
			clReleaseMemObject(Mesh_Element_Neighborindex4);
			
			clReleaseMemObject(r);
			clReleaseMemObject(s);
			clReleaseMemObject(t);
			clReleaseMemObject(eid);
			
		}

		// Cleaning Velocity variables
		
			free(velocity.u0);
			free(velocity.v0);
			free(velocity.w0);
			free(velocity.u1);
			free(velocity.v1);
			free(velocity.w1);
			free(velocity.time0);
			free(velocity.time1);
			
			free(Tracer.x);
			Tracer.x = NULL;
			free(Tracer.y);
			Tracer.y = NULL;
			
		
			free(Tracer.z);
			Tracer.z = NULL;
		
			free(Tracer.ElementIndex);
			Tracer.ElementIndex = NULL;
			free(Tracer.Start_time);
			Tracer.Start_time = NULL;
			free(Tracer.Stop_time);
			Tracer.Stop_time = NULL;
			free(Tracer.LeftDomain);
			Tracer.LeftDomain = NULL;
			
			if (Trace_ReleaseStrategy == 1) {
				free(Tracer1.x);
				Tracer1.x = NULL;
				free(Tracer1.y);
				Tracer1.y = NULL;
			
		
				free(Tracer1.z);
				Tracer1.z = NULL;
		
				free(Tracer1.ElementIndex);
				Tracer1.ElementIndex = NULL;
				free(Tracer1.Start_time);
				Tracer1.Start_time = NULL;
				free(Tracer1.Stop_time);
				Tracer1.Stop_time = NULL;
				free(Tracer1.LeftDomain);
				Tracer1.LeftDomain = NULL;
			
				free(index1);
				index1 = NULL;
				
				free(Tracer.Status);
				Tracer.Status = NULL;
			}
			free(DataTime1);
			free(Output_time);
			free(Launch_time);
			
			
		clReleaseMemObject(Vel_U0);
		clReleaseMemObject(Vel_U1);
		clReleaseMemObject(Vel_V0);
		clReleaseMemObject(Vel_V1);
		clReleaseMemObject(Vel_W0);
		clReleaseMemObject(Vel_W1);
		
		clReleaseMemObject(x_dev);
		clReleaseMemObject(y_dev);
		
		clReleaseMemObject(posx);
		clReleaseMemObject(posy);
		clReleaseMemObject(xn0);
		clReleaseMemObject(xn1);
		clReleaseMemObject(integrate);
		
		if (Dimensions == 3) {
	
			clReleaseMemObject(z_dev);
			clReleaseMemObject(posz);
			clReleaseMemObject(xn2);
		}
	
		clReleaseMemObject(Start_time_dev);
		clReleaseMemObject(Stop_time_dev);
		
		clReleaseMemObject(ElementIndex_dev);
		clReleaseMemObject(LeftDomain_dev);
		
		// Remove Temp file containing tracer release information
		if (!Keep_Tempfile) {
			char BinFile[LONGSTRING];
			sprintf(BinFile, "%s%s.bin", Path_Output, Temp_OutFilePrefix);
			if(remove(BinFile))
					fprintf(stderr, "Warning: Could not delete file %s\n", BinFile);
		}
		
		CL_CHECK(clReleaseKernel(kernel1));
		CL_CHECK(clReleaseKernel(kernel2));
		CL_CHECK(clReleaseKernel(kernel3));
		CL_CHECK(clReleaseKernel(kernel4));
		CL_CHECK(clReleaseKernel(kernel5));
	
   		CL_CHECK(clReleaseProgram(program));
		printf("Cleaning Successfull \n\n");

}
