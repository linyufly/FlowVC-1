
# include "clean.h"

void tempclean(Point temp) {

	free(temp.x);
	free(temp.y);
	free(temp.z);
	free(temp.ElementIndex);
	free(temp.Start_time);
	free(temp.Stop_time);
	free(temp.LeftDomain); 


}

void clean_gs(void) {

	cudaFree(MeshNodeArray_device.x);
	cudaFree(MeshNodeArray_device.y);
	cudaFree(MeshNodeArray_device.z);

	if (Memory_Usage_GS == Use_Pinned_Memory) {
	
		cudaFreeHost(MeshNodeArray.x);
		cudaFreeHost(MeshNodeArray.y);
		cudaFreeHost(MeshNodeArray.z);
	} else {
		
		free(MeshNodeArray.x);
		free(MeshNodeArray.y);
		free(MeshNodeArray.z);
	
	}
	
}


void clean_all(void) {


		
		// Release all memory allocated

		if (Memory_Usage_GS == Use_Pinned_Memory) {

			cudaFreeHost(MeshElementArray.Node1);
			cudaFreeHost(MeshElementArray.Node2);
			cudaFreeHost(MeshElementArray.Node3);
			cudaFreeHost(MeshElementArray.Node4);

		//	cudaFreeHost(MeshNodeArray.x);
		//	cudaFreeHost(MeshNodeArray.y);
		//	cudaFreeHost(MeshNodeArray.z);
			
			cudaFreeHost(MeshNodeArray_double.x);
			cudaFreeHost(MeshNodeArray_double.y);
			cudaFreeHost(MeshNodeArray_double.z);

		} else {

			free(MeshElementArray.Node1);
			free(MeshElementArray.Node2);
			free(MeshElementArray.Node3);
			free(MeshElementArray.Node4);

		//	free(MeshNodeArray.x);
		//	free(MeshNodeArray.y);
		//	free(MeshNodeArray.z);
			
			free(MeshNodeArray_double.x);
			free(MeshNodeArray_double.y);
			free(MeshNodeArray_double.z);
			
			
		}
		
		cudaFree(MeshElementArray_device.Node1);
		cudaFree(MeshElementArray_device.Node2);
		cudaFree(MeshElementArray_device.Node3);
		cudaFree(MeshElementArray_device.Node4);

//		cudaFree(MeshNodeArray_device.x);
//		cudaFree(MeshNodeArray_device.y);
//		cudaFree(MeshNodeArray_device.z);
		
		cudaFree(MeshNodeArray_double_device.x);
		cudaFree(MeshNodeArray_double_device.y);
		cudaFree(MeshNodeArray_double_device.z);




		free(MeshElementArray.Neighborindex1);
		free(MeshElementArray.Neighborindex2);
		free(MeshElementArray.Neighborindex3);
		free(MeshElementArray.Neighborindex4);
		

		// Cleaning Velocity variables
		
		if (Memory_Usage_Tracer == Use_Pinned_Memory) {
	
			cudaFreeHost(velocity.u0);
			cudaFreeHost(velocity.v0);
			cudaFreeHost(velocity.w0);
			cudaFreeHost(velocity.u1);
			cudaFreeHost(velocity.v1);
			cudaFreeHost(velocity.w1);
			cudaFreeHost(velocity.time0);
			cudaFreeHost(velocity.time1);
			
			cudaFreeHost(Tracer.x);
			cudaFreeHost(Tracer.y);
			
			//if(Dimensions == 3)
				cudaFreeHost(Tracer.z);
				
			cudaFreeHost(Tracer.ElementIndex);
			cudaFreeHost(Tracer.LeftDomain);
			cudaFreeHost(Tracer.Start_time);
			cudaFreeHost(Tracer.Stop_time);
			if (Trace_ReleaseStrategy == 1) {
				cudaFreeHost(Tracer1.x);
				cudaFreeHost(Tracer1.y);
			
				//if(Dimensions == 3)
					cudaFreeHost(Tracer1.z);
				
				cudaFreeHost(Tracer1.ElementIndex);
				cudaFreeHost(Tracer1.LeftDomain);
				cudaFreeHost(Tracer1.Start_time);
				cudaFreeHost(Tracer1.Stop_time);
			
				cudaFreeHost(index1);
				cudaFreeHost(Tracer.Status);
			}
			cudaFreeHost(DataTime1);
			cudaFreeHost(Output_time);
			cudaFreeHost(Launch_time);
		
		} else {
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
		}
		
		cudaFree(velocity_dev.u0);
		cudaFree(velocity_dev.v0);
		cudaFree(velocity_dev.w0);
		cudaFree(velocity_dev.u1);
		cudaFree(velocity_dev.v1);
		cudaFree(velocity_dev.w1);
		

		cudaFree(Tracer_dev.x);
		cudaFree(Tracer_dev.y);
		cudaFree(Tracer_dev.ElementIndex);
		cudaFree(Tracer_dev.Start_time);
		cudaFree(Tracer_dev.Stop_time);
		
		
		free(x_host);
		free(y_host);
		
		if(Dimensions == 3) {
			cudaFree(Tracer_dev.z);
			free(z_host);
		}
		
	
		
		//cudaFree(tempx);
		//cudaFree(tempy);
		cudaFree(posx);
		cudaFree(posy);
		
		
		cudaFree(xn0);
		cudaFree(xn1);
		
		if (Data_MeshType == UNSTRUCTURED) {
		//	cudaFree(eid);
			cudaFree(r);
			cudaFree(s);
			cudaFree(t);
		
		}
		if (Dimensions == 3) {
			cudaFree(posz);
			cudaFree(xn2);
		
		}
		
		cudaFree(integrate);
		
		// Remove Temp file containing tracer release information
		if (!Keep_Tempfile) {
			char BinFile[LONGSTRING];
			sprintf(BinFile, "%s%s.bin", Path_Output, Temp_OutFilePrefix);
			if(remove(BinFile))
					fprintf(stderr, "Warning: Could not delete file %s\n", BinFile);
		}

}
