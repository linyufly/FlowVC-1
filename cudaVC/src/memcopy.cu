
# include "memcopy.h"



void host2device_const(void) {

if (Data_MeshType == UNSTRUCTURED) {

err = cudaMemcpy( MeshElementArray_device.Node1, MeshElementArray.Node1 , Vel_MeshNumElements * sizeof(int) , cudaMemcpyHostToDevice);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Host structure array MeshElementArray.Node1 to Device structure array MeshElementArray_device failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
exit(1);
}	

err = cudaMemcpy( MeshElementArray_device.Node2 , MeshElementArray.Node2, Vel_MeshNumElements * sizeof(int) , cudaMemcpyHostToDevice);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Host structure array MeshElementArray.Node2 to Device structure array MeshElementArray_device failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
exit(1);
}	

err = cudaMemcpy( MeshElementArray_device.Node3 , MeshElementArray.Node3, Vel_MeshNumElements * sizeof(int) , cudaMemcpyHostToDevice);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Host structure array MeshElementArray.Node3 to Device structure array MeshElementArray_device failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
exit(1);
}	

err = cudaMemcpy( MeshElementArray_device.Node4 , MeshElementArray.Node4, Vel_MeshNumElements * sizeof(int) , cudaMemcpyHostToDevice);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Host structure array MeshElementArray.Node4 to Device structure array MeshElementArray_device failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
exit(1);
}	

err = cudaMemcpy( MeshElementArray_device.Neighborindex1 , MeshElementArray.Neighborindex1, Vel_MeshNumElements * sizeof(int) , cudaMemcpyHostToDevice);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Host structure array MeshElementArray.Neighborindex1 to Device structure array MeshElementArray_device failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
exit(1);
}	

err = cudaMemcpy(MeshElementArray_device.Neighborindex2 ,MeshElementArray.Neighborindex2, Vel_MeshNumElements * sizeof(int) , cudaMemcpyHostToDevice);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Host structure array MeshElementArray.Neighborindex2 to Device structure array MeshElementArray_device failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
exit(1);
}	

err = cudaMemcpy( MeshElementArray_device.Neighborindex3 , MeshElementArray.Neighborindex3, Vel_MeshNumElements * sizeof(int) , cudaMemcpyHostToDevice);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Host structure array MeshElementArray.Neighborindex3 to Device structure array MeshElementArray_device failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
exit(1);
}	

err = cudaMemcpy( MeshElementArray_device.Neighborindex4 , MeshElementArray.Neighborindex4, Vel_MeshNumElements * sizeof(int) , cudaMemcpyHostToDevice);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Host structure array MeshElementArray.Neighborindex4 to Device structure array MeshElementArray_device failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
exit(1);
}	


err = cudaMemcpy( MeshNodeArray_double_device.x , MeshNodeArray_double.x, Vel_MeshNumNodes * sizeof(double) , cudaMemcpyHostToDevice);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Host structure array MeshNodeArray_double.x to Device structure array MeshNodeArray_double_device failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
exit(1);
}	

err = cudaMemcpy( MeshNodeArray_double_device.y , MeshNodeArray_double.y, Vel_MeshNumNodes * sizeof(double) , cudaMemcpyHostToDevice);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Host structure array MeshNodeArray_double.y to Device structure array MeshNodeArray_double_device failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
exit(1);
}	

err = cudaMemcpy( MeshNodeArray_double_device.z , MeshNodeArray_double.z, Vel_MeshNumNodes * sizeof(double) , cudaMemcpyHostToDevice);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Host structure array MeshNodeArray_double.z to Device structure array MeshNodeArray_double_device failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
exit(1);
}


//printf("host2device_const completed ......\n\n\n");

} else {
double XX_host[3];

// Memcopy Velocity mesh data to constant memory on device..
XX_host[0] = Vel_CartMesh.XMin;
XX_host[1] = Vel_CartMesh.XMax;
XX_host[2] = (1/Vel_CartMesh.XDelta);

// Memcpy of constant memory
err = cudaMemcpyToSymbol( XX, &XX_host, 3 * sizeof(double), 0, cudaMemcpyHostToDevice);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Host to Device constant XX failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
fprintf(stderr, " \n\n");
exit(1);
}

XX_host[0] = Vel_CartMesh.YMin;
XX_host[1] = Vel_CartMesh.YMax;
XX_host[2] = (1/Vel_CartMesh.YDelta);

// Memcpy of constant memory
err = cudaMemcpyToSymbol( YY, &XX_host, 3 * sizeof(double), 0, cudaMemcpyHostToDevice);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Host to Device constant YY failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
fprintf(stderr, " \n\n");
exit(1);
}

XX_host[0] = Vel_CartMesh.ZMin;
XX_host[1] = Vel_CartMesh.ZMax;
XX_host[2] = (1/Vel_CartMesh.ZDelta);

// Memcpy of constant memory
err = cudaMemcpyToSymbol( ZZ, &XX_host, 3 * sizeof(double), 0, cudaMemcpyHostToDevice);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Host to Device constant ZZ failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
fprintf(stderr, " \n\n");
exit(1);
}

int Res_host[3];

Res_host[0] = Vel_CartMesh.XRes;
Res_host[1] = Vel_CartMesh.YRes;
Res_host[2] = Vel_CartMesh.ZRes;

// Memcpy of constant memory
err = cudaMemcpyToSymbol( RES, &Res_host, 3 * sizeof(int), 0, cudaMemcpyHostToDevice);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Host to Device constant RES failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
fprintf(stderr, " \n\n");
exit(1);
}

// Memcopy Data mesh data to constant memory on device..
XX_host[0] = Data_MeshBounds.XMin ;
XX_host[1] = Data_MeshBounds.XMax;
XX_host[2] = Data_MeshBounds.XDelta;

// Memcpy of constant memory
err = cudaMemcpyToSymbol( XX_Data, &XX_host, 3 * sizeof(double), 0, cudaMemcpyHostToDevice);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Host to Device constant XX_Data failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
fprintf(stderr, " \n\n");
exit(1);
}

XX_host[0] = Data_MeshBounds.YMin ;
XX_host[1] = Data_MeshBounds.YMax;
XX_host[2] = Data_MeshBounds.YDelta;

// Memcpy of constant memory
err = cudaMemcpyToSymbol( YY_Data, &XX_host, 3 * sizeof(double), 0, cudaMemcpyHostToDevice);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Host to Device constant YY_Data failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
fprintf(stderr, " \n\n");
exit(1);
}
XX_host[0] = Data_MeshBounds.ZMin ;
XX_host[1] = Data_MeshBounds.ZMax;
XX_host[2] = Data_MeshBounds.ZDelta;


// Memcpy of constant memory
err = cudaMemcpyToSymbol( ZZ_Data, &XX_host, 3 * sizeof(double), 0, cudaMemcpyHostToDevice);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Host to Device constant ZZ_Data failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
fprintf(stderr, " \n\n");
exit(1);
}



Res_host[0] = Data_MeshBounds.XRes;
Res_host[1] = Data_MeshBounds.YRes;
Res_host[2] = Data_MeshBounds.ZRes;

// Memcpy of constant memory
err = cudaMemcpyToSymbol( RES_Data, &Res_host, 3 * sizeof(int), 0, cudaMemcpyHostToDevice);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Host to Device constant RES_Data failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
fprintf(stderr, " \n\n");
exit(1);
}




}
}

void host2device_gs(int ss) {


if (Data_MeshType == UNSTRUCTURED) {
err = cudaMemcpy(outputid_dev, outputid , ss * sizeof(int) , cudaMemcpyHostToDevice);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Host structure array Vel_MeshNodeArray to Device structure array Vel_MeshNodeArray_device failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
exit(1);
}	

err = cudaMemcpy( MeshNodeArray_device.x , MeshNodeArray.x, Vel_MeshNumNodes * sizeof(float) , cudaMemcpyHostToDevice);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Host structure array MeshNodeArray.x to Device structure array MeshNodeArray_device failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
exit(1);
}	

err = cudaMemcpy( MeshNodeArray_device.y , MeshNodeArray.y, Vel_MeshNumNodes * sizeof(float) , cudaMemcpyHostToDevice);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Host structure array MeshNodeArray.y to Device structure array MeshNodeArray_device failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
exit(1);
}	

err = cudaMemcpy( MeshNodeArray_device.z , MeshNodeArray.z, Vel_MeshNumNodes * sizeof(float) , cudaMemcpyHostToDevice);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Host structure array MeshNodeArray.z to Device structure array MeshNodeArray_device failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
exit(1);
}



//printf("Cuda Memcopy Successful .... \n\n");
}
}

void device2host_gs( int ss) {

if (Data_MeshType == UNSTRUCTURED) {
err = cudaMemcpy(outputid, outputid_dev , (ss) * sizeof(int) , cudaMemcpyDeviceToHost);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Device structure array outputid_dev to Host structure array failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
exit(1);
}
}

}

void host2const_mem_gs(void) {

const size_t ssize = sizeof(float) * size_t(CONSTANT_MEMORY);


// Memcpy of constant memory
err = cudaMemcpyToSymbol( x_dev, &x_host, ssize, 0, cudaMemcpyHostToDevice);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Host to Device constant x_dev failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
fprintf(stderr, " \n\n");
exit(1);
}

err = cudaMemcpyToSymbol( y_dev, &y_host, ssize,0, cudaMemcpyHostToDevice);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Host to Device constant y_dev failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
exit(1);
}


if (Dimensions == 3) {

err = cudaMemcpyToSymbol( z_dev, &z_host, ssize, 0 ,cudaMemcpyHostToDevice );
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Host to Device constant z_dev failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
exit(1);
}
}

//printf("Constant memory loaded  .... \n\n ");


}


void host2device(int ss) {


//	if (Trace_ReleaseStrategy == 0) {
err = cudaMemcpy(Tracer_dev.x, Tracer.x , ss * sizeof(double) , cudaMemcpyHostToDevice);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Host structure array Tracer to Device structure array Tracer_dev (x values) failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
exit(1);
}	

err = cudaMemcpy(Tracer_dev.y, Tracer.y , ss * sizeof(double) , cudaMemcpyHostToDevice);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Host structure array Tracer to Device structure array Tracer_dev (y values) failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
exit(1);
}	

if (Dimensions == 3) {
err = cudaMemcpy(Tracer_dev.z, Tracer.z , ss * sizeof(double) , cudaMemcpyHostToDevice);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Host structure array Tracer to Device structure array Tracer_dev (z values) failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
exit(1);
}	
}

if (Data_MeshType == UNSTRUCTURED) {

err = cudaMemcpy(Tracer_dev.ElementIndex, Tracer.ElementIndex , ss * sizeof(int) , cudaMemcpyHostToDevice);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Host structure array Tracer1 to Device structure array Tracer1_dev (Element ID's) failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
exit(1);
}	

}

err = cudaMemcpy(Tracer_dev.LeftDomain, Tracer.LeftDomain , ss * sizeof(int) , cudaMemcpyHostToDevice);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Host structure array Tracer1 to Device structure array Tracer1_dev (LeftDomain) failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
exit(1);
}	

err = cudaMemcpy(Tracer_dev.Start_time, Tracer.Start_time , ss * sizeof(double) , cudaMemcpyHostToDevice);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Host structure array Tracer1 to Device structure array Tracer1_dev (Start Times) failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
exit(1);
}	

err = cudaMemcpy(Tracer_dev.Stop_time, Tracer.Stop_time , ss * sizeof(double) , cudaMemcpyHostToDevice);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Host structure array Tracer1 to Device structure array Tracer1_dev (Stop times) failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
exit(1);
}	

/*	} else if (Trace_ReleaseStrategy == 1) { // Staggered release

err = cudaMemcpy(Tracer_dev.x, Tracer1.x , ss * sizeof(double) , cudaMemcpyHostToDevice);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Host structure array Tracer to Device structure array Tracer_dev (x values) failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
exit(1);
}	

err = cudaMemcpy(Tracer_dev.y, Tracer1.y , ss * sizeof(double) , cudaMemcpyHostToDevice);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Host structure array Tracer to Device structure array Tracer_dev (y values) failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
exit(1);
}	

if (Dimensions == 3) {
err = cudaMemcpy(Tracer_dev.z, Tracer1.z , ss * sizeof(double) , cudaMemcpyHostToDevice);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Host structure array Tracer to Device structure array Tracer_dev (z values) failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
exit(1);
}	
}

if (Data_MeshType == UNSTRUCTURED) {

err = cudaMemcpy(Tracer_dev.ElementIndex, Tracer1.ElementIndex , ss * sizeof(int) , cudaMemcpyHostToDevice);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Host structure array Tracer1 to Device structure array Tracer1_dev (Element ID's) failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
exit(1);
}	

}

err = cudaMemcpy(Tracer_dev.LeftDomain, Tracer1.LeftDomain , ss * sizeof(int) , cudaMemcpyHostToDevice);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Host structure array Tracer1 to Device structure array Tracer1_dev (LeftDomain) failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
exit(1);
}	

err = cudaMemcpy(Tracer_dev.Start_time, Tracer1.Start_time , ss * sizeof(double) , cudaMemcpyHostToDevice);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Host structure array Tracer1 to Device structure array Tracer1_dev (Start Times) failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
exit(1);
}	

err = cudaMemcpy(Tracer_dev.Stop_time, Tracer1.Stop_time , ss * sizeof(double) , cudaMemcpyHostToDevice);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Host structure array Tracer1 to Device structure array Tracer1_dev (Stop times) failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
exit(1);
}	
}*/

int num;

if (Data_MeshType == CARTESIAN)
num = N_Vel;
else
num = Vel_MeshNumNodes;


err = cudaMemcpy(velocity_dev.u0, velocity.u0 , num * sizeof(double) , cudaMemcpyHostToDevice);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Host structure array velocity to Device structure array velocity_dev (u0) failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
exit(1);
}	

err = cudaMemcpy(velocity_dev.u1, velocity.u1 , num * sizeof(double) , cudaMemcpyHostToDevice);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Host structure array velocity to Device structure array velocity_dev (u1) failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
exit(1);
}	
err = cudaMemcpy(velocity_dev.v0, velocity.v0 , num * sizeof(double) , cudaMemcpyHostToDevice);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Host structure array velocity to Device structure array velocity_dev (v0) failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
exit(1);
}	
err = cudaMemcpy(velocity_dev.v1, velocity.v1 , num * sizeof(double) , cudaMemcpyHostToDevice);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Host structure array velocity to Device structure array velocity_dev (v1) failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
exit(1);
}	

err = cudaMemcpy(velocity_dev.w0, velocity.w0 , num * sizeof(double) , cudaMemcpyHostToDevice);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Host structure array velocity to Device structure array velocity_dev (w0) failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
exit(1);
}	
err = cudaMemcpy(velocity_dev.w1, velocity.w1 , num * sizeof(double) , cudaMemcpyHostToDevice);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Host structure array velocity to Device structure array velocity_dev (w1) failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
exit(1);
}	

err = cudaMemcpy(velocity_dev.time0, velocity.time0 , num * sizeof(double) , cudaMemcpyHostToDevice);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Host structure array velocity to Device structure array velocity_dev (time0) failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
exit(1);
}	

err = cudaMemcpy(velocity_dev.time1, velocity.time1 , num * sizeof(double) , cudaMemcpyHostToDevice);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Host structure array velocity to Device structure array velocity_dev (time1) failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
exit(1);
}	

}

void device2host(int ss) {

//	if (Trace_ReleaseStrategy == 0) {
err = cudaMemcpy(Tracer.x, Tracer_dev.x , (ss) * sizeof(double) , cudaMemcpyDeviceToHost);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Device structure array Tracer1_dev to Host structure array (x value) failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
exit(1);
}

err = cudaMemcpy(Tracer.y, Tracer_dev.y , (ss) * sizeof(double) , cudaMemcpyDeviceToHost);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Device structure array Tracer1_dev to Host structure array (y value) failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
exit(1);
}

if (Dimensions == 3) {
err = cudaMemcpy(Tracer.z, Tracer_dev.z , (ss) * sizeof(double) , cudaMemcpyDeviceToHost);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Device structure array Tracer1_dev to Host structure array (z value) failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
exit(1);
}
}

if (Data_MeshType == UNSTRUCTURED) {


err = cudaMemcpy(Tracer.ElementIndex, Tracer_dev.ElementIndex , (ss) * sizeof(int) , cudaMemcpyDeviceToHost);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Device structure array Tracer1_dev to Host structure array (Element ID's) failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
exit(1);
}


}

err = cudaMemcpy(Tracer.LeftDomain, Tracer_dev.LeftDomain , (ss) * sizeof(int) , cudaMemcpyDeviceToHost);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Device structure array Tracer1_dev to Host structure array (LeftDomain) failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
exit(1);
}

/*	} else if (Trace_ReleaseStrategy == 1) {

err = cudaMemcpy(Tracer1.x, Tracer_dev.x , (ss) * sizeof(double) , cudaMemcpyDeviceToHost);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Device structure array Tracer1_dev to Host structure array (x value) failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
exit(1);
}

err = cudaMemcpy(Tracer1.y, Tracer_dev.y , (ss) * sizeof(double) , cudaMemcpyDeviceToHost);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Device structure array Tracer1_dev to Host structure array (y value) failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
exit(1);
}

if (Dimensions == 3) {
err = cudaMemcpy(Tracer1.z, Tracer_dev.z , (ss) * sizeof(double) , cudaMemcpyDeviceToHost);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Device structure array Tracer1_dev to Host structure array (z value) failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
exit(1);
}
}

if (Data_MeshType == UNSTRUCTURED) {


err = cudaMemcpy(Tracer1.ElementIndex, Tracer_dev.ElementIndex , (ss) * sizeof(int) , cudaMemcpyDeviceToHost);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Device structure array Tracer1_dev to Host structure array (Element ID's) failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
exit(1);
}


}

err = cudaMemcpy(Tracer1.LeftDomain, Tracer_dev.LeftDomain , (ss) * sizeof(int) , cudaMemcpyDeviceToHost);
if(err != cudaSuccess) {
fprintf(stderr, "Memory copy from Device structure array Tracer1_dev to Host structure array (LeftDomain) failed\n");
printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
exit(1);
}





}*/

}



