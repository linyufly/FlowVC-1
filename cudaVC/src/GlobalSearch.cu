

# include "GlobalSearch.h"
# include "velocity.h"
# include "clean.h"

// This file defines all functions related to searching algorithm on Host and on device. Device Global searching is not used as it takes a lot of time. 

__global__ void GlobalSearch2D(Element MeshElementArray_device, Node MeshNodeArray_device, int Num_Frame, int *outputid_dev, int Vel_MeshNumElements, int ss, int N_Frame) {


	int tid = (blockIdx.x*blockDim.x)+threadIdx.x; 
	
	if (tid >= Vel_MeshNumElements)
		return;
	
	float V;
	float r,s,d;

	float x0,y0,x1,y1,x2,y2;

	
	

	// Determinant of mapping from natural to physical coordinates of test element 

	

	x0 = MeshNodeArray_device.x[MeshElementArray_device.Node1[tid]];
	y0 = MeshNodeArray_device.y[MeshElementArray_device.Node1[tid]];

	x1 = MeshNodeArray_device.x[MeshElementArray_device.Node2[tid]];
	y1 = MeshNodeArray_device.y[MeshElementArray_device.Node2[tid]];

	x2 = MeshNodeArray_device.x[MeshElementArray_device.Node3[tid]];
	y2 = MeshNodeArray_device.y[MeshElementArray_device.Node3[tid]];	


	

// Entries for mapping of physical to natural coordinates of test element

//		a11 = y2 - y0
//		a12 = x0 - x2
//		a21 = y0 - y1
//		a22 = x1 - x0

// Determinant of mapping from natural to physical coordinates of test element 
	V = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);

	


	// Main Loop over the number of points passed .......
	for (int N1 = 0; N1 < Num_Frame; N1++ ) {

		if(((Num_Frame*ss) + N1) > N_Frame)
			return;
	
		if (outputid_dev[N1] == -1) { // Search only if the point is not found ...
			// Natural coordinates of point to be interpolated ....
	
			r = ( (y2 - y0) * (x_dev[N1] - x0) + (x0 - x2) * ( y_dev[N1]- y0)) / V;
			s = ( (y0 - y1) * (x_dev[N1] - x0) + (x1 - x0) * ( y_dev[N1]- y0)) / V;	

			d = fmin(r, fmin(s, 1 - r - s));

			if(d > -TINY_GS) { // Point inside test element

				outputid_dev[N1] = tid;
		
			}
		
		}
	}
	
}


// For  improved performance in case of 3D kernel, Try to reduce register usage to 20 .....

// Using pinned memory for Element1 and Node1 ...




__global__ void  GlobalSearch3D(Element MeshElementArray_device, Node MeshNodeArray_device, int Num_Frame, int *outputid_dev, int Vel_MeshNumElements, int ss, int N_Frame) {


	int tid;
	tid = (blockIdx.x*blockDim.x)+threadIdx.x; 
	
	if (tid >= Vel_MeshNumElements)
		return;
	
	float V;
	float r,s,t,d;

	float x0,y0,x1,y1,x2,y2,z0,z1,z2,z3,x3,y3;

	
	

	// Determinant of mapping from natural to physical coordinates of test element 

	

	x0 = MeshNodeArray_device.x[MeshElementArray_device.Node1[tid]];
	y0 = MeshNodeArray_device.y[MeshElementArray_device.Node1[tid]];
	z0 = MeshNodeArray_device.z[MeshElementArray_device.Node1[tid]];

	x1 = MeshNodeArray_device.x[MeshElementArray_device.Node2[tid]];
	y1 = MeshNodeArray_device.y[MeshElementArray_device.Node2[tid]];
	z1 = MeshNodeArray_device.z[MeshElementArray_device.Node2[tid]];

	x2 = MeshNodeArray_device.x[MeshElementArray_device.Node3[tid]];
	y2 = MeshNodeArray_device.y[MeshElementArray_device.Node3[tid]];
	z2 = MeshNodeArray_device.z[MeshElementArray_device.Node3[tid]];

	x3 = MeshNodeArray_device.x[MeshElementArray_device.Node4[tid]];
	y3 = MeshNodeArray_device.y[MeshElementArray_device.Node4[tid]];
	z3 = MeshNodeArray_device.z[MeshElementArray_device.Node4[tid]];



	


/*	// Entries for mapping of physical to natural coordinates of test element
			a11 = (z3 - z0) * (y2 - y3) - (z2 - z3) * (y3 - y0);
			a21 = (z3 - z0) * (y0 - y1) - (z0 - z1) * (y3 - y0);
			a31 = (z1 - z2) * (y0 - y1) - (z0 - z1) * (y1 - y2);
			a12 = (x3 - x0) * (z2 - z3) - (x2 - x3) * (z3 - z0);
			a22 = (x3 - x0) * (z0 - z1) - (x0 - x1) * (z3 - z0);
			a32 = (x1 - x2) * (z0 - z1) - (x0 - x1) * (z1 - z2);
			a13 = (y3 - y0) * (x2 - x3) - (y2 - y3) * (x3 - x0);
			a23 = (y3 - y0) * (x0 - x1) - (y0 - y1) * (x3 - x0);
			a33 = (y1 - y2) * (x0 - x1) - (y0 - y1) * (x1 - x2);	
*/

// Determinant of mapping from natural to physical coordinates of test element 
	V = (x1 - x0) * ((y2 - y0) * (z3 - z0) - (z2 - z0) * (y3 - y0)) + (x2 - x0) * ((y0 - y1) * (z3 - z0) - (z0 - z1) * (y3 - y0)) + (x3 - x0) * ((y1 - y0) * (z2 - z0) - (z1 - z0) * (y2 - y0));


	


	// Main Loop over the number of points passed .......

	for (int N1 = 0; N1 < Num_Frame; N1++ ) {

		if(((Num_Frame*ss) + N1) > N_Frame) {
			return;
		} else {
		if (outputid_dev[N1] == -1) {
	
	
			// Natural coordinates of point to be interpolated ...

			r = (((z3 - z0) * (y2 - y3) - (z2 - z3) * (y3 - y0)) * (x_dev[N1] - x0) + ((x3 - x0) * (z2 - z3) - (x2 - x3) * (z3 - z0)) * (y_dev[N1] - y0) + ((y3 - y0) * (x2 - x3) - (y2 - y3) * (x3 - x0)) * (z_dev[N1] - z0)) / V;
			
			s = (((z3 - z0) * (y0 - y1) - (z0 - z1) * (y3 - y0)) * (x_dev[N1] - x0) +  ((x3 - x0) * (z2 - z3) - (x2 - x3) * (z3 - z0)) * (y_dev[N1] - y0) +  ((y3 - y0) * (x0 - x1) - (y0 - y1) * (x3 - x0)) * (z_dev[N1] - z0)) / V;

			t = ( ((z1 - z2) * (y0 - y1) - (z0 - z1) * (y1 - y2)) * (x_dev[N1] - x0) +  ((x1 - x2) * (z0 - z1) - (x0 - x1) * (z1 - z2)) * (y_dev[N1] - y0) +  ((y1 - y2) * (x0 - x1) - (y0 - y1) * (x1 - x2)) * (z_dev[N1] - z0)) / V;

			d = fmin(r, fmin(s, fmin(t, 1 - r - s - t)));

			if(d > -TINY_GS) { // Point inside test element

				outputid_dev[N1] = tid;
		
			}
		
		}
		
		}
		
	}
	
}
		

void globalsearch_device(void) {

		cudaEvent_t start_gs, stop_gs;
		cudaEventCreate( &start_gs );
		cudaEventCreate( &stop_gs );
		cudaEventRecord( start_gs, 0 );
		
		int i;

		printf("Global searching on Device......... \n\n");
		
		// Since we are computing in chunks of data loaded on memory, we will have a loop for all this chunks.
		int Number_loop;
		Number_loop = (int) (N_Frame/CONSTANT_MEMORY) + 1;
		
		//printf("Number_Points_loaded = %d and number of loops are %d \n",Number_Points_loaded, Number_loop);
	
		// Memcpy some constant parameters ....
		host2device_const();

		// Number of Threads per block.
		nthreads = POINTS_BLOCKSIZE_GS;

		// Calculation of Block size for Global search kernel ...
		if ((Vel_MeshNumElements) % nthreads == 0 ) {
			nblocks = (int)((Vel_MeshNumElements)/nthreads);
		}
		else {
			nblocks = (int)(((Vel_MeshNumElements)/nthreads) + 1);
		}

		//printf("nblocks = %d and nthreads = %d \n\n",nblocks, nthreads);

		int k;

		// Main Loop for Global Search
		//#pragma unroll 5
		for (i = 0; i < Number_loop; i++ ) {

			// Initializing the temp buffer outputid1, x_host, y_host;
			//#pragma unroll 200
			for(k = 0; k < Number_Points_loaded; k++) {

				// For number of tracers in the current frame only ...
				if (((i*Number_Points_loaded) + k) < N_Frame) {
					
					outputid[k] = Tracer.ElementIndex[(i*Number_Points_loaded) + k];
					// Initializing host constant array 
					x_host[k] = (float)Tracer.x[(i*Number_Points_loaded) + k];
					y_host[k] = (float)Tracer.y[(i*Number_Points_loaded) + k];

					if(Dimensions == 3)
						z_host[k] = (float)Tracer.z[(i*Number_Points_loaded) + k];
				
				}

			}
			// Copying outputid1 to device.....
			host2device_gs(Number_Points_loaded);


			host2const_mem_gs();
		
			cudaPrintfInit(); 

			if (Dimensions == 3 ) {
				// Global search function for 3D
				GlobalSearch3D<<< nblocks, nthreads >>>(MeshElementArray_device, MeshNodeArray_device,Number_Points_loaded , outputid_dev, Vel_MeshNumElements , i, N_Frame);
			} else {
				GlobalSearch2D<<< nblocks, nthreads >>>(MeshElementArray_device, MeshNodeArray_device,Number_Points_loaded , outputid_dev, Vel_MeshNumElements , i, N_Frame);
			
			}
			err = cudaThreadSynchronize();
			if(err != cudaSuccess) {
				fprintf(stderr, "Something went terribly wrong in Global_Search.  %i\n", nblocks1);
				printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
				exit(1);
			}
			cudaPrintfDisplay(stdout, true);
			

			// Copying variables back to host from device ....
			device2host_gs(Number_Points_loaded );

			// copy element index calculated from kernel into our temp buffer which will be used for initialization of new tracer releases.
			for (k = 0; k < Number_Points_loaded; k++) {
				
				if (((i*Number_Points_loaded) + k) < N_Frame) {
					
					// Set element index of Tracer
					Tracer.ElementIndex[(i*Number_Points_loaded) + k] = outputid[k];
				
					if(outputid[k] > -1) {
						// Set leftdomain number = 0 if they are in domain. By default this number is 1.
						Tracer.LeftDomain[k] = 0;
					}
				}
			}

		//	printf("Loop number %d for global search completed .....\n\n", i);


		}
		printf("Global Search Successful ...... \n");
		int sum = 0;
		int sum_nf = 0;
	
	
		// Calculate how many elements are inside the domain and how many left
		for (i = 0; i < N_Frame; i++) {

			//printf("Tracer.ElementIndex[i] = %d \t",Tracer.ElementIndex[i]);
			if(Tracer.ElementIndex[i] > -1)
				sum++;
			if(Tracer.ElementIndex[i] == -1)
				sum_nf++;

		}
		N_Present = sum;
		printf("%d points out of %d found inside the domain \n %d points out of %d left domain \n\n\n",sum,N_Frame,sum_nf,N_Frame);


		// Calculate Performance
		cudaEventRecord( stop_gs, 0 );
		cudaEventSynchronize( stop_gs );
		float elapsedTime;
		cudaEventElapsedTime( &elapsedTime, start_gs, stop_gs );
		printf( "Time to complete global search: %3.2f ms\n\n", elapsedTime );
		cudaEventDestroy( start_gs ) ;
		cudaEventDestroy( stop_gs ) ;
		
		// Clean constant memory ...
		clean_gs();



}


int Get_Element_Global_Search(const double *X) {
	 
	// Returns index of Vel_MeshElementArray for element that contains point X
	// Returns -1 if point does not belong to any element 
	
	
	int i;
	double x, y, z, x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3; 
	double a11, a12, a13, a21, a22, a23, a31, a32, a33;
	double r, s, t, d;
	double V;
	
	x = X[0];
	y = X[1];
	z = X[2];
	
	i = 0;
	while(i < Vel_MeshNumElements) { // Search over all elements
		if(Dimensions == 3) {
			// Physical coordinates of nodes of the test element 
			x0 = MeshNodeArray.x[MeshElementArray.Node1[i]];
			y0 = MeshNodeArray.y[MeshElementArray.Node1[i]];
			z0 = MeshNodeArray.z[MeshElementArray.Node1[i]];
			x1 = MeshNodeArray.x[MeshElementArray.Node2[i]];
			y1 = MeshNodeArray.y[MeshElementArray.Node2[i]];
			z1 = MeshNodeArray.z[MeshElementArray.Node2[i]];
			x2 = MeshNodeArray.x[MeshElementArray.Node3[i]];
			y2 = MeshNodeArray.y[MeshElementArray.Node3[i]];
			z2 = MeshNodeArray.z[MeshElementArray.Node3[i]];
			x3 = MeshNodeArray.x[MeshElementArray.Node4[i]];
			y3 = MeshNodeArray.y[MeshElementArray.Node4[i]];
			z3 = MeshNodeArray.z[MeshElementArray.Node4[i]];
			
			// Entries for mapping of physical to natural coordinates of test element 
			a11 = (z3 - z0) * (y2 - y3) - (z2 - z3) * (y3 - y0);
			a21 = (z3 - z0) * (y0 - y1) - (z0 - z1) * (y3 - y0);
			a31 = (z1 - z2) * (y0 - y1) - (z0 - z1) * (y1 - y2);
			a12 = (x3 - x0) * (z2 - z3) - (x2 - x3) * (z3 - z0);
			a22 = (x3 - x0) * (z0 - z1) - (x0 - x1) * (z3 - z0);
			a32 = (x1 - x2) * (z0 - z1) - (x0 - x1) * (z1 - z2);
			a13 = (y3 - y0) * (x2 - x3) - (y2 - y3) * (x3 - x0);
			a23 = (y3 - y0) * (x0 - x1) - (y0 - y1) * (x3 - x0);
			a33 = (y1 - y2) * (x0 - x1) - (y0 - y1) * (x1 - x2);			 
			
			// Determinant of mapping from natural to physical coordinates of test element
			V = (x1 - x0) * ((y2 - y0) * (z3 - z0) - (z2 - z0) * (y3 - y0)) + 
			(x2 - x0) * ((y0 - y1) * (z3 - z0) - (z0 - z1) * (y3 - y0)) +
			(x3 - x0) * ((y1 - y0) * (z2 - z0) - (z1 - z0) * (y2 - y0));
			
			// Natural coordinates of point to be interpolated 
			r = (a11 * (x - x0) + a12 * (y - y0) + a13 * (z - z0)) / V;
			s = (a21 * (x - x0) + a22 * (y - y0) + a23 * (z - z0)) / V;
			t = (a31 * (x - x0) + a32 * (y - y0) + a33 * (z - z0)) / V;
			
			d = fmin(r, fmin(s, fmin(t, 1 - r - s - t)));
			
			if(d > -TINY) // Point inside test element 
				return(i);
			else // Check next element 
				i++;
		} 
		else { // Dimensions == 2 
			// Physical coordinates of nodes of the test element
			x0 = MeshNodeArray_double.x[MeshElementArray.Node1[i]];
			y0 = MeshNodeArray_double.y[MeshElementArray.Node1[i]];
			x1 = MeshNodeArray_double.x[MeshElementArray.Node2[i]];
			y1 = MeshNodeArray_double.y[MeshElementArray.Node2[i]];
			x2 = MeshNodeArray_double.x[MeshElementArray.Node3[i]];
			y2 = MeshNodeArray_double.y[MeshElementArray.Node3[i]];
			
			// Entries for mapping of physical to natural coordinates of test element
			a11 = y2 - y0;
			a12 = x0 - x2;
			a21 = y0 - y1;
			a22 = x1 - x0;
			
			// Determinant of mapping from natural to physical coordinates of test element 
			V = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
			
			// Natural coordinates of point to be interpolated
			r = (a11 * (x - x0) + a12 * (y - y0)) / V;
			s = (a21 * (x - x0) + a22 * (y - y0)) / V;
			
			d = fmin(r, fmin(s, 1 - r - s));
			
			if(d > -TINY) // Point inside test element
				return(i);
			else // Check next element
				i++;
		}
	}
	
	return -1;
	
}



int Get_Element_Local_Search(const double *X, int guess) {
	
//	 Returns index to Vel_MeshElementArray for element that contains point X
//	 Returns -1 if element not located
	 
	double x, y, z, x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3;
	double a11, a12, a13, a21, a22, a23, a31, a32, a33;
	double r, s, t, d;
	double V;
	
	x = X[0];
	y = X[1];
	z = X[2];
	
	
	
	while(1) {
		if(Dimensions == 3) {
			// Physical coordinates of nodes of the test element 
			x0 = MeshNodeArray_double.x[MeshElementArray.Node1[guess]];
			y0 = MeshNodeArray_double.y[MeshElementArray.Node1[guess]];
			z0 = MeshNodeArray_double.z[MeshElementArray.Node1[guess]];
			x1 = MeshNodeArray_double.x[MeshElementArray.Node2[guess]];
			y1 = MeshNodeArray_double.y[MeshElementArray.Node2[guess]];
			z1 = MeshNodeArray_double.z[MeshElementArray.Node2[guess]];
			x2 = MeshNodeArray_double.x[MeshElementArray.Node3[guess]];
			y2 = MeshNodeArray_double.y[MeshElementArray.Node3[guess]];
			z2 = MeshNodeArray_double.z[MeshElementArray.Node3[guess]];
			x3 = MeshNodeArray_double.x[MeshElementArray.Node4[guess]];
			y3 = MeshNodeArray_double.y[MeshElementArray.Node4[guess]];
			z3 = MeshNodeArray_double.z[MeshElementArray.Node4[guess]];
			
			
			//printf("guess = %d \n", guess);
			
			// Entries for mapping of physical to natural coordinates of test element 
			a11 = (z3 - z0) * (y2 - y3) - (z2 - z3) * (y3 - y0);
			a21 = (z3 - z0) * (y0 - y1) - (z0 - z1) * (y3 - y0);
			a31 = (z1 - z2) * (y0 - y1) - (z0 - z1) * (y1 - y2);
			a12 = (x3 - x0) * (z2 - z3) - (x2 - x3) * (z3 - z0);
			a22 = (x3 - x0) * (z0 - z1) - (x0 - x1) * (z3 - z0);
			a32 = (x1 - x2) * (z0 - z1) - (x0 - x1) * (z1 - z2);
			a13 = (y3 - y0) * (x2 - x3) - (y2 - y3) * (x3 - x0);
			a23 = (y3 - y0) * (x0 - x1) - (y0 - y1) * (x3 - x0);
			a33 = (y1 - y2) * (x0 - x1) - (y0 - y1) * (x1 - x2);			 
			
			// Determinant of mapping from natural to physical coordinates of test element 
			V = (x1 - x0) * ((y2 - y0) * (z3 - z0) - (z2 - z0) * (y3 - y0)) + 
			(x2 - x0) * ((y0 - y1) * (z3 - z0) - (z0 - z1) * (y3 - y0)) +
			(x3 - x0) * ((y1 - y0) * (z2 - z0) - (z1 - z0) * (y2 - y0));
			
			// Natural coordinates of point to be interpolated 
			r = (a11 * (x - x0) + a12 * (y - y0) + a13 * (z - z0)) / V;
			s = (a21 * (x - x0) + a22 * (y - y0) + a23 * (z - z0)) / V;
			t = (a31 * (x - x0) + a32 * (y - y0) + a33 * (z - z0)) / V;
			
			d = fmin(r, fmin(s, fmin(t, 1 - r - s - t)));
			//
			if(d > -TINY) { // Point inside test element 
				
				//printf("ls\t");
				return(guess);
			}
			else { // Reset test element to neighbor 
				if(fabs(r - d) < TINY) 
					guess = MeshElementArray.Neighborindex1[guess]; 
				else if(fabs(s - d) < TINY)  
					guess = MeshElementArray.Neighborindex2[guess]; 
				else if(fabs(t - d) < TINY) 
					guess = MeshElementArray.Neighborindex3[guess];  
				else if(fabs(d - (1 - r - s - t)) < TINY) 
					guess = MeshElementArray.Neighborindex4[guess];
				else 
					FatalError("Indeterminate neighbor in function Get_Element_Local_Search");
					
				if(guess == -1) // Neighbor not present, could not find element from local search (likely point left domain) 
					return(-1);
			}
		}
		else { // Dimensions == 2
			// Physical coordinates of nodes of the test element 
			x0 = MeshNodeArray_double.x[MeshElementArray.Node1[guess]];
			y0 = MeshNodeArray_double.y[MeshElementArray.Node1[guess]];
			x1 = MeshNodeArray_double.x[MeshElementArray.Node2[guess]];
			y1 = MeshNodeArray_double.y[MeshElementArray.Node2[guess]];
			x2 = MeshNodeArray_double.x[MeshElementArray.Node3[guess]];
			y2 = MeshNodeArray_double.y[MeshElementArray.Node3[guess]];
			
			// printf("x = %f; y =  %f; in element %d x0 = %f; y0 = %f;  x1 = %f;  y1 = %f;  x2 = %f y2 = %f\n", x, y, guess, x0, y0, x1, y1, x2, y2);
			
			// Entries for mapping of physical to natural coordinates of test element 
			a11 = y2 - y0;
			a12 = x0 - x2;
			a21 = y0 - y1;
			a22 = x1 - x0;
			
			// Determinant of mapping from natural to physical coordinates of test element 
			V = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
			
			// Natural coordinates of point to be interpolated 
			r = (a11 * (x - x0) + a12 * (y - y0)) / V;
			s = (a21 * (x - x0) + a22 * (y - y0)) / V;
			
			d = fmin(r, fmin(s, 1 - r - s));
			
			if(d > -TINY) // Point inside test element 
				return(guess);
			else { // Reset test element to neighbor 
				if(fabs(r - d) < TINY) 
					guess = MeshElementArray.Neighborindex1[guess];
				else if(fabs(s - d) < TINY) 
					guess = MeshElementArray.Neighborindex2[guess];
				else if(fabs(d - (1 - r - s)) < TINY) 
					guess = MeshElementArray.Neighborindex3[guess];
				else {
					printf("Indeterminate neighbor in function Get_Element_Local_Search");
					exit(0);
				}
				if(guess == -1) // Neighbor not present, could not find element from local search (point may have exited domain)
					return(-1);
			}
		}
	}
	
}

void search_host(void) {

	cudaEvent_t start_gs, stop_gs;
	cudaEventCreate( &start_gs );
	cudaEventCreate( &stop_gs );
	cudaEventRecord( start_gs, 0 );
	
	printf("Searching on host ... \n");
	
	
	int ii, seed = -1, index = -1, found = 0, foundGS = 0, guess = -1;

	double Xseed[3], Xpoint[3];
	
	Xseed[0] = Trace_CartMesh.XMin + ((int)(Trace_CartMesh.XRes/2)) * Trace_CartMesh.XDelta;
	Xseed[1] = Trace_CartMesh.YMin + ((int)(Trace_CartMesh.YRes/2)) * Trace_CartMesh.YDelta;
	Xseed[2] = Trace_CartMesh.ZMin + ((int)(Trace_CartMesh.ZRes/2)) * Trace_CartMesh.ZDelta;
	
	//printf("Xseed[0] = %f Xseed[1] = %f and Xseed[2] = %f ..\n\n",Xseed[0], Xseed[1],Xseed[2]);
	seed = Get_Element_Global_Search(Xseed);
	
	if(seed >= 0) {
	
		printf("  Using center element seed.\n");
		fflush(stdout);
		guess  = seed;
	
	}
	
	for (ii = 0; ii < N_Frame; ii++) {
	
		// Set coordinates of point
		Xpoint[0] = Tracer.x[ii];
		Xpoint[1] = Tracer.y[ii];
		Xpoint[2] = Tracer.z[ii];
		
		// See if coordinates outside of data bounding box
		Tracer.LeftDomain[ii] = TestOutsideDomain_host(Xpoint);

	
		// If inside bounding box and using a tetmesh, find element the point is in
		if(!Tracer.LeftDomain[ii]) {
		
			
			if (seed < 0) { // Use global search until a point has been located
			
				seed = Get_Element_Global_Search(Xpoint);
				
				if (seed < 0) { // Point outside mesh boundary 
					Tracer.LeftDomain[ii] = 1;
					Tracer.ElementIndex[ii] = -1;

				}
				else {
					printf("  Using first found element seed.\n");
					Tracer.ElementIndex[ii] = seed;
					
					found++;
					guess  = seed;
				}
			}
			else { // A point has been located; use local searching
			
				// Find location of point using local search
				index = Get_Element_Local_Search(Xpoint, guess);
			
				if (index < 0) { // Local search from guess not working
				
					// Try local search from seed 
					index = Get_Element_Local_Search(Xpoint, seed);
					
					if (index < 0) { // Local search from seed not working
						// Do global search if requested.
						
						if (LocalSearchChecking) {
						
							index = Get_Element_Global_Search(Xpoint);
						
							if(index < 0) { // Point Outside the domain  
							
								Tracer.LeftDomain[ii] = 1;
								Tracer.ElementIndex[ii] = -1;
						
							} else { // Global search succeeded 
						
								Tracer.ElementIndex[ii] = index;
								Tracer.LeftDomain[ii] = 0;
							
								guess = index;
								found++;
								foundGS++;
						
							}
						
						
						} else { // Local search failed and not attempting global search
					
							Tracer.LeftDomain[ii] = 1;
							Tracer.ElementIndex[ii] = -1;
					
						}
				
					} else { // Local search from seed succeeded
				
						Tracer.ElementIndex[ii] = index;
						Tracer.LeftDomain[ii] = 0;
							
						guess = index;
						found++;
					
					}
					
				} else { // Local search from guess succeeded 
				
					Tracer.ElementIndex[ii] = index;
					Tracer.LeftDomain[ii] = 0;
							
					guess = index;
					
					found++;
				}
				
			}
		
		}
	}
	
	if(LocalSearchChecking) 
		printf("  %d of %d located in domain (%d caught by global search checking)\n", found, Trace_CartMesh.XRes * Trace_CartMesh.YRes * Trace_CartMesh.ZRes, foundGS);
	else
		printf("  %d of %d located in domain\n", found, Trace_CartMesh.XRes * Trace_CartMesh.YRes * Trace_CartMesh.ZRes);  
	
	fflush(stdout);
	
	// Calculate Performance
	cudaEventRecord( stop_gs, 0 );
	cudaEventSynchronize( stop_gs );
	float elapsedTime;
	cudaEventElapsedTime( &elapsedTime, start_gs, stop_gs );
	printf( "Time to complete global search: %3.2f ms\n\n", elapsedTime );
	cudaEventDestroy( start_gs ) ;
	cudaEventDestroy( stop_gs ) ;
	
}








