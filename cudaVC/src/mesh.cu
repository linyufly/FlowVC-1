
#include "mesh.h"

// This file will load mesh related data.


void LoadMeshData(void) {

	printf("Loading mesh data... \n");

		if (Data_MeshType == CARTESIAN) 
			LoadCartMeshData(); // Cartesian mesh
		else // Unstructured mesh
			LoadUnstructMeshData();

	printf("Mesh Data Loaded .... \n\n");

}

void LoadCartMeshData(void) {

	char Data_BinFilePath[LONGSTRING];
	FILE *Data_BinFileID;
	
	// Open binary mesh file for reading 
	sprintf(Data_BinFilePath, "%s%s_Cartesian.bin", Path_Data, Data_InFilePrefix);
	if((Data_BinFileID = fopen(Data_BinFilePath, "rb")) == NULL) 
		FatalError("Could not open file %s", Data_BinFilePath);
	// Read grid parameters 
	if(fread(&Vel_CartMesh.XMin, sizeof(double), 1, Data_BinFileID) < 1) 
		FatalError("Could not read Vel_CartMesh.XMin from file %s", Data_BinFilePath);
	if(fread(&Vel_CartMesh.XMax, sizeof(double), 1, Data_BinFileID) < 1) 
		FatalError("Could not read Vel_CartMesh.XMax from file %s", Data_BinFilePath);
	if(fread(&Vel_CartMesh.XRes, sizeof(int), 1, Data_BinFileID) < 1) 
		FatalError("Could not read Vel_CartMesh.XRes from file %s", Data_BinFilePath);
	if(fread(&Vel_CartMesh.YMin, sizeof(double), 1, Data_BinFileID) < 1) 
		FatalError("Could not read Vel_CartMesh.YMin from file %s", Data_BinFilePath);
	if(fread(&Vel_CartMesh.YMax, sizeof(double), 1, Data_BinFileID) < 1) 
		FatalError("Could not read Vel_CartMesh.YMax from file %s", Data_BinFilePath);
	if(fread(&Vel_CartMesh.YRes, sizeof(int), 1, Data_BinFileID) < 1) 
		FatalError("Could not read Vel_CartMesh.YRes from file %s", Data_BinFilePath);
	if(fread(&Vel_CartMesh.ZMin, sizeof(double), 1, Data_BinFileID) < 1) 
		FatalError("Could not read Vel_CartMesh.ZMin from file %s", Data_BinFilePath);
	if(fread(&Vel_CartMesh.ZMax, sizeof(double), 1, Data_BinFileID) < 1) 
		FatalError("Could not read Vel_CartMesh.ZMax from file %s", Data_BinFilePath);
	if(fread(&Vel_CartMesh.ZRes, sizeof(int), 1, Data_BinFileID) < 1) 
		FatalError("Could not read Vel_CartMesh.ZRes from file %s", Data_BinFilePath);
		
	printf("XRes = %d, YRes = %d and ZRES = %d \n\n",Vel_CartMesh.XRes, Vel_CartMesh.YRes,Vel_CartMesh.ZRes);
	// Close file 
	fclose(Data_BinFileID);
	// Derived parameters 
	if(Vel_CartMesh.XRes < 2 || Vel_CartMesh.XMin >= Vel_CartMesh.XMax)
		FatalError("Insufficient resolution of velocity data in x-direction");
	else
		Vel_CartMesh.XDelta = (Vel_CartMesh.XMax - Vel_CartMesh.XMin) / (Vel_CartMesh.XRes - 1);
	if(Vel_CartMesh.YRes < 2 || Vel_CartMesh.YMin >= Vel_CartMesh.YMax)
		FatalError("Insufficient resolution of velocity data in y-direction");
	else
		Vel_CartMesh.YDelta = (Vel_CartMesh.YMax - Vel_CartMesh.YMin) / (Vel_CartMesh.YRes - 1);
	if(Dimensions == 3) {
		if(Vel_CartMesh.ZRes < 2 || Vel_CartMesh.ZMin >= Vel_CartMesh.ZMax)
			FatalError("Insufficient resolution of velocity data in z-direction");
		else
			Vel_CartMesh.ZDelta = (Vel_CartMesh.ZMax - Vel_CartMesh.ZMin) / (Vel_CartMesh.ZRes - 1);
	}
	else { // Dimensions == 2
		Vel_CartMesh.ZMin = 0;
		Vel_CartMesh.ZMax = 0;
		Vel_CartMesh.ZRes = 1;
		Vel_CartMesh.ZDelta = 0;
	}
	
	N_Vel = Vel_CartMesh.XRes * Vel_CartMesh.YRes * Vel_CartMesh.ZRes;
		
	printf("Total number of points in velocity data frame = %d \n\n",N_Vel);

}



void LoadUnstructMeshData(void) {
	/*
	 *	Reads in unstructured mesh information into Vel_MeshElementArray and Vel_MeshNodeArray 
	 */
	
	int i, elements;
	char Mesh_BinFilePath[LONGSTRING], FileName[SHORTSTRING];
	FILE *Mesh_BinFileID;
	
	//****************************** NODES *********************************
	// Open file to read node coordinates 
	sprintf(FileName, "%s_coordinates.bin", Data_InFilePrefix);
	sprintf(Mesh_BinFilePath, "%s%s", Path_Data, FileName);
	if((Mesh_BinFileID = fopen(Mesh_BinFilePath, "rb")) == NULL) 
		FatalError("Could not open file %s", Mesh_BinFilePath);
	
	// Read number of nodes
	if(fread(&Vel_MeshNumNodes, sizeof(int), 1, Mesh_BinFileID) < 1) 
		FatalError("Could not read number of mesh nodes from %s", Mesh_BinFilePath);
	
	printf("Loading coordinate data for %d nodes... \n", Vel_MeshNumNodes);  
	fflush(stdout);
	
/*	// Allocate memory for Vel_MeshNodeArray 
	if((Vel_MeshNodeArray = (double **)malloc(Vel_MeshNumNodes * sizeof(double *))) == NULL)
		FatalError("Malloc failed for Vel_MeshNodeArray");
	for(i = 0; i < Vel_MeshNumNodes; i++) 
		if((Vel_MeshNodeArray[i] = (double *)malloc(3 * sizeof(double))) == NULL)
			FatalError("Malloc failed for Vel_MeshNodeArray[%d]", i);
	
	// Read Vel_MeshNodeArray information from file
	for(i = 0; i < Vel_MeshNumNodes; i++) 
		if(fread(Vel_MeshNodeArray[i], sizeof(double), 3, Mesh_BinFileID) < 3)
			FatalError("Could not read nodal coordinates completely from %s", Mesh_BinFilePath);
*/	

	// Allocate arrays of structure for NODE struct....


	if(Memory_Usage_GS == Use_Pinned_Memory) { // Using Pinned memory

		// Allocate for MeshNodeArray which has float values (Used for global search)
		err = cudaHostAlloc( (void**)&MeshNodeArray.x, Vel_MeshNumNodes * sizeof( float ), cudaHostAllocDefault );
		if(err != cudaSuccess) {
			fprintf(stderr, "Something went terribly wrong in allocating Host Memory for Mesh Node Array(x position).\n");
			printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
			exit(1);
		}

		err = cudaHostAlloc( (void**)&MeshNodeArray.y, Vel_MeshNumNodes * sizeof( float ), cudaHostAllocDefault );
		if(err != cudaSuccess) {
			fprintf(stderr, "Something went terribly wrong in allocating Host Memory for Mesh Node Array(y position).\n");
			printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
			exit(1);
		}

		err = cudaHostAlloc( (void**)&MeshNodeArray.z, Vel_MeshNumNodes * sizeof( float ), cudaHostAllocDefault );
		if(err != cudaSuccess) {
			fprintf(stderr, "Something went terribly wrong in allocating Host Memory for Mesh Node Array(z position).\n");
			printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
			exit(1);
		}
		  
		
		// Allocate for MeshNodeArray_double which has double values (Used for integration)
		err = cudaHostAlloc( (void**)&MeshNodeArray_double.x, Vel_MeshNumNodes * sizeof( double ), cudaHostAllocDefault );
		if(err != cudaSuccess) {
			fprintf(stderr, "Something went terribly wrong in allocating Host Memory for Mesh Node Array_double(x position).\n");
			printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
			exit(1);
		}

		err = cudaHostAlloc( (void**)&MeshNodeArray_double.y, Vel_MeshNumNodes * sizeof( double ), cudaHostAllocDefault );
		if(err != cudaSuccess) {
			fprintf(stderr, "Something went terribly wrong in allocating Host Memory for Mesh Node Array_double(y position).\n");
			printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
			exit(1);
		}

		err = cudaHostAlloc( (void**)&MeshNodeArray_double.z, Vel_MeshNumNodes * sizeof( double ), cudaHostAllocDefault );
		if(err != cudaSuccess) {
			fprintf(stderr, "Something went terribly wrong in allocating Host Memory for Mesh Node Array_double(z position).\n");
			printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
			exit(1);
		}
		  

	} else { // Using DRAM in CPU

		// Allocate for MeshNodeArray which has float values (Used for global search)
		if(( MeshNodeArray.x = (float *)malloc(Vel_MeshNumNodes * sizeof(float))) == NULL)
			FatalError("Malloc failed for MeshNodeArray.x failed");
		if(( MeshNodeArray.y = (float *)malloc(Vel_MeshNumNodes * sizeof(float))) == NULL)
			FatalError("Malloc failed for MeshNodeArray.y failed");
		if(( MeshNodeArray.z = (float *)malloc(Vel_MeshNumNodes * sizeof(float))) == NULL)
			FatalError("Malloc failed for MeshNodeArray.z failed");
		
		// Allocate for MeshNodeArray_double which has double values (Used for integration)
		if(( MeshNodeArray_double.x = (double *)malloc(Vel_MeshNumNodes * sizeof(double))) == NULL)
			FatalError("Malloc failed for MeshNodeArray_double.x failed");
		if(( MeshNodeArray_double.y = (double *)malloc(Vel_MeshNumNodes * sizeof(double))) == NULL)
			FatalError("Malloc failed for MeshNodeArray_double.y failed");
		if(( MeshNodeArray_double.z = (double *)malloc(Vel_MeshNumNodes * sizeof(double))) == NULL)
			FatalError("Malloc failed for MeshNodeArray_double.z failed");

	}



	for (i = 0; i < Vel_MeshNumNodes; i++) {
	
		if(fread(&MeshNodeArray_double.x[i], sizeof(double), 1, Mesh_BinFileID) < 1)
			FatalError("Could not read x nodal coordinates completely from %s", Mesh_BinFilePath);
		
		if(fread(&MeshNodeArray_double.y[i], sizeof(double), 1, Mesh_BinFileID) < 1)
			FatalError("Could not read y nodal coordinates completely from %s", Mesh_BinFilePath);
		
		if(fread(&MeshNodeArray_double.z[i], sizeof(double), 1, Mesh_BinFileID) < 1)
			FatalError("Could not read z nodal coordinates completely from %s", Mesh_BinFilePath);
		
		
		MeshNodeArray.x[i] = (float)MeshNodeArray_double.x[i];
		MeshNodeArray.y[i] = (float)MeshNodeArray_double.y[i];
		MeshNodeArray.z[i] = (float)MeshNodeArray_double.z[i];
	}
	
	
	
		
	
	
	
	fclose(Mesh_BinFileID);
	printf("OK!\n"); 
	fflush(stdout);
	
	// ***************** CONNECTIVITY AND ADJACENCY ***********************
	// Open connectivity file to read binary 
	sprintf(FileName, "%s_connectivity.bin", Data_InFilePrefix);
	sprintf(Mesh_BinFilePath, "%s%s", Path_Data, FileName);
	if((Mesh_BinFileID = fopen(Mesh_BinFilePath, "rb")) == NULL) 
		FatalError("Could not open file %s", Mesh_BinFilePath);
	
	// Read number of elements 
	if(fread(&Vel_MeshNumElements, sizeof(int), 1, Mesh_BinFileID) < 1)
		FatalError("Could not read number of mesh elements from %s", Mesh_BinFilePath);
	
	printf("Loading connectivity and adjacency data for %d elements... \n", Vel_MeshNumElements); 
	fflush(stdout);
	
	// Allocate memory for arrays of MeshElementArray

	if (Memory_Usage_GS == Use_Pinned_Memory) {

		err = cudaHostAlloc( (void**)&MeshElementArray.Node1, Vel_MeshNumElements * sizeof( int ), cudaHostAllocDefault );
		if(err != cudaSuccess) {
			fprintf(stderr, "Something went terribly wrong in allocating Host Memory for Mesh Element Array(Node 1).\n");
			printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
			exit(1);
		  }

		err = cudaHostAlloc( (void**)&MeshElementArray.Node2, Vel_MeshNumElements * sizeof( int ), cudaHostAllocDefault );
		if(err != cudaSuccess) {
			fprintf(stderr, "Something went terribly wrong in allocating Host Memory for Mesh Element Array(Node 2).\n");
			printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
			exit(1);
		  }
		
		err = cudaHostAlloc( (void**)&MeshElementArray.Node3, Vel_MeshNumElements * sizeof( int ), cudaHostAllocDefault );
		if(err != cudaSuccess) {
			fprintf(stderr, "Something went terribly wrong in allocating Host Memory for Mesh Element Array(Node 3).\n");
			printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
			exit(1);
		  }

		err = cudaHostAlloc( (void**)&MeshElementArray.Node4, Vel_MeshNumElements * sizeof( int ), cudaHostAllocDefault );
		if(err != cudaSuccess) {
			fprintf(stderr, "Something went terribly wrong in allocating Host Memory for Mesh Element Array(Node 4).\n");
			printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
			exit(1);
		  }
		



	} else { 
		if((MeshElementArray.Node1 = (int *)malloc(Vel_MeshNumElements * sizeof(int))) == NULL)
			FatalError("Malloc failed for MeshElementArray.Node1");

		if((MeshElementArray.Node2 = (int *)malloc(Vel_MeshNumElements * sizeof(int))) == NULL)
			FatalError("Malloc failed for MeshElementArray.Node2");

		if((MeshElementArray.Node3 = (int *)malloc(Vel_MeshNumElements * sizeof(int))) == NULL)
			FatalError("Malloc failed for MeshElementArray.Node3");

		if((MeshElementArray.Node4 = (int *)malloc(Vel_MeshNumElements * sizeof(int))) == NULL)
			FatalError("Malloc failed for MeshElementArray.Node4");
	}
	


	if((MeshElementArray.Neighborindex1 = (int *)malloc(Vel_MeshNumElements * sizeof(int))) == NULL)
		FatalError("Malloc failed for MeshElementArray.Neighborindex1");

	if((MeshElementArray.Neighborindex2 = (int *)malloc(Vel_MeshNumElements * sizeof(int))) == NULL)
		FatalError("Malloc failed for MeshElementArray.Neighborindex2");

	if((MeshElementArray.Neighborindex3 = (int *)malloc(Vel_MeshNumElements * sizeof(int))) == NULL)
		FatalError("Malloc failed for MeshElementArray.Neighborindex3");

	if((MeshElementArray.Neighborindex4 = (int *)malloc(Vel_MeshNumElements * sizeof(int))) == NULL)
		FatalError("Malloc failed for MeshElementArray.Neighborindex4");

	
	// Read connectivity information from file
	for(i = 0; i < Vel_MeshNumElements; i++) {
		

		if(fread(&MeshElementArray.Node1[i], sizeof(int), 1, Mesh_BinFileID) < 1)
			FatalError("Could not read connectivity(Node1) for element %d from %s", i, Mesh_BinFilePath);

		if(fread(&MeshElementArray.Node2[i], sizeof(int), 1, Mesh_BinFileID) < 1)
			FatalError("Could not read connectivity(Node2) for element %d from %s", i, Mesh_BinFilePath);

		if(fread(&MeshElementArray.Node3[i], sizeof(int), 1, Mesh_BinFileID) < 1)
			FatalError("Could not read connectivity(Node3) for element %d from %s", i, Mesh_BinFilePath);

		if(fread(&MeshElementArray.Node4[i], sizeof(int), 1, Mesh_BinFileID) < 1)
			FatalError("Could not read connectivity(Node4) for element %d from %s", i, Mesh_BinFilePath);	

	}
	
	// Close connectivity file
	fclose(Mesh_BinFileID);
	
	// Open adjacency file to read binary
	sprintf(FileName, "%s_adjacency.bin", Data_InFilePrefix);
	sprintf(Mesh_BinFilePath, "%s%s", Path_Data, FileName);
	if((Mesh_BinFileID = fopen(Mesh_BinFilePath, "rb")) == NULL) 
		FatalError("Could not open file %s", Mesh_BinFilePath);
	
	// Read number of elements 
	if(fread(&elements, sizeof(int), 1, Mesh_BinFileID) < 1)
		FatalError("Could not read number of mesh elements from %s", Mesh_BinFilePath);
	if(elements != Vel_MeshNumElements)
		FatalError("Incompatible number of elements listed in connectivity and adjacency files");
	
	// Read adjacency information from file 
	for(i = 0; i < Vel_MeshNumElements; i++) {


		if(fread(&MeshElementArray.Neighborindex1[i], sizeof(int), 1, Mesh_BinFileID) < 1)
			FatalError("Could not read adjacency(NeighborIndex1) for element %d from %s", i, Mesh_BinFilePath);

		if(fread(&MeshElementArray.Neighborindex2[i], sizeof(int), 1, Mesh_BinFileID) < 1)
			FatalError("Could not read adjacency(NeighborIndex1) for element %d from %s", i, Mesh_BinFilePath);

		if(fread(&MeshElementArray.Neighborindex3[i], sizeof(int), 1, Mesh_BinFileID) < 1)
			FatalError("Could not read adjacency(NeighborIndex1) for element %d from %s", i, Mesh_BinFilePath);

		if(fread(&MeshElementArray.Neighborindex4[i], sizeof(int), 1, Mesh_BinFileID) < 1)
			FatalError("Could not read adjacency(NeighborIndex1) for element %d from %s", i, Mesh_BinFilePath);



}
	
	
	// Close adjacency file
	fclose(Mesh_BinFileID);
	
	printf("OK!\n");  
	fflush(stdout);
	
}





