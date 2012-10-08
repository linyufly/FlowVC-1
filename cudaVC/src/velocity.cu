#include "velocity.h"



int TestOutsideDomain_host(double point[3]) {
	/*  
	 inside data bounding box --> return 0
	 outside data bounding box --> return 1
	 */ 
	
	if(point[0] < (Data_MeshBounds.XMin - TINY) || point[0] > (Data_MeshBounds.XMax + TINY)) 
		return(1);
	else if(point[1] < (Data_MeshBounds.YMin - TINY) || point[1] > (Data_MeshBounds.YMax + TINY)) 
		return(1);
	else if(point[2] < (Data_MeshBounds.ZMin - TINY) || point[2] > (Data_MeshBounds.ZMax + TINY)) 
		return(1);
	else
		return(0);
	
}

// Velocity cartesian mesh data in constant memory.

__device__ void GetVel_cartesian2D(double posx, double posy, double tloc, double *dXdt, VelData_double  vel, int ld ) {

	
	if (ld == 1) {
		dXdt[0] = 0.00;
		dXdt[1] = 0.00;
	} else {
	
		double xloc, yloc;
		int    i, j;

		i = (int)floor((posx - XX[0]) * XX[2]);
		if(i == (RES[0] - 1)) {
			i = RES[0] - 2;
			xloc = 1.0;
		}
		else
			xloc = ((posx - XX[0])* XX[2]) - i;
			
		j = (int)floor((posy - YY[0]) * YY[2]);
		if(j == (RES[1] - 1)) {
			j = RES[1] - 2;
			yloc = 1.0;
		}
		else
			yloc = ((posy - YY[0])* YY[2]) - j;
		
		/////////////////////////////////////////////
		// Linear Interpolation coefficients 
		//V1 = (tloc)*(xloc)*(yloc)*(zloc) = 0
		//V2 = (tloc)*(xloc)*(yloc)*(1-zloc);
		//V3 = (tloc)*(xloc)*(1-yloc)*(zloc) = 0
		//V4 = (tloc)*(xloc)*(1-yloc)*(1-zloc);
		//V5 = (tloc)*(1-xloc)*(yloc)*(zloc) = 0
		//V6 = (tloc)*(1-xloc)*(yloc)*(1-zloc);
		//V7 = (tloc)*(1-xloc)*(1-yloc)*(zloc) = 0
		//V8 = (tloc)*(1-xloc)*(1-yloc)*(1-zloc);
		//V9 = (1-tloc)*(xloc)*(yloc)*(zloc) = 0
		//V10 = (1-tloc)*(xloc)*(yloc)*(1-zloc);
		//V11 = (1-tloc)*(xloc)*(1-yloc)*(zloc) = 0
		//V12 = (1-tloc)*(xloc)*(1-yloc)*(1-zloc);
		//V13 = (1-tloc)*(1-xloc)*(yloc)*(zloc) = 0
		//V14 = (1-tloc)*(1-xloc)*(yloc)*(1-zloc);
		//V15 = (1-tloc)*(1-xloc)*(1-yloc)*(zloc) = 0
		//V16 = (1-tloc)*(1-xloc)*(1-yloc)*(1-zloc);
		///////////////////////////////////////////////
		
		double f000, f001, f010, f011;
		double f100, f101, f110, f111;
		
		// Vertices of space-time voxel for U
		f000 = vel.u0[getindex(i,j,0)];
		f001 = vel.u0[getindex(i,j+1,0)];
		f010 = vel.u0[getindex(i+1,j,0)];
		f011 = vel.u0[getindex(i+1,j+1,0)];
		
		f100 = vel.u1[getindex(i,j,0)];
		f101 = vel.u1[getindex(i,j+1,0)];
		f110 = vel.u1[getindex(i+1,j,0)];
		f111 = vel.u1[getindex(i+1,j+1,0)];
		
		
		// Linear interpolation for U
		
		dXdt[0] =  ( (1-tloc)*(1-xloc)*(1-yloc)*f000 ) + ( (1-tloc)*(1-xloc)*(yloc)*f001 ) + ( (1-tloc)*(xloc)*(1-yloc)*f010 ) + ( (1-tloc)*(xloc)*(yloc)*f011 ) 					+  ( (tloc)*(1-xloc)*(1-yloc)*f100 ) + ( (tloc)*(1-xloc)*(yloc)*f101 ) + ( (tloc)*(xloc)*(1-yloc)*f110 ) + ( (tloc)*(xloc)*(yloc)*f111 );
	
		
		// Vertices of space-time voxel for V
		f000 = vel.v0[getindex(i,j,0)];
		f001 = vel.v0[getindex(i,j+1,0)];
		f010 = vel.v0[getindex(i+1,j,0)];
		f011 = vel.v0[getindex(i+1,j+1,0)];
		
		f100 = vel.v1[getindex(i,j,0)];
		f101 = vel.v1[getindex(i,j+1,0)];
		f110 = vel.v1[getindex(i+1,j,0)];
		f111 = vel.v1[getindex(i+1,j+1,0)];
		
		// Linear interpolation for V
		
		dXdt[1] =  ( (1-tloc)*(1-xloc)*(1-yloc)*f000 ) + ( (1-tloc)*(1-xloc)*(yloc)*f001 ) + ( (1-tloc)*(xloc)*(1-yloc)*f010 ) + ( (1-tloc)*(xloc)*(yloc)*f011 ) 					+  ( (tloc)*(1-xloc)*(1-yloc)*f100 ) + ( (tloc)*(1-xloc)*(yloc)*f101 ) + ( (tloc)*(xloc)*(1-yloc)*f110 ) + ( (tloc)*(xloc)*(yloc)*f111 );
	
	
	}
	

}


__device__ void GetVel_cartesian3D(double posx, double posy, double posz, double tloc, double *dXdt, VelData_double  vel, int ld ) {

	
	if (ld == 1) {
		dXdt[0] = 0.00;
		dXdt[1] = 0.00;
	} else {
	
		double xloc, yloc, zloc;
		int    i, j, k;
		
		
		
		i = (int)floor((posx - XX[0]) * XX[2]);
		if(i == (RES[0] - 1)) {
			i = RES[0] - 2;
			xloc = 1.0;
		}
		else
			
			xloc = ((posx - XX[0])* XX[2]) - i;
	
		
		j = (int)floor((posy - YY[0]) * YY[2]);
		if(j == (RES[1] - 1)) {
			j = RES[1] - 2;
			yloc = 1.0;
		}
		else
			
			yloc = ((posy - YY[0])* YY[2]) - j;
			
			
		
		k = (int)floor((posz - ZZ[0]) * ZZ[2]);
		if(k == (RES[2] - 1)) {
			k = RES[2] - 2;
			zloc = 1.0;
		}
		else
			
			zloc = ((posz - ZZ[0]) * ZZ[2]) - k;
		
		
		dXdt[0] = 0.00;
		dXdt[1] = 0.00;
		dXdt[2] = 0.00;	
		
		/////////////////////////////////////////////
		// Linear Interpolation coefficients 
		//V1 = (tloc)*(xloc)*(yloc)*(zloc) 
		//V2 = (tloc)*(xloc)*(yloc)*(1-zloc)
		//V3 = (tloc)*(xloc)*(1-yloc)*(zloc)
		//V4 = (tloc)*(xloc)*(1-yloc)*(1-zloc)
		//V5 = (tloc)*(1-xloc)*(yloc)*(zloc) 
		//V6 = (tloc)*(1-xloc)*(yloc)*(1-zloc)
		//V7 = (tloc)*(1-xloc)*(1-yloc)*(zloc) 
		//V8 = (tloc)*(1-xloc)*(1-yloc)*(1-zloc)
		//V9 = (1-tloc)*(xloc)*(yloc)*(zloc) 
		//V10 = (1-tloc)*(xloc)*(yloc)*(1-zloc)
		//V11 = (1-tloc)*(xloc)*(1-yloc)*(zloc) 
		//V12 = (1-tloc)*(xloc)*(1-yloc)*(1-zloc)
		//V13 = (1-tloc)*(1-xloc)*(yloc)*(zloc) 
		//V14 = (1-tloc)*(1-xloc)*(yloc)*(1-zloc)
		//V15 = (1-tloc)*(1-xloc)*(1-yloc)*(zloc) 
		//V16 = (1-tloc)*(1-xloc)*(1-yloc)*(1-zloc)
		///////////////////////////////////////////////
		
		double f000, f001, f010, f011;
		double f100, f101, f110, f111;
		
		// Vertices of space-time voxel for U with x and y axis only
		f000 = vel.u0[getindex(i,j,k)];
		f001 = vel.u0[getindex(i,j+1,k)];
		f010 = vel.u0[getindex(i+1,j,k)];
		f011 = vel.u0[getindex(i+1,j+1,k)];
		
		f100 = vel.u1[getindex(i,j,k)];
		f101 = vel.u1[getindex(i,j+1,k)];
		f110 = vel.u1[getindex(i+1,j,k)];
		f111 = vel.u1[getindex(i+1,j+1,k)];
		
		
		// Linear interpolation for U with x and y axis only
		
		dXdt[0] =  ( (1-tloc)*(1-xloc)*(1-yloc)*(1-zloc)*f000 ) + ( (1-tloc)*(1-xloc)*(yloc)*(1-zloc)*f001 ) + ( (1-tloc)*(xloc)*(1-yloc)*(1-zloc)*f010 ) + ( (1-tloc)*(xloc)*(yloc)*(1-zloc)*f011 ) +  ( (tloc)*(1-xloc)*(1-yloc)*(1-zloc)*f100 ) + ( (tloc)*(1-xloc)*(yloc)*(1-zloc)*f101 ) + ( (tloc)*(xloc)*(1-yloc)*(1-zloc)*f110 ) + ( (tloc)*(xloc)*(yloc)*(1-zloc)*f111 );
		
		
		// Vertices of space-time voxel for U
		f000 = vel.u0[getindex(i,j,k+1)];
		f001 = vel.u0[getindex(i,j+1,k+1)];
		f010 = vel.u0[getindex(i+1,j,k+1)];
		f011 = vel.u0[getindex(i+1,j+1,k+1)];
		
		f100 = vel.u1[getindex(i,j,k+1)];
		f101 = vel.u1[getindex(i,j+1,k+1)];
		f110 = vel.u1[getindex(i+1,j,k+1)];
		f111 = vel.u1[getindex(i+1,j+1,k+1)];
		
		// Linear interpolation for U
		dXdt[0] = dXdt[0] +  ( (1-tloc)*(1-xloc)*(1-yloc)*(zloc)*f000 ) + ( (1-tloc)*(1-xloc)*(yloc)*(zloc)*f001 ) + ( (1-tloc)*(xloc)*(1-yloc)*(zloc)*f010 ) + ( (1-tloc)*(xloc)*(yloc)*(zloc)*f011 ) +  ( (tloc)*(1-xloc)*(1-yloc)*(zloc)*f100 ) + ( (tloc)*(1-xloc)*(yloc)*(zloc)*f101 ) + ( (tloc)*(xloc)*(1-yloc)*(zloc)*f110 ) + ( (tloc)*(xloc)*(yloc)*(zloc)*f111 );
		
		
		
		
		
		
		// Vertices of space-time voxel for V with x and y axis only
		f000 = vel.v0[getindex(i,j,k)];
		f001 = vel.v0[getindex(i,j+1,k)];
		f010 = vel.v0[getindex(i+1,j,k)];
		f011 = vel.v0[getindex(i+1,j+1,k)];
		
		f100 = vel.v1[getindex(i,j,k)];
		f101 = vel.v1[getindex(i,j+1,k)];
		f110 = vel.v1[getindex(i+1,j,k)];
		f111 = vel.v1[getindex(i+1,j+1,k)];
		
		// Linear interpolation for V with x and y axis only ..
		
		dXdt[1] =  ( (1-tloc)*(1-xloc)*(1-yloc)*(1-zloc)*f000 ) + ( (1-tloc)*(1-xloc)*(yloc)*(1-zloc)*f001 ) + ( (1-tloc)*(xloc)*(1-yloc)*(1-zloc)*f010 ) + ( (1-tloc)*(xloc)*(yloc)*(1-zloc)*f011 ) +  ( (tloc)*(1-xloc)*(1-yloc)*(1-zloc)*f100 ) + ( (tloc)*(1-xloc)*(yloc)*(1-zloc)*f101 ) + ( (tloc)*(xloc)*(1-yloc)*(1-zloc)*f110 ) + ( (tloc)*(xloc)*(yloc)*(1-zloc)*f111 );
		
		// Vertices of space-time voxel for V
		f000 = vel.v0[getindex(i,j,k+1)];
		f001 = vel.v0[getindex(i,j+1,k+1)];
		f010 = vel.v0[getindex(i+1,j,k+1)];
		f011 = vel.v0[getindex(i+1,j+1,k+1)];
		
		f100 = vel.v1[getindex(i,j,k+1)];
		f101 = vel.v1[getindex(i,j+1,k+1)];
		f110 = vel.v1[getindex(i+1,j,k+1)];
		f111 = vel.v1[getindex(i+1,j+1,k+1)];
		
		// Linear interpolation for V
		
		dXdt[1] = dXdt[1] +  ( (1-tloc)*(1-xloc)*(1-yloc)*(zloc)*f000 ) + ( (1-tloc)*(1-xloc)*(yloc)*(zloc)*f001 ) + ( (1-tloc)*(xloc)*(1-yloc)*(zloc)*f010 ) + ( (1-tloc)*(xloc)*(yloc)*(zloc)*f011 ) +  ( (tloc)*(1-xloc)*(1-yloc)*(zloc)*f100 ) + ( (tloc)*(1-xloc)*(yloc)*(zloc)*f101 ) + ( (tloc)*(xloc)*(1-yloc)*(zloc)*f110 ) + ( (tloc)*(xloc)*(yloc)*(zloc)*f111 );
		

	// Vertices of space-time voxel for W with x and y axis only
		f000 = vel.w0[getindex(i,j,k)];
		f001 = vel.w0[getindex(i,j+1,k)];
		f010 = vel.w0[getindex(i+1,j,k)];
		f011 = vel.w0[getindex(i+1,j+1,k)];
		
		f100 = vel.w1[getindex(i,j,k)];
		f101 = vel.w1[getindex(i,j+1,k)];
		f110 = vel.w1[getindex(i+1,j,k)];
		f111 = vel.w1[getindex(i+1,j+1,k)];
		
		// Linear interpolation for W with x and y axis only ..
		
		dXdt[2] =  ( (1-tloc)*(1-xloc)*(1-yloc)*(1-zloc)*f000 ) + ( (1-tloc)*(1-xloc)*(yloc)*(1-zloc)*f001 ) + ( (1-tloc)*(xloc)*(1-yloc)*(1-zloc)*f010 ) + ( (1-tloc)*(xloc)*(yloc)*(1-zloc)*f011 ) +  ( (tloc)*(1-xloc)*(1-yloc)*(1-zloc)*f100 ) + ( (tloc)*(1-xloc)*(yloc)*(1-zloc)*f101 ) + ( (tloc)*(xloc)*(1-yloc)*(1-zloc)*f110 ) + ( (tloc)*(xloc)*(yloc)*(1-zloc)*f111 );
		
		// Vertices of space-time voxel for W
		f000 = vel.w0[getindex(i,j,k+1)];
		f001 = vel.w0[getindex(i,j+1,k+1)];
		f010 = vel.w0[getindex(i+1,j,k+1)];
		f011 = vel.w0[getindex(i+1,j+1,k+1)];
		
		f100 = vel.w1[getindex(i,j,k+1)];
		f101 = vel.w1[getindex(i,j+1,k+1)];
		f110 = vel.w1[getindex(i+1,j,k+1)];
		f111 = vel.w1[getindex(i+1,j+1,k+1)];
		
		// Linear interpolation for W
		
		dXdt[2] = dXdt[2] +  ( (1-tloc)*(1-xloc)*(1-yloc)*(zloc)*f000 ) + ( (1-tloc)*(1-xloc)*(yloc)*(zloc)*f001 ) + ( (1-tloc)*(xloc)*(1-yloc)*(zloc)*f010 ) + ( (1-tloc)*(xloc)*(yloc)*(zloc)*f011 ) +  ( (tloc)*(1-xloc)*(1-yloc)*(zloc)*f100 ) + ( (tloc)*(1-xloc)*(yloc)*(zloc)*f101 ) + ( (tloc)*(xloc)*(1-yloc)*(zloc)*f110 ) + ( (tloc)*(xloc)*(yloc)*(zloc)*f111 );
	}
	

}

__device__ int TestOutDomain_dev(double x, double y) {

		if(x < (XX_Data[0] - TINY) || x > (XX_Data[1] + TINY)) {

			return 1;
						
		} else if(y < (YY_Data[0] - TINY) || y > (YY_Data[1] + TINY)) {
			
			return 1;
		
		} else {
		
			return  0;
		
		}

}

__device__ int TestOutDomain3D_dev(double x, double y, double z) {

		if(x < (XX_Data[0] - TINY) || x > (XX_Data[1] + TINY)) {

			return 1;
						
		} else if(y < (YY_Data[0] - TINY) || y > (YY_Data[1] + TINY)) {
			
			return 1;
		
		} else if(z < (ZZ_Data[0] - TINY) || z > (ZZ_Data[1] + TINY))  {
		
			return 1;
		
		} else {
		
			return  0;
		
		}

}

void readveldata(int ss) {

	if (Data_MeshType == CARTESIAN) {

		char Data_BinFilePath[LONGSTRING];
		FILE *Data_BinFileID;

		int i,j,k;
		int a;
		
		ss = (ss % (Data_TRes - 1));
		
		if(ss < 0) {
			ss += Data_TRes - 1;
		}
		
		a = Data_SuffixTMin + (ss*Data_SuffixTDelta);
		
		
	
		sprintf(Data_BinFilePath, "%s%s_vel.%d.bin",Path_Data, Data_InFilePrefix , a);
		if((Data_BinFileID = fopen(Data_BinFilePath, "rb")) == NULL) {
			fprintf(stderr, "Could not open file %s \n", Data_BinFilePath);
			exit(1);
		}
	
		 //Read time stamp 
		double time1;
		if(fread(&time1, sizeof(double), 1, Data_BinFileID) < 1) {
			fprintf(stderr,"Could not read time stamp from file %s \n", Data_BinFilePath);
			exit(1);
		}
		for (i = 0; i < N_Vel; i++) {	

			velocity.time0[i] = time1;

		}
		
		printf("Loading velocity data from %s_vel.%d.bin (time stamp = %f)\n", Data_InFilePrefix , a, (double)time1);
		fflush(stdout);
		
		double velu_d, velv_d, velw_d;
		for(k = 0; k < Vel_CartMesh.ZRes; k++) {
			for(j = 0; j < Vel_CartMesh.YRes; j++) {
				for(i = 0; i < Vel_CartMesh.XRes; i++) {
					if((fread(&velu_d, sizeof(double), 1, Data_BinFileID) < 1) ||
					   (fread(&velv_d, sizeof(double), 1, Data_BinFileID) < 1) ||
					   (fread(&velw_d, sizeof(double), 1, Data_BinFileID) < 1)) { 
						fprintf(stderr,"Could not load velocity data for index [%d][%d][%d] from the file \n\n", i, j, k);
						exit(1);
					}
					else {
						velocity.u0[getindexhost(i,j,k)] = velu_d;
						velocity.v0[getindexhost(i,j,k)] = velv_d;
						velocity.w0[getindexhost(i,j,k)] = velw_d;
					
					}
				}
			}
		}

		/* Close file */
		fclose(Data_BinFileID);


		//load our next file

		sprintf(Data_BinFilePath, "%s%s_vel.%d.bin",Path_Data, Data_InFilePrefix , a + Data_SuffixTDelta);
		if((Data_BinFileID = fopen(Data_BinFilePath, "rb")) == NULL) {
			fprintf(stderr, "Could not open file %s \n", Data_BinFilePath);
			exit(1);
		}
	
		/* Read time stamp */
	
		if(fread(&time1, sizeof(double), 1, Data_BinFileID) < 1) {
			fprintf(stderr,"Could not read time stamp from file %s \n", Data_BinFilePath);
			exit(1);
		}
	
		for (i = 0; i < N_Vel ; i++) {
			velocity.time1[i] = time1;
		}	

		printf("Loading velocity data from %s_vel.%d.bin (time stamp = %f)\n", Data_InFilePrefix , a+ Data_SuffixTDelta, (double)time1);
		fflush(stdout);
		


		for(k = 0; k < Vel_CartMesh.ZRes; k++) {
			for(j = 0; j < Vel_CartMesh.YRes; j++) {
				for(i = 0; i < Vel_CartMesh.XRes; i++) {
					if((fread(&velu_d, sizeof(double), 1, Data_BinFileID) < 1) ||
					   (fread(&velv_d, sizeof(double), 1, Data_BinFileID) < 1) ||
					   (fread(&velw_d, sizeof(double), 1, Data_BinFileID) < 1)) { 
						fprintf(stderr,"Could not load velocity data for index [%d][%d][%d] from the file \n\n", i, j, k);
						exit(1);
					} else {
							velocity.u1[getindexhost(i,j,k)] = velu_d;
							velocity.v1[getindexhost(i,j,k)] = velv_d;
							velocity.w1[getindexhost(i,j,k)] = velw_d;


					}
				}
			}
		}
					

		/* Close file */
		fclose(Data_BinFileID);

		//printf("Data Loaded from files into variables ...... \n\n");

	} else {	 // Unstructured Grid 
		char Data_BinFilePath[LONGSTRING];
		FILE *Data_BinFileID;

		int i;
		int a;
		
		ss = (ss % (Data_TRes - 1));
		
		if(ss < 0) {
			ss += (Data_TRes - 1);
		}
		
		a = Data_SuffixTMin + (ss*Data_SuffixTDelta);
	
		sprintf(Data_BinFilePath, "%s%s_vel.%d.bin", Path_Data, Data_InFilePrefix , a);
		if((Data_BinFileID = fopen(Data_BinFilePath, "rb")) == NULL) {
			fprintf(stderr, "Could not open file %s \n", Data_BinFilePath);
			exit(1);
		}
	
		 //Read time stamp 
		double time1;
		if(fread(&time1, sizeof(double), 1, Data_BinFileID) < 1) {
			fprintf(stderr,"Could not read time stamp from file %s \n", Data_BinFilePath);
			exit(1);
		}
		for (i = 0; i < Vel_MeshNumNodes; i++) {	

			velocity.time0[i] = time1;

		}
		printf("Loading velocity data from %s_vel.%d.bin (time stamp = %f)\n", Data_InFilePrefix , a, (double)time1);
		fflush(stdout);

		double velu_d, velv_d, velw_d;
		

		for(i = 0; i < Vel_MeshNumNodes; i++) {
			if((fread(&velu_d, sizeof(double), 1, Data_BinFileID) < 1) || (fread(&velv_d, sizeof(double), 1, Data_BinFileID) < 1) || (fread(&velw_d, sizeof(double), 1, Data_BinFileID) < 1)) { 
				fprintf(stderr,"Could not load velocity data for index [%d] from the file \n\n", i);
				exit(1);
			}
			else {
			
				
				velocity.u0[i] = velu_d;
				velocity.v0[i] = velv_d;
				velocity.w0[i] = velw_d;
			
			}
		}

		/* Close file */
		fclose(Data_BinFileID);


		//load our next file

		sprintf(Data_BinFilePath, "%s%s_vel.%d.bin",Path_Data, Data_InFilePrefix , a + Data_SuffixTDelta);
		if((Data_BinFileID = fopen(Data_BinFilePath, "rb")) == NULL) {
			fprintf(stderr, "Could not open file %s \n", Data_BinFilePath);
			exit(1);
		}
	
		/* Read time stamp */
	
		if(fread(&time1, sizeof(double), 1, Data_BinFileID) < 1) {
			fprintf(stderr,"Could not read time stamp from file %s \n", Data_BinFilePath);
			exit(1);
		}
	
		for (i = 0; i <  Vel_MeshNumNodes; i++) {
			velocity.time1[i] = time1;
		}	


		printf("Loading velocity data from %s_vel.%d.bin (time stamp = %f)\n",Data_InFilePrefix , a + Data_SuffixTDelta , (double)time1);
		
		fflush(stdout);

		for(i = 0; i < Vel_MeshNumNodes; i++) {
			if((fread(&velu_d, sizeof(double), 1, Data_BinFileID) < 1) ||  (fread(&velv_d, sizeof(double), 1, Data_BinFileID) < 1) || (fread(&velw_d, sizeof(double), 1, Data_BinFileID) < 1)) { 
				fprintf(stderr,"Could not load velocity data for index [%d] from the file \n\n", i);
				exit(1);
			} else {
			
				
			
				velocity.u1[i] = velu_d;
				velocity.v1[i] = velv_d;
				velocity.w1[i] = velw_d;
				
				//printf("u1 = %f \n", velu_d);

			}
		}
					

		/* Close file */
		fclose(Data_BinFileID);

		//printf("Data Loaded from files into variables ...... \n\n");
	
	}



}

__global__ void LocalSearch2D(double *posx, double *posy, Element MeshElementArray_device, Node_double MeshNodeArray_double_device, int *eid, int ss, int offset, double *r, double *s, int *integrate) {

	int tid;
	
	// Get thread index
	tid=(blockIdx.x*blockDim.x)+threadIdx.x; 
	if( tid >= ss) { // Redundant thread ...
		return; 
	} else {
	
		tid = tid + offset;
	
		if(integrate[tid] == 1) {
		
		double x0,y0,x1,y1,x2,y2;
		double d,V;
		int ID = eid[tid];
		
		if(ID == -1) { // Point already left.
		
			return;
		
		} else {
			while(1) {
				
				// Physical coordinates of nodes of the test element
				x0 = MeshNodeArray_double_device.x[MeshElementArray_device.Node1[ID]];
				
				y0 = MeshNodeArray_double_device.y[MeshElementArray_device.Node1[ID]];

				x1 = MeshNodeArray_double_device.x[MeshElementArray_device.Node2[ID]];
				y1 = MeshNodeArray_double_device.y[MeshElementArray_device.Node2[ID]];

				x2 = MeshNodeArray_double_device.x[MeshElementArray_device.Node3[ID]];
				y2 = MeshNodeArray_double_device.y[MeshElementArray_device.Node3[ID]];	

			
				// Determinant of mapping from natural to physical coordinates of test element 
				V = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
			
				// Natural coordinates of point to be interpolated
				r[tid] = ( (y2 - y0) * (posx[tid] - x0) + (x0 - x2) * (posy[tid] - y0)) / V;
				s[tid] = ( (y0 - y1) * (posx[tid] - x0) + (x1 - x0) * (posy[tid] - y0)) / V;
			
			
				d = fmin(r[tid], fmin(s[tid], 1 - r[tid] - s[tid]));

				if((d + TINY) >= 0) { // Point inside the element .....
					eid[tid] = ID;
					return;

				} else { // Point is outside ....

					if(fabs(r[tid] - d) < TINY) 
						ID = MeshElementArray_device.Neighborindex1[ID];
					else if(fabs(s[tid] - d) < TINY) 
						ID = MeshElementArray_device.Neighborindex2[ID];
					else if(fabs(d - (1 - r[tid] - s[tid])) < TINY) 
						ID = MeshElementArray_device.Neighborindex3[ID];

					if (ID == -1) {
						eid[tid] = ID;
						return;
					}
				}

			}
			
		}
		
		
		} else
			return;
	
	}
}

// 3D local search ...

__global__ void LocalSearch3D(double *posx, double *posy, double *posz, Element MeshElementArray_device, Node_double MeshNodeArray_double_device, int *eid, int ss, int offset, double *r, double *s, double *t, int *integrate) {

	int tid;
	
	// Get thread index
	tid=(blockIdx.x*blockDim.x)+threadIdx.x; 
	
	
	if( tid >= ss) { // Redundant thread ...
		return; 
	} else {
		
		tid = tid + offset;
	
		if (integrate[tid] == 1) { // Only Do Local Search If a point has to be integrated over the interval.
		
			double x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3;
			double d,V;
			int ID = eid[tid];
			if(ID == -1) { // Point already left.
		
				return;
		
			} else {
		
			
				while(1) {
				
					// Physical coordinates of nodes of the test element
					x0 = MeshNodeArray_double_device.x[MeshElementArray_device.Node1[ID]];
					y0 = MeshNodeArray_double_device.y[MeshElementArray_device.Node1[ID]];
					z0 = MeshNodeArray_double_device.z[MeshElementArray_device.Node1[ID]];

					x1 = MeshNodeArray_double_device.x[MeshElementArray_device.Node2[ID]];
					y1 = MeshNodeArray_double_device.y[MeshElementArray_device.Node2[ID]];
					z1 = MeshNodeArray_double_device.z[MeshElementArray_device.Node2[ID]];

					x2 = MeshNodeArray_double_device.x[MeshElementArray_device.Node3[ID]];
					y2 = MeshNodeArray_double_device.y[MeshElementArray_device.Node3[ID]];
					z2 = MeshNodeArray_double_device.z[MeshElementArray_device.Node3[ID]];
				
					x3 = MeshNodeArray_double_device.x[MeshElementArray_device.Node4[ID]];
					y3 = MeshNodeArray_double_device.y[MeshElementArray_device.Node4[ID]];
					z3 = MeshNodeArray_double_device.z[MeshElementArray_device.Node4[ID]];
			
					// Determinant of mapping from natural to physical coordinates of test element
					V = (x1 - x0) * ((y2 - y0) * (z3 - z0) - (z2 - z0) * (y3 - y0)) + (x2 - x0) * ((y0 - y1) * (z3 - z0) - (z0 - z1) * (y3 - y0)) +
					(x3 - x0) * ((y1 - y0) * (z2 - z0) - (z1 - z0) * (y2 - y0));
			
			
					// Natural coordinates of point to be interpolated
					r[tid] = ( (((z3 - z0) * (y2 - y3) - (z2 - z3) * (y3 - y0)) * (posx[tid] - x0)) + ( ((x3 - x0) * (z2 - z3) - (x2 - x3) * (z3 - z0)) * 							(posy[tid] - y0)) + ( ((y3 - y0) * (x2 - x3) - (y2 - y3) * (x3 - x0)) *(posz[tid] - z0) ))/ V;
			
					s[tid] = ( (((z3 - z0) * (y0 - y1) - (z0 - z1) * (y3 - y0)) * (posx[tid] - x0)) + (((x3 - x0) * (z0 - z1) - (x0 - x1) * (z3 - z0)) * 							(posy[tid] - y0)) + (((y3 - y0) * (x0 - x1) - (y0 - y1) * (x3 - x0)) * (posz[tid] - z0)) ) / V;
			
			
					t[tid] = ( (((z1 - z2) * (y0 - y1) - (z0 - z1) * (y1 - y2)) * (posx[tid] - x0)) + (((x1 - x2) * (z0 - z1) - (x0 - x1) * (z1 - z2)) * 							(posy[tid] - y0)) + (((y1 - y2) * (x0 - x1) - (y0 - y1) * (x1 - x2)) * (posz[tid] - z0)) ) / V;
			
			
					d = fmin(r[tid], fmin(s[tid], fmin(t[tid], 1 - r[tid] - s[tid] - t[tid])));

				
					if((d + TINY) >= 0) { // Point inside the element .....
						eid[tid] = ID;
						return;

					} else { // Point is outside ....

						if(fabs(r[tid] - d) < TINY) 
							ID = MeshElementArray_device.Neighborindex1[ID];
						else if(fabs(s[tid] - d) < TINY) 
							ID = MeshElementArray_device.Neighborindex2[ID];
						else if(fabs(t[tid] - d) < TINY) 
							ID = MeshElementArray_device.Neighborindex3[ID];
						else if(fabs(d - (1 - r[tid] - s[tid] - t[tid])) < TINY) 
							ID = MeshElementArray_device.Neighborindex4[ID];
				
						if (ID == -1) {
					
							eid[tid] = ID;
							return;
						
						}
						
					}

				}
			
			}
		
		
		}
		else
			return;
	
	}
}

__device__ inline int get_local_search_2D(double x, double y, Element MeshElementArray_device, Node_double MeshNodeArray_double_device, int eid, double r, double s) {

		double x0,y0,x1,y1,x2,y2;
		double d;

		if (eid == -1) { // Point is outside the domain

			return eid;
		}

		while(1) {

			// Physical coordinates of nodes of the test element
			x0 = MeshNodeArray_double_device.x[MeshElementArray_device.Node1[eid]];
			y0 = MeshNodeArray_double_device.y[MeshElementArray_device.Node1[eid]];

			x1 = MeshNodeArray_double_device.x[MeshElementArray_device.Node2[eid]];
			y1 = MeshNodeArray_double_device.y[MeshElementArray_device.Node2[eid]];

			x2 = MeshNodeArray_double_device.x[MeshElementArray_device.Node3[eid]];
			y2 = MeshNodeArray_double_device.y[MeshElementArray_device.Node3[eid]];	

			
			// Determinant of mapping from natural to physical coordinates of test element 
			//	V = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
			
			// Natural coordinates of point to be interpolated
			r = ( (y2 - y0) * (x - x0) + (x0 - x2) * (y - y0)) / ( (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0) );
			s = ( (y0 - y1) * (x - x0) + (x1 - x0) * (y - y0)) / ( (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0) );
			
			
			d = fmin(r, fmin(s, 1 - r - s));

			if((d + TINY) >= 0) { // Point inside the element .....
	
				return eid;

			} else { // Point is outside ....

				if(fabs(r - d) < TINY) 
					eid = MeshElementArray_device.Neighborindex1[eid];
				else if(fabs(s - d) < TINY) 
					eid = MeshElementArray_device.Neighborindex2[eid];
				else if(fabs(d - (1 - r - s)) < TINY) 
					eid = MeshElementArray_device.Neighborindex3[eid];
					
				if (eid == -1)
					return eid;

			}

		}

}

__device__ inline int get_local_search_3D(double x, double y, double z, Element MeshElementArray_device, Node_double MeshNodeArray_double_device, int eid, double r, double s, double t) {

		double x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3;
		double d;

		if (eid == -1) { // Point is outside the domain

			return eid;
		}

		while(1) {
			
			// Physical coordinates of nodes of the test element
			x0 = MeshNodeArray_double_device.x[MeshElementArray_device.Node1[eid]];
			y0 = MeshNodeArray_double_device.y[MeshElementArray_device.Node1[eid]];
			z0 = MeshNodeArray_double_device.z[MeshElementArray_device.Node1[eid]];

			x1 = MeshNodeArray_double_device.x[MeshElementArray_device.Node2[eid]];
			y1 = MeshNodeArray_double_device.y[MeshElementArray_device.Node2[eid]];
			z1 = MeshNodeArray_double_device.z[MeshElementArray_device.Node2[eid]];

			x2 = MeshNodeArray_double_device.x[MeshElementArray_device.Node3[eid]];
			y2 = MeshNodeArray_double_device.y[MeshElementArray_device.Node3[eid]];
			z2 = MeshNodeArray_double_device.z[MeshElementArray_device.Node3[eid]];
				
			x3 = MeshNodeArray_double_device.x[MeshElementArray_device.Node4[eid]];
			y3 = MeshNodeArray_double_device.y[MeshElementArray_device.Node4[eid]];
			z3 = MeshNodeArray_double_device.z[MeshElementArray_device.Node4[eid]];
			
			// Determinant of mapping from natural to physical coordinates of test element
			//V = (x1 - x0) * ((y2 - y0) * (z3 - z0) - (z2 - z0) * (y3 - y0)) + (x2 - x0) * ((y0 - y1) * (z3 - z0) - (z0 - z1) * (y3 - y0)) + (x3 - x0) * ((y1 - y0) * (z2 - z0) - (z1 - z0) * (y2 - y0));
			
			
			r = ( (((z3 - z0) * (y2 - y3) - (z2 - z3) * (y3 - y0)) * (x - x0)) + ( ((x3 - x0) * (z2 - z3) - (x2 - x3) * (z3 - z0)) * (y - y0)) + ( ((y3 - y0) * (x2 - x3) - (y2 - y3) * (x3 - x0)) *(z - z0) ))/ ( (x1 - x0) * ((y2 - y0) * (z3 - z0) - (z2 - z0) * (y3 - y0)) + (x2 - x0) * ((y0 - y1) * (z3 - z0) - (z0 - z1) * (y3 - y0)) + (x3 - x0) * ((y1 - y0) * (z2 - z0) - (z1 - z0) * (y2 - y0)) );
			
			s = ( (((z3 - z0) * (y0 - y1) - (z0 - z1) * (y3 - y0)) * (x - x0)) + (((x3 - x0) * (z0 - z1) - (x0 - x1) * (z3 - z0)) * (y - y0)) + (((y3 - y0) * (x0 - x1) - (y0 - y1) * (x3 - x0)) * (z - z0)) ) / ( (x1 - x0) * ((y2 - y0) * (z3 - z0) - (z2 - z0) * (y3 - y0)) + (x2 - x0) * ((y0 - y1) * (z3 - z0) - (z0 - z1) * (y3 - y0)) + (x3 - x0) * ((y1 - y0) * (z2 - z0) - (z1 - z0) * (y2 - y0)) );
			
			
			t = ( (((z1 - z2) * (y0 - y1) - (z0 - z1) * (y1 - y2)) * (x - x0)) + (((x1 - x2) * (z0 - z1) - (x0 - x1) * (z1 - z2)) * (y - y0)) + (((y1 - y2) * (x0 - x1) - (y0 - y1) * (x1 - x2)) * (z - z0)) ) / ( (x1 - x0) * ((y2 - y0) * (z3 - z0) - (z2 - z0) * (y3 - y0)) + (x2 - x0) * ((y0 - y1) * (z3 - z0) - (z0 - z1) * (y3 - y0)) + (x3 - x0) * ((y1 - y0) * (z2 - z0) - (z1 - z0) * (y2 - y0)) );
			
			
			d = fmin(r, fmin(s, fmin(t, 1 - r - s - t)));

			if((d + TINY) >= 0) { // Point inside the element .....
	
				return eid;

			} else { // Point is outside ....

				if(fabs(r - d) < TINY) 
					eid = MeshElementArray_device.Neighborindex1[eid];
				else if(fabs(s - d) < TINY) 
					eid = MeshElementArray_device.Neighborindex2[eid];
				else if(fabs(t - d) < TINY) 
					eid = MeshElementArray_device.Neighborindex3[eid];
				else if(fabs(d - (1 - r - s - t)) < TINY) 
					eid = MeshElementArray_device.Neighborindex4[eid];	
			
			
				if (eid == -1) // Point left domain
					return eid;

			}

		}

}





// Save r and s while doing local search .....

__device__ inline void GetVel_unstruct2D(double tloc, double *dXdt, VelData_double vel, Element MeshElementArray_device, int eid, double r, double s ) {

	if (eid == -1) {  // Point out of Domain ....

		dXdt[0] = 0.0;
		dXdt[1] = 0.0;
		return;

	} else {
	
		double U1t, U2t, U3t, V1t, V2t, V3t;
		// Interpolate velocity at nodes in time */
		U1t = (1 - tloc) * vel.u0[MeshElementArray_device.Node1[eid]] + tloc * vel.u1[MeshElementArray_device.Node1[eid]];
		V1t = (1 - tloc) * vel.v0[MeshElementArray_device.Node1[eid]] + tloc * vel.v1[MeshElementArray_device.Node1[eid]];

		U2t = (1 - tloc) * vel.u0[MeshElementArray_device.Node2[eid]] + tloc * vel.u1[MeshElementArray_device.Node2[eid]];
		V2t = (1 - tloc) * vel.v0[MeshElementArray_device.Node2[eid]] + tloc * vel.v1[MeshElementArray_device.Node2[eid]];
				
		U3t = (1 - tloc) * vel.u0[MeshElementArray_device.Node3[eid]] + tloc * vel.u1[MeshElementArray_device.Node3[eid]];
		V3t = (1 - tloc) * vel.v0[MeshElementArray_device.Node3[eid]] + tloc * vel.v1[MeshElementArray_device.Node3[eid]];
			
		// Get interpolate velocity at point in space using time-interpolated velocity at nodes
		dXdt[0] = U1t + (U2t - U1t) * fabs(r) + (U3t - U1t) * fabs(s);
		dXdt[1] = V1t + (V2t - V1t) * fabs(r) + (V3t - V1t) * fabs(s); 

		return;
	}
}


__device__  inline void GetVel_unstruct3D(double tloc, double *dXdt, VelData_double vel, Element MeshElementArray_device, int eid, double r, double s, double t ) {

	if (eid == -1) {  // Point out of Domain ....

		dXdt[0] = 0.0;
		dXdt[1] = 0.0;
		dXdt[2] = 0.0;
		
		return;

	} else {
	
		double U1t, U2t, U3t, U4t, V1t, V2t, V3t, V4t, W1t, W2t, W3t, W4t;
		
		// Interpolate velocity at nodes in time */
		U1t = (1 - tloc) * vel.u0[MeshElementArray_device.Node1[eid]] + tloc * vel.u1[MeshElementArray_device.Node1[eid]];
		V1t = (1 - tloc) * vel.v0[MeshElementArray_device.Node1[eid]] + tloc * vel.v1[MeshElementArray_device.Node1[eid]];
		W1t = (1 - tloc) * vel.w0[MeshElementArray_device.Node1[eid]] + tloc * vel.w1[MeshElementArray_device.Node1[eid]];

		U2t = (1 - tloc) * vel.u0[MeshElementArray_device.Node2[eid]] + tloc * vel.u1[MeshElementArray_device.Node2[eid]];
		V2t = (1 - tloc) * vel.v0[MeshElementArray_device.Node2[eid]] + tloc * vel.v1[MeshElementArray_device.Node2[eid]];
		W2t = (1 - tloc) * vel.w0[MeshElementArray_device.Node2[eid]] + tloc * vel.w1[MeshElementArray_device.Node2[eid]];
				
		U3t = (1 - tloc) * vel.u0[MeshElementArray_device.Node3[eid]] + tloc * vel.u1[MeshElementArray_device.Node3[eid]];
		V3t = (1 - tloc) * vel.v0[MeshElementArray_device.Node3[eid]] + tloc * vel.v1[MeshElementArray_device.Node3[eid]];
		W3t = (1 - tloc) * vel.w0[MeshElementArray_device.Node3[eid]] + tloc * vel.w1[MeshElementArray_device.Node3[eid]];
		
		U4t = (1 - tloc) * vel.u0[MeshElementArray_device.Node4[eid]] + tloc * vel.u1[MeshElementArray_device.Node4[eid]];
		V4t = (1 - tloc) * vel.v0[MeshElementArray_device.Node4[eid]] + tloc * vel.v1[MeshElementArray_device.Node4[eid]];
		W4t = (1 - tloc) * vel.w0[MeshElementArray_device.Node4[eid]] + tloc * vel.w1[MeshElementArray_device.Node4[eid]];
	
		// Get interpolate velocity at point in space using time-interpolated velocity at nodes
		dXdt[0] = U1t + (U2t - U1t) * fabs(r) + (U3t - U1t) * fabs(s) + (U4t - U1t) * fabs(t);
		dXdt[1] = V1t + (V2t - V1t) * fabs(r) + (V3t - V1t) * fabs(s) + (V4t - V1t) * fabs(t);
		dXdt[2] = W1t + (W2t - W1t) * fabs(r) + (W3t - W1t) * fabs(s) + (W4t - W1t) * fabs(t);
		

		return;
	}
}



__device__ inline void getvel_unstruct_2D(double *output, double t, double *dXdt, VelData_double vel, double TMin, double TMax, Element MeshElementArray_device, Node_double MeshNodeArray_double_device, int eid) {

	
		double x0,y0,x1,y1,x2,y2;
		double r,s,d,V;
		double U1t,U2t,V1t,V2t,U3t,V3t;
			


		if (eid == -1) {  // Point out of Domain ....

			dXdt[0] = 0.0;
			dXdt[1] = 0.0;
			return;

		} else {

		double tloc = (t - TMin) / (TMax - TMin); // I am not checking whether it is between tmin and tmax.

		while(1) {

			// Physical coordinates of nodes of the test element

			x0 = MeshNodeArray_double_device.x[MeshElementArray_device.Node1[eid]];
			y0 = MeshNodeArray_double_device.y[MeshElementArray_device.Node1[eid]];

			x1 = MeshNodeArray_double_device.x[MeshElementArray_device.Node2[eid]];
			y1 = MeshNodeArray_double_device.y[MeshElementArray_device.Node2[eid]];

			x2 = MeshNodeArray_double_device.x[MeshElementArray_device.Node3[eid]];
			y2 = MeshNodeArray_double_device.y[MeshElementArray_device.Node3[eid]];	

			
			// Determinant of mapping from natural to physical coordinates of test element 
			V = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
			
			// Natural coordinates of point to be interpolated
			r = ( (y2 - y0) * (output[0] - x0) + (x0 - x2) * (output[1] - y0)) / V;
			s = ( (y0 - y1) * (output[0] - x0) + (x1 - x0) * (output[1] - y0)) / V;
			
			
			d = fmin(r, fmin(s, 1 - r - s));

			if((d + TINY) >= 0) {   // Point inside test element */
				// Interpolate velocity at nodes in time */
				U1t = (1 - tloc) * vel.u0[MeshElementArray_device.Node1[eid]] + tloc * vel.u1[MeshElementArray_device.Node1[eid]];
				V1t = (1 - tloc) * vel.v0[MeshElementArray_device.Node1[eid]] + tloc * vel.v1[MeshElementArray_device.Node1[eid]];

				U2t = (1 - tloc) * vel.u0[MeshElementArray_device.Node2[eid]] + tloc * vel.u1[MeshElementArray_device.Node2[eid]];
				V2t = (1 - tloc) * vel.v0[MeshElementArray_device.Node2[eid]] + tloc * vel.v1[MeshElementArray_device.Node2[eid]];
				
				U3t = (1 - tloc) * vel.u0[MeshElementArray_device.Node3[eid]] + tloc * vel.u1[MeshElementArray_device.Node3[eid]];
				V3t = (1 - tloc) * vel.v0[MeshElementArray_device.Node3[eid]] + tloc * vel.v1[MeshElementArray_device.Node3[eid]];

				// Get interpolate velocity at point in space using time-interpolated velocity at nodes*/
				dXdt[0] = U1t + (U2t - U1t) * fabs(r) + (U3t - U1t) * fabs(s);
				dXdt[1] = V1t + (V2t - V1t) * fabs(r) + (V3t - V1t) * fabs(s);

				return;
		

			} else { // Point not inside the test element ......


				if(fabs(r - d) < TINY) 
					eid = MeshElementArray_device.Neighborindex1[eid];
				else if(fabs(s - d) < TINY) 
					eid = MeshElementArray_device.Neighborindex2[eid];
				else if(fabs(d - (1 - r - s)) < TINY) 
					eid = MeshElementArray_device.Neighborindex3[eid];
				
				if (eid == -1)
					return;
			}

		}
	}
			
}


int TestOutsideCartVelDomain(double X[3]) {
	/***  
	 inside data bounding box --> return 0
	 outside data bounding box --> return 1
	 ***/ 
	
	int imin, jmin, kmin, imax, jmax, kmax;
	
	imin = (int)floor((X[0] - Vel_CartMesh.XMin) / Vel_CartMesh.XDelta);
	jmin = (int)floor((X[1] - Vel_CartMesh.YMin) / Vel_CartMesh.YDelta);
	if(Dimensions == 3) 
		kmin = (int)floor((X[2] - Vel_CartMesh.ZMin) / Vel_CartMesh.ZDelta);
	else 
		kmin = 0;
	imax = (int)ceil((X[0] - Vel_CartMesh.XMin) / Vel_CartMesh.XDelta);
	jmax = (int)ceil((X[1] - Vel_CartMesh.YMin) / Vel_CartMesh.YDelta);
	if(Dimensions == 3) 
		kmax = (int)ceil((X[2] - Vel_CartMesh.ZMin) / Vel_CartMesh.ZDelta);
	else 
		kmax = 0;
	
	if(imin < 0 || imax > (Vel_CartMesh.XRes - 1) || jmin < 0 || jmax > (Vel_CartMesh.YRes - 1) 
	   || kmin < 0 || kmax > (Vel_CartMesh.ZRes - 1)) {
		return 1;
	}
	else 
		return 0;
}
	






void GetVelocity_Cartesian(double tq, LagrangianPoint *pt, double *dXdt) {
	
	double xloc, yloc, zloc, tloc;
	int    i, j, k;
	double V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12, V13, V14, V15, V16;
	double f0000, f0001, f0010, f0011, f0100, f0101, f0110, f0111, f1000, f1001, f1010, f1011, f1100, f1101, f1110, f1111;
	

	/* check outside domain */
	if(TestOutsideCartVelDomain(pt->X)) {
		dXdt[0] = 0.0;
		dXdt[1] = 0.0;
		dXdt[2] = 0.0;
		pt->LeftDomain = 1;
	}
	else {  
		/* Set "local" coordinates (relative to space-time voxel) */
		/* t */
		tloc = (tq - Data_LoadedTMin) / (Data_LoadedTMax - Data_LoadedTMin);
		if((tloc > 1 + TINY) || (tloc < 0 - TINY))
			FatalError("tloc must be between 0 and 1");
		
		/* x */
		i = (int)floor((pt->X[0] - Vel_CartMesh.XMin) / Vel_CartMesh.XDelta);
		if(i == (Vel_CartMesh.XRes - 1)) {
			i = Vel_CartMesh.XRes - 2;
			xloc = 1.0;
		}
		else
			xloc = (pt->X[0] - Vel_CartMesh.XMin - i * Vel_CartMesh.XDelta) / Vel_CartMesh.XDelta;
			
			
		//printf("velxmin = %0.9f and velxdelta = %0.9f \n",Vel_CartMesh.XMin, Vel_CartMesh.XDelta);
		
		/* y */
		j = (int)floor((pt->X[1] - Vel_CartMesh.YMin) / Vel_CartMesh.YDelta);
		if(j == (Vel_CartMesh.YRes - 1)) {
			j = Vel_CartMesh.YRes - 2;
			yloc = 1.0;
		}
		else
			yloc = (pt->X[1] - Vel_CartMesh.YMin - j * Vel_CartMesh.YDelta) / Vel_CartMesh.YDelta;
		
		//printf("yloc = %0.9f \n", yloc);
		/* z */
		if(Dimensions == 3) {
			k = (int)floor((pt->X[2] - Vel_CartMesh.ZMin) / Vel_CartMesh.ZDelta);
			if(k == (Vel_CartMesh.ZRes - 1)) {
				k = Vel_CartMesh.ZRes - 2;
				zloc = 1.0;
			}
			else 
				zloc = (pt->X[2] - Vel_CartMesh.ZMin - k * Vel_CartMesh.ZDelta) / Vel_CartMesh.ZDelta;
		}
		else {
			k = 0;
			zloc = 0.0;
		}
		
		/* Linear Interpolation coefficients */
		V1 = (tloc)*(xloc)*(yloc)*(zloc);
		V2 = (tloc)*(xloc)*(yloc)*(1-zloc);
		V3 = (tloc)*(xloc)*(1-yloc)*(zloc);
		V4 = (tloc)*(xloc)*(1-yloc)*(1-zloc);
		V5 = (tloc)*(1-xloc)*(yloc)*(zloc);
		V6 = (tloc)*(1-xloc)*(yloc)*(1-zloc);
		V7 = (tloc)*(1-xloc)*(1-yloc)*(zloc);
		V8 = (tloc)*(1-xloc)*(1-yloc)*(1-zloc);
		V9 = (1-tloc)*(xloc)*(yloc)*(zloc);
		V10 = (1-tloc)*(xloc)*(yloc)*(1-zloc);
		V11 = (1-tloc)*(xloc)*(1-yloc)*(zloc);
		V12 = (1-tloc)*(xloc)*(1-yloc)*(1-zloc);
		V13 = (1-tloc)*(1-xloc)*(yloc)*(zloc);
		V14 = (1-tloc)*(1-xloc)*(yloc)*(1-zloc);
		V15 = (1-tloc)*(1-xloc)*(1-yloc)*(zloc);
		V16 = (1-tloc)*(1-xloc)*(1-yloc)*(1-zloc);
		
		/* Vertices of space-time voxel for U */
		if(Dimensions == 3) {
			f0000 = velocity.u0[getindexhost(i,j,k)];
			f0001 = velocity.u0[getindexhost(i,j,k+1)];
			f0010 = velocity.u0[getindexhost(i,j+1,k)];
			f0011 = velocity.u0[getindexhost(i,j+1,k+1)];
			f0100 = velocity.u0[getindexhost(i+1,j,k)];
			f0101 = velocity.u0[getindexhost(i+1,j,k+1)];
			f0110 = velocity.u0[getindexhost(i+1,j+1,k)];
			f0111 = velocity.u0[getindexhost(i+1,j+1,k+1)];
			
			f1000 = velocity.u1[getindexhost(i,j,k)];
			f1001 = velocity.u1[getindexhost(i,j,k+1)];
			f1010 = velocity.u1[getindexhost(i,j+1,k)];
			f1011 = velocity.u1[getindexhost(i,j+1,k+1)];
			f1100 = velocity.u1[getindexhost(i+1,j,k)];
			f1101 = velocity.u1[getindexhost(i+1,j,k+1)];
			f1110 = velocity.u1[getindexhost(i+1,j+1,k)];
			f1111 = velocity.u1[getindexhost(i+1,j+1,k+1)];
			
			
		}
		else {
		
			f0000 = velocity.u0[getindexhost(i,j,k)];
			f0001 = velocity.u0[getindexhost(i,j,k)];
			f0010 = velocity.u0[getindexhost(i,j+1,k)];
			f0011 = velocity.u0[getindexhost(i,j+1,k)];
			f0100 = velocity.u0[getindexhost(i+1,j,k)];
			f0101 = velocity.u0[getindexhost(i+1,j,k)];
			f0110 = velocity.u0[getindexhost(i+1,j+1,k)];
			f0111 = velocity.u0[getindexhost(i+1,j+1,k)];
			
			f1000 = velocity.u1[getindexhost(i,j,k)];
			f1001 = velocity.u1[getindexhost(i,j,k)];
			f1010 = velocity.u1[getindexhost(i,j+1,k)];
			f1011 = velocity.u1[getindexhost(i,j+1,k)];
			f1100 = velocity.u1[getindexhost(i+1,j,k)];
			f1101 = velocity.u1[getindexhost(i+1,j,k)];
			f1110 = velocity.u1[getindexhost(i+1,j+1,k)];
			f1111 = velocity.u1[getindexhost(i+1,j+1,k)];
		
			
		}
		
		/* Linear interpolation for U */
		dXdt[0] = f0000*V16 + f0001*V15 + f0010*V14 + f0011*V13 + f0100*V12 + f0101*V11 + f0110*V10
		+ f0111*V9 + f1000*V8 + f1001*V7 + f1010*V6 + f1011*V5 + f1100*V4 + f1101*V3 + f1110*V2 + f1111*V1; 
		
		/* Vertices of space-time voxel for V */
		if(Dimensions == 3) {
		
			f0000 = velocity.v0[getindexhost(i,j,k)];
			f0001 = velocity.v0[getindexhost(i,j,k+1)];
			f0010 = velocity.v0[getindexhost(i,j+1,k)];
			f0011 = velocity.v0[getindexhost(i,j+1,k+1)];
			f0100 = velocity.v0[getindexhost(i+1,j,k)];
			f0101 = velocity.v0[getindexhost(i+1,j,k+1)];
			f0110 = velocity.v0[getindexhost(i+1,j+1,k)];
			f0111 = velocity.v0[getindexhost(i+1,j+1,k+1)];
			
			f1000 = velocity.v1[getindexhost(i,j,k)];
			f1001 = velocity.v1[getindexhost(i,j,k+1)];
			f1010 = velocity.v1[getindexhost(i,j+1,k)];
			f1011 = velocity.v1[getindexhost(i,j+1,k+1)];
			f1100 = velocity.v1[getindexhost(i+1,j,k)];
			f1101 = velocity.v1[getindexhost(i+1,j,k+1)];
			f1110 = velocity.v1[getindexhost(i+1,j+1,k)];
			f1111 = velocity.v1[getindexhost(i+1,j+1,k+1)];
		
		}
		else {
		
			
			f0000 = velocity.v0[getindexhost(i,j,k)];
			f0001 = velocity.v0[getindexhost(i,j,k)];
			f0010 = velocity.v0[getindexhost(i,j+1,k)];
			f0011 = velocity.v0[getindexhost(i,j+1,k)];
			f0100 = velocity.v0[getindexhost(i+1,j,k)];
			f0101 = velocity.v0[getindexhost(i+1,j,k)];
			f0110 = velocity.v0[getindexhost(i+1,j+1,k)];
			f0111 = velocity.v0[getindexhost(i+1,j+1,k)];
			
			f1000 = velocity.v1[getindexhost(i,j,k)];
			f1001 = velocity.v1[getindexhost(i,j,k)];
			f1010 = velocity.v1[getindexhost(i,j+1,k)];
			f1011 = velocity.v1[getindexhost(i,j+1,k)];
			f1100 = velocity.v1[getindexhost(i+1,j,k)];
			f1101 = velocity.v1[getindexhost(i+1,j,k)];
			f1110 = velocity.v1[getindexhost(i+1,j+1,k)];
			f1111 = velocity.v1[getindexhost(i+1,j+1,k)];
		
		}
		
		/* Linear interpolation for V */
		dXdt[1] = f0000*V16 + f0001*V15 + f0010*V14 + f0011*V13 + f0100*V12 + f0101*V11 + f0110*V10
		+ f0111*V9 + f1000*V8 + f1001*V7 + f1010*V6 + f1011*V5 + f1100*V4 + f1101*V3 + f1110*V2 + f1111*V1; 

		
		/* Vertices of space-time voxel for W */
		if(Dimensions == 3) {
		
			f0000 = velocity.w0[getindexhost(i,j,k)];
			f0001 = velocity.w0[getindexhost(i,j,k+1)];
			f0010 = velocity.w0[getindexhost(i,j+1,k)];
			f0011 = velocity.w0[getindexhost(i,j+1,k+1)];
			f0100 = velocity.w0[getindexhost(i+1,j,k)];
			f0101 = velocity.w0[getindexhost(i+1,j,k+1)];
			f0110 = velocity.w0[getindexhost(i+1,j+1,k)];
			f0111 = velocity.w0[getindexhost(i+1,j+1,k+1)];
			
			f1000 = velocity.w1[getindexhost(i,j,k)];
			f1001 = velocity.w1[getindexhost(i,j,k+1)];
			f1010 = velocity.w1[getindexhost(i,j+1,k)];
			f1011 = velocity.w1[getindexhost(i,j+1,k+1)];
			f1100 = velocity.w1[getindexhost(i+1,j,k)];
			f1101 = velocity.w1[getindexhost(i+1,j,k+1)];
			f1110 = velocity.w1[getindexhost(i+1,j+1,k)];
			f1111 = velocity.w1[getindexhost(i+1,j+1,k+1)];
		
		}
		else {
		
			f0000 = velocity.w0[getindexhost(i,j,k)];
			f0001 = velocity.w0[getindexhost(i,j,k)];
			f0010 = velocity.w0[getindexhost(i,j+1,k)];
			f0011 = velocity.w0[getindexhost(i,j+1,k)];
			f0100 = velocity.w0[getindexhost(i+1,j,k)];
			f0101 = velocity.w0[getindexhost(i+1,j,k)];
			f0110 = velocity.w0[getindexhost(i+1,j+1,k)];
			f0111 = velocity.w0[getindexhost(i+1,j+1,k)];
			
			f1000 = velocity.w1[getindexhost(i,j,k)];
			f1001 = velocity.w1[getindexhost(i,j,k)];
			f1010 = velocity.w1[getindexhost(i,j+1,k)];
			f1011 = velocity.w1[getindexhost(i,j+1,k)];
			f1100 = velocity.w1[getindexhost(i+1,j,k)];
			f1101 = velocity.w1[getindexhost(i+1,j,k)];
			f1110 = velocity.w1[getindexhost(i+1,j+1,k)];
			f1111 = velocity.w1[getindexhost(i+1,j+1,k)];
			
		
		}
		
		/* Linear interpolation for W */
		dXdt[2] = f0000*V16 + f0001*V15 + f0010*V14 + f0011*V13 + f0100*V12 + f0101*V11 + f0110*V10
		+ f0111*V9 + f1000*V8 + f1001*V7 + f1010*V6 + f1011*V5 + f1100*V4 + f1101*V3 + f1110*V2 + f1111*V1; 
	}
}

void GetVelocity_Unstructured(const double tq, LagrangianPoint *pt, double *dXdt) {
	
	double x, y, z, x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3;
	double a11, a12, a13, a21, a22, a23, a31, a32, a33;
	double U1t, U2t, U3t, U4t, V1t, V2t, V3t, V4t, W1t, W2t, W3t, W4t;
	double r, s, t, d;
	double V;
	double tloc;
	
	if(pt->ElementIndex == -1) 
		FatalError("Attempting to interpolate velocity at point with Element_Index = -1");
	
	x = pt->X[0];
	y = pt->X[1];
	z = pt->X[2];
	
	/* Set tloc defining where in between loaded data to interpolate in time */
	tloc = (tq - Data_LoadedTMin) / (Data_LoadedTMax - Data_LoadedTMin);
	if((tloc > 1 + TINY) || (tloc < 0 - TINY)) 
		FatalError("tloc = %f and must be between 0 and 1", tloc);
	
	while(1) {   
		if(Dimensions == 3) {
			/* Physical coordinates of nodes of the test element */
			x0 = MeshNodeArray_double.x[MeshElementArray.Node1[pt->ElementIndex]];
			y0 = MeshNodeArray_double.y[MeshElementArray.Node1[pt->ElementIndex]];
			z0 = MeshNodeArray_double.z[MeshElementArray.Node1[pt->ElementIndex]];
			x1 = MeshNodeArray_double.x[MeshElementArray.Node2[pt->ElementIndex]];
			y1 = MeshNodeArray_double.y[MeshElementArray.Node2[pt->ElementIndex]];
			z1 = MeshNodeArray_double.z[MeshElementArray.Node2[pt->ElementIndex]];
			x2 = MeshNodeArray_double.x[MeshElementArray.Node3[pt->ElementIndex]];
			y2 = MeshNodeArray_double.y[MeshElementArray.Node3[pt->ElementIndex]];
			z2 = MeshNodeArray_double.z[MeshElementArray.Node3[pt->ElementIndex]];
			x3 = MeshNodeArray_double.x[MeshElementArray.Node4[pt->ElementIndex]];
			y3 = MeshNodeArray_double.y[MeshElementArray.Node4[pt->ElementIndex]];
			z3 = MeshNodeArray_double.z[MeshElementArray.Node4[pt->ElementIndex]];
			
			/* Entries for mapping of physical to natural coordinates of test element */
			a11 = (z3 - z0) * (y2 - y3) - (z2 - z3) * (y3 - y0);
			a21 = (z3 - z0) * (y0 - y1) - (z0 - z1) * (y3 - y0);
			a31 = (z1 - z2) * (y0 - y1) - (z0 - z1) * (y1 - y2);
			a12 = (x3 - x0) * (z2 - z3) - (x2 - x3) * (z3 - z0);
			a22 = (x3 - x0) * (z0 - z1) - (x0 - x1) * (z3 - z0);
			a32 = (x1 - x2) * (z0 - z1) - (x0 - x1) * (z1 - z2);
			a13 = (y3 - y0) * (x2 - x3) - (y2 - y3) * (x3 - x0);
			a23 = (y3 - y0) * (x0 - x1) - (y0 - y1) * (x3 - x0);
			a33 = (y1 - y2) * (x0 - x1) - (y0 - y1) * (x1 - x2);			 
			
			/* Determinant of mapping from natural to physical coordinates of test element */
			V = (x1 - x0) * ((y2 - y0) * (z3 - z0) - (z2 - z0) * (y3 - y0)) + 
			(x2 - x0) * ((y0 - y1) * (z3 - z0) - (z0 - z1) * (y3 - y0)) +
			(x3 - x0) * ((y1 - y0) * (z2 - z0) - (z1 - z0) * (y2 - y0));
			
			/* Natural coordinates of point to be interpolated */
			r = (a11 * (x - x0) + a12 * (y - y0) + a13 * (z - z0)) / V;
			s = (a21 * (x - x0) + a22 * (y - y0) + a23 * (z - z0)) / V;
			t = (a31 * (x - x0) + a32 * (y - y0) + a33 * (z - z0)) / V;
			
			d = fmin(r, fmin(s, fmin(t, 1 - r - s - t)));
		}
		else { /* Dimensions == 2 */
			/* Physical coordinates of nodes of the test element */
			
			x0 = MeshNodeArray_double.x[MeshElementArray.Node1[pt->ElementIndex]];
			y0 = MeshNodeArray_double.y[MeshElementArray.Node1[pt->ElementIndex]];
			
			x1 = MeshNodeArray_double.x[MeshElementArray.Node2[pt->ElementIndex]];
			y1 = MeshNodeArray_double.y[MeshElementArray.Node2[pt->ElementIndex]];
			
			x2 = MeshNodeArray_double.x[MeshElementArray.Node3[pt->ElementIndex]];
			y2 = MeshNodeArray_double.y[MeshElementArray.Node3[pt->ElementIndex]];
			
			/* Entries for mapping of physical to natural coordinates of test element */
			a11 = y2 - y0;
			a12 = x0 - x2;
			a21 = y0 - y1;
			a22 = x1 - x0;
			
			/* Determinant of mapping from natural to physical coordinates of test element */
			V = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
			
			/* Natural coordinates of point to be interpolated */
			r = (a11 * (x - x0) + a12 * (y - y0)) / V;
			s = (a21 * (x - x0) + a22 * (y - y0)) / V;
			t = 0;
			
			d = fmin(r, fmin(s, 1 - r - s));
		}
		
		if((d + TINY) >= 0) {   /* Point inside test element */
			/* Interpolate velocity at nodes in time */
			U1t = (1 - tloc) * velocity.u0[MeshElementArray.Node1[pt->ElementIndex]] + 
			tloc * velocity.u1[MeshElementArray.Node1[pt->ElementIndex]];
			V1t = (1 - tloc) * velocity.v0[MeshElementArray.Node1[pt->ElementIndex]]+ 
			tloc * velocity.v1[MeshElementArray.Node1[pt->ElementIndex]];
			U2t = (1 - tloc) * velocity.u0[MeshElementArray.Node2[pt->ElementIndex]] + 
			tloc * velocity.u1[MeshElementArray.Node2[pt->ElementIndex]];
			V2t = (1 - tloc) * velocity.v0[MeshElementArray.Node2[pt->ElementIndex]] + 
			tloc * velocity.v1[MeshElementArray.Node2[pt->ElementIndex]];
			U3t = (1 - tloc) * velocity.u0[MeshElementArray.Node3[pt->ElementIndex]] + 
			tloc * velocity.u1[MeshElementArray.Node3[pt->ElementIndex]];
			V3t = (1 - tloc) * velocity.v0[MeshElementArray.Node3[pt->ElementIndex]] + 
			tloc * velocity.v1[MeshElementArray.Node3[pt->ElementIndex]];
			if(Dimensions == 3) {
				W1t = (1 - tloc) * velocity.w0[MeshElementArray.Node1[pt->ElementIndex]] + 
				tloc * velocity.w1[MeshElementArray.Node1[pt->ElementIndex]];
				W2t = (1 - tloc) * velocity.w0[MeshElementArray.Node2[pt->ElementIndex]] + 
				tloc * velocity.w1[MeshElementArray.Node2[pt->ElementIndex]];       
				W3t = (1 - tloc) * velocity.w0[MeshElementArray.Node3[pt->ElementIndex]] + 
				tloc * velocity.w1[MeshElementArray.Node3[pt->ElementIndex]]; 
				U4t = (1 - tloc) * velocity.u0[MeshElementArray.Node4[pt->ElementIndex]] + 
				tloc * velocity.u1[MeshElementArray.Node4[pt->ElementIndex]];
				V4t = (1 - tloc) * velocity.v0[MeshElementArray.Node4[pt->ElementIndex]] + 
				tloc * velocity.v1[MeshElementArray.Node4[pt->ElementIndex]];
				W4t = (1 - tloc) * velocity.w0[MeshElementArray.Node4[pt->ElementIndex]] + 
				tloc * velocity.w1[MeshElementArray.Node4[pt->ElementIndex]];
				
				/* Get interpolate velocity at point in space using time-interpolated velocity at nodes*/
				dXdt[0] = U1t + (U2t - U1t) * fabs(r) + (U3t - U1t) * fabs(s) + (U4t - U1t) * fabs(t);
				dXdt[1] = V1t + (V2t - V1t) * fabs(r) + (V3t - V1t) * fabs(s) + (V4t - V1t) * fabs(t);
				dXdt[2] = W1t + (W2t - W1t) * fabs(r) + (W3t - W1t) * fabs(s) + (W4t - W1t) * fabs(t);
				
				
			}
			else { /* Dimensions = 2 */
				/* Get interpolate velocity at point in space using time-interpolated velocity at nodes*/
				dXdt[0] = U1t + (U2t - U1t) * fabs(r) + (U3t - U1t) * fabs(s);
				dXdt[1] = V1t + (V2t - V1t) * fabs(r) + (V3t - V1t) * fabs(s);
				dXdt[2] = 0.0;
				
				
			}
			return;
		}
		else { /* Point not inside test element */
			/* Reset test element to one of neighbors */
			if(Dimensions == 3) {
				if(fabs(r - d) < TINY) 
					pt->ElementIndex = MeshElementArray.Neighborindex1[pt->ElementIndex];
				else if(fabs(s - d)  < TINY) 
					pt->ElementIndex = MeshElementArray.Neighborindex2[pt->ElementIndex];
				else if(fabs(t - d) < TINY) 
					pt->ElementIndex = MeshElementArray.Neighborindex3[pt->ElementIndex];
				else if(fabs((1 - r - s - t) - d) < TINY) 
					pt->ElementIndex = MeshElementArray.Neighborindex4[pt->ElementIndex];
				else 
					FatalError("Indeterminate neighbor in function GetVelocity_Unstructured()");
			}
			else {
				if(fabs(r - d) < TINY) 
					pt->ElementIndex = MeshElementArray.Neighborindex1[pt->ElementIndex];
				else if(fabs(s - d) < TINY) 
					pt->ElementIndex = MeshElementArray.Neighborindex2[pt->ElementIndex];
				else if(fabs(d - (1 - r - s)) < TINY) 
					pt->ElementIndex = MeshElementArray.Neighborindex3[pt->ElementIndex];
				else 
					FatalError("Indeterminate neighbor in function GetVelocity_Unstructured()");
			}
			/* If point has left domain, flag it and return with velocity set to zero */
			if(pt->ElementIndex == -1) {
				pt->LeftDomain = 1;
				dXdt[0] = 0.0;
				dXdt[1] = 0.0;
				dXdt[2] = 0.0;
				return;
			}
		}
	}
	
}



