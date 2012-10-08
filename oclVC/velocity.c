# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <CL/cl.h>


# include "structs.h"
# include "globals.h"
# include "settings.h"
# include "velocity.h"
# include "index.h"



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
						//printf("vel u = %f \n",velocity[getindexwosshost(i,j,k)].u0 );
					
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
				//printf("vel u = %f \n",velocity[getindexwosshost(i,j,k)].u0 );
				
				
					
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
		
		//printf("i = %d and j = %d \n",i,j);
		//printf("f000 = %0.9f f001 = %0.9f f010 = %0.9f and f011 = %0.9f \n", f0000, f0001, f0010, f0011);
		
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
	//printf("xloc = %0.9f and yloc = %0.9f \n", xloc , yloc);
	//printf("\n xvel = %0.9f, yvel = %0.9f \n for x = %0.9f, y = %0.9f tloc = %0.9f  \n", dXdt[0],dXdt[1],pt->X[0],pt->X[1], tloc);
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
				
				//printf("U1 = %0.9f \n", U1t);
				
				//printf("vel.u0 = %0.9f and vel.u1 = %0.9f tloc = %0.9f \n",Vel_UnstructVelArray_U[0][Vel_MeshElementArray[pt->ElementIndex].Nodes[0]], Vel_UnstructVelArray_U[1][Vel_MeshElementArray[pt->ElementIndex].Nodes[0]] , tloc);
				
				//printf("eid = %d \n", pt->ElementIndex);
				//printf("\n xvel = %0.9f, yvel = %0.9f and zvel = %0.9f \n for x = %0.9f, y = %0.9f and z = %0.9f  \n", dXdt[0],dXdt[1], dXdt[2],x,y,z);
			}
			else { /* Dimensions = 2 */
				/* Get interpolate velocity at point in space using time-interpolated velocity at nodes*/
				dXdt[0] = U1t + (U2t - U1t) * fabs(r) + (U3t - U1t) * fabs(s);
				dXdt[1] = V1t + (V2t - V1t) * fabs(r) + (V3t - V1t) * fabs(s);
				dXdt[2] = 0.0;
				
				//printf("\n U1t = %0.9f, U2t = %0.9f, and U3t = %0.9f \n", V1t, V2t, V3t);
				//printf("u1 = %f \n", Vel_UnstructVelArray_V[1][Vel_MeshElementArray[pt->ElementIndex].Nodes[0]]);
				//printf("U1t = %0.9f \n", U1t);
			//	printf("Mesh Element = %d \n", Vel_MeshElementArray[pt->ElementIndex].Nodes[0]);
			//	printf("eid = %d \n", pt->ElementIndex);
			//	printf("tloc = %0.9f, vel.u0 = %0.9f,  vel.u1 = %0.9f \n", tloc, Vel_UnstructVelArray_U[0][Vel_MeshElementArray[pt->ElementIndex].Nodes[0]], Vel_UnstructVelArray_U[1][Vel_MeshElementArray[pt->ElementIndex].Nodes[0]] );
				//printf("\n xvel = %0.9f and yvel = %0.9f \n for x = %0.9f and y = %0.9f where r = %0.9f and s = %0.9f  \n", dXdt[0],dXdt[1],x,y,r,s);
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



