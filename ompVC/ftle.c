/*
 *  ftle.c
 *  flowVC
 *
 *  Created by Shawn Shadden.
 *  Copyright 2010 Flow Physics Group. All rights reserved.
 *
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "ftle.h"
#include "globals.h"
#include "integration.h"
#include "io.h"
#include "macros.h"
#include "mesh.h"
#include "mymath.h"
#include "structs.h"
#include "velocity.h"


void ComputeFTLE(int df) {

	double t1, t2;
	LagrangianPoint TempPoint;
	int i,j,k,ss,num_integrated;
	
	if (UPPERLOOP == 1) {
		// Upper loop optimization ...
		
		//int dyn_num = floor(Output_TRes / num_processors);
		//int dyn_num1 = floor((FTLE_CartMesh.XRes)/num_processors);
		
		// main parallel omp loop 
		#pragma omp parallel private(ss,t1,t2, TempPoint, num_integrated,i,j,k, FTLE_MeshPt, FTLE_NewArray) 
		{ 
		
			// Allocate memory for FTLE_MeshPt
			if((FTLE_MeshPt = (FTLEPoint ***)calloc(FTLE_CartMesh.XRes , sizeof(FTLEPoint **))) == NULL)
				FatalError("Malloc failed for FTLE_MeshPt");
			//#pragma omp for schedule(dynamic, dyn_num1)
			for(i = 0; i < FTLE_CartMesh.XRes; i++) {
				if((FTLE_MeshPt[i] = (FTLEPoint **)calloc(FTLE_CartMesh.YRes , sizeof(FTLEPoint *))) == NULL) 
					FatalError("Malloc failed for FTLE_MeshPt[%d]", i);
				for(j = 0; j < FTLE_CartMesh.YRes; j++) 
					if((FTLE_MeshPt[i][j] = (FTLEPoint *)calloc(FTLE_CartMesh.ZRes , sizeof(FTLEPoint))) == NULL) 
						FatalError("Malloc failed for FTLE_MeshPt[%d][%d]", i, j);
			} 
		
			// Allocate memory for FTLE_NewArray
			if((FTLE_NewArray = (double ****) calloc(FTLE_CartMesh.XRes, sizeof(double ***))) == NULL)
				FatalError("Malloc failed for FTLE_NewArray failed ...");
			//#pragma omp for schedule(dynamic, dyn_num1)
			for(i = 0; i < FTLE_CartMesh.XRes; i++) {
				if((FTLE_NewArray[i] = (double ***)malloc(FTLE_CartMesh.YRes * sizeof(double **))) == NULL) 
					FatalError("Malloc failed for FTLE_NewArray[%d]", i);
				for(j = 0; j < FTLE_CartMesh.YRes; j++) {
					if((FTLE_NewArray[i][j] = (double **)malloc(FTLE_CartMesh.ZRes * sizeof(double *))) == NULL) 
						FatalError("Malloc failed for FTLE_NewArray[%d][%d]", i, j);
					for(k = 0; k < FTLE_CartMesh.ZRes; k++) 
						if((FTLE_NewArray[i][j][k] = (double *)malloc(3 * sizeof(double))) == NULL) 
							FatalError("Malloc failed for FTLE_NewArray[%d][%d][%d]", i, j, k);
				}
			}
			
			#pragma omp for schedule(dynamic) 
			for(ss = 0; ss < Output_TRes; ss++) {
			
				/* Consider only output times that are not complete AND, have been launched OR need to be launched */ 
				if(FTLE_Launches[ss].Status != COMPLETE && (FTLE_Launches[ss].Status == LAUNCHED || 
							(FTLE_Launches[ss].StartTime <= Data_LoadedTMax + TINY && FTLE_Launches[ss].StartTime >= Data_LoadedTMin - TINY))) {
				
					ReadInFTLELaunch(ss, FTLE_MeshPt); // Load data for release ss from bin file to FTLE_MeshPt
					// Read in from the file .....
					//printf("FTLE_Launches[ss].StartTime = %f \n", FTLE_MeshPt[14][14][0].Pt.X[1]);
							   
					/* Initialize slide if it hasn't been launched yet and set integration bounds */
					if(FTLE_Launches[ss].Status == UNLAUNCHED) {
						FTLE_Launches[ss].Status = LAUNCHED;
						t1 = FTLE_Launches[ss].StartTime;
						/* Loop over all the points */
						for(i = 0; i < FTLE_CartMesh.XRes; i++) 
							for(j = 0; j < FTLE_CartMesh.YRes; j++) 
								for(k = 0; k < FTLE_CartMesh.ZRes; k++) 
									if(!FTLE_MeshPt[i][j][k].Pt.LeftDomain) {
										if(Particle_Radius > TINY && Particle_ICType == 0) {
											FTLE_MeshPt[i][j][k].Pt.V[0] = 0.0;
											FTLE_MeshPt[i][j][k].Pt.V[1] = 0.0;
											FTLE_MeshPt[i][j][k].Pt.V[2] = 0.0;
										}
										else {
											GetVelocity(FTLE_Launches[ss].StartTime, &(FTLE_MeshPt[i][j][k].Pt), FTLE_MeshPt[i][j][k].Pt.V);
										}
									}
						printf("Initiated FTLE output %d at time %g\n", ss, t1); 
						fflush(stdout);
					}
					else 
						t1 = (Int_TimeDirection > 0) ? Data_LoadedTMin : Data_LoadedTMax;
					t2 = (Int_TimeDirection > 0) ? fmin(Data_LoadedTMax, FTLE_Launches[ss].StopTime) : fmax(Data_LoadedTMin, FTLE_Launches[ss].StopTime);
							   
					/* Integrate release ss from time t1 to time t2 */
					//printf("Integrating FTLE grid initiated at %g from %g to %g...\n", FTLE_Launches[ss].StartTime, t1, t2); 
					//fflush(stdout);
					num_integrated = 0; 
				
					for(i = 0; i < FTLE_CartMesh.XRes; i++) {
						for(j = 0; j < FTLE_CartMesh.YRes; j++) {
							for(k = 0; k < FTLE_CartMesh.ZRes; k++) { 
								/* Make sure point not initially outside of domain AND: it hasn't since left, unless extrapolating position */ 
								if(FTLE_MeshPt[i][j][k].Pt.LeftDomainTime > -TINY && (!FTLE_MeshPt[i][j][k].Pt.LeftDomain || Int_Extrapolate)) {
									/* Create TempPoint to integrate in case point leaves the domain */
							
						
									TempPoint.X[0] = FTLE_MeshPt[i][j][k].Pt.X[0];
									TempPoint.X[1] = FTLE_MeshPt[i][j][k].Pt.X[1];
									TempPoint.X[2] = FTLE_MeshPt[i][j][k].Pt.X[2];
									TempPoint.LeftDomain = FTLE_MeshPt[i][j][k].Pt.LeftDomain;
									TempPoint.LeftDomainTime = FTLE_MeshPt[i][j][k].Pt.LeftDomainTime;
									TempPoint.V[0] = FTLE_MeshPt[i][j][k].Pt.V[0];
									TempPoint.V[1] = FTLE_MeshPt[i][j][k].Pt.V[1];
									TempPoint.V[2] = FTLE_MeshPt[i][j][k].Pt.V[2];
									if(Data_MeshType == UNSTRUCTURED)
										TempPoint.ElementIndex = FTLE_MeshPt[i][j][k].Pt.ElementIndex;
										   
									Advect(&TempPoint, t1, t2);
									num_integrated++;
									if(TempPoint.LeftDomain && !Int_Extrapolate) {
										/* Compute FTLE at point (early) */
										if(!FTLE_MeshPt[i][j][k].HaveFTLE) 
											GetFTLEForPointEarly(ss, i, j, k, t1, TempPoint.LeftDomainTime, FTLE_MeshPt);
										/* Compute FTLE for neighbors (early) */
										if(i > 0)  /* Get FTLE for i-1 neighbor */
											if(!FTLE_MeshPt[i-1][j][k].HaveFTLE) 
												GetFTLEForPointEarly(ss, i-1, j, k, t1, TempPoint.LeftDomainTime, FTLE_MeshPt);
										if(i < FTLE_CartMesh.XRes - 1)  /* Get FTLE for i+1 neighbor */
											if(!FTLE_MeshPt[i+1][j][k].HaveFTLE) 
												GetFTLEForPointEarly(ss, i+1, j, k, t1, TempPoint.LeftDomainTime, FTLE_MeshPt);
										if(j > 0)  /* Get FTLE for j-1 neighbor */
											if(!FTLE_MeshPt[i][j-1][k].HaveFTLE)
												GetFTLEForPointEarly(ss, i, j-1, k, t1, TempPoint.LeftDomainTime, FTLE_MeshPt);
										if(j < FTLE_CartMesh.YRes - 1) /* Get FTLE for j+1 neighbor */
											if(!FTLE_MeshPt[i][j+1][k].HaveFTLE) 
												GetFTLEForPointEarly(ss, i, j+1, k, t1, TempPoint.LeftDomainTime, FTLE_MeshPt);
										if(k > 0)  /* Get FTLE for k-1 neighbor */
											if(!FTLE_MeshPt[i][j][k-1].HaveFTLE) 
												GetFTLEForPointEarly(ss, i, j, k-1, t1, TempPoint.LeftDomainTime, FTLE_MeshPt);
										if(k < FTLE_CartMesh.ZRes - 1)  /* Get FTLE for k+1 neighbor */
											if(!FTLE_MeshPt[i][j][k+1].HaveFTLE) 
												GetFTLEForPointEarly(ss, i, j, k+1, t1, TempPoint.LeftDomainTime, FTLE_MeshPt);
												   
										/* Set LeftDomain flag last, as GetFTLEForPointEarly calls require integration of FTLE_MeshPt[i][j][k].Pt, which will fail if this flag is true */
										FTLE_MeshPt[i][j][k].Pt.LeftDomain = 1; 
									}
									else {
										/* Point still in domain, save position for update of FTLE_MeshPt once all points integrated over current time interval */
										FTLE_NewArray[i][j][k][0] = TempPoint.X[0];
										FTLE_NewArray[i][j][k][1] = TempPoint.X[1];
										FTLE_NewArray[i][j][k][2] = TempPoint.X[2];
										FTLE_MeshPt[i][j][k].Pt.V[0] = TempPoint.V[0];
										FTLE_MeshPt[i][j][k].Pt.V[1] = TempPoint.V[1];
										FTLE_MeshPt[i][j][k].Pt.V[2] = TempPoint.V[2];
										if(Int_Extrapolate)
											FTLE_MeshPt[i][j][k].Pt.LeftDomain = TempPoint.LeftDomain;
										if(Data_MeshType == UNSTRUCTURED) 
											FTLE_MeshPt[i][j][k].Pt.ElementIndex = TempPoint.ElementIndex;
									}
								}
							}
						}
					}
					
					UpdateFTLELocations(FTLE_MeshPt, FTLE_NewArray); /* Update position of points still in domain, i.e. copy FTLE_NewArray values to FTLE_MeshPt */  
					WriteOutFTLELaunch(ss, FTLE_MeshPt); /* Store FTLE_MeshPt in a bin file for slide ss */
			
					printf(" %d points from the FTLE grid integrated initiated at %g from %g to %g...\n",num_integrated, FTLE_Launches[ss].StartTime, t1, t2); 
					fflush(stdout);
							   
					/* Write FTLE to output file if slide ss is complete  */
					if(((Int_TimeDirection > 0) ? (FTLE_Launches[ss].StopTime - Data_LoadedTMax) : (Data_LoadedTMin - FTLE_Launches[ss].StopTime)) < TINY) {
						printf("Stopping FTLE output slide %d and writing to output file.\n", ss); 
						fflush(stdout);
						OutputFTLE(ss, FTLE_IntTLength, FTLE_MeshPt);
						FTLE_Launches[ss].Status = COMPLETE;
					}
					/* Or write FTLE to output file if computing variation of FTLE with integration time */
					else if(FTLE_ComputeVariation && !fmod(df - Data_FirstFrame, FTLE_VariationOutFreq)) {
						printf("Exporting FTLE values at integration length %f.\n", fabs(t2 - FTLE_Launches[ss].StartTime)); fflush(stdout);
						OutputFTLE(ss, fabs(t2 - FTLE_Launches[ss].StartTime), FTLE_MeshPt); 
					}
				
					   
				} // End of IF
			}   // End of omp for ... 
			
			// Free memory for FTLE_NewArray and FTLE_MeshPt
			for(i = 0; i < FTLE_CartMesh.XRes; i++) {
				for(j = 0; j < FTLE_CartMesh.YRes; j++) {
					for(k = 0; k < FTLE_CartMesh.ZRes; k++) {
						free(FTLE_NewArray[i][j][k]);
						FTLE_NewArray[i][j][k] = NULL;
					}
					free(FTLE_NewArray[i][j]);
					FTLE_NewArray[i][j] = NULL;
					free(FTLE_MeshPt[i][j]);
					FTLE_MeshPt[i][j] = NULL;
				}
				free(FTLE_NewArray[i]);
				FTLE_NewArray[i] = NULL;
				free(FTLE_MeshPt[i]);
				FTLE_MeshPt[i] = NULL;
			}
			free(FTLE_NewArray);
			FTLE_NewArray = NULL;
			free(FTLE_MeshPt);
			FTLE_MeshPt = NULL;
			
		} // End of OMP parallel section ...
		
		
		
	} else {
		// Inner loop optimization ... 
		
		for(ss = 0; ss < Output_TRes; ss++) {
			
				/* Consider only output times that are not complete AND, have been launched OR need to be launched */ 
				if(FTLE_Launches[ss].Status != COMPLETE && (FTLE_Launches[ss].Status == LAUNCHED || 
							(FTLE_Launches[ss].StartTime <= Data_LoadedTMax + TINY && FTLE_Launches[ss].StartTime >= Data_LoadedTMin - TINY))) {
				
					ReadInFTLELaunch(ss, FTLE_MeshPt); // Load data for release ss from bin file to FTLE_MeshPt
					// Read in from the file .....
					//printf("FTLE_Launches[ss].StartTime = %f \n", FTLE_MeshPt[14][14][0].Pt.X[1]);
							   
					/* Initialize slide if it hasn't been launched yet and set integration bounds */
					if(FTLE_Launches[ss].Status == UNLAUNCHED) {
						FTLE_Launches[ss].Status = LAUNCHED;
						t1 = FTLE_Launches[ss].StartTime;
						/* Loop over all the points */
						for(i = 0; i < FTLE_CartMesh.XRes; i++) 
							for(j = 0; j < FTLE_CartMesh.YRes; j++) 
								for(k = 0; k < FTLE_CartMesh.ZRes; k++) 
									if(!FTLE_MeshPt[i][j][k].Pt.LeftDomain) {
										if(Particle_Radius > TINY && Particle_ICType == 0) {
											FTLE_MeshPt[i][j][k].Pt.V[0] = 0.0;
											FTLE_MeshPt[i][j][k].Pt.V[1] = 0.0;
											FTLE_MeshPt[i][j][k].Pt.V[2] = 0.0;
										}
										else {
											GetVelocity(FTLE_Launches[ss].StartTime, &(FTLE_MeshPt[i][j][k].Pt), FTLE_MeshPt[i][j][k].Pt.V);
										}
									}
						printf("Initiated FTLE output %d at time %g\n", ss, t1); 
						fflush(stdout);
					}
					else 
						t1 = (Int_TimeDirection > 0) ? Data_LoadedTMin : Data_LoadedTMax;
					t2 = (Int_TimeDirection > 0) ? fmin(Data_LoadedTMax, FTLE_Launches[ss].StopTime) : fmax(Data_LoadedTMin, FTLE_Launches[ss].StopTime);
							   
					/* Integrate release ss from time t1 to time t2 */
					//printf("Integrating FTLE grid initiated at %g from %g to %g...\n", FTLE_Launches[ss].StartTime, t1, t2); 
					//fflush(stdout);
					num_integrated = 0; 
					int num_thread;
					int dyn_num = floor(FTLE_CartMesh.XRes/num_processors);
					#pragma omp parallel private(num_thread, TempPoint)
					{
					num_thread = 0;
					
					#pragma omp for private(i,j,k) schedule(dynamic,dyn_num) collapse(3) nowait
					for(i = 0; i < FTLE_CartMesh.XRes; i++) {
						for(j = 0; j < FTLE_CartMesh.YRes; j++) {
							for(k = 0; k < FTLE_CartMesh.ZRes; k++) { 
								/* Make sure point not initially outside of domain AND: it hasn't since left, unless extrapolating position */ 
								if(FTLE_MeshPt[i][j][k].Pt.LeftDomainTime > -TINY && (!FTLE_MeshPt[i][j][k].Pt.LeftDomain || Int_Extrapolate)) {
									/* Create TempPoint to integrate in case point leaves the domain */
							
						
									TempPoint.X[0] = FTLE_MeshPt[i][j][k].Pt.X[0];
									TempPoint.X[1] = FTLE_MeshPt[i][j][k].Pt.X[1];
									TempPoint.X[2] = FTLE_MeshPt[i][j][k].Pt.X[2];
									TempPoint.LeftDomain = FTLE_MeshPt[i][j][k].Pt.LeftDomain;
									TempPoint.LeftDomainTime = FTLE_MeshPt[i][j][k].Pt.LeftDomainTime;
									TempPoint.V[0] = FTLE_MeshPt[i][j][k].Pt.V[0];
									TempPoint.V[1] = FTLE_MeshPt[i][j][k].Pt.V[1];
									TempPoint.V[2] = FTLE_MeshPt[i][j][k].Pt.V[2];
									if(Data_MeshType == UNSTRUCTURED)
										TempPoint.ElementIndex = FTLE_MeshPt[i][j][k].Pt.ElementIndex;
										   
									Advect(&TempPoint, t1, t2);
									//num_integrated++;
									num_thread++;
									if(TempPoint.LeftDomain && !Int_Extrapolate) {
										/* Compute FTLE at point (early) */
										if(!FTLE_MeshPt[i][j][k].HaveFTLE) 
											GetFTLEForPointEarly(ss, i, j, k, t1, TempPoint.LeftDomainTime, FTLE_MeshPt);
										/* Compute FTLE for neighbors (early) */
										if(i > 0)  /* Get FTLE for i-1 neighbor */
											if(!FTLE_MeshPt[i-1][j][k].HaveFTLE) 
												GetFTLEForPointEarly(ss, i-1, j, k, t1, TempPoint.LeftDomainTime, FTLE_MeshPt);
										if(i < FTLE_CartMesh.XRes - 1)  /* Get FTLE for i+1 neighbor */
											if(!FTLE_MeshPt[i+1][j][k].HaveFTLE) 
												GetFTLEForPointEarly(ss, i+1, j, k, t1, TempPoint.LeftDomainTime, FTLE_MeshPt);
										if(j > 0)  /* Get FTLE for j-1 neighbor */
											if(!FTLE_MeshPt[i][j-1][k].HaveFTLE)
												GetFTLEForPointEarly(ss, i, j-1, k, t1, TempPoint.LeftDomainTime, FTLE_MeshPt);
										if(j < FTLE_CartMesh.YRes - 1) /* Get FTLE for j+1 neighbor */
											if(!FTLE_MeshPt[i][j+1][k].HaveFTLE) 
												GetFTLEForPointEarly(ss, i, j+1, k, t1, TempPoint.LeftDomainTime, FTLE_MeshPt);
										if(k > 0)  /* Get FTLE for k-1 neighbor */
											if(!FTLE_MeshPt[i][j][k-1].HaveFTLE) 
												GetFTLEForPointEarly(ss, i, j, k-1, t1, TempPoint.LeftDomainTime, FTLE_MeshPt);
										if(k < FTLE_CartMesh.ZRes - 1)  /* Get FTLE for k+1 neighbor */
											if(!FTLE_MeshPt[i][j][k+1].HaveFTLE) 
												GetFTLEForPointEarly(ss, i, j, k+1, t1, TempPoint.LeftDomainTime, FTLE_MeshPt);
												   
										/* Set LeftDomain flag last, as GetFTLEForPointEarly calls require integration of FTLE_MeshPt[i][j][k].Pt, which will fail if this flag is true */
										FTLE_MeshPt[i][j][k].Pt.LeftDomain = 1; 
									}
									else {
										/* Point still in domain, save position for update of FTLE_MeshPt once all points integrated over current time interval */
										FTLE_NewArray[i][j][k][0] = TempPoint.X[0];
										FTLE_NewArray[i][j][k][1] = TempPoint.X[1];
										FTLE_NewArray[i][j][k][2] = TempPoint.X[2];
										FTLE_MeshPt[i][j][k].Pt.V[0] = TempPoint.V[0];
										FTLE_MeshPt[i][j][k].Pt.V[1] = TempPoint.V[1];
										FTLE_MeshPt[i][j][k].Pt.V[2] = TempPoint.V[2];
										if(Int_Extrapolate)
											FTLE_MeshPt[i][j][k].Pt.LeftDomain = TempPoint.LeftDomain;
										if(Data_MeshType == UNSTRUCTURED) 
											FTLE_MeshPt[i][j][k].Pt.ElementIndex = TempPoint.ElementIndex;
									}
								}
							}
						}
					} // End of omp for
					
					#pragma omp critical
					{
						num_integrated = num_integrated + num_thread;
					}// End of critical ... 
					
					}   // End of omp parallel ...
					
					UpdateFTLELocations(FTLE_MeshPt, FTLE_NewArray); /* Update position of points still in domain, i.e. copy FTLE_NewArray values to FTLE_MeshPt */  
					WriteOutFTLELaunch(ss, FTLE_MeshPt); /* Store FTLE_MeshPt in a bin file for slide ss */
			
					printf(" %d points from the FTLE grid integrated initiated at %g from %g to %g...\n",num_integrated, FTLE_Launches[ss].StartTime, t1, t2); 
					fflush(stdout);
							   
					/* Write FTLE to output file if slide ss is complete  */
					if(((Int_TimeDirection > 0) ? (FTLE_Launches[ss].StopTime - Data_LoadedTMax) : (Data_LoadedTMin - FTLE_Launches[ss].StopTime)) < TINY) {
						printf("Stopping FTLE output slide %d and writing to output file.\n", ss); 
						fflush(stdout);
						OutputFTLE(ss, FTLE_IntTLength, FTLE_MeshPt);
						FTLE_Launches[ss].Status = COMPLETE;
					}
					/* Or write FTLE to output file if computing variation of FTLE with integration time */
					else if(FTLE_ComputeVariation && !fmod(df - Data_FirstFrame, FTLE_VariationOutFreq)) {
						printf("Exporting FTLE values at integration length %f.\n", fabs(t2 - FTLE_Launches[ss].StartTime)); fflush(stdout);
						OutputFTLE(ss, fabs(t2 - FTLE_Launches[ss].StartTime), FTLE_MeshPt); 
					}
				
					   
				} // End of IF
			}
		
		
		
		
		
		
	}

}





void InitializeFTLEArray(void) {
	
	int ss, i, j, k, seed = -1, index = -1, found = 0, guess = -1, iguess = -1, jguess = -1, count = 0, foundGS = 0, percentage = 0;
	double Xseed[3];
	double FTLE_outside = -1.0; /* Used to mask FTLE values outside of domain */
	FILE *FTLE_BinFileID;
	char FTLE_BinFilePath[LONGSTRING];
	
	FTLEPoint ***FTLE_MeshPt_temp; // temp pointer for storing output slide information ...
	
	if(FTLE_GenerateMesh) {
		
		printf("\nGenerating FTLE mesh...\n");
		fflush(stdout);
		
		/* Allocate memory */
		if((FTLE_MeshPt_temp = (FTLEPoint ***)malloc(FTLE_CartMesh.XRes * sizeof(FTLEPoint **))) == NULL)
			FatalError("Malloc failed for FTLE_temp");
		for(i = 0; i < FTLE_CartMesh.XRes; i++) {
			if((FTLE_MeshPt_temp[i] = (FTLEPoint **)malloc(FTLE_CartMesh.YRes * sizeof(FTLEPoint *))) == NULL) 
				FatalError("Malloc failed for FTLE_temp[%d]", i);
			for(j = 0; j < FTLE_CartMesh.YRes; j++) 
				if((FTLE_MeshPt_temp[i][j] = (FTLEPoint *)malloc(FTLE_CartMesh.ZRes * sizeof(FTLEPoint))) == NULL)
					FatalError("Malloc failed for FTLE_temp[%d][%d]", i, j);
		}
		
		/* Initialize FTLE_MeshPt */
		for(i = 0; i < FTLE_CartMesh.XRes; i++) {
			for(j = 0; j < FTLE_CartMesh.YRes; j++) {
				for(k = 0; k < FTLE_CartMesh.ZRes; k++) {
					FTLE_MeshPt_temp[i][j][k].Pt.X[0] = FTLE_CartMesh.XMin + i * FTLE_CartMesh.XDelta;
					FTLE_MeshPt_temp[i][j][k].Pt.X[1] = FTLE_CartMesh.YMin + j * FTLE_CartMesh.YDelta;
					FTLE_MeshPt_temp[i][j][k].Pt.X[2] = FTLE_CartMesh.ZMin + k * FTLE_CartMesh.ZDelta;
					FTLE_MeshPt_temp[i][j][k].FTLEwT = 0.0;
					FTLE_MeshPt_temp[i][j][k].FTLEwoT = 0.0;
					FTLE_MeshPt_temp[i][j][k].HaveFTLE = 0;
					FTLE_MeshPt_temp[i][j][k].Pt.LeftDomain = 0;
					FTLE_MeshPt_temp[i][j][k].Pt.LeftDomainTime = 0.0;
				}
			}
		}    
		
		/* Check if any points outside bounding box of velocity data */
		for(i = 0; i < FTLE_CartMesh.XRes; i++) {
			for(j = 0; j < FTLE_CartMesh.YRes; j++) {
				for(k = 0; k < FTLE_CartMesh.ZRes; k++) {
					if(TestOutsideDomain(FTLE_MeshPt_temp[i][j][k].Pt.X)) {
						FTLE_MeshPt_temp[i][j][k].Pt.LeftDomain = 1;
						FTLE_MeshPt_temp[i][j][k].Pt.LeftDomainTime = -1.0;
						FTLE_MeshPt_temp[i][j][k].FTLEwT = FTLE_outside;
						FTLE_MeshPt_temp[i][j][k].FTLEwoT = FTLE_outside;
						/* No need to calculate FTLE for this point */
						FTLE_MeshPt_temp[i][j][k].HaveFTLE = 1;	    
						/* Or neighboring points */
						if(i > 0) 
							FTLE_MeshPt_temp[i-1][j][k].HaveFTLE = 1;
						if(i < FTLE_CartMesh.XRes - 1)
							FTLE_MeshPt_temp[i+1][j][k].HaveFTLE = 1;
						if(j > 0) 
							FTLE_MeshPt_temp[i][j-1][k].HaveFTLE = 1;
						if(j < FTLE_CartMesh.YRes - 1)
							FTLE_MeshPt_temp[i][j+1][k].HaveFTLE = 1;
						if(k > 0) 
							FTLE_MeshPt_temp[i][j][k-1].HaveFTLE = 1;
						if(k < FTLE_CartMesh.ZRes - 1)
							FTLE_MeshPt_temp[i][j][k+1].HaveFTLE = 1;
					}
				}
			}
		}
		
		/* If velocity field specified on an unstructured mesh, determine which element each FTLE_MeshPt point is located in  */
		if(Data_MeshType == UNSTRUCTURED) {
			/* First, try to "seed" a local element search by doing a global search for the center FTLE_MeshPt grid point */
			Xseed[0] = FTLE_CartMesh.XMin + ((int)(FTLE_CartMesh.XRes/2)) * FTLE_CartMesh.XDelta;
			Xseed[1] = FTLE_CartMesh.YMin + ((int)(FTLE_CartMesh.YRes/2)) * FTLE_CartMesh.YDelta;
			Xseed[2] = FTLE_CartMesh.ZMin + ((int)(FTLE_CartMesh.ZRes/2)) * FTLE_CartMesh.ZDelta;
			seed = Get_Element_Global_Search(Xseed);
			if(seed >= 0) {
				printf("  Using center element seed.\n");
				fflush(stdout);
				guess = seed;
			}
			
			/* Search over all the points */
			for(i = 0; i < FTLE_CartMesh.XRes; i++) {
				iguess = -1;
				for(j = 0; j < FTLE_CartMesh.YRes; j++) {
					jguess = -1;
					for(k = 0; k < FTLE_CartMesh.ZRes; k++) {
						if(!FTLE_MeshPt_temp[i][j][k].Pt.LeftDomain) { 
							if(seed < 0) {
								seed = Get_Element_Global_Search(FTLE_MeshPt_temp[i][j][k].Pt.X);
								if(seed < 0) {
									FTLE_MeshPt_temp[i][j][k].Pt.LeftDomain = 1;
									FTLE_MeshPt_temp[i][j][k].Pt.LeftDomainTime = -1.0;
									/* No need to calculate FTLE for this point */
									FTLE_MeshPt_temp[i][j][k].HaveFTLE = 1;
									FTLE_MeshPt_temp[i][j][k].FTLEwT = FTLE_outside;
									FTLE_MeshPt_temp[i][j][k].FTLEwoT = FTLE_outside;
									/* Or its neighbors */
									if(i > 0) 
										FTLE_MeshPt_temp[i-1][j][k].HaveFTLE = 1;
									if(i < FTLE_CartMesh.XRes - 1)
										FTLE_MeshPt_temp[i+1][j][k].HaveFTLE = 1;
									if(j > 0) 
										FTLE_MeshPt_temp[i][j-1][k].HaveFTLE = 1;
									if(j < FTLE_CartMesh.YRes - 1)
										FTLE_MeshPt_temp[i][j+1][k].HaveFTLE = 1;
									if(k > 0) 
										FTLE_MeshPt_temp[i][j][k-1].HaveFTLE = 1;
									if(k < FTLE_CartMesh.ZRes - 1)
										FTLE_MeshPt_temp[i][j][k+1].HaveFTLE = 1;
								}
								else {
									printf("  Using first located element seed.\n");
									fflush(stdout);
									FTLE_MeshPt_temp[i][j][k].Pt.ElementIndex = seed;
									found++;
									guess  = seed;
									jguess = seed;
									iguess = seed;
								}
							}
							else {
								count++;
								if(100 * count / (FTLE_CartMesh.XRes * FTLE_CartMesh.YRes * FTLE_CartMesh.ZRes) > percentage) {
									printf("  %d%%", percentage);
									fflush(stdout);
									percentage = percentage + 10;
								}  
								index = Get_Element_Local_Search(FTLE_MeshPt_temp[i][j][k].Pt.X, guess);
								if(index < 0) { 
									/* An element not found, try searching from seed location */ 
									index = Get_Element_Local_Search(FTLE_MeshPt_temp[i][j][k].Pt.X, seed);
									if(index < 0) {
										/* If an element still not found, do global search if local search checking requested */
										if(LocalSearchChecking) {
											index = Get_Element_Global_Search(FTLE_MeshPt_temp[i][j][k].Pt.X);
											if(index < 0) { /* Point is definitely not in any element */
												FTLE_MeshPt_temp[i][j][k].Pt.LeftDomain = 1;
												FTLE_MeshPt_temp[i][j][k].Pt.LeftDomainTime = -1.0;
												/* No need to calculate FTLE for this point or its neighbors */
												FTLE_MeshPt_temp[i][j][k].HaveFTLE = 1;
												FTLE_MeshPt_temp[i][j][k].FTLEwT = FTLE_outside;
												FTLE_MeshPt_temp[i][j][k].FTLEwoT = FTLE_outside;
												if(i > 0) 
													FTLE_MeshPt_temp[i-1][j][k].HaveFTLE = 1;
												if(i < FTLE_CartMesh.XRes - 1)
													FTLE_MeshPt_temp[i+1][j][k].HaveFTLE = 1;
												if(j > 0) 
													FTLE_MeshPt_temp[i][j-1][k].HaveFTLE = 1;
												if(j < FTLE_CartMesh.YRes - 1)
													FTLE_MeshPt_temp[i][j+1][k].HaveFTLE = 1;
												if(k > 0) 
													FTLE_MeshPt_temp[i][j][k-1].HaveFTLE = 1;
												if(k < FTLE_CartMesh.ZRes - 1)
													FTLE_MeshPt_temp[i][j][k+1].HaveFTLE = 1;
											}
											else {
												FTLE_MeshPt_temp[i][j][k].Pt.ElementIndex = index;	
												foundGS++;
												found++;
												guess = index;
												if(jguess < 0)
													jguess = index;
												if(iguess < 0)
													iguess = index;
											}
										}
										else {
											FTLE_MeshPt_temp[i][j][k].Pt.LeftDomain = 1;
											FTLE_MeshPt_temp[i][j][k].Pt.LeftDomainTime = -1.0;
											/* No need to calculate FTLE for this point or its neighbors */
											FTLE_MeshPt_temp[i][j][k].HaveFTLE = 1;
											FTLE_MeshPt_temp[i][j][k].FTLEwT = FTLE_outside;
											FTLE_MeshPt_temp[i][j][k].FTLEwoT = FTLE_outside;
											if(i > 0) 
												FTLE_MeshPt_temp[i-1][j][k].HaveFTLE = 1;
											if(i < FTLE_CartMesh.XRes - 1)
												FTLE_MeshPt_temp[i+1][j][k].HaveFTLE = 1;
											if(j > 0) 
												FTLE_MeshPt_temp[i][j-1][k].HaveFTLE = 1;
											if(j < FTLE_CartMesh.YRes - 1)
												FTLE_MeshPt_temp[i][j+1][k].HaveFTLE = 1;
											if(k > 0) 
												FTLE_MeshPt_temp[i][j][k-1].HaveFTLE = 1;
											if(k < FTLE_CartMesh.ZRes - 1)
												FTLE_MeshPt_temp[i][j][k+1].HaveFTLE = 1;
										}
									}
									else {
										FTLE_MeshPt_temp[i][j][k].Pt.ElementIndex = index;	
										found++;
										guess = index;
										if(jguess < 0)
											jguess = index;
										if(iguess < 0)
											iguess = index;
									}
								}
								else {
									FTLE_MeshPt_temp[i][j][k].Pt.ElementIndex = index;
									found++;
									guess = index;
									if(jguess < 0)
										jguess = index;
									if(iguess < 0)
										iguess = index;
								}
							}  
						}
					}
				}
				if(jguess >= 0)
					guess = jguess;
			}
			if(iguess >= 0)
				guess = iguess;
		}
		
		if(LocalSearchChecking) 
			printf("  %d of %d located in domain (%d caught by local search checking)\n", found, FTLE_CartMesh.XRes * FTLE_CartMesh.YRes * FTLE_CartMesh.ZRes, foundGS);
		else
			printf("  %d of %d located in domain\n", found, FTLE_CartMesh.XRes * FTLE_CartMesh.YRes * FTLE_CartMesh.ZRes);
		fflush(stdout);
		
		
		/* Open bin file to store FTLE_MeshPt_temp initialization data */
		sprintf(FTLE_BinFilePath, "%s%s", Path_Data, FTLE_ICFile);
		if((FTLE_BinFileID = fopen(FTLE_BinFilePath, "wb")) == NULL) 
			FatalError("Could not open file %s", FTLE_BinFilePath);
		
		/* Write FTLE grid dimensions and FTLE_MeshPt_temp information to binary file */
		if(fwrite(&FTLE_CartMesh.XRes, sizeof(int), 1, FTLE_BinFileID) < 1 || 
		   fwrite(&FTLE_CartMesh.YRes, sizeof(int), 1, FTLE_BinFileID) < 1 || 
		   fwrite(&FTLE_CartMesh.ZRes, sizeof(int), 1, FTLE_BinFileID) < 1) 
			FatalError("Could not completely write FTLE grid resolution data to %s", FTLE_BinFilePath);
		for(i = 0; i < FTLE_CartMesh.XRes; i++)
			for(j = 0; j < FTLE_CartMesh.YRes; j++) 
				if((int)fwrite(FTLE_MeshPt_temp[i][j], sizeof(FTLEPoint), FTLE_CartMesh.ZRes, FTLE_BinFileID) < FTLE_CartMesh.ZRes) 
					FatalError("Could not completely write FTLE_MeshPt_temp data to %s", FTLE_BinFilePath);
		
		/* Close binary file */
		fclose(FTLE_BinFileID);
	}
	else { /* Don't generate mesh, read it from IC file */
		printf("\nLoading FTLE array data...");
		fflush(stdout);
		
		/* Open file */
		sprintf(FTLE_BinFilePath, "%s%s", Path_Data, FTLE_ICFile);
		if((FTLE_BinFileID = fopen(FTLE_BinFilePath, "rb")) == NULL) {
			fprintf(stderr, "Error: Could not open file %s (file must be generated by flowVC with FTLE_GenerateMesh = 0)\n", FTLE_BinFilePath);
			exit(1);
		}
		
		/* Read grid dimensions */
		if(fread(&FTLE_CartMesh.XRes, sizeof(int), 1, FTLE_BinFileID) < 1 || 
		   fread(&FTLE_CartMesh.YRes, sizeof(int), 1, FTLE_BinFileID) < 1 || 
		   fread(&FTLE_CartMesh.ZRes, sizeof(int), 1, FTLE_BinFileID) < 1)
			FatalError("Could not completely read FTLE grid resolution from %s", FTLE_BinFilePath);
		
		/* Allocate memory */
		if((FTLE_MeshPt_temp = (FTLEPoint ***)malloc(FTLE_CartMesh.XRes * sizeof(FTLEPoint **))) == NULL)
			FatalError("Malloc failed for FTLE_MeshPt_temp");
		for(i = 0; i < FTLE_CartMesh.XRes; i++) {
			if((FTLE_MeshPt_temp[i] = (FTLEPoint **)malloc(FTLE_CartMesh.YRes * sizeof(FTLEPoint *))) == NULL) 
				FatalError("Malloc failed for FTLE_MeshPt_temp[%d]", i);
			for(j = 0; j < FTLE_CartMesh.YRes; j++) 
				if((FTLE_MeshPt_temp[i][j] = (FTLEPoint *)malloc(FTLE_CartMesh.ZRes * sizeof(FTLEPoint))) == NULL) 
					FatalError("Malloc failed for FTLE_MeshPt_temp[%d][%d]", i, j);
		}
		
		/* Read FTLE_MeshPt_temp information from file */
		for(i = 0; i < FTLE_CartMesh.XRes; i++)
			for(j = 0; j < FTLE_CartMesh.YRes; j++)
				if((int)fread(FTLE_MeshPt_temp[i][j], sizeof(FTLEPoint), FTLE_CartMesh.ZRes, FTLE_BinFileID) < FTLE_CartMesh.ZRes) 
					FatalError("Could not completely read FTLE_MeshPt_temp data from %s", FTLE_BinFilePath);
		
		/* Close file */
		fclose(FTLE_BinFileID);
	}
	
	
	if((FTLE_Launches = (Launch *)malloc(Output_TRes * sizeof(Launch))) == NULL) 
		FatalError("Malloc failed for FTLE");
	
	for(ss = 0; ss < Output_TRes; ss++) {
		/* Initiate slide information for release */
		FTLE_Launches[ss].Status = UNLAUNCHED;
		FTLE_Launches[ss].StartTime = Output_TStart + ss * Int_TimeDirection * Output_TDelta;
		if(!Data_TPeriodic) {
			FTLE_Launches[ss].StopTime = fmin(Output_TStart + ss * Int_TimeDirection * Output_TDelta + (double)Int_TimeDirection * FTLE_IntTLength, Data_TMax);
			FTLE_Launches[ss].StopTime = fmax(FTLE_Launches[ss].StopTime, Data_TMin);
			if(fabs(FTLE_Launches[ss].StopTime - FTLE_Launches[ss].StartTime) < FTLE_IntTLength - TINY)
				fprintf(stderr, "Warning: FTLE at time %f will have a maximum integration length of %f\n", FTLE_Launches[ss].StartTime, fabs(FTLE_Launches[ss].StopTime - FTLE_Launches[ss].StartTime));
		}
		else 
			FTLE_Launches[ss].StopTime = Output_TStart + ss * Int_TimeDirection * Output_TDelta + (double)Int_TimeDirection * FTLE_IntTLength;
		
		/* Open binary file to store FTLE_MeshPt_temp for given slide */
		sprintf(FTLE_BinFilePath, "%s%s.%d.bin", Path_Data, FTLE_OutFilePrefix, ss);
		if((FTLE_BinFileID = fopen(FTLE_BinFilePath, "wb")) == NULL) 
			FatalError("Could not open file %s", FTLE_BinFilePath);
		
		/* Write FTLE_MeshPt_temp information to binary file */
		for(i = 0; i < FTLE_CartMesh.XRes; i++)
			for(j = 0; j < FTLE_CartMesh.YRes; j++)
				if((int)fwrite(FTLE_MeshPt_temp[i][j], sizeof(FTLEPoint), FTLE_CartMesh.ZRes, FTLE_BinFileID) < FTLE_CartMesh.ZRes)
					FatalError("Could not completely write FTLE_MeshPt_temp to file %s", FTLE_BinFilePath);
		
		/* Close binary file */
		fclose(FTLE_BinFileID);
	}
	
	// Free memory for FTLE_MeshPt_temp ...
	for(i = 0; i < FTLE_CartMesh.XRes; i++) {
		for(j = 0; j < FTLE_CartMesh.YRes; j++) {
			free(FTLE_MeshPt_temp[i][j]);
			FTLE_MeshPt_temp[i][j] = NULL;
		}
		free(FTLE_MeshPt_temp[i]);
		FTLE_MeshPt_temp[i] = NULL;
	}
	free(FTLE_MeshPt_temp);
	FTLE_MeshPt_temp = NULL;

	if (LOWERLOOP == 1) {
	
		printf("memory allocating .. \n\n\n\n\n");
		/* Allocate memory */
		if((FTLE_MeshPt = (FTLEPoint ***)malloc(FTLE_CartMesh.XRes * sizeof(FTLEPoint **))) == NULL)
			FatalError("Malloc failed for FTLE");
		for(i = 0; i < FTLE_CartMesh.XRes; i++) {
			if((FTLE_MeshPt[i] = (FTLEPoint **)malloc(FTLE_CartMesh.YRes * sizeof(FTLEPoint *))) == NULL) 
				FatalError("Malloc failed for FTLE[%d]", i);
			for(j = 0; j < FTLE_CartMesh.YRes; j++) 
				if((FTLE_MeshPt[i][j] = (FTLEPoint *)malloc(FTLE_CartMesh.ZRes * sizeof(FTLEPoint))) == NULL)
					FatalError("Malloc failed for FTLE[%d][%d]", i, j);
		}
	
		/* Allocate memory for an array to hold updated FTLE grid positions */
		if((FTLE_NewArray = (double ****)malloc(FTLE_CartMesh.XRes * sizeof(double ***))) == NULL) 
			FatalError("Malloc failed for FTLE_NewArray");
		for(i = 0; i < FTLE_CartMesh.XRes; i++) {
			if((FTLE_NewArray[i] = (double ***)malloc(FTLE_CartMesh.YRes * sizeof(double **))) == NULL) 
				FatalError("Malloc failed for FTLE_NewArray[%d]", i);
			for(j = 0; j < FTLE_CartMesh.YRes; j++) {
				if((FTLE_NewArray[i][j] = (double **)malloc(FTLE_CartMesh.ZRes * sizeof(double *))) == NULL) 
					FatalError("Malloc failed for FTLE_NewArray[%d][%d]", i, j);
				for(k = 0; k < FTLE_CartMesh.ZRes; k++) 
					if((FTLE_NewArray[i][j][k] = (double *)malloc(3 * sizeof(double))) == NULL) 
						FatalError("Malloc failed for FTLE_NewArray[%d][%d][%d]", i, j, k);
			}
		}
	
	}
	
	printf("OK!\n");
	fflush(stdout);
}

void ReadInFTLELaunch(int ss, FTLEPoint ***FTLE_MeshPt) {
	/*** Function reads in FTLE data structures for slide specified by ss into FTLE_MeshPt ***/
	
	int i, j;
	FILE *FTLE_BinFileID;
	char FTLE_BinFilePath[LONGSTRING];
	
	/* Open file to read binary */
	sprintf(FTLE_BinFilePath, "%s%s.%d.bin", Path_Data, FTLE_OutFilePrefix, ss);
	if((FTLE_BinFileID = fopen(FTLE_BinFilePath, "rb")) == NULL) 
		FatalError("Could not open file %s", FTLE_BinFilePath);
	
	/* Read FTLE_MeshPt information from file */
	for(i = 0; i < FTLE_CartMesh.XRes; i++) 
		for(j = 0; j < FTLE_CartMesh.YRes; j++) 
			if((int)fread(FTLE_MeshPt[i][j], sizeof(FTLEPoint), FTLE_CartMesh.ZRes, FTLE_BinFileID) < FTLE_CartMesh.ZRes)
				FatalError("Could not completely read FTLE_MeshPt[%d][%d] from %s", i, j, FTLE_BinFilePath);
	
	fclose(FTLE_BinFileID);
	
}

void WriteOutFTLELaunch(int ss, FTLEPoint ***FTLE_MeshPt) {
	/* Writes contents of FTLE_MeshPt to binary file for slide ss */
	
	int i, j;
	FILE *FTLE_BinFileID;
	char FTLE_BinFilePath[LONGSTRING];
	
	sprintf(FTLE_BinFilePath, "%s%s.%d.bin", Path_Data, FTLE_OutFilePrefix, ss);
	if((FTLE_BinFileID = fopen(FTLE_BinFilePath, "wb")) == NULL) 
		FatalError("Could not open file %s", FTLE_BinFilePath);
	
	for(i = 0; i < FTLE_CartMesh.XRes; i++) 
		for(j = 0; j < FTLE_CartMesh.YRes; j++) 
			if((int)fwrite(FTLE_MeshPt[i][j], sizeof(FTLEPoint), FTLE_CartMesh.ZRes, FTLE_BinFileID) < FTLE_CartMesh.ZRes) 
				FatalError("Could not completely write FTLE_MeshPt data to %s", FTLE_BinFilePath);
	
	fclose(FTLE_BinFileID);
	
}

void GetFTLEForPointEarly(int ss, int i, int j, int k, double t1, double te, FTLEPoint ***FTLE_MeshPt) {
	
	double ptLeftX, ptLeftY, ptLeftZ;
	double ptRightX, ptRightY, ptRightZ;
	double DiffXFact;
	
	double ptBottomX, ptBottomY, ptBottomZ;
	double ptTopX, ptTopY, ptTopZ;
	double DiffYFact;
	
	double ptBackX, ptBackY, ptBackZ;
	double ptFrontX, ptFrontY, ptFrontZ;
	double DiffZFact;
	
	double Dphi[3][3];
	double Delta[3][3];
	double lambda;
	
	LagrangianPoint TempPoint, PassedPoint; 
	double LeftDomainTime[6], MinExitTime;
	
	int ii, jj, mm;
	
	/* Advect all neighboring FTLE nodes and see who leaves first */
	/* Left Neighbor */
	if(i > 0) {
		TempPoint = Advect_FTLEPoint(i-1, j, k, t1, te, FTLE_MeshPt);
		if(TempPoint.LeftDomain)
			LeftDomainTime[0] = TempPoint.LeftDomainTime;
		else
			LeftDomainTime[0] = te;
	}
	else
		LeftDomainTime[0] = te;
	
	/* Right Neighbor */
	if(i < FTLE_CartMesh.XRes - 1) {
		TempPoint = Advect_FTLEPoint(i+1, j, k, t1, te, FTLE_MeshPt);
		if(TempPoint.LeftDomain)
			LeftDomainTime[1] = TempPoint.LeftDomainTime;
		else
			LeftDomainTime[1] = te;
	}
	else
		LeftDomainTime[1] = te;
	
	/* Bottom Neighbor */
	if(j > 0) {
		TempPoint = Advect_FTLEPoint(i, j-1, k, t1, te, FTLE_MeshPt);
		if(TempPoint.LeftDomain)
			LeftDomainTime[2] = TempPoint.LeftDomainTime;
		else
			LeftDomainTime[2] = te;
	}
	else
		LeftDomainTime[2] = te;
	
	/* Top Neighbor */
	if(j < FTLE_CartMesh.YRes - 1) {
		TempPoint = Advect_FTLEPoint(i, j+1, k, t1, te, FTLE_MeshPt);
		if(TempPoint.LeftDomain)
			LeftDomainTime[3] = TempPoint.LeftDomainTime;
		else
			LeftDomainTime[3] = te;
	}
	else
		LeftDomainTime[3] = te;
	
	/* Back Neighbor */
	if(k > 0) {
		TempPoint = Advect_FTLEPoint(i, j, k-1, t1, te, FTLE_MeshPt);
		if(TempPoint.LeftDomain)
			LeftDomainTime[4] = TempPoint.LeftDomainTime;
		else
			LeftDomainTime[4] = te;
	}
	else
		LeftDomainTime[4] = te;
	
	/* Front Neighbor */
	if(k < FTLE_CartMesh.ZRes - 1) {
		TempPoint = Advect_FTLEPoint(i, j, k+1, t1, te, FTLE_MeshPt);
		if(TempPoint.LeftDomain)
			LeftDomainTime[5] = TempPoint.LeftDomainTime;
		else
			LeftDomainTime[5] = te;
	}
	else
		LeftDomainTime[5] = te;
	
	MinExitTime = fmin(te, fmin(LeftDomainTime[5], fmin(LeftDomainTime[4], fmin(LeftDomainTime[3], fmin(LeftDomainTime[2], fmin(LeftDomainTime[1], LeftDomainTime[0]))))));
	
	PassedPoint = Advect_FTLEPoint(i, j, k, t1, MinExitTime, FTLE_MeshPt);
	
	/* Differencing in X-Dir */
	if(i == 0) {
		/* Left Point */
		ptLeftX  = PassedPoint.X[0];
		ptLeftY  = PassedPoint.X[1];
		ptLeftZ  = PassedPoint.X[2];
		/* Right Point */
		TempPoint = Advect_FTLEPoint(i+1, j, k, t1, MinExitTime, FTLE_MeshPt);
		ptRightX = TempPoint.X[0];
		ptRightY = TempPoint.X[1];
		ptRightZ = TempPoint.X[2];
		DiffXFact = 1.0;
	}
	else if(i == FTLE_CartMesh.XRes - 1) {
		/* Left Point */
		TempPoint = Advect_FTLEPoint(i-1, j, k, t1, MinExitTime, FTLE_MeshPt);
		ptLeftX = TempPoint.X[0];
		ptLeftY = TempPoint.X[1];
		ptLeftZ = TempPoint.X[2];
		/* Right Point */
		ptRightX = PassedPoint.X[0];
		ptRightY = PassedPoint.X[1];
		ptRightZ = PassedPoint.X[2];
		DiffXFact = 1.0;
	}
	else {
		/* Left Point */
		TempPoint = Advect_FTLEPoint(i-1, j, k, t1, MinExitTime, FTLE_MeshPt);
		ptLeftX = TempPoint.X[0];
		ptLeftY = TempPoint.X[1];
		ptLeftZ = TempPoint.X[2];
		/* Right Point */
		TempPoint = Advect_FTLEPoint(i+1, j, k, t1, MinExitTime, FTLE_MeshPt);
		ptRightX = TempPoint.X[0];
		ptRightY = TempPoint.X[1];
		ptRightZ = TempPoint.X[2];
		DiffXFact = 2.0;
	}
	
	/* Differencing in Y-Dir */
	if(j == 0) {
		/* Bottom Point */
		ptBottomX  = PassedPoint.X[0];
		ptBottomY  = PassedPoint.X[1];
		ptBottomZ  = PassedPoint.X[2];
		/* Top Point */
		TempPoint = Advect_FTLEPoint(i, j+1, k, t1, MinExitTime, FTLE_MeshPt);
		ptTopX = TempPoint.X[0];
		ptTopY = TempPoint.X[1];
		ptTopZ = TempPoint.X[2];
		DiffYFact = 1.0;
	}
	else if(j == FTLE_CartMesh.YRes - 1) {
		/* Bottom Point */
		TempPoint = Advect_FTLEPoint(i, j-1, k, t1, MinExitTime, FTLE_MeshPt);
		ptBottomX = TempPoint.X[0];
		ptBottomY = TempPoint.X[1];
		ptBottomZ = TempPoint.X[2];
		/* Top Point */
		ptTopX = PassedPoint.X[0];
		ptTopY = PassedPoint.X[1];
		ptTopZ = PassedPoint.X[2];
		DiffYFact = 1.0;
	}
	else {
		/* Bottom Point */
		TempPoint = Advect_FTLEPoint(i, j-1, k, t1, MinExitTime, FTLE_MeshPt);
		ptBottomX = TempPoint.X[0];
		ptBottomY = TempPoint.X[1];
		ptBottomZ = TempPoint.X[2];
		/* Top Point */
		TempPoint = Advect_FTLEPoint(i, j+1, k, t1, MinExitTime, FTLE_MeshPt);
		ptTopX = TempPoint.X[0];
		ptTopY = TempPoint.X[1];
		ptTopZ = TempPoint.X[2];
		DiffYFact = 2.0;
	}
	
	if(Dimensions == 3) {
		/* Differencing in Z-Dir */
		if(k == 0) {
			/* Back Point */
			ptBackX  = PassedPoint.X[0];
			ptBackY  = PassedPoint.X[1];
			ptBackZ  = PassedPoint.X[2];
			/* Front Point */
			TempPoint = Advect_FTLEPoint(i, j, k+1, t1, MinExitTime, FTLE_MeshPt);
			ptFrontX = TempPoint.X[0];
			ptFrontY = TempPoint.X[1];
			ptFrontZ = TempPoint.X[2];
			DiffZFact = 1.0;
		}
		else if(k == FTLE_CartMesh.ZRes - 1) {
			/* Back Point */
			TempPoint = Advect_FTLEPoint(i, j, k-1, t1, MinExitTime, FTLE_MeshPt);
			ptBackX = TempPoint.X[0];
			ptBackY = TempPoint.X[1];
			ptBackZ = TempPoint.X[2];
			/* Front Point */
			ptFrontX = PassedPoint.X[0];
			ptFrontY = PassedPoint.X[1];
			ptFrontZ = PassedPoint.X[2];
			DiffZFact = 1.0;
		}
		else {
			/* Back Point */
			TempPoint = Advect_FTLEPoint(i, j, k-1, t1, MinExitTime, FTLE_MeshPt);
			ptBackX = TempPoint.X[0];
			ptBackY = TempPoint.X[1];
			ptBackZ = TempPoint.X[2];
			/* Front Point */
			TempPoint = Advect_FTLEPoint(i, j, k+1, t1, MinExitTime, FTLE_MeshPt);
			ptFrontX = TempPoint.X[0];
			ptFrontY = TempPoint.X[1];
			ptFrontZ = TempPoint.X[2];
			DiffZFact = 2.0;
		}
		
		/* Linearized flow map */
		Dphi[0][2] = (ptFrontX - ptBackX) / (DiffZFact * FTLE_CartMesh.ZDelta);
		Dphi[1][2] = (ptFrontY - ptBackY) / (DiffZFact * FTLE_CartMesh.ZDelta);
		Dphi[2][0] = (ptRightZ - ptLeftZ) / (DiffXFact * FTLE_CartMesh.XDelta);
		Dphi[2][1] = (ptTopZ - ptBottomZ) / (DiffYFact * FTLE_CartMesh.YDelta);
		Dphi[2][2] = (ptFrontZ - ptBackZ) / (DiffZFact * FTLE_CartMesh.ZDelta);  
	}
	else { /* Dimensions == 2 */
		/* Linearized flow map */
		Dphi[0][2] = 0;
		Dphi[1][2] = 0;
		Dphi[2][0] = 0;
		Dphi[2][1] = 0;
		Dphi[2][2] = 0;
	}
	/* Linearized flow map */
	Dphi[0][0] = (ptRightX - ptLeftX) / (DiffXFact * FTLE_CartMesh.XDelta);
	Dphi[0][1] = (ptTopX - ptBottomX) / (DiffYFact * FTLE_CartMesh.YDelta);
	Dphi[1][0] = (ptRightY - ptLeftY) / (DiffXFact * FTLE_CartMesh.XDelta);
	Dphi[1][1] = (ptTopY - ptBottomY) / (DiffYFact * FTLE_CartMesh.YDelta);
	
	/* Deformation tensor */
	for(ii = 0; ii < 3; ii++) 
		for(jj = 0; jj < 3; jj++) {
			Delta[ii][jj] = 0;
			for(mm = 0; mm < 3; mm++)
				Delta[ii][jj] += Dphi[mm][ii] * Dphi[mm][jj];
		}
	
    lambda = GetMaxEigenvalue(Delta);  
    if(lambda > 1.0 - TINY) {
		if(fabs(MinExitTime - FTLE_Launches[ss].StartTime) < TINY) {
			FTLE_MeshPt[i][j][k].FTLEwT = 0.0;
			FTLE_MeshPt[i][j][k].FTLEwoT = 0.0;
		}
		else {
			FTLE_MeshPt[i][j][k].FTLEwT = (0.5 * log(lambda)) / fabs(MinExitTime - FTLE_Launches[ss].StartTime);
			FTLE_MeshPt[i][j][k].FTLEwoT = 0.5 * log(lambda);
		}
		
    }
    else {
		FTLE_MeshPt[i][j][k].FTLEwT = 0.0;
		FTLE_MeshPt[i][j][k].FTLEwoT = 0.0;
    }
	
    FTLE_MeshPt[i][j][k].HaveFTLE = 1;
	
}

void GetFTLEForPoint(int i, int j, int k, double IntTime, FTLEPoint ***FTLE_MeshPt) {
	
	double ptLeftX, ptLeftY, ptLeftZ;
	double ptRighX, ptRighY, ptRighZ;
	double DiffXFact;
	
	double ptBotX, ptBotY, ptBotZ;
	double ptTopX, ptTopY, ptTopZ;
	double DiffYFact;
	
	double ptBackX, ptBackY, ptBackZ;
	double ptFronX, ptFronY, ptFronZ;
	double DiffZFact;
	
	double Dphi[3][3];
	double Delta[3][3];
	double lambda;
	
	int  ii, jj, mm;
	
	if(IntTime < TINY) {
		FTLE_MeshPt[i][j][k].FTLEwT = 0.0;
		FTLE_MeshPt[i][j][k].FTLEwoT = 0.0;
	}
	else {
		
		if(i == 0) { /* At the minimum x value in grid */
			ptLeftX = FTLE_MeshPt[i][j][k].Pt.X[0];
			ptLeftY = FTLE_MeshPt[i][j][k].Pt.X[1];
			ptLeftZ = FTLE_MeshPt[i][j][k].Pt.X[2];
			ptRighX = FTLE_MeshPt[i+1][j][k].Pt.X[0];
			ptRighY = FTLE_MeshPt[i+1][j][k].Pt.X[1];
			ptRighZ = FTLE_MeshPt[i+1][j][k].Pt.X[2];
			DiffXFact   = 1.0;
		}
		else if(i == FTLE_CartMesh.XRes - 1) { /* At the maximum x value in grid */
			ptLeftX = FTLE_MeshPt[i-1][j][k].Pt.X[0];
			ptLeftY = FTLE_MeshPt[i-1][j][k].Pt.X[1];
			ptLeftZ = FTLE_MeshPt[i-1][j][k].Pt.X[2];
			ptRighX = FTLE_MeshPt[i][j][k].Pt.X[0];
			ptRighY = FTLE_MeshPt[i][j][k].Pt.X[1];
			ptRighZ = FTLE_MeshPt[i][j][k].Pt.X[2];
			DiffXFact = 1.0;
		} else { /* Central difference in x */ 
			ptLeftX = FTLE_MeshPt[i-1][j][k].Pt.X[0];
			ptLeftY = FTLE_MeshPt[i-1][j][k].Pt.X[1];
			ptLeftZ = FTLE_MeshPt[i-1][j][k].Pt.X[2];
			ptRighX = FTLE_MeshPt[i+1][j][k].Pt.X[0];
			ptRighY = FTLE_MeshPt[i+1][j][k].Pt.X[1];
			ptRighZ = FTLE_MeshPt[i+1][j][k].Pt.X[2];
			DiffXFact = 2.0;
		}
		
		if(j == 0) { /* At the minimum y value in grid */
			ptBotX = FTLE_MeshPt[i][j][k].Pt.X[0];
			ptBotY = FTLE_MeshPt[i][j][k].Pt.X[1];
			ptBotZ = FTLE_MeshPt[i][j][k].Pt.X[2];
			ptTopX  = FTLE_MeshPt[i][j+1][k].Pt.X[0];
			ptTopY  = FTLE_MeshPt[i][j+1][k].Pt.X[1];
			ptTopZ  = FTLE_MeshPt[i][j+1][k].Pt.X[2];
			DiffYFact = 1.0;
		}
		else if(j == FTLE_CartMesh.YRes - 1) { /* At the maximum y value in grid */
			ptBotX = FTLE_MeshPt[i][j-1][k].Pt.X[0];
			ptBotY = FTLE_MeshPt[i][j-1][k].Pt.X[1];
			ptBotZ = FTLE_MeshPt[i][j-1][k].Pt.X[2];
			ptTopX = FTLE_MeshPt[i][j][k].Pt.X[0];
			ptTopY = FTLE_MeshPt[i][j][k].Pt.X[1];
			ptTopZ = FTLE_MeshPt[i][j][k].Pt.X[2];
			DiffYFact = 1.0;
		}
		else { /* central difference in y */
			ptBotX = FTLE_MeshPt[i][j-1][k].Pt.X[0];
			ptBotY = FTLE_MeshPt[i][j-1][k].Pt.X[1];
			ptBotZ = FTLE_MeshPt[i][j-1][k].Pt.X[2];
			ptTopX = FTLE_MeshPt[i][j+1][k].Pt.X[0];
			ptTopY = FTLE_MeshPt[i][j+1][k].Pt.X[1];
			ptTopZ = FTLE_MeshPt[i][j+1][k].Pt.X[2];
			DiffYFact = 2.0; 
		}
		
		if(Dimensions == 3) {
			if(k == 0) { /* At the minimum z value in grid */
				ptBackX = FTLE_MeshPt[i][j][k].Pt.X[0];
				ptBackY = FTLE_MeshPt[i][j][k].Pt.X[1];
				ptBackZ = FTLE_MeshPt[i][j][k].Pt.X[2];
				ptFronX = FTLE_MeshPt[i][j][k+1].Pt.X[0];
				ptFronY = FTLE_MeshPt[i][j][k+1].Pt.X[1];
				ptFronZ = FTLE_MeshPt[i][j][k+1].Pt.X[2];
				DiffZFact = 1.0;
			}
			else if(k == FTLE_CartMesh.ZRes - 1) { /* At the maximum z value in grid */
				ptBackX = FTLE_MeshPt[i][j][k-1].Pt.X[0];
				ptBackY = FTLE_MeshPt[i][j][k-1].Pt.X[1];
				ptBackZ = FTLE_MeshPt[i][j][k-1].Pt.X[2];
				ptFronX = FTLE_MeshPt[i][j][k].Pt.X[0];
				ptFronY = FTLE_MeshPt[i][j][k].Pt.X[1];
				ptFronZ = FTLE_MeshPt[i][j][k].Pt.X[2];
				DiffZFact = 1.0;
			}
			else { /* central difference in z */
				ptBackX = FTLE_MeshPt[i][j][k-1].Pt.X[0];
				ptBackY = FTLE_MeshPt[i][j][k-1].Pt.X[1];
				ptBackZ = FTLE_MeshPt[i][j][k-1].Pt.X[2];
				ptFronX = FTLE_MeshPt[i][j][k+1].Pt.X[0];
				ptFronY = FTLE_MeshPt[i][j][k+1].Pt.X[1];
				ptFronZ = FTLE_MeshPt[i][j][k+1].Pt.X[2];
				DiffZFact = 2.0; 
			}
			/* Linearized flow map */
			Dphi[0][2] = (ptFronX - ptBackX) / (DiffZFact * FTLE_CartMesh.ZDelta);
			Dphi[1][2] = (ptFronY - ptBackY) / (DiffZFact * FTLE_CartMesh.ZDelta);
			Dphi[2][0] = (ptRighZ - ptLeftZ) / (DiffXFact * FTLE_CartMesh.XDelta);
			Dphi[2][1] = (ptTopZ - ptBotZ) / (DiffYFact * FTLE_CartMesh.YDelta);
			Dphi[2][2] = (ptFronZ - ptBackZ) / (DiffZFact * FTLE_CartMesh.ZDelta);  
		}
		else { /* Dimensions == 2 */
			/* Linearized flow map */
			Dphi[0][2] = 0;
			Dphi[1][2] = 0;
			Dphi[2][0] = 0;
			Dphi[2][1] = 0;
			Dphi[2][2] = 0;
		}
		/* Linearized flow map */
		Dphi[0][0] = (ptRighX - ptLeftX) / (DiffXFact * FTLE_CartMesh.XDelta);
		Dphi[0][1] = (ptTopX - ptBotX) / (DiffYFact * FTLE_CartMesh.YDelta);
		Dphi[1][0] = (ptRighY - ptLeftY) / (DiffXFact * FTLE_CartMesh.XDelta);
		Dphi[1][1] = (ptTopY - ptBotY) / (DiffYFact * FTLE_CartMesh.YDelta);
		
		/* Deformation tensor */
		for(ii = 0; ii < 3; ii++) {
			for(jj = 0; jj < 3; jj++) {
				Delta[ii][jj] = 0;
				for(mm = 0; mm < 3; mm++)
					Delta[ii][jj] += Dphi[mm][ii] * Dphi[mm][jj];
			}
		}
		
		lambda = GetMaxEigenvalue(Delta);
		if(lambda > 1.0 - TINY) {
			FTLE_MeshPt[i][j][k].FTLEwT = (0.5 * log(lambda)) / IntTime;
			FTLE_MeshPt[i][j][k].FTLEwoT = 0.5 * log(lambda);
		}
		else {
			FTLE_MeshPt[i][j][k].FTLEwT = 0.0;
			FTLE_MeshPt[i][j][k].FTLEwoT = 0.0;
		}
		
	}
	
	/* FTLE_MeshPt[i][j][k].HaveFTLE = 1; */
	
}

LagrangianPoint Advect_FTLEPoint(int i, int j, int k, double t1, double t2, FTLEPoint ***FTLE_MeshPt) {
	
	LagrangianPoint TempPoint;
	
	/* Initialize Structure */
	TempPoint.X[0] = FTLE_MeshPt[i][j][k].Pt.X[0];
	TempPoint.X[1] = FTLE_MeshPt[i][j][k].Pt.X[1];
	TempPoint.X[2] = FTLE_MeshPt[i][j][k].Pt.X[2];	
	TempPoint.LeftDomain = FTLE_MeshPt[i][j][k].Pt.LeftDomain;
	if(Data_MeshType == UNSTRUCTURED)
		TempPoint.ElementIndex = FTLE_MeshPt[i][j][k].Pt.ElementIndex;
	if(Particle_Radius > TINY) {
		TempPoint.V[0] = FTLE_MeshPt[i][j][k].Pt.V[0];
		TempPoint.V[1] = FTLE_MeshPt[i][j][k].Pt.V[1];
		TempPoint.V[2] = FTLE_MeshPt[i][j][k].Pt.V[2];
	}
	
	/* Integrate */
	Advect(&TempPoint, t1, t2); /* Default */
	
	return TempPoint;
	
}

void UpdateFTLELocations(FTLEPoint ***FTLE_MeshPt, double ****FTLE_NewArray) {     
	/* Copies contents of FTLE_NewArray to FTLE_MeshPt */
	int i, j, k;
	
	for(i = 0; i < FTLE_CartMesh.XRes; i++) 
		for(j = 0; j < FTLE_CartMesh.YRes; j++) 
			for( k = 0; k < FTLE_CartMesh.ZRes; k++) 
				if(!FTLE_MeshPt[i][j][k].Pt.LeftDomain || Int_Extrapolate) {
					FTLE_MeshPt[i][j][k].Pt.X[0] = FTLE_NewArray[i][j][k][0];
					FTLE_MeshPt[i][j][k].Pt.X[1] = FTLE_NewArray[i][j][k][1];
					FTLE_MeshPt[i][j][k].Pt.X[2] = FTLE_NewArray[i][j][k][2];
				}
	
}

void OutputFTLE(int ss, double IntT, FTLEPoint ***FTLE_MeshPt) {
	/* Write FTLE values to output file */
	
	int    i, j, k;
	char   OutFilePath[LONGSTRING], OutFilePath1[LONGSTRING], OutFilePath2[LONGSTRING];
	FILE   *OutFileID, *OutFileID1, *OutFileID2;
	char   BinFile[LONGSTRING];
	
	/* Output parameters */
	if((ss == 0 && !FTLE_ComputeVariation) || (FTLE_ComputeVariation && !FTLE_NumOutput)) {
		/* Create a file that lists key simulation parameters */
		sprintf(OutFilePath, "%s%s.info", Path_Output, FTLE_OutFilePrefix);
		if((OutFileID = fopen(OutFilePath, "w")) == NULL) 
			FatalError("Could not open file %s", OutFilePath);
		fprintf(OutFileID,"Int_Type\t\t=\t%10d\n",Int_Type);
		fprintf(OutFileID,"Int_TimeStep\t\t=\t%10g\n",Int_TimeStep);
		fprintf(OutFileID,"Int_Accuracy\t\t=\t%10g\n",Int_Accuracy);
		fprintf(OutFileID,"Int_MinTimeStep\t\t=\t%10g\n",Int_MinTimeStep);
		fprintf(OutFileID,"Int_MaxTimeStep\t\t=\t%10g\n",Int_MaxTimeStep);
		fprintf(OutFileID,"Int_TimeDirection\t=\t%10d\n",Int_TimeDirection);
		fprintf(OutFileID,"Int_NormalFlow\t\t=\t%10d\n",Int_NormalFlow);
		fprintf(OutFileID,"Int_NormalFlowScaling\t=\t%10g\n\n",Int_NormalFlowScaling);
		fprintf(OutFileID,"LocalSearchChecking\t=\t%10d\n\n",LocalSearchChecking);
		fprintf(OutFileID,"FTLE_GenerateMesh\t=\t%10d\n",FTLE_GenerateMesh);
		fprintf(OutFileID,"FTLE_ICFile\t\t=\t%s\n",FTLE_ICFile);
		fprintf(OutFileID,"FTLE_IntTLength\t\t=\t%10g\n",FTLE_IntTLength);
		fclose(OutFileID);
		
		/* Creat a file that lists mesh data */
		sprintf(OutFilePath, "%s%s_Cartesian.bin", Path_Output, FTLE_OutFilePrefix);
		if((OutFileID = fopen(OutFilePath, "wb")) == NULL) 
			FatalError("Could not open file %s", OutFilePath);
		if(fwrite(&FTLE_CartMesh.XMin, sizeof(double), 1, OutFileID) < 1 ||
		   fwrite(&FTLE_CartMesh.XMax, sizeof(double), 1, OutFileID) < 1 ||
		   fwrite(&FTLE_CartMesh.XRes, sizeof(int),    1, OutFileID) < 1 ||
		   fwrite(&FTLE_CartMesh.YMin, sizeof(double), 1, OutFileID) < 1 ||
		   fwrite(&FTLE_CartMesh.YMax, sizeof(double), 1, OutFileID) < 1 ||
		   fwrite(&FTLE_CartMesh.YRes, sizeof(int),    1, OutFileID) < 1 ||
		   fwrite(&FTLE_CartMesh.ZMin, sizeof(double), 1, OutFileID) < 1 ||
		   fwrite(&FTLE_CartMesh.ZMax, sizeof(double), 1, OutFileID) < 1 ||
		   fwrite(&FTLE_CartMesh.ZRes, sizeof(int),    1, OutFileID) < 1)
			FatalError("Could not write FTLE_CartMesh data to %s", OutFilePath);
		fclose(OutFileID);
	}
	
	/* Open output files for output number ss */
	sprintf(OutFilePath1, "%s%s.%d.bin", Path_Output, FTLE_OutFilePrefix, FTLE_ComputeVariation? FTLE_NumOutput:(Int_TimeDirection > 0? ss:(Output_TRes - 1 - ss)));
	if((OutFileID1 = fopen(OutFilePath1, "wb")) == NULL) 
		FatalError("Could not open file %s", OutFilePath1);
	sprintf(OutFilePath2, "%s%s_noT.%d.bin", Path_Output, FTLE_OutFilePrefix, FTLE_ComputeVariation? FTLE_NumOutput:(Int_TimeDirection > 0? ss:(Output_TRes - 1 - ss)));
	if((OutFileID2 = fopen(OutFilePath2, "wb")) == NULL) 
		FatalError("Could not open file %s", OutFilePath2);
	
	/* Write time stamp */ 
	if(fwrite(&FTLE_Launches[ss].StartTime, sizeof(double), 1, OutFileID1) < 1)
		FatalError("Could not write FTLE_Launches[%d].StartTime to %s", ss, OutFilePath1);
	if(fwrite(&FTLE_Launches[ss].StartTime, sizeof(double), 1, OutFileID2) < 1)
		FatalError("Could not write FTLE_Launches[%d].StartTime to %s", ss, OutFilePath2);
	
	/* Write data */    
	for(k = 0; k < FTLE_CartMesh.ZRes; k++) {
		for(j = 0; j < FTLE_CartMesh.YRes; j++) {
			for(i = 0; i < FTLE_CartMesh.XRes; i++) {
				if(!FTLE_MeshPt[i][j][k].HaveFTLE) 
					GetFTLEForPoint(i, j, k, IntT, FTLE_MeshPt);
				if(fwrite(&FTLE_MeshPt[i][j][k].FTLEwT, sizeof(double), 1, OutFileID1) < 1)
					FatalError("Could not write FTLE_MeshPt[%d][%d][%d].FTLEwT to %s", i, j, k, OutFilePath1);
				if(fwrite(&FTLE_MeshPt[i][j][k].FTLEwoT, sizeof(double), 1, OutFileID2) < 1)
					FatalError("Could not write FTLE_MeshPt[%d][%d][%d].FTLEwoT to %s", i, j, k, OutFilePath2);	
			}
		}
	}
	
	fclose(OutFileID1);
	fclose(OutFileID2);
	
	/* Unless computing variation of FTLE with integration time, delete binary file for slide ss */
	if(FTLE_ComputeVariation) 
		FTLE_NumOutput++; 
	else {
		sprintf(BinFile, "%s%s.%d.bin", Path_Data, FTLE_OutFilePrefix, ss);
		if(remove(BinFile))
			fprintf(stderr, "Warning: Could not delete file %s\n", BinFile);
	}
	
}

void FreeFTLEData(void) {
	
	int i, j, k;
	
	free(FTLE_Launches);
	FTLE_Launches = NULL;
	
	if (Output_TRes < opt_num*num_processors) {
		for(i = 0; i < FTLE_CartMesh.XRes; i++) {
			for(j = 0; j < FTLE_CartMesh.YRes; j++) {
				for(k = 0; k < FTLE_CartMesh.ZRes; k++) {
					free(FTLE_NewArray[i][j][k]);
					FTLE_NewArray[i][j][k] = NULL;
				}
				free(FTLE_NewArray[i][j]);
				FTLE_NewArray[i][j] = NULL;
				free(FTLE_MeshPt[i][j]);
				FTLE_MeshPt[i][j] = NULL;
			}
			free(FTLE_NewArray[i]);
			FTLE_NewArray[i] = NULL;
			free(FTLE_MeshPt[i]);
			FTLE_MeshPt[i] = NULL;
		}
		free(FTLE_NewArray);
		FTLE_NewArray = NULL;
		free(FTLE_MeshPt);
		FTLE_MeshPt = NULL;
	}
}

