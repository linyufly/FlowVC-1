/*
 *  tracers.c
 *  flowVC
 *
 *  Created by Shawn Shadden.
 *  Copyright 2010 Flow Physics Group. All rights reserved.
 *  
 *
 */


#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "globals.h"
#include "io.h"
#include "macros.h"
#include "memory.h"
#include "mesh.h"
#include "mymath.h"
#include "structs.h"
#include "tracers.h"
#include "velocity.h"

double Compute_Tracer (int *outputframe, double nextoutputtime, FILE *Trace_BinFileID) {



	double t1,t2;
	int ss, ii, num_integrated;
	
	
	// Set start of integration interval 
	t1 = (Int_TimeDirection > 0) ? Data_LoadedTMin : Data_LoadedTMax;
	
	
	
	// Integrate until end of loaded reached 
	while(((Int_TimeDirection > 0) ? (t1 - Data_LoadedTMax) : (Data_LoadedTMin - t1)) < 0) { 
	
		// Set end of integration interval to next output time, or end of loaded data 
		t2 = (Int_TimeDirection > 0) ? fmin(nextoutputtime, Data_LoadedTMax) : fmax(Data_LoadedTMin, nextoutputtime);
		
		if (TRACER_UPPERLOOP) { // Upperloop optimization ...
		
		
		// OMP implementation .....
		
		#pragma omp parallel private(ss, num_integrated, Trace_MeshPt, ii)
		{
		
			FILE *OutFileID;
			char OutFilePath[LONGSTRING];
			int tid = omp_get_thread_num();
			
			// Each thread will write data to separate output temp files ...
			// Open temp output files for the current thread
			//sprintf(OutFilePath, "%s%s_temp.%d.bin", Path_Data, Trace_OutFilePrefix, tid);
			//printf("Output temp file path = %s \n",OutFilePath);
			//if((OutFileID = fopen(OutFilePath, "ab+")) == NULL) 
			//	FatalError("Could not open file %s", OutFilePath);
		
			// Allocate memory for Trace_MeshPt
			if((Trace_MeshPt = (LagrangianPoint *)malloc(Trace_NumTracers * sizeof(LagrangianPoint))) == NULL) 
				FatalError("Malloc failed for Trace_MeshPt");
		
			// Loop over each release 
			#pragma omp for schedule(dynamic) nowait
			for(ss = 0; ss < Trace_NumLaunchTimes; ss++) { 
				if(Trace_Launches[ss].Status != COMPLETE && (Trace_Launches[ss].Status == LAUNCHED || (Trace_Launches[ss].StartTime <= fmax(t1, t2) && Trace_Launches[ss].StartTime >= fmin(t1, t2)))) {
				
					// Load tracer data from file to Trace_MeshPt 
					Trace_MeshPt = ReadInTraceLaunch(ss, Trace_MeshPt);
				
				
					// Initiate release if needed
					if(Trace_Launches[ss].Status == UNLAUNCHED) {
					
						Trace_Launches[ss].Status = LAUNCHED;
					
						// Update encounter counter for elements containing each release point 
						if(Trace_CETCompute)
							for(ii = 0; ii < Trace_NumTracers; ii++) 
								if(!Trace_MeshPt[ii].LeftDomain)
									Trace_CETArray[Trace_CETAuxillaryMesh? Trace_MeshPt[ii].AuxElementIndex : Trace_MeshPt[ii].ElementIndex].Encounters++;	
				
						// Reset start of integration to release time
						t1 = Trace_Launches[ss].StartTime;
					
						// Initialize velocity of each particle
						for(ii = 0; ii < Trace_NumTracers; ii++) {
							if(!Trace_MeshPt[ii].LeftDomain) {
								if(Particle_Radius > TINY && Particle_ICType == 0) {
									Trace_MeshPt[ii].V[0] = 0.0;
									Trace_MeshPt[ii].V[1] = 0.0;
									Trace_MeshPt[ii].V[2] = 0.0;
								}
								else
									GetVelocity(Trace_Launches[ss].StartTime, &Trace_MeshPt[ii], Trace_MeshPt[ii].V);
							}
						}
						printf("Initiated release %d for %d tracers with a launch time at %g\n", ss, Trace_NumTracers, t1); fflush(stdout);
					}	
				
					//printf("Integrating release %d from %g to %g...", ss, t1, t2); fflush(stdout);
					num_integrated = 0;
				
					for(ii = 0; ii < Trace_NumTracers; ii++) {
						if(Trace_MeshPt[ii].LeftDomainTime > -TINY && (!Trace_MeshPt[ii].LeftDomain || Int_Extrapolate)) { 
							Advect(&Trace_MeshPt[ii], t1, t2);
							num_integrated++;
						}
					}
					printf("Release %d from %g to %g integrated ... %d tracers integrated.\n",ss,t1,t2, num_integrated); 
					fflush(stdout);
					// If nextoutputtime reached, output tracer positions to file
					if(fabs(nextoutputtime - t2) < TINY) {
						//printf("Outputing for output frame = %d and tid = %d nextoutput time = %f and t2 = %f \n", *outputframe, tid, nextoutputtime,t2);
						// Loop over all tracers in the current release	

					

						// Open Output file
						sprintf(OutFilePath, "%s%s_release%d.%d.bin", Path_Data, Trace_OutFilePrefix, ss, *outputframe);
						//printf("Output temp file path = %s \n",OutFilePath);
						if((OutFileID = fopen(OutFilePath, "wb")) == NULL) 
							FatalError("Could not open file %s", OutFilePath);
					
						for(ii = 0; ii < Trace_NumTracers; ii++) {
				
							if(Trace_MeshPt[ii].LeftDomainTime > -TINY && (!Trace_MeshPt[ii].LeftDomain || Trace_AlwaysOutput)) {
				
								if(fwrite(Trace_MeshPt[ii].X, sizeof(double), 3, OutFileID) < 3) 
									FatalError("Could not completely write Trace_MeshPt[%d].X data to %s", ii, TraceOUT_BinFilePath);
				
								if(Trace_APCompute)  
									if(fwrite(&Trace_MeshPt[ii].Scalar, sizeof(double), 1, OutFileID) < 1)
										FatalError("Could not write Trace_MeshPt[%d].Scalar value to %s", ii, TraceOUT_BinFilePath);
							
								Trace_NumOutput_OMP[*outputframe][ss]++;
							}
						}	
					
						fclose(OutFileID);					
					}
					
				
					WriteOutTraceLaunch(ss,Trace_MeshPt); /* Saves Trace_MeshPt array to binary file for current slide */
						   
				}
			
				if(((Int_TimeDirection > 0) ? (Trace_Launches[ss].StopTime - t2) : (t2 - Trace_Launches[ss].StopTime)) < TINY) 
					Trace_Launches[ss].Status = COMPLETE;
			
			} // End of omp for
		
		
		
			// Free Trace_MeshPt 
			free(Trace_MeshPt);  
			Trace_MeshPt = NULL;
		
		} // end of omp parallel ...
		
		} else if (TRACER_LOWERLOOP){ // LOWER LOOP OPTIMIZATION ...
		
			
		//	FILE   *Trace_BinFileID;
			
			
			/* Loop over each release */
			for(ss = 0; ss < Trace_NumLaunchTimes; ss++) { 
				/* Consider only releases that have been launched, or need to be launched */
				if(Trace_Launches[ss].Status != COMPLETE && 
				   (Trace_Launches[ss].Status == LAUNCHED || (Trace_Launches[ss].StartTime <= fmax(t1, t2) && Trace_Launches[ss].StartTime >= fmin(t1, t2)))) {
							
					/* Load tracer data from file to Trace_MeshPt */
					Trace_MeshPt = ReadInTraceLaunch(ss, Trace_MeshPt);
					//printf("%f \n", Trace_MeshPt[23].X[0]);
					/* Initiate release if needed */
					if(Trace_Launches[ss].Status == UNLAUNCHED) { 
						Trace_Launches[ss].Status = LAUNCHED;
						/* Update encounter counter for elements containing each release point */
						if(Trace_CETCompute)
							for(ii = 0; ii < Trace_NumTracers; ii++) 
								if(!Trace_MeshPt[ii].LeftDomain)
									Trace_CETArray[Trace_CETAuxillaryMesh? Trace_MeshPt[ii].AuxElementIndex : Trace_MeshPt[ii].ElementIndex].Encounters++;
						/* Reset start of integration to release time */
						t1 = Trace_Launches[ss].StartTime;
						/* Initialize velocity of each particle */
						// OpenMP Implementation # 2	
						
						int dynamicnum = floor(Trace_NumTracers/num_processors);
						# pragma omp parallel for private(ii) schedule(dynamic, dynamicnum)
						for(ii = 0; ii < Trace_NumTracers; ii++) {
						
							if(!Trace_MeshPt[ii].LeftDomain) {
								if(Particle_Radius > TINY && Particle_ICType == 0) {
									Trace_MeshPt[ii].V[0] = 0.0;
									Trace_MeshPt[ii].V[1] = 0.0;
									Trace_MeshPt[ii].V[2] = 0.0;
								}
								else
									GetVelocity(Trace_Launches[ss].StartTime, &Trace_MeshPt[ii], Trace_MeshPt[ii].V);
							}
						
						}
						printf("Initiated release %d for %d tracers with a launch time at %g\n", ss, Trace_NumTracers, t1); fflush(stdout);
					} 
							
					printf("Integrating release %d from %g to %g...", ss, t1, t2); fflush(stdout);
					num_integrated = 0;
					/* Loop over tracers */
					int num_int_thread;
							
					// OpenMP Implementation # 1
					# pragma omp parallel private(num_int_thread)
					{
						//threadid = omp_get_thread_num();
						num_int_thread = 0;
						//printf("%d thread forked .. \n\n",omp_get_thread_num());
						int dynamicnum = floor(Trace_NumTracers/num_processors);
						# pragma omp for private(ii) schedule(dynamic, dynamicnum) nowait
						for(ii = 0; ii < Trace_NumTracers; ii++) {
							if(Trace_MeshPt[ii].LeftDomainTime > -TINY && (!Trace_MeshPt[ii].LeftDomain || Int_Extrapolate)) 										{ 
								Advect(&Trace_MeshPt[ii], t1, t2);
								num_int_thread++;
							}
						}


					// Adding the num_int_thread counter to the main num_integrated counter. 
						#pragma omp critical
						{
							num_integrated = num_integrated + num_int_thread;
							
									

						} // End of omp critical block	
					} // End of parallel section


					for(ii = 0; ii < Trace_NumTracers; ii++) {
			
						/* If nextoutputtime reached, output tracer positions to file */
						if(fabs(nextoutputtime - t2) < TINY) {
							if(Trace_MeshPt[ii].LeftDomainTime > -TINY && (!Trace_MeshPt[ii].LeftDomain || Trace_AlwaysOutput)) {
								if(fwrite(Trace_MeshPt[ii].X, sizeof(double), 3, Trace_BinFileID) < 3) 
									FatalError("Could not completely write Trace_MeshPt[%d].X data to %s", ii, TraceOUT_BinFilePath);
								if(Trace_APCompute)  
									if(fwrite(&Trace_MeshPt[ii].Scalar, sizeof(double), 1, Trace_BinFileID) < 1)
										FatalError("Could not write Trace_MeshPt[%d].Scalar value to %s", ii, TraceOUT_BinFilePath);
								Trace_NumOutput[*outputframe]++;
							}		    
						}
					}
					printf(" %d tracers integrated.\n", num_integrated); 
					fflush(stdout);
					WriteOutTraceLaunch(ss, Trace_MeshPt); /* Saves Trace_MeshPt array to binary file for current slide */
				}
				
				if(((Int_TimeDirection > 0) ? (Trace_Launches[ss].StopTime - t2) : (t2 - Trace_Launches[ss].StopTime)) < TINY) 
					Trace_Launches[ss].Status = COMPLETE;
			}
			
			
		
		} else {
		
			FatalError("ERROR: Compiler did not recognize type of optimized loop for Tracer calculation. \n");
		
		}
		
		if(fabs(t2 - nextoutputtime) < TINY) {
			nextoutputtime += Int_TimeDirection * Output_TDelta;
			//printf("nextoutputtime = %f \n",nextoutputtime);
			*outputframe= *outputframe+1;
			//printf("next output frame = %d \n\n", *outputframe);
			if(*outputframe >= Output_TRes)
				break;
		}
		t1 = t2;

	} // End of while loop ...
	
	
	return nextoutputtime;
	
}



void Tracer_StaggeredRelease(FILE *Trace_BinFileID, FileData *Trace_BinFileArray) {
	
				int num_integrated = 0; 
				int ii,outputcount, ss;
				double t1,t2, nextoutputtime;
				ReleaseLocation *rnode;
				 
				 
				
				/* Loop over all release locations */
				int num_int_thread;
				int dyn_num = floor(Trace_NumReleaseLocations/ omp_get_num_procs());
				#pragma omp parallel private(num_int_thread)
				{
		
					num_int_thread = 0;

				# pragma omp for private(ss,t1,t2,ii,rnode, outputcount,nextoutputtime) schedule(dynamic,2) nowait
				for(ss = 0; ss < Trace_NumReleaseLocations;  ss++) {
					// Loop over all release times 
					for(rnode = Trace_ReleaseList[ss]; rnode != NULL; rnode = rnode->next) {
						// Make sure point launched, or needs to be launched and not complete 
						if(rnode->slide.Status != COMPLETE && 
						   (rnode->slide.Status == LAUNCHED || (rnode->slide.StartTime >= Data_LoadedTMin && rnode->slide.StartTime <= Data_LoadedTMax))) { 
							// Make sure point still in domain  
							if(!rnode->pt.LeftDomain || Int_Extrapolate) {
								// Check if point hasn't been launched yet
								if(rnode->slide.Status == UNLAUNCHED) {
									rnode->slide.Status = LAUNCHED;
									// Increment number of encounters in element where point is launched
									if(Trace_CETCompute)
										Trace_CETArray[Trace_CETAuxillaryMesh? rnode->pt.AuxElementIndex : rnode->pt.ElementIndex].Encounters++;
									// If launchtime is Output_TStart, then write out locations
									if(fabs(rnode->slide.StartTime - Output_TStart) < TINY) {
										if(fwrite(rnode->pt.X, sizeof(double), 3, Trace_BinFileArray[0].FileID) < 3) 
											FatalError("Could not completely write rnode->pt.X to %s", Trace_BinFileArray[0].FilePath); 
										if(Trace_APCompute)  
											if(fwrite(&(rnode->pt.Scalar), sizeof(double), 1, Trace_BinFileArray[0].FileID) < 1)
												FatalError("Could not write rnode->pt.Scalar value to %s", Trace_BinFileArray[0].FilePath);
										Trace_NumOutput[0]++;
									}
								} 
								// Integrate point forward
								t1 = fmax(rnode->slide.StartTime, Data_LoadedTMin);
								t2 = fmin(rnode->slide.StopTime, Data_LoadedTMax);
								// See how many times we need to output data over time interval from t1 to t2
								outputcount = 0;
								for(ii = 0; ii < Output_TRes; ii++) {
									if((Output_TStart + ii*Output_TDelta > t1 + TINY) && (Output_TStart + ii*Output_TDelta <= t2 + TINY))
										outputcount++;
										
								}
															
								// Break up integration if needed
								if(outputcount == 0)  // not needed
									Advect(&rnode->pt, t1, t2);
								else { //break up integration 
									// Determine next output time 
									for(ii = 0; ii < Output_TRes; ii++) {
										nextoutputtime = Output_TStart + ii*Output_TDelta;
										if(nextoutputtime > t1 + TINY) 
											break;
									}
									t2 = nextoutputtime;
									Advect(&rnode->pt, t1, t2);
									// Output to file 
									if(!rnode->pt.LeftDomain || Int_Extrapolate) {
										if(fwrite(rnode->pt.X, sizeof(double), 3, Trace_BinFileArray[ii].FileID) < 3) 
											FatalError("Could not completely write rnode->pt.X to %s", Trace_BinFileArray[ii].FilePath); 
										if(Trace_APCompute)  
											if(fwrite(&(rnode->pt.Scalar), sizeof(double), 1, Trace_BinFileArray[ii].FileID) < 1)
												FatalError("Could not write rnode->pt.Scalar value to %s", Trace_BinFileArray[ii].FilePath);
										Trace_NumOutput[ii]++;
									}
									// Break up integration for further output times (if needed)
									for(ii = 1; ii < outputcount; ii++) {
										if(!rnode->pt.LeftDomain || Int_Extrapolate) {
											t1 = t2;
											t2 = t1 + Output_TDelta;
											Advect(&rnode->pt, t1, t2);
											// Output to file
											if(!rnode->pt.LeftDomain || Int_Extrapolate) {
												if(fwrite(rnode->pt.X, sizeof(double), 3, Trace_BinFileArray[(int)((t2 + TINY)/Output_TDelta)].FileID) < 3) 
													FatalError("Could not completely write rnode->pt.X to %s", Trace_BinFileArray[(int)((t2 + TINY)/Output_TDelta)].FilePath); 
												if(Trace_APCompute)  
													if(fwrite(&(rnode->pt.Scalar), sizeof(double), 1, Trace_BinFileArray[(int)((t2 + TINY)/Output_TDelta)].FileID) < 1)
														FatalError("Could not write rnode->pt.Scalar value to %s", Trace_BinFileArray[(int)((t2 + TINY)/Output_TDelta)].FilePath);
												Trace_NumOutput[(int)((t2 + TINY)/Output_TDelta)]++;
											} 
										}
									}
									// Continue until end of interval or StopTime
									if(!rnode->pt.LeftDomain || Int_Extrapolate) {
										t1 = t2;
										t2 = fmin(rnode->slide.StopTime, Data_LoadedTMax);
										if(t2 - t1 > 0) {
											Advect(&rnode->pt, t1, t2);
										}
										else if(t2 - t1 < -TINY)
											FatalError("Integration of point beyond requested stop time or range of loaded data");
									}
								}	     
								if(rnode->slide.StopTime <= Data_LoadedTMax)
									rnode->slide.Status = COMPLETE; 
								num_int_thread++; 
							} // Particle inside domain 
						} // If particle should be advected during this time interval
					} // For each particle released from this location 
				} // For each release location // End of omp for section

				#pragma omp critical
				{
				//found = found + found_thread;	
					//printf("Thread %d adding its iterations (%d) to the sum (%d)...\n",omp_get_thread_num(), num_int_thread, num_integrated);
					num_integrated = num_integrated + num_int_thread;
									
									
				} // End of omp critical block	


				} // End of omp parallel section .....

				printf(" %d tracers integrated.\n", num_integrated); fflush(stdout);	
}

void GenerateStaggeredRelease(void) {
	
	double M, N, ti, tii, a, b, c, d, voxdim, *voxdim0, usbar, uebar, uibar, us[3], ue[3], **u0;
	int ss, ii, df, nff, *ptdone;
	ReleaseLocation *nnode, **pnode;
	
	LoadReleasePoints(&voxdim);  /* Sets voxdim, Trace_ReleasePoints array and Trace_NumReleaseLocations */ 
	
	printf("\nDetermining release times...");
	fflush(stdout);
	
	/* Allocate memory for global linked list that will store release information */
	if((Trace_ReleaseList = (ReleaseLocation **)malloc(Trace_NumReleaseLocations * sizeof(ReleaseLocation *))) == NULL) 
		FatalError("Malloc failed for Trace_ReleaseList");
	
	/* Allocate memory for pnode, which is a local linked list used to connect nodes of Trace_ReleaseList */
	if((pnode = (ReleaseLocation **)malloc(Trace_NumReleaseLocations * sizeof(ReleaseLocation *))) == NULL) 
		FatalError("Malloc failed for pnode");
	
	/* Allocate memory for voxdim0, which keeps track of distance release points integrated from previous cycle(s) */
	if((voxdim0 = (double *)calloc(Trace_NumReleaseLocations, sizeof(double))) == NULL)
		FatalError("Calloc failed for voxdim0");
	
	/* Allocate memory for u0 */
	if((u0 = (double **)malloc(Trace_NumReleaseLocations * sizeof(double *))) == NULL)
		FatalError("Malloc failed for u0");
	for(ss = 0; ss < Trace_NumReleaseLocations; ss++) 
		if((u0[ss] = (double *)malloc(3 * sizeof(double))) == NULL)
			FatalError("Malloc failed for u0[%d]", ss);
	
	/* Allocate memory for ptdone */
	if((ptdone = (int *)calloc(Trace_NumReleaseLocations, sizeof(int))) == NULL)
		FatalError("Malloc failed for ptdone");
	
	Trace_NumTracers = 0; /* number of points that get released */
	
	/* Initialize first release to Output_TStart */
	for(ss = 0; ss < Trace_NumReleaseLocations; ss++) {
		if(Trace_ReleasePoints[ss].LeftDomain) 
			Trace_ReleaseList[ss] = NULL;
		else {
			/* Alocate memory for what will be the first node in the global linked list Trace_ReleaseList */
			if((nnode = (ReleaseLocation *)malloc(sizeof(ReleaseLocation))) == NULL)
				FatalError("Malloc failed for nnode");
			
			/* Initialize values */
			nnode->slide.StartTime = Output_TStart;
			nnode->slide.StopTime = fmin(Output_TEnd, nnode->slide.StartTime + Trace_IntTLength);
			nnode->slide.Status = UNLAUNCHED;
			for(ii = 0; ii < 3; ii++) 
				nnode->pt.X[ii] = Trace_ReleasePoints[ss].X[ii];
			nnode->pt.LeftDomain = 0;
			nnode->pt.ElementIndex = Trace_ReleasePoints[ss].ElementIndex;
			nnode->pt.Scalar = 0.0;
			if(Particle_Radius > TINY) {
				if(Particle_ICType == 1) /* Set initial velocity of particle to fluid velocity */
					GetVelocity(nnode->slide.StartTime, &(nnode->pt), nnode->pt.V);
				else /* Set initial velocity to zero */
					for(ii = 0; ii < 3; ii++) 
						nnode->pt.V[ii] = 0.0;
			}
			if(Trace_CETCompute && Trace_CETAuxillaryMesh)
				nnode->pt.AuxElementIndex = Trace_ReleasePoints[ss].AuxElementIndex;
			/* Link node to lists */
			Trace_ReleaseList[ss] = nnode;
			pnode[ss] = nnode;  
			Trace_NumTracers++;
		}
	}
	
	/* If normal flow on, temporarily turn it off */
	if(Int_NormalFlow) {
		Int_NormalFlow = 0;
		nff = 1;
	}    
	else
		nff = 0;
	
	/* Load first velocity data frame */
	if(Data_MeshType == CARTESIAN) 
		LoadCartVelDataFrame(Data_FirstFrame);
	else if(Data_MeshType == UNSTRUCTURED) 
		LoadUnstructVelDataFrame(Data_FirstFrame);
	
	/* Step through the velocity data and determine release times */
	for(df = Data_FirstFrame; df < Data_LastFrame; df++) {
		
		/* Load next velocity data slide */
		if(Data_MeshType == CARTESIAN) 
			LoadCartVelDataFrame(df + 1);
		else if(Data_MeshType == UNSTRUCTURED) 
			LoadUnstructVelDataFrame(df + 1);
		
		/* Determine time interval of loaded data */
		Data_LoadedTMin = Data_TMin + df * Data_TDelta;
		Data_LoadedTMax = Data_LoadedTMin + Data_TDelta;
		printf("Data loaded in memory spans t = %g to %g\n", Data_LoadedTMin, Data_LoadedTMax);
		fflush(stdout);
		
		/* Store velocity from Data_FirstFrame for each release point--it will be used as the "positive" flow direction */
		if(df == Data_FirstFrame) {
			/* Loop over release points */
			for(ss = 0; ss < Trace_NumReleaseLocations; ss++) {
				/* Only consider release points located inside the velocity domain */
				if(!Trace_ReleasePoints[ss].LeftDomain) {
					/* Interpolate velocity at release point at start of interval */
					if(Data_MeshType == CARTESIAN) {
						GetVelocity_Cartesian(Data_LoadedTMin, &Trace_ReleasePoints[ss], u0[ss]);
					}     
					else if(Data_MeshType == UNSTRUCTURED) { 
						GetVelocity_Unstructured(Data_LoadedTMin, &Trace_ReleasePoints[ss], u0[ss]);
					}
				}
			}
		}
		
		if(Data_LoadedTMin >= Trace_ReleaseTMax) 
			break;
		
		/* Loop over release points and determine release times during interval of currently loaded data */
		for(ss = 0; ss < Trace_NumReleaseLocations; ss++) {
			/* Only consider release points located inside the velocity domain */
			if(!Trace_ReleasePoints[ss].LeftDomain && !ptdone[ss]) {
				/* Interpolate velocity at release point at bounds of loaded data */
				if(Data_MeshType == CARTESIAN) {
					GetVelocity_Cartesian(Data_LoadedTMin, &Trace_ReleasePoints[ss], us);
					GetVelocity_Cartesian(Data_LoadedTMax, &Trace_ReleasePoints[ss], ue);
				}     
				else if(Data_MeshType == UNSTRUCTURED) { 
					GetVelocity_Unstructured(Data_LoadedTMin, &Trace_ReleasePoints[ss], us);
					GetVelocity_Unstructured(Data_LoadedTMax, &Trace_ReleasePoints[ss], ue);
				}	
				/* Calculate "signed magnitude" of velocity */
				usbar = SIGN(sqrt(vdot(us, us, 3)), vdot(us, u0[ss], 3));
				uebar = SIGN(sqrt(vdot(ue, ue, 3)), vdot(ue, u0[ss], 3));
				
				/*printf("Point %d: usbar = %.9f, uebar = %.9f\n", ss, usbar, uebar);*/
				
				if(usbar < TINY && uebar < TINY) { /* Reverse or stagnant flow, do nothing */
					/*printf("Reverse flow, continuing\n");*/
					continue;
				}
				else if(fabs(usbar - uebar) < 1.0e-6) { /* Flow steady, solve linear equation for next release */
					/*printf("Steady flow\n");*/
					
					if(voxdim0[ss] > TINY)
						ti = Data_LoadedTMin;
					else if(voxdim0[ss] < -TINY) {
						ti = Data_LoadedTMin;
						voxdim0[ss] = 0;
					}
					else
						ti = pnode[ss]->slide.StartTime;
					
					while(1) {
						/*printf("ti = %g, %s", ti);*/
						d = voxdim - voxdim0[ss];
						if(d < 0) 
							FatalError("d < 0"); 
						
						tii = ti + d / usbar; /* Linear equation */
						
						/*printf("d = %g, tii = %g, %s", d, tii);*/
						
						if(tii <= ti) 
							FatalError("Next output time not after previous");
						else {
							if(tii >= Trace_ReleaseTMax) { /* Release time greater than requested, break */
								/*printf("Terminating point %d\n", ss);*/
								ptdone[ss] = 1;
								break; 
							}
							else if(tii > Data_LoadedTMax) {
								/* Save distance point will travel during the remainder of current interval and break */
								voxdim0[ss] = voxdim0[ss] + uebar * (Data_LoadedTMax - ti);
								/*printf("tii > Data_LoadedTMax, setting voxdim0[%d] = %.9f, %.9f\n", ss, voxdim0[ss], 0.5 * (uebar + uibar) * (Data_LoadedTMax - ti)); */
								if(voxdim0[ss] < 0) 
									FatalError("voxdim0 < 0");
								break;
							}
							else {
								/*printf("creating new release time at %g\n", tii);*/
								/* Create next entry in global linked list Trace_ReleaseList */ 
								if((nnode = (ReleaseLocation *)malloc(sizeof(ReleaseLocation))) == NULL)
									FatalError("Malloc failed for nnode");
								nnode->slide.StartTime = tii;
								nnode->slide.StopTime = fmin(Output_TEnd, nnode->slide.StartTime + Trace_IntTLength);
								nnode->slide.Status = UNLAUNCHED;
								for(ii = 0; ii < 3; ii++) 
									nnode->pt.X[ii] = Trace_ReleasePoints[ss].X[ii];
								nnode->pt.LeftDomain = 0;
								nnode->pt.ElementIndex = Trace_ReleasePoints[ss].ElementIndex;
								nnode->pt.Scalar = 0.0;
								if(Particle_Radius > TINY) {
									if(Particle_ICType == 1) /* Set initial velocity of particle to fluid velocity */
										GetVelocity(nnode->slide.StartTime, &(nnode->pt), nnode->pt.V);
									else /* Set initial velocity to zero */
										for(ii = 0; ii < 3; ii++) 
											nnode->pt.V[ii] = 0.0;
								}
								if(Trace_CETCompute && Trace_CETAuxillaryMesh)
									nnode->pt.AuxElementIndex = Trace_ReleasePoints[ss].AuxElementIndex;
								/* Link node to lists */
								pnode[ss]->next = nnode;	  
								pnode[ss] = nnode;
								/* Reset voxdim0[ss] */
								voxdim0[ss] = 0;
								/* Update ti */
								ti = tii;
								Trace_NumTracers++;
							}
						}
					}	  
				}
				else { /* Flow changing, solve quadratic equation for next release */
					/* Velocity varies linearly in time, e.g. u = M*t + N, solve for M and N */
					M = (usbar - uebar) / (Data_LoadedTMin - Data_LoadedTMax);
					N = usbar - M * Data_LoadedTMin;
					
					if(usbar < 0 && uebar > 0) {
						/* Initiate a release when flow becomes positive */
						if((nnode = (ReleaseLocation *)malloc(sizeof(ReleaseLocation))) == NULL)
							FatalError("Malloc failed for nnode");
						nnode->slide.StartTime = -N/M;
						nnode->slide.StopTime = fmin(Output_TEnd, nnode->slide.StartTime + Trace_IntTLength);
						nnode->slide.Status = UNLAUNCHED;
						for(ii = 0; ii < 3; ii++) 
							nnode->pt.X[ii] = Trace_ReleasePoints[ss].X[ii];
						nnode->pt.LeftDomain = 0;
						nnode->pt.ElementIndex = Trace_ReleasePoints[ss].ElementIndex;
						nnode->pt.Scalar = 0.0;
						if(Particle_Radius > TINY) {
							if(Particle_ICType == 1) /* Set initial velocity of particle to fluid velocity */
								GetVelocity(nnode->slide.StartTime, &(nnode->pt), nnode->pt.V);
							else /* Set initial velocity to zero */
								for(ii = 0; ii < 3; ii++) 
									nnode->pt.V[ii] = 0.0;
						}
						if(Trace_CETCompute && Trace_CETAuxillaryMesh)
							nnode->pt.AuxElementIndex = Trace_ReleasePoints[ss].AuxElementIndex;
						/* Link node to lists */
						pnode[ss]->next = nnode;	  
						pnode[ss] = nnode;
						/* Reset voxdim0[ss] */
						voxdim0[ss] = 0;
						ti = -N/M;
						Trace_NumTracers++;
					}
					else {
						if(voxdim0[ss] > TINY)
							ti = Data_LoadedTMin;
						else if(voxdim0[ss] < -TINY) {
							ti = Data_LoadedTMin;
							voxdim0[ss] = 0;
						}
						else
							ti = pnode[ss]->slide.StartTime;
					}
					
					
					while(1) { 
						/*printf("ti = %g, %s", ti);*/
						
						/* Interpolate velocity at time ti */
						uibar = M*ti + N;
						
						/* Set target distance to advect point */
						d = voxdim - voxdim0[ss];
						if(d < 0) {fprintf(stderr, "Error: d < 0\n"); exit(1);}
						
						/* Set up quadratic equation for next release time, tii = (-b + sqrt(b * b - 4 * a * c)) / (2 * a) */ 
						a = M;
						b = N + uibar - M * ti;
						c = -N * ti - uibar * ti - 2*d;
						
						/* 
						 Since velocity is linear in time, distance a point travels from time ti to time tii is given by 
						 d = 0.5 * (uiibar + uibar) * (tii - ti), where uiibar is the velocity magnitude at time tii. 
						 However uiibar = M * tii + N, so we solve d = 0.5 * (M * tii + N + uibar) * (tii - ti), or rewriting, 
						 M * tii^s + (N + uibar - M * ti) * tii + (-N * ti - uibar * ti - 2d) = 0 
						 */
						
						/* Check if imaginary roots */
						if(b * b - 4 * a * c < 0) { 
							if(usbar < uebar)  /* Flow rate increasing, solution should exist */
								FatalError("Next output time has imaginary roots even with increasing flow");
							else { /* Flow rate must be decreasing too fast, break */
								/*printf("flow rate decreasing too fast, breaking\n");*/
								voxdim0[ss] = -1; /* Used so next consideration of point ss starts from Data_LoadedTMin */
								break;
							}
						}
						else {
							/* Compute next release time */
							tii = (-b + sqrt(b * b - 4 * a * c)) / (2 * a);
							
							/*printf("d = %g, tii = %g, ", d, tii); */
							
							if(tii <= ti) 
								FatalError("Next output time not after previous");
							else {
								if(tii >= Trace_ReleaseTMax) { /* Release time greater than requested, break */
									/*printf("Terminating point %d\n", ss);*/
									ptdone[ss] = 1;
									break; 
								}
								else if(tii > Data_LoadedTMax) {
									/* Save distance point will travel during the remainder of current interval and break */
									voxdim0[ss] = voxdim0[ss] + 0.5 * (uebar + uibar) * (Data_LoadedTMax - ti);
									/*printf("tii > Data_LoadedTMax, setting voxdim0[%d] = %.9f, %.9f\n", ss, voxdim0[ss], 0.5 * (uebar + uibar) * (Data_LoadedTMax - ti)); */
									if(voxdim0[ss] < 0) 
										FatalError("voxdim0 < 0");
									break;
								}
								else {
									/*printf("creating new release time at %g\n", tii);*/
									/* Create next entry in global linked list Trace_ReleaseList */ 
									if((nnode = (ReleaseLocation *)malloc(sizeof(ReleaseLocation))) == NULL)
										FatalError("Malloc failed for nnode");
									nnode->slide.StartTime = tii;
									nnode->slide.StopTime = fmin(Output_TEnd, nnode->slide.StartTime + Trace_IntTLength);
									nnode->slide.Status = UNLAUNCHED;
									for(ii = 0; ii < 3; ii++) 
										nnode->pt.X[ii] = Trace_ReleasePoints[ss].X[ii];
									nnode->pt.LeftDomain = 0;
									nnode->pt.ElementIndex = Trace_ReleasePoints[ss].ElementIndex;
									nnode->pt.Scalar = 0.0;
									if(Particle_Radius > TINY) {
										if(Particle_ICType == 1) /* Set initial velocity of particle to fluid velocity */
											GetVelocity(nnode->slide.StartTime, &(nnode->pt), nnode->pt.V);
										else /* Set initial velocity to zero */
											for(ii = 0; ii < 3; ii++) 
												nnode->pt.V[ii] = 0.0;
									}
									if(Trace_CETCompute && Trace_CETAuxillaryMesh)
										nnode->pt.AuxElementIndex = Trace_ReleasePoints[ss].AuxElementIndex;
									/* Link node to lists */
									pnode[ss]->next = nnode;	  
									pnode[ss] = nnode;
									/* Reset voxdim0[ss] */
									voxdim0[ss] = 0;
									/* Update ti */
									ti = tii;
									Trace_NumTracers++;
								}
							}
						}
					}
				}
			}
		}
	}
	
	
	printf("OK!\n");
	printf("Total number of points to be released is %d\n", Trace_NumTracers);
	fflush(stdout);
	
	/* Add terminator to last node of global linked list Trace_ReleaseList */
	for(ss = 0; ss < Trace_NumReleaseLocations; ss++) 
		if(!Trace_ReleasePoints[ss].LeftDomain) 
			pnode[ss]->next = NULL;
	
	/* Free local linked list used to connect nodes of Trace_ReleaseList and Trace_ReleasePoints and voxdim0 array */
	free(pnode);
	pnode = NULL;
	free(Trace_ReleasePoints);
	Trace_ReleasePoints = NULL;
	free(voxdim0);
	voxdim0 = NULL;
	free(ptdone);
	ptdone = NULL;
	
	/* Turn back on normal flow if it was turned off */
	if(nff)
		Int_NormalFlow = 1;
	
	
}

void LoadReleasePoints(double *voxdim) {
	/***
	 Reads in number of release points and release point locations. Determines whether 
	 each release point is outside velocity data domain and, if unstructured velocity data
	 is being used, the tetrahedral element the release point is located in. 
	 ***/
	
	char BinFilePath[LONGSTRING];
	FILE *BinFileID;
	int i, index = -1, found = 0, incompat = 0, guess = -1, auxguess = -1;
	ssize_t err;
	
	/* Open binary file for reading */
	sprintf(BinFilePath, "%s%s", Path_Data, Trace_InFile);
	if((BinFileID = fopen(BinFilePath, "rb")) == NULL)
		FatalError("Could not open %s", BinFilePath);
	
	/* Read in the number of release locations */
	if(fread(&Trace_NumReleaseLocations, sizeof(int), 1, BinFileID) < 1)
		FatalError("Could not read Trace_NumReleaseLocations from %s", BinFilePath);
	
	printf("Loading %d release locations and determining if any outside of domain...", Trace_NumReleaseLocations);
	fflush(stdout);
	
	/* Allocate memory to store material point information for the release points */
	if((Trace_ReleasePoints = (LagrangianPoint *)malloc(Trace_NumReleaseLocations * sizeof(LagrangianPoint))) == NULL) 
		FatalError("Malloc failed for Trace_ReleasePoints in function LoadReleasePoints()");
	
	/* Read in coordinates */
	for(i = 0; i < Trace_NumReleaseLocations; i++) 
		if(fread(&Trace_ReleasePoints[i].X[0], sizeof(double), 1, BinFileID) < 1 || 
		   fread(&Trace_ReleasePoints[i].X[1], sizeof(double), 1, BinFileID) < 1 || 
		   fread(&Trace_ReleasePoints[i].X[2], sizeof(double), 1, BinFileID) < 1) 
			FatalError("Could not read complete coordinate list from %s", BinFilePath);
	
	/* Compute spacing between points (assumed release points are spaced uniformly) */
	if(Trace_NumReleaseLocations > 1) {
		*voxdim = dist(Trace_ReleasePoints[0].X, Trace_ReleasePoints[1].X, 3);
		for(i = 1; i < Trace_NumReleaseLocations - 1; i++) 
			*voxdim = fmin(*voxdim, dist(Trace_ReleasePoints[i].X, Trace_ReleasePoints[i+1].X, 3));
		if(*voxdim < TINY) 
			FatalError("Spacing between release points becomes nonpositive");
	}
	else if(Trace_NumReleaseLocations > 0) {
		fprintf(stderr, "\nWarning: Only 1 release location\n");
		fprintf(stderr, "Enter desired spacing:");
		err = fscanf(stdin, "%lf", voxdim); 
	} 
	else
		FatalError("Number of release locations cannot be zero");
	
	/* Close file */
	fclose(BinFileID);
	
	
	for(i = 0; i < Trace_NumReleaseLocations; i++) {
		
		/* Check if release point outside bounding box for data */
		Trace_ReleasePoints[i].LeftDomain = TestOutsideDomain(Trace_ReleasePoints[i].X);
		
		/* If point in bounding box and velocity data is unstrutured, find element that release point is in */
		if(Data_MeshType == UNSTRUCTURED && !Trace_ReleasePoints[i].LeftDomain) { 
			if(!found) {
				/* First point located using global search, and result serves as starting point for further (local) searching */
				index = Get_Element_Global_Search(Trace_ReleasePoints[i].X);
				if(index < 0) {
					Trace_ReleasePoints[i].LeftDomain = 1;
					Trace_ReleasePoints[i].ElementIndex = -1;
					if(Trace_CETCompute && Trace_CETAuxillaryMesh)
						Trace_ReleasePoints[i].AuxElementIndex = -1;
				}
				else { /* Global search succeeded */
					if(Trace_CETCompute && Trace_CETAuxillaryMesh) {
						Trace_ReleasePoints[i].AuxElementIndex = Get_Element_Global_Search_Aux(Trace_ReleasePoints[i].X);
						if(Trace_ReleasePoints[i].AuxElementIndex < 0) { /* Inside vel mesh, but outside aux mesh */
							Trace_ReleasePoints[i].LeftDomain = 1;
							Trace_ReleasePoints[i].ElementIndex = -1;
							incompat++;
						}
						else {
							Trace_ReleasePoints[i].ElementIndex = index;
							Trace_ReleasePoints[i].LeftDomain = 0;
							guess = index;
							auxguess = Trace_ReleasePoints[i].AuxElementIndex;
							found++;
						}
					}
					else { /* Not using aux mesh */
						Trace_ReleasePoints[i].ElementIndex = index;
						Trace_ReleasePoints[i].LeftDomain = 0;
						guess = index;
						found++;
					}
				}
			}
			else { /* A point has been located */ 
				/* Find location of point using local search */
				index = Get_Element_Local_Search(Trace_ReleasePoints[i].X, guess);
				if(index < 0) {
					/* If local seeach fails, always double-check with a global search */
					index = Get_Element_Global_Search(Trace_ReleasePoints[i].X);
					if(index < 0) {
						Trace_ReleasePoints[i].LeftDomain = 1;
						Trace_ReleasePoints[i].ElementIndex = -1;
						if(Trace_CETCompute && Trace_CETAuxillaryMesh) 
							Trace_ReleasePoints[i].AuxElementIndex = -1;
					}
					else {
						if(Trace_CETCompute && Trace_CETAuxillaryMesh) {
							Trace_ReleasePoints[i].AuxElementIndex = Get_Element_Global_Search_Aux(Trace_ReleasePoints[i].X);
							if(Trace_ReleasePoints[i].AuxElementIndex < 0) {
								Trace_ReleasePoints[i].LeftDomain = 1;
								Trace_ReleasePoints[i].ElementIndex = -1;
								incompat++;
							}
							else {
								Trace_ReleasePoints[i].ElementIndex = index;
								Trace_ReleasePoints[i].LeftDomain = 0;
								guess = index;
								auxguess = Trace_ReleasePoints[i].AuxElementIndex;
								found++;
							}
						}
						else { /* Not using aux mesh */
							Trace_ReleasePoints[i].ElementIndex = index;
							Trace_ReleasePoints[i].LeftDomain = 0;
							guess = index;
							found++;
						}
					}
				}
				else { /* Local search succeeded */
					if(Trace_CETCompute && Trace_CETAuxillaryMesh) {
						Trace_ReleasePoints[i].AuxElementIndex = Get_Element_Local_Search_Aux(Trace_ReleasePoints[i].X, auxguess);
						if(Trace_ReleasePoints[i].AuxElementIndex < 0) {
							/* Check with global search */
							Trace_ReleasePoints[i].AuxElementIndex = Get_Element_Global_Search_Aux(Trace_ReleasePoints[i].X);
							if(Trace_ReleasePoints[i].AuxElementIndex < 0) {
								printf("Local search succeeded in velocity mesh, but global search failed in auxillary mesh.\n");
								Trace_ReleasePoints[i].LeftDomain = 1;
								Trace_ReleasePoints[i].ElementIndex = -1;
							}
							else {
								Trace_ReleasePoints[i].ElementIndex = index;
								Trace_ReleasePoints[i].LeftDomain = 0;
								guess = index;
								auxguess = Trace_ReleasePoints[i].AuxElementIndex;
								found++;
							}
						}
						else { /* Local Auxillary search succeeded */
							Trace_ReleasePoints[i].ElementIndex = index;
							Trace_ReleasePoints[i].LeftDomain = 0;
							guess = index;
							auxguess = Trace_ReleasePoints[i].AuxElementIndex;
							found++;
						}
					}
					else { /* Not using aux mesh */
						Trace_ReleasePoints[i].ElementIndex = index;
						Trace_ReleasePoints[i].LeftDomain = 0;
						guess = index;
						found++;
					}
				} 
			}
		}
		else if(!Trace_ReleasePoints[i].LeftDomain) 
			found++; 
	}
	
	printf("OK!\n");
	fflush(stdout);
	
	/* Rewrite release point file with only points inside domain if needed */
	if(found < Trace_NumReleaseLocations) {
		printf("%d located outside of velocity domain, ", Trace_NumReleaseLocations - found);
		if(incompat)
			printf("%d in interior of velocity mesh, but exterior of auxillary mesh\n", incompat);
		printf("Rewriting %s to contain interior points only.\n", Trace_InFile);
		fflush(stdout);
		/* Open file to (re)write data */
		if((BinFileID = fopen(BinFilePath, "wb")) == NULL) 
			FatalError("Could not open file %s", BinFilePath);
		/* Write out the number of release locations */
		if(fwrite(&found, sizeof(int), 1, BinFileID) < 1) 
			FatalError("Could not write number of points to file %s", BinFilePath);
		/* Write out coordinates */
		for(i = 0; i < Trace_NumReleaseLocations; i++) 
			if(!Trace_ReleasePoints[i].LeftDomain) 
				if(fwrite(&Trace_ReleasePoints[i].X[0], sizeof(double), 1, BinFileID) < 1 ||
				   fwrite(&Trace_ReleasePoints[i].X[1], sizeof(double), 1, BinFileID) < 1 ||
				   fwrite(&Trace_ReleasePoints[i].X[2], sizeof(double), 1, BinFileID) < 1) 
					FatalError("Could not write complete coordinate list to file %s", BinFilePath);
		/* Close file */
		fclose(BinFileID);   
	}
	
}

void ReadInTraceMesh(void) {
	
	int ii, ff, index = -1, found = 0, foundGS = 0, guess = -1, auxguess = -1, percentage = 0, numfiles = 1;
	FILE *TraceIC_BinFileID, *Trace_InFileID;
	char Trace_InFilePath[LONGSTRING], buf[LONGSTRING], Trace_ICFile[LONGSTRING];
	LagrangianPoint *TraceIC_Array;
	ssize_t err;
	
	if(Trace_InFileFormat != INITIALIZED) {  
		if(Trace_MultipleInFiles)
			numfiles = Trace_NumLaunchTimes;
		for(ff = 0; ff < numfiles; ff++) {
			if(Trace_InFileFormat == BINARY_LIST) {
				/* Open input file */
				if(Trace_MultipleInFiles)
					sprintf(Trace_InFilePath, "%s%s.%d.bin", Path_Data, Trace_InFile, ff);
				else
					sprintf(Trace_InFilePath, "%s%s", Path_Data, Trace_InFile);
				if((Trace_InFileID = fopen(Trace_InFilePath, "rb")) == NULL)  
					FatalError("Could not open file %s", Trace_InFilePath);
				
				/* First line should contain number of nodal data points to read in */
				if(fread(&Trace_NumTracers, sizeof(int), 1, Trace_InFileID) < 1)
					FatalError("Could not read Trace_NumTracers from %s", Trace_InFilePath);
				printf("\nReading, initializing, and writing data for %d tracers...", Trace_NumTracers);
				fflush(stdout);
				
				/* Allocate memory */
				if((TraceIC_Array = (LagrangianPoint *)malloc(Trace_NumTracers * sizeof(LagrangianPoint))) == NULL)
					FatalError("Malloc failed for TraceIC_Array in ReadInTraceMesh()");	
				
				/* Read in nodal coordinates */
				for(ii = 0; ii < Trace_NumTracers; ii++) { 
					if((int)fread(&TraceIC_Array[ii].X, sizeof(double), Dimensions,  Trace_InFileID) < Dimensions)
						FatalError("Could not read TraceIC_Array[%d].X from %s", ii, Trace_InFilePath);	
					if(Dimensions == 2) 
						TraceIC_Array[ii].X[2] = 0;
				}
				fclose(Trace_InFileID);
			}
			else { 
				if(Trace_InFileFormat == ASCII_LIST) { /* ASCII list of coordinates */
					/* Open input file */
					if(Trace_MultipleInFiles)
						sprintf(Trace_InFilePath, "%s%s.%d.dat", Path_Data, Trace_InFile, ff);
					else 
						sprintf(Trace_InFilePath, "%s%s", Path_Data, Trace_InFile);
					if((Trace_InFileID = fopen(Trace_InFilePath, "r")) == NULL)  
						FatalError("Could not open file %s", Trace_InFilePath);
					/* First line should contain number of nodal data points to read in */
					err = fscanf(Trace_InFileID, "%d\n", &Trace_NumTracers);
					printf("\nReading, initializing, and writing data for %d tracers...", Trace_NumTracers);
					fflush(stdout);
					
					/* Allocate memory */
					if((TraceIC_Array = (LagrangianPoint *)malloc(Trace_NumTracers * sizeof(LagrangianPoint))) == NULL)
						FatalError("Malloc failed for TraceIC_Array in ReadInTraceMesh()");	
					
					/* Read in nodal coordinates */
					for(ii = 0; ii < Trace_NumTracers; ii++) 
						if(Dimensions == 3)
							err = fscanf(Trace_InFileID, "%lf %lf %lf\n", &TraceIC_Array[ii].X[0], &TraceIC_Array[ii].X[1], &TraceIC_Array[ii].X[2]);
						else { /* Dimensions == 2 */
							err = fscanf(Trace_InFileID, "%lf %lf\n", &TraceIC_Array[ii].X[0], &TraceIC_Array[ii].X[1]);
							TraceIC_Array[ii].X[2] = 0;
						}
					
					fclose(Trace_InFileID);
				}
				else if(Trace_InFileFormat == VTK_UNSTRUCTURED || Trace_InFileFormat == VTK_POLYDATA) { /* Legacy vtk ASCII format (polydata or unstructured grid) */
					if(Dimensions == 2) {
						TraceIC_Array = NULL;
						FatalError("The specified Trace_InFileFormat is not currently supported for 2D analysis");
					}
					else {
						/* Open input file */
						if(Trace_MultipleInFiles)
							sprintf(Trace_InFilePath, "%s%s.%d.vtk", Path_Data, Trace_InFile, ff);
						else 
							sprintf(Trace_InFilePath, "%s%s", Path_Data, Trace_InFile);
						if((Trace_InFileID = fopen(Trace_InFilePath, "r")) == NULL)  
							FatalError("Could not open file %s", Trace_InFilePath);
						
						/* Skip first 4 header lines */
						for(ii = 0; ii < 4; ii++)
							if (fgets(buf, LONGSTRING, Trace_InFileID) == NULL);
								FatalError("Skipping of first four header lines failed.");
						/* Get number of nodal data points to read in */
						err = fscanf(Trace_InFileID, "%*s %d %*s\n", &Trace_NumTracers); 
						printf("\nReading, initializing, and writing data for %d tracers...", Trace_NumTracers); 
						fflush(stdout);
						
						/* Allocate memory */
						if((TraceIC_Array = (LagrangianPoint *)malloc(Trace_NumTracers * sizeof(LagrangianPoint))) == NULL) 
							FatalError("Malloc failed for TraceIC_Array in ReadInTraceMesh()");	
						
						/* Read in nodal coordinates */
						for(ii = 0; ii < Trace_NumTracers; ii++)
							err = fscanf(Trace_InFileID, "%lf %lf %lf", &TraceIC_Array[ii].X[0], &TraceIC_Array[ii].X[1], &TraceIC_Array[ii].X[2]);
						
						fclose(Trace_InFileID);
					}
				}
				else {
					TraceIC_Array = NULL;
					FatalError("Unsupported value for Trace_InFileFormat in ReadInTraceMesh()");	
				}
			}
			
			/* Initialize Data */
			percentage = 0;
			found = 0;
			foundGS = 0;
			for(ii = 0; ii < Trace_NumTracers; ii++) {
				if(100 * ii / Trace_NumTracers > percentage) {
					printf("  %d%%", percentage);
					fflush(stdout);
					percentage = percentage + 10;
				}
				
				if(Trace_APCompute) 
					TraceIC_Array[ii].Scalar = 0.0;
				
				/* Check if point outside bounding box for data */
				TraceIC_Array[ii].LeftDomain = TestOutsideDomain(TraceIC_Array[ii].X);
				
				/* If point in bounding box and velocity data is unstrutured, find element that point is in */
				if(Data_MeshType == UNSTRUCTURED && !TraceIC_Array[ii].LeftDomain) { /* Unstructured mesh */
					if(!found) {
						/* First point located using global search, and result serves as starting point for further (local) searching */
						index = Get_Element_Global_Search(TraceIC_Array[ii].X);
						if(index < 0) {
							TraceIC_Array[ii].LeftDomain = 1;
							TraceIC_Array[ii].ElementIndex = -1;
							if(Trace_CETCompute && Trace_CETAuxillaryMesh)
								TraceIC_Array[ii].AuxElementIndex = -1;
						}
						else { /* Global search succeeded */
							if(Trace_CETCompute && Trace_CETAuxillaryMesh) {
								TraceIC_Array[ii].AuxElementIndex = Get_Element_Global_Search_Aux(TraceIC_Array[ii].X);
								if(TraceIC_Array[ii].AuxElementIndex < 0) {
									TraceIC_Array[ii].LeftDomain = 1;
									TraceIC_Array[ii].ElementIndex = -1;
								}
								else {
									TraceIC_Array[ii].ElementIndex = index;
									guess = index;
									auxguess = TraceIC_Array[ii].AuxElementIndex;
									found++;
								}
							}
							else { /* Not using aux mesh */
								TraceIC_Array[ii].ElementIndex = index;
								guess = index;
								found++; 
							}
						}
					}
					else { /* A point has been located, use local search */ 
						/* Find location of point using local search */
						index = Get_Element_Local_Search(TraceIC_Array[ii].X, guess);
						if(index < 0) { /* Local search failed */
							/* If local seeach fails, double-check with a global search if requested */
							if(LocalSearchChecking) {
								index = Get_Element_Global_Search(TraceIC_Array[ii].X);
								if(index < 0) { 
									TraceIC_Array[ii].LeftDomain = 1;
									TraceIC_Array[ii].ElementIndex = -1;
									if(Trace_CETCompute && Trace_CETAuxillaryMesh)
										TraceIC_Array[ii].AuxElementIndex = -1;
								}
								else { /* Global search succeeded */
									if(Trace_CETCompute && Trace_CETAuxillaryMesh) {
										TraceIC_Array[ii].AuxElementIndex = Get_Element_Global_Search_Aux(TraceIC_Array[ii].X);
										if(TraceIC_Array[ii].AuxElementIndex < 0) {
											TraceIC_Array[ii].LeftDomain = 1;
											TraceIC_Array[ii].ElementIndex = -1;
										}
										else { /* Global Auxillary search succeeded */
											TraceIC_Array[ii].ElementIndex = index;
											TraceIC_Array[ii].LeftDomain = 0;
											guess = index;
											auxguess = TraceIC_Array[ii].AuxElementIndex;
											found++;
											foundGS++;
										}
									}
									else { /* Not using aux mesh */
										TraceIC_Array[ii].ElementIndex = index;
										TraceIC_Array[ii].LeftDomain = 0;
										guess = index;
										found++;
										foundGS++;
									}
								}
							}
							else {
								TraceIC_Array[ii].LeftDomain = 1;
								TraceIC_Array[ii].ElementIndex = -1;
								if(Trace_CETCompute && Trace_CETAuxillaryMesh)
									TraceIC_Array[ii].AuxElementIndex = -1;
							}
						}
						else { /* Local search succeeded */
							if(Trace_CETCompute && Trace_CETAuxillaryMesh) {
								TraceIC_Array[ii].AuxElementIndex = Get_Element_Local_Search_Aux(TraceIC_Array[ii].X, auxguess);
								if(TraceIC_Array[ii].AuxElementIndex < 0) { /* Local Auxillary search failed */
									if(LocalSearchChecking) {
										TraceIC_Array[ii].AuxElementIndex = Get_Element_Global_Search_Aux(TraceIC_Array[ii].X);
										if(TraceIC_Array[ii].AuxElementIndex < 0) {
											TraceIC_Array[ii].LeftDomain = 1;
											TraceIC_Array[ii].ElementIndex = -1;
										}
										else { /* Global Auxillary search succeeded */
											TraceIC_Array[ii].ElementIndex = index;
											TraceIC_Array[ii].LeftDomain = 0;
											guess = index;
											auxguess = TraceIC_Array[ii].AuxElementIndex;
											found++;
											foundGS++;
										}
									}
									else {
										TraceIC_Array[ii].LeftDomain = 1;
										TraceIC_Array[ii].ElementIndex = -1;
										TraceIC_Array[ii].AuxElementIndex = -1;
									}
								}
								else { /* Local Auxillary search succeeded */ 
									TraceIC_Array[ii].ElementIndex = index;
									TraceIC_Array[ii].LeftDomain = 0;
									guess = index;
									auxguess = TraceIC_Array[ii].AuxElementIndex;
									found++;
								}
							}
							else { /* Not using aux mesh */
								TraceIC_Array[ii].ElementIndex = index;
								TraceIC_Array[ii].AuxElementIndex = TraceIC_Array[ii].ElementIndex;
								TraceIC_Array[ii].LeftDomain = 0;
								guess = index;
								auxguess = TraceIC_Array[ii].AuxElementIndex;
								found++;
							}
						}
					} /* A point has been located, use local search */
				} /* if(Data_MeshType == UNSTRUCTURED && !TraceIC_Array[ii].LeftDomain) */
				
				if(TraceIC_Array[ii].LeftDomain)
					TraceIC_Array[ii].LeftDomainTime = -1.0; /* Flag for points initially outside of domain  */
				else
					TraceIC_Array[ii].LeftDomainTime = 0.0; /* Initialize to 0 */
				
			}
			
			if(LocalSearchChecking) 
				printf("  %d of %d located in domain (%d caught by local search checking)\n", found, Trace_NumTracers, foundGS);
			else
				printf("  %d of %d located in domain\n", found, Trace_NumTracers);       
			
			/* Store intitial condition data in binary file to be used to generate new slides */
			if(Trace_MultipleInFiles)
				sprintf(Trace_ICFile, "%s.%d", TraceIC_BinFilePath, ff);
			else
				sprintf(Trace_ICFile, "%s", TraceIC_BinFilePath);	
			if((TraceIC_BinFileID = fopen(Trace_ICFile, "wb")) == NULL) 
				FatalError("Could not open file %s", Trace_ICFile);
			
			/* Write Number of Tracers to file */
			if(fwrite(&Trace_NumTracers, sizeof(int), 1, TraceIC_BinFileID) < 1)
				FatalError("Could not write Trace_NumTracers to %s", TraceIC_BinFilePath);
			
			/* Write TraceIC_Array information to binary file */
			if((int)fwrite(TraceIC_Array, sizeof(LagrangianPoint), Trace_NumTracers, TraceIC_BinFileID) < Trace_NumTracers)
				FatalError("Could not completely write TraceIC_Array to %s", TraceIC_BinFilePath);
			
			/* Free TraceIC_Array */
			free(TraceIC_Array);
			TraceIC_Array = NULL;
			
			fclose(TraceIC_BinFileID);
			
		} /* for(ff = 0; ff < numfiles; ff++) */
	} /* if(Trace_InFileFormat != INITIALIZED) */
	
	/* Initialize launch time information */
	if((Trace_Launches = (Launch *)malloc(Trace_NumLaunchTimes * sizeof(Launch))) == NULL)
		FatalError("Malloc failed for Trace_Launches");
	
	for(ii = 0; ii < Trace_NumLaunchTimes; ii++) {
		Trace_Launches[ii].StartTime = Output_TStart + Int_TimeDirection * ii * Trace_LaunchTimeSpacing;
		Trace_Launches[ii].StopTime = (Int_TimeDirection > 0) ? fmin(Trace_Launches[ii].StartTime + Trace_IntTLength, Output_TEnd) : fmax(Trace_Launches[ii].StartTime - Trace_IntTLength, Output_TEnd);
		Trace_Launches[ii].Status = UNLAUNCHED;
	}
	
	printf("OK!\n");
	fflush(stdout);
	
}

void GenerateTracerMesh(void) {
	
	int ii, jj, kk, seed = -1, auxseed = -1, index = -1, count = 0, found = 0, foundGS = 0, guess = -1, auxguess = -1;
	double Xseed[3];
	FILE *TraceIC_BinFileID;
	LagrangianPoint *TraceIC_Array;
	
	printf("\nInitializing tracer mesh...\n");
	fflush(stdout);
	
	Trace_NumTracers = Trace_CartMesh.XRes * Trace_CartMesh.YRes * Trace_CartMesh.ZRes;
	
	/* Allocate memory for TraceIC_Array */
	if((TraceIC_Array = (LagrangianPoint *)malloc(Trace_NumTracers * sizeof(LagrangianPoint))) == NULL) 
		FatalError("Malloc failed for TraceIC_Array in GenerateTracerMesh()");
	
	/* Try to "seed" local element search by doing global search for center point */
	if(Data_MeshType == UNSTRUCTURED) {
		Xseed[0] = Trace_CartMesh.XMin + ((int)(Trace_CartMesh.XRes/2)) * Trace_CartMesh.XDelta;
		Xseed[1] = Trace_CartMesh.YMin + ((int)(Trace_CartMesh.YRes/2)) * Trace_CartMesh.YDelta;
		Xseed[2] = Trace_CartMesh.ZMin + ((int)(Trace_CartMesh.ZRes/2)) * Trace_CartMesh.ZDelta;
		seed = Get_Element_Global_Search(Xseed);
		
		if(Trace_CETCompute && Trace_CETAuxillaryMesh)
			auxseed = Get_Element_Global_Search_Aux(Xseed);
		if(seed >= 0) {
			if(Trace_CETCompute && Trace_CETAuxillaryMesh) {
				if(auxseed >= 0) {     
					printf("  Using center element seed.\n");
					fflush(stdout);
					guess  = seed;
					auxguess = auxseed;
				}
			}
			else {
				printf("  Using center element seed.\n");
				fflush(stdout);
				guess  = seed;
			}
		}
	}
	
	for(kk = 0; kk < Trace_CartMesh.ZRes; kk++) {
		for(jj = 0; jj < Trace_CartMesh.YRes; jj++) {
			for(ii = 0; ii < Trace_CartMesh.XRes; ii++) {
				
				/* Set coordinates of point */
				TraceIC_Array[count].X[0] = Trace_CartMesh.XMin + ii*Trace_CartMesh.XDelta;
				TraceIC_Array[count].X[1] = Trace_CartMesh.YMin + jj*Trace_CartMesh.YDelta;
				TraceIC_Array[count].X[2] = Trace_CartMesh.ZMin + kk*Trace_CartMesh.ZDelta;
				
				if(Trace_APCompute) 
					TraceIC_Array[count].Scalar = 0.0;
				
				/* See if coordinates outside of data bounding box */
				TraceIC_Array[count].LeftDomain = TestOutsideDomain(TraceIC_Array[count].X);
				
				/* If inside bounding box and using a tetmesh, find element the point is in */
				if(!TraceIC_Array[count].LeftDomain) {
					if(Data_MeshType == UNSTRUCTURED) { /* Unstructured mesh */
						if((seed < 0) || ((Trace_CETCompute && Trace_CETAuxillaryMesh)? auxseed < 0 : seed < 0)) { /* Use global search until a point has been located */
							seed = Get_Element_Global_Search(TraceIC_Array[count].X);
							if(Trace_CETCompute && Trace_CETAuxillaryMesh)
								auxseed = Get_Element_Global_Search_Aux(TraceIC_Array[count].X);
							if((seed < 0) || ((Trace_CETCompute && Trace_CETAuxillaryMesh)? auxseed < 0 : seed < 0)) { /* Point outside mesh boundary */
								TraceIC_Array[count].LeftDomain = 1;
								TraceIC_Array[count].LeftDomainTime = -1.0;
								TraceIC_Array[count].ElementIndex = -1;
								TraceIC_Array[count].AuxElementIndex = -1;
							}
							else {
								printf("  Using first found element seed.\n");
								TraceIC_Array[count].ElementIndex = seed;
								if(Trace_CETCompute && Trace_CETAuxillaryMesh) {
									TraceIC_Array[count].AuxElementIndex = auxseed;
									auxguess = auxseed;
								}
								found++;
								guess  = seed;
							}
						}
						else /* A point has been located; use local searching */ {
							/* Find location of point using local search */
							index = Get_Element_Local_Search(TraceIC_Array[count].X, guess);
							if(index < 0)  /* Local search from guess failed */ {
								/* Try local search from seed */
								index = Get_Element_Local_Search(TraceIC_Array[count].X, seed);
								if(index < 0) /* Local search from seed also not working */ {
									/* Use a global search if requested */
									if(LocalSearchChecking) {
										index = Get_Element_Global_Search(TraceIC_Array[count].X);
										if(index < 0) {  
											TraceIC_Array[count].LeftDomain = 1;
											TraceIC_Array[count].LeftDomainTime = -1.0;
											TraceIC_Array[count].ElementIndex = -1;
											if(Trace_CETCompute && Trace_CETAuxillaryMesh)
												TraceIC_Array[count].AuxElementIndex = -1;
										}
										else { /* Global search succeeded */
											if(Trace_CETCompute && Trace_CETAuxillaryMesh) {
												TraceIC_Array[count].AuxElementIndex = Get_Element_Global_Search_Aux(TraceIC_Array[count].X);
												if(TraceIC_Array[count].AuxElementIndex < 0) {
													TraceIC_Array[count].LeftDomain = 1;
													TraceIC_Array[count].LeftDomainTime = -1.0;
													TraceIC_Array[count].ElementIndex = -1;
												}
												else { /* Global Auxillary search succeeded */
													TraceIC_Array[count].ElementIndex = index;
													TraceIC_Array[count].LeftDomain = 0;
													TraceIC_Array[count].LeftDomainTime = 0.0;
													guess = index;
													auxguess = TraceIC_Array[count].AuxElementIndex;
													found++;
													foundGS++;
												}
											}
											else { /* Not using aux mesh */
												TraceIC_Array[count].ElementIndex = index;
												TraceIC_Array[count].LeftDomain = 0;
												TraceIC_Array[count].LeftDomainTime = 0.0;
												guess = index;
												found++;
												foundGS++;
											}
										}
									}
									else { /* Local search failed and not attempting global search */
										TraceIC_Array[count].LeftDomain = 1;
										TraceIC_Array[count].LeftDomainTime = -1.0;
										TraceIC_Array[count].ElementIndex = -1;
										if(Trace_CETCompute && Trace_CETAuxillaryMesh)
											TraceIC_Array[count].AuxElementIndex = -1;
									}
								}
								else { /* Local search from seed succeeded */
									if(Trace_CETCompute && Trace_CETAuxillaryMesh) {
										TraceIC_Array[count].AuxElementIndex = Get_Element_Local_Search_Aux(TraceIC_Array[count].X, auxseed);
										if(TraceIC_Array[count].AuxElementIndex < 0) { /* Local Auxillary search from auxseed failed */
											if(LocalSearchChecking) {
												TraceIC_Array[count].AuxElementIndex = Get_Element_Global_Search_Aux(TraceIC_Array[count].X);
												if(TraceIC_Array[count].AuxElementIndex < 0) {
													/* Global search failed */
													TraceIC_Array[count].LeftDomain = 1;
													TraceIC_Array[count].LeftDomainTime = -1.0;
													TraceIC_Array[count].ElementIndex = -1;
												}
												else { /* Global Auxillary search succeeded */
													TraceIC_Array[count].ElementIndex = index;
													TraceIC_Array[count].LeftDomain = 0;
													TraceIC_Array[count].LeftDomainTime = 0.0;
													guess = index;
													auxguess = TraceIC_Array[count].AuxElementIndex;
													found++;
													foundGS++;
												}
											}
											else { /* No local search checking */
												TraceIC_Array[count].LeftDomain = 1;
												TraceIC_Array[count].LeftDomainTime = -1.0;
												TraceIC_Array[count].ElementIndex = -1;
												TraceIC_Array[count].AuxElementIndex = -1;
											}
										}
										else { /* Local Auxillary search from auxseed succeeded */ 
											TraceIC_Array[count].ElementIndex = index;
											TraceIC_Array[count].LeftDomain = 0;
											TraceIC_Array[count].LeftDomainTime = 0.0;
											guess = index;
											auxguess = TraceIC_Array[count].AuxElementIndex;
											found++;
										}  
									} /* Not using aux mesh */
									else {
										TraceIC_Array[count].ElementIndex = index;
										TraceIC_Array[count].LeftDomain = 0;
										TraceIC_Array[count].LeftDomainTime = 0.0;
										guess = index;
										found++;
									}
								}
							}
							else { /* Local search from guess succeeded */
								if(Trace_CETCompute && Trace_CETAuxillaryMesh) {
									TraceIC_Array[count].AuxElementIndex = Get_Element_Local_Search_Aux(TraceIC_Array[count].X, auxguess);
									if(TraceIC_Array[count].AuxElementIndex < 0) { /* Local Auxillary search from auxguess failed */
										if(LocalSearchChecking) {
											TraceIC_Array[count].AuxElementIndex = Get_Element_Global_Search_Aux(TraceIC_Array[count].X);
											if(TraceIC_Array[count].AuxElementIndex < 0) {
												TraceIC_Array[count].LeftDomain = 1;
												TraceIC_Array[count].LeftDomainTime = -1.0;
												TraceIC_Array[count].ElementIndex = -1;
											}
											else { /* Global Auxillary search succeeded */
												TraceIC_Array[count].ElementIndex = index;
												TraceIC_Array[count].LeftDomain = 0;
												TraceIC_Array[count].LeftDomainTime = 0.0;
												guess = index;
												auxguess = TraceIC_Array[count].AuxElementIndex;
												found++;
												foundGS++;
											}
										}
										else { /* No local search checking */
											TraceIC_Array[count].LeftDomain = 1;
											TraceIC_Array[count].LeftDomainTime = -1.0;
											TraceIC_Array[count].ElementIndex = -1;
											TraceIC_Array[count].AuxElementIndex = -1;
										}
									}
									else { /* Local Auxillary search from auxguess succeeded */ 
										TraceIC_Array[count].ElementIndex = index;
										TraceIC_Array[count].LeftDomain = 0;
										TraceIC_Array[count].LeftDomainTime = 0.0;
										guess = index;
										auxguess = TraceIC_Array[count].AuxElementIndex;
										found++;
									}
								} 
								else { /* Not using aux mesh */
									TraceIC_Array[count].ElementIndex = index;
									TraceIC_Array[count].LeftDomain = 0;
									TraceIC_Array[count].LeftDomainTime = 0.0;
									guess = index;
									found++;
								}
							}
						}
					}
					else
						found++;
				}
				count++;
			}
		}
	}
	
	
	if(LocalSearchChecking) 
		printf("  %d of %d located in domain (%d caught by local search checking)\n", found, Trace_CartMesh.XRes * Trace_CartMesh.YRes * Trace_CartMesh.ZRes, foundGS);
	else
		printf("  %d of %d located in domain\n", found, Trace_CartMesh.XRes * Trace_CartMesh.YRes * Trace_CartMesh.ZRes);  
	fflush(stdout);
	
	/* Store Tracer intitial conditions in binary file to be used to generate new slides */
	if((TraceIC_BinFileID = fopen(TraceIC_BinFilePath, "wb")) == NULL) 
		FatalError("Could not open file %s", TraceIC_BinFilePath);
	
	/* 
	 if(fwrite(&found, sizeof(int), 1, TraceIC_BinFileID) < 1) 
	 FatalError("Could not write tracer count to %s", TraceIC_BinFilePath);
	 */
	Trace_NumTracers = Trace_CartMesh.XRes * Trace_CartMesh.YRes * Trace_CartMesh.ZRes ;
	if(fwrite(&Trace_NumTracers, sizeof(int), 1, TraceIC_BinFileID) < 1) 
		FatalError("Could not write Trace_NumTracers to %s", TraceIC_BinFilePath);
	
	/* Write TraceIC_Array information to binary file for tracers initially inside domain */
	for(count = 0; count < Trace_NumTracers; count++) 
    /* if(!TraceIC_Array[count].LeftDomain) */
		if(fwrite(&TraceIC_Array[count], sizeof(LagrangianPoint), 1, TraceIC_BinFileID) < 1) 
			FatalError("Could not write TraceIC_Array[%d] to %s", count, TraceIC_BinFilePath);
	
	/* Close binary file */
	fclose(TraceIC_BinFileID);
	
	/* Free TraceIC_Array */
	free(TraceIC_Array);
	TraceIC_Array = NULL;
	
	if((Trace_Launches = (Launch *)malloc(Trace_NumLaunchTimes * sizeof(Launch))) == NULL)
		FatalError("Malloc failed for Trace_Launches");
	
	for(ii = 0; ii < Trace_NumLaunchTimes; ii++) {
		Trace_Launches[ii].StartTime = Output_TStart + Int_TimeDirection * ii * Trace_LaunchTimeSpacing;
		Trace_Launches[ii].StopTime = (Int_TimeDirection > 0) ? fmin(Trace_Launches[ii].StartTime + Trace_IntTLength, Output_TEnd) : fmax(Trace_Launches[ii].StartTime - Trace_IntTLength, Output_TEnd);
		Trace_Launches[ii].Status = UNLAUNCHED;
	}
	
	printf("OK!\n");
	fflush(stdout);
}

void CreateNewTraceLaunch(int ss) { 
	/*** 
	 Copies contents of TraceIC_BinFile to TraceN_BinFile for release ss 
	 ***/
	
	char Trace_BinFilePath[LONGSTRING], Trace_ICFile[LONGSTRING];
	FILE *TraceIC_BinFileID, *Trace_BinFileID;
	LagrangianPoint *TraceIC_Array;
	
	/* Open input file */
	if(!Trace_GenerateMesh && Trace_MultipleInFiles)
		sprintf(Trace_ICFile, "%s.%d", TraceIC_BinFilePath, ss);
	else
		sprintf(Trace_ICFile, "%s", TraceIC_BinFilePath);
	
	if((TraceIC_BinFileID = fopen(Trace_ICFile, "rb")) == NULL) 
		FatalError("Could not open file %s", Trace_ICFile);
	
	/* Read in number of tracers from file */ 
	if(fread(&Trace_NumTracers, sizeof(int), 1, TraceIC_BinFileID) < 1)
		FatalError("Could not read Trace_NumTracers from %s", TraceIC_BinFilePath);
	
	/* Allocate memory for TraceIC_Array */
	if((TraceIC_Array = (LagrangianPoint *)malloc(Trace_NumTracers * sizeof(LagrangianPoint))) == NULL)  
		FatalError("Malloc failed for TraceIC_Array in CreateNewTraceLaunch()");
	
	/* Read TraceIC_Array information from file */
	if((int)fread(TraceIC_Array, sizeof(LagrangianPoint), Trace_NumTracers, TraceIC_BinFileID) < Trace_NumTracers)
		FatalError("Could not completely read TraceIC_Array data from %s", TraceIC_BinFilePath);
	
	/* Close IC file */
	fclose(TraceIC_BinFileID);
	
	/* Open output file */
	sprintf(Trace_BinFilePath, "%s.%d.bin", TraceN_BinFilePathPrefix, ss);
	if((Trace_BinFileID = fopen(Trace_BinFilePath, "wb")) == NULL) 
		FatalError("Could not open file %s", Trace_BinFilePath);
	
	/* Write Number of Tracers to file */
	if(fwrite(&Trace_NumTracers, sizeof(int), 1, Trace_BinFileID) < 1)
		FatalError("Could not write Trace_NumTracers to %s", Trace_BinFilePath);
	
	/* Write TraceIC_Array information to file */
	if((int)fwrite(TraceIC_Array, sizeof(LagrangianPoint), Trace_NumTracers, Trace_BinFileID) < Trace_NumTracers)
		FatalError("Could not completely write TraceIC_Array data to %s", Trace_BinFilePath);
	
	/* Close output file */
	fclose(Trace_BinFileID);
	
	/* Free TraceIC_Array */
	free(TraceIC_Array);
	TraceIC_Array = NULL;
	
}

LagrangianPoint* ReadInTraceLaunch(int ss, LagrangianPoint *Trace_MeshPt) {
	/* Function reads in data structures for slide specified by ss into Trace_MeshPt */
	
	FILE *Trace_BinFileID1;
	char Trace_BinFilePath1[LONGSTRING];
	
	/* Open file to read binary */
	sprintf(Trace_BinFilePath1, "%s.%d.bin", TraceN_BinFilePathPrefix, ss);
	if((Trace_BinFileID1 = fopen(Trace_BinFilePath1, "rb")) == NULL) 
		FatalError("Could not open file %s", Trace_BinFilePath1);
	
	/* Read in number of tracers from file */ 
	if(fread(&Trace_NumTracers, sizeof(int), 1, Trace_BinFileID1) < 1)
		FatalError("Could not read Trace_NumTracers from %s", Trace_BinFilePath1);
	
	if (TRACER_LOWERLOOP) {
		// Allocate memory for Trace_MeshPt when parallizing according to number of points to be integrated. 
		if((Trace_MeshPt = (LagrangianPoint *)malloc(Trace_NumTracers * sizeof(LagrangianPoint))) == NULL) 
			FatalError("Malloc failed for Trace_MeshPt");
			
		//printf("Memory allocated .... \n\n");
	
	}
	
	
	/* Read Trace_MeshPt information from file */
	if((int)fread(Trace_MeshPt, sizeof(LagrangianPoint), Trace_NumTracers, Trace_BinFileID1) < Trace_NumTracers)
		FatalError("Could not completely read Trace_MeshPt data from %s", Trace_BinFilePath1);
	
	//printf("%f \n", Trace_MeshPt[23].X[0]);
	fclose(Trace_BinFileID1); 
	
	return Trace_MeshPt;
	
}

void WriteOutTraceLaunch(int ss, LagrangianPoint *Trace_MeshPt) {
	/* Writes contents of Trace_MeshPt to binary file for slide ss */
	
	FILE *Trace_BinFileID1;
	char Trace_BinFilePath1[LONGSTRING];
	
	/* Open file to write binary */
	sprintf(Trace_BinFilePath1, "%s.%d.bin", TraceN_BinFilePathPrefix, ss);
	if((Trace_BinFileID1 = fopen(Trace_BinFilePath1, "wb")) == NULL) 
		FatalError("Could not open file %s", Trace_BinFilePath1);
	
	/* Write Number of Tracers to file */
	if(fwrite(&Trace_NumTracers, sizeof(int), 1, Trace_BinFileID1) < 1)
		FatalError("Could not write Trace_NumTracers to %s", Trace_BinFilePath1);
	
	/* Write Trace_MeshPt information to file */
	if((int)fwrite(Trace_MeshPt, sizeof(LagrangianPoint), Trace_NumTracers, Trace_BinFileID1) < Trace_NumTracers)
		FatalError("Could not completely write Trace_MeshPt data to %s", Trace_BinFilePath1);
	
	if (TRACER_LOWERLOOP) {
		
		/* Free Trace_MeshPt for LOWERLOOP optimization only */
		free(Trace_MeshPt);  
		Trace_MeshPt = NULL;
	
	}
	fclose(Trace_BinFileID1);
	
}

void OutputTracers_OMP(void) {

	int i,j,k;
	double X[3], time, ap;
	char Trace_OutFilePath[LONGSTRING], Trace_BinFilePath[LONGSTRING];
	FILE *Trace_BinFileID, *Trace_OutFileID;

	printf("\nOutputting tracer data to file... \n");
	fflush(stdout);

	// Club all releases into one single output file
	
	for (i = 0; i < Output_TRes; i++) {

		// Open output file to write data
		sprintf(Trace_OutFilePath, "%s%s.%d.bin", Path_Output, Trace_OutFilePrefix, i);
		if((Trace_OutFileID = fopen(Trace_OutFilePath , "wb")) == NULL) 
			FatalError("Could not open file %s", Trace_OutFilePath);
		
		// Print time stamp 
		time = Output_TStart + i * Int_TimeDirection * Output_TDelta;
			
		if(fwrite(&time, sizeof(double), 1, Trace_OutFileID) < 1)
			FatalError("Could not write time stamp to %s", Trace_OutFilePath);

		for (j = 0; j < Trace_NumLaunchTimes; j++) {
			// Open file to write only if release has any points for the current output
			if(Trace_NumOutput_OMP[i][j] > 0) {

				sprintf(Trace_BinFilePath, "%s%s_release%d.%d.bin", Path_Data, Trace_OutFilePrefix, j, i);
				//printf("Output temp file path = %s \n",OutFilePath);
				if((Trace_BinFileID = fopen(Trace_BinFilePath, "rb")) == NULL) 
					FatalError("Could not open file %s", Trace_BinFilePath);
					
				for (k = 0; k < Trace_NumOutput_OMP[i][j]; k++) {

					if(fread(X, sizeof(double), 3, Trace_BinFileID) < 3) 
						FatalError("Could not read tracer location from %s", Trace_BinFilePath);
					if(Trace_APCompute)
						if(fread(&ap, sizeof(double), 1, Trace_BinFileID) < 1) 
							FatalError("Could not read activation potential for tracer %d from %s", k, Trace_BinFilePath);
					if(fwrite(X, sizeof(double), 3, Trace_OutFileID) < 3) 
						FatalError("Could not write tracer location to %s", Trace_OutFilePath);
					if(Trace_APCompute)
						if(fwrite(&ap, sizeof(double), 1, Trace_OutFileID) < 1) 
							FatalError("Could not write activation potential for tracer %d from %s", k, Trace_OutFilePath);

				}
				
				fclose(Trace_BinFileID);
				
				if(remove(Trace_BinFilePath))
					fprintf(stderr, "Warning: Could not delete file %s\n", Trace_BinFilePath);
			}
		}
		
		fclose(Trace_OutFileID);
	}
	
	printf("OK!\n");
	fflush(stdout);
}

void OutputTracers(void) {
	
	int i, j;
	double X[3], time, ap;
	char Trace_OutFilePath[LONGSTRING], Trace_BinFilePath[LONGSTRING];
	FILE *Trace_BinFileID, *Trace_OutFileID;
	
	printf("\nOutputting tracer data to file...");
	fflush(stdout);
	
	/* Open bin file containing all output locations */
	if((Trace_BinFileID = fopen(TraceOUT_BinFilePath, "rb")) == NULL) 
		FatalError("Could not open %s", TraceOUT_BinFilePath);
	
	for(i = 0; i < Output_TRes; i++) {
		if(Trace_NumOutput[i] > 0) {
			/* Open bin output file for locations at current output time */
			sprintf(Trace_OutFilePath, "%s%s.%d.bin", Path_Output, Trace_OutFilePrefix, i);
			if((Trace_OutFileID = fopen(Trace_OutFilePath , "wb")) == NULL) 
				FatalError("Could not open file %s", Trace_OutFilePath);
			
			/* Print time stamp */
			time = Output_TStart + i * Int_TimeDirection * Output_TDelta;
			if(fwrite(&time, sizeof(double), 1, Trace_OutFileID) < 1)
				FatalError("Could not write time stamp to %s", Trace_OutFilePath);
			
			/* Read/Write data */
			for(j = 0; j < Trace_NumOutput[i]; j++) {
				if(fread(X, sizeof(double), 3, Trace_BinFileID) < 3) 
					FatalError("Could not read tracer location from %s", Trace_BinFilePath);
				if(Trace_APCompute)
					if(fread(&ap, sizeof(double), 1, Trace_BinFileID) < 1) 
						FatalError("Could not read activation potential for tracer %d from %s", j, Trace_BinFilePath);
				if(fwrite(X, sizeof(double), 3, Trace_OutFileID) < 3) 
					FatalError("Could not write tracer location to %s", Trace_OutFilePath);
				if(Trace_APCompute)
					if(fwrite(&ap, sizeof(double), 1, Trace_OutFileID) < 1) 
						FatalError("Could not write activation potential for tracer %d from %s", j, Trace_OutFilePath);
			}
			/* Close output file */
			fclose(Trace_OutFileID);  
		}
	}
	
	/* Close and delete file */
	if(fclose(Trace_BinFileID))
		fprintf(stderr, "Warning: Could not close file %s\n", Trace_BinFilePath);
	if(remove(TraceOUT_BinFilePath))
		fprintf(stderr, "Warning: Could not delete file %s\n", TraceOUT_BinFilePath);
	
	printf("OK!\n");
	fflush(stdout);
	
}


void FreeTracerData(void) {
	
	int i;
	ReleaseLocation *n, *p;
	
	free(Trace_NumOutput);
	Trace_NumOutput = NULL;
	
	for (i = 0; i < Output_TRes; i++) {
		free(Trace_NumOutput_OMP[i]);
		Trace_NumOutput_OMP[i] = NULL;
		
	}
	free(Trace_NumOutput_OMP);
	Trace_NumOutput_OMP = NULL;

	
	if(Trace_ReleaseStrategy == STAGGERED) {
		free(Trace_ReleasePoints);
		for(i = 0; i < Trace_NumReleaseLocations; i++) {
			n = Trace_ReleaseList[i];
			while(n != NULL) {
				p = n->next;
				free(n);
				n = p;
			}
		}
		free(Trace_ReleaseList);
		Trace_ReleaseList = NULL;
	}
	else {
		free(Trace_Launches);
		Trace_Launches = NULL;
	}
	if(Trace_CETCompute) {
		free(Trace_CETArray);
		Trace_CETArray = NULL;
	}
}
