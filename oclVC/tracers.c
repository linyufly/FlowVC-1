#include <stdlib.h>
#include <stdio.h>
#include <CL/cl.h>

# include "globals.h"
# include "tracers.h"
# include "structs.h"
# include "clean.h"
# include "GlobalSearch.h"
# include "velocity.h"
# include "fileoutput.h"
# include "mymath.h"



void GenerateStaggeredRelease(void) {

	double M, N, ti, tii, a, b, c, d, voxdim, *voxdim0, usbar, uebar, uibar, us[3], ue[3], **u0;
	int ss, ii, df, *ptdone;
	ReleaseLocation *nnode, **pnode;
	ReleaseLocation *rnode;
	
	LoadReleasePoints(&voxdim);  // Sets voxdim, Trace_ReleasePoints array and Trace_NumReleaseLocations 

	printf("\nDetermining release times...\n\n");
	fflush(stdout);
	
	// Allocate memory for global linked list that will store release information 
	if((Trace_ReleaseList = (ReleaseLocation **)malloc(Trace_NumReleaseLocations * sizeof(ReleaseLocation *))) == NULL) 
		FatalError("Malloc failed for Trace_ReleaseList");
	
	// Allocate memory for pnode, which is a local linked list used to connect nodes of Trace_ReleaseList 
	if((pnode = (ReleaseLocation **)malloc(Trace_NumReleaseLocations * sizeof(ReleaseLocation *))) == NULL) 
		FatalError("Malloc failed for pnode");
	
	// Allocate memory for voxdim0, which keeps track of distance release points integrated from previous cycle(s) 
	if((voxdim0 = (double *)calloc(Trace_NumReleaseLocations, sizeof(double))) == NULL)
		FatalError("Calloc failed for voxdim0");
	
	// Allocate memory for u0 
	if((u0 = (double **)malloc(Trace_NumReleaseLocations * sizeof(double *))) == NULL)
		FatalError("Malloc failed for u0");
	for(ss = 0; ss < Trace_NumReleaseLocations; ss++) 
		if((u0[ss] = (double *)malloc(3 * sizeof(double))) == NULL)
			FatalError("Malloc failed for u0[%d]", ss);
	
	// Allocate memory for ptdone 
	if((ptdone = (int *)calloc(Trace_NumReleaseLocations, sizeof(int))) == NULL)
		FatalError("Malloc failed for ptdone");
	
	Trace_NumTracers = 0; // number of points that get released 
	
	// Initialize first release to Output_TStart 
	for(ss = 0; ss < Trace_NumReleaseLocations; ss++) {
		if(Trace_ReleasePoints[ss].LeftDomain) 
			Trace_ReleaseList[ss] = NULL;
		else {
			// Alocate memory for what will be the first node in the global linked list Trace_ReleaseList 
			if((nnode = (ReleaseLocation *)malloc(sizeof(ReleaseLocation))) == NULL)
				FatalError("Malloc failed for nnode");
			
			// Initialize values 
			nnode->slide.Start_time = Output_TStart;
			nnode->slide.Stop_time = fmin(Output_TEnd, nnode->slide.Start_time + Trace_IntTLength);
			nnode->slide.Status = UNLAUNCHED;
			for(ii = 0; ii < 3; ii++) 
				nnode->pt.X[ii] = Trace_ReleasePoints[ss].X[ii];
			nnode->pt.LeftDomain = 0;
			nnode->pt.ElementIndex = Trace_ReleasePoints[ss].ElementIndex;
			nnode->pt.Scalar = 0.0;
			
			// Link node to lists 
			Trace_ReleaseList[ss] = nnode;
			pnode[ss] = nnode;  
			Trace_NumTracers++;
		}
	}
	
	
	// Step through the velocity data and determine release times
	for(df = Data_FirstFrame; df < Data_LastFrame; df++) {
	
		// Load velocity data slide 
		readveldata(df);
		
		// Determine time interval of loaded data
		Data_LoadedTMin = Data_TMin + df * Data_TDelta;
		Data_LoadedTMax = Data_LoadedTMin + Data_TDelta;
		printf("Data loaded in memory spans t = %g to %g\n", Data_LoadedTMin, Data_LoadedTMax);
		fflush(stdout);
		
		// Store velocity from Data_FirstFrame for each release point--it will be used as the "positive" flow direction 
		if(df == Data_FirstFrame) {
			// Loop over release points 
			for(ss = 0; ss < Trace_NumReleaseLocations; ss++) {
				// Only consider release points located inside the velocity domain 
				if(!Trace_ReleasePoints[ss].LeftDomain) {
					// Interpolate velocity at release point at start of interval 
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
			
		// Loop over release points and determine release times during interval of currently loaded data
		for(ss = 0; ss < Trace_NumReleaseLocations; ss++) {
			
			// Only consider release points located inside the velocity domain
			if(!Trace_ReleasePoints[ss].LeftDomain && !ptdone[ss]) {
			
				// Interpolate velocity at release point at bounds of loaded data 
				if(Data_MeshType == CARTESIAN) {
					GetVelocity_Cartesian(Data_LoadedTMin, &Trace_ReleasePoints[ss], us);
					GetVelocity_Cartesian(Data_LoadedTMax, &Trace_ReleasePoints[ss], ue);
				}     
				else if(Data_MeshType == UNSTRUCTURED) { 
					GetVelocity_Unstructured(Data_LoadedTMin, &Trace_ReleasePoints[ss], us);
					GetVelocity_Unstructured(Data_LoadedTMax, &Trace_ReleasePoints[ss], ue);
				}	
				// Calculate "signed magnitude" of velocity 
				usbar = SIGN(sqrt(vdot(us, us, 3)), vdot(us, u0[ss], 3));
				uebar = SIGN(sqrt(vdot(ue, ue, 3)), vdot(ue, u0[ss], 3));
				
				//printf("Point %d: usbar = %.9f, uebar = %.9f\n", ss, usbar, uebar);
				
				if(usbar < TINY && uebar < TINY) { // Reverse or stagnant flow, do nothing 
					//printf("Reverse flow, continuing\n");
					continue;
				}
				else if(fabs(usbar - uebar) < 1.0e-6) { // Flow steady, solve linear equation for next release 
				
					if(voxdim0[ss] > TINY)
						ti = Data_LoadedTMin;
					else if(voxdim0[ss] < -TINY) {
						ti = Data_LoadedTMin;
						voxdim0[ss] = 0;
					}
					else
						ti = pnode[ss]->slide.Start_time;
					
					while(1) {
						//printf("ti = %g, %s", ti);
						d = voxdim - voxdim0[ss];
						if(d < 0) 
							FatalError("d < 0"); 
						
						tii = ti + d / usbar; // Linear equation 
						
						//printf("d = %g, tii = %g, %s", d, tii);
						if (tii >= Output_TEnd) {
							
							ptdone[ss] = 1;
							break;
						} else {
						
							if(tii <= ti) 
								FatalError("Next output time not after previous");
							else {
								if(tii >= Trace_ReleaseTMax) { // Release time greater than requested, break 
									//printf("Terminating point %d\n", ss);
									ptdone[ss] = 1;
									break; 
								}
								else if(tii > Data_LoadedTMax) {
									// Save distance point will travel during the remainder of current interval and break 
									voxdim0[ss] = voxdim0[ss] + uebar * (Data_LoadedTMax - ti);
									//printf("tii > Data_LoadedTMax, setting voxdim0[%d] = %.9f, %.9f\n", ss, voxdim0[ss], 0.5 * (uebar + uibar) * (Data_LoadedTMax - ti)); 
									if(voxdim0[ss] < 0) 
										FatalError("voxdim0 < 0");
									break;
								}
								else {
									//printf("creating new release time at %g\n", tii);
									// Create next entry in global linked list Trace_ReleaseList 
									if((nnode = (ReleaseLocation *)malloc(sizeof(ReleaseLocation))) == NULL)
										FatalError("Malloc failed for nnode");
									nnode->slide.Start_time = tii;
									nnode->slide.Stop_time = fmin(Output_TEnd, nnode->slide.Start_time + Trace_IntTLength);
									nnode->slide.Status = UNLAUNCHED;
									for(ii = 0; ii < 3; ii++) 
										nnode->pt.X[ii] = Trace_ReleasePoints[ss].X[ii];
									nnode->pt.LeftDomain = 0;
									nnode->pt.ElementIndex = Trace_ReleasePoints[ss].ElementIndex;
									nnode->pt.Scalar = 0.0;
								
									// Link node to lists 
									pnode[ss]->next = nnode;	  
									pnode[ss] = nnode;
									// Reset voxdim0[ss] 
									voxdim0[ss] = 0;
									// Update ti 
									ti = tii;
									Trace_NumTracers++;
								}
							}
						
						}
					}	

				
				} else { // Flow changing, solve quadratic equation for next release 
					// Velocity varies linearly in time, e.g. u = M*t + N, solve for M and N 
					
					M = (usbar - uebar) / (Data_LoadedTMin - Data_LoadedTMax);
					N = usbar - M * Data_LoadedTMin;
					
					if(usbar < 0 && uebar > 0) {
						// Initiate a release when flow becomes positive 
						if((nnode = (ReleaseLocation *)malloc(sizeof(ReleaseLocation))) == NULL)
							FatalError("Malloc failed for nnode");
						nnode->slide.Start_time = -N/M;
						nnode->slide.Stop_time = fmin(Output_TEnd, nnode->slide.Start_time + Trace_IntTLength);
						nnode->slide.Status = UNLAUNCHED;
						for(ii = 0; ii < 3; ii++) 
							nnode->pt.X[ii] = Trace_ReleasePoints[ss].X[ii];
						nnode->pt.LeftDomain = 0;
						nnode->pt.ElementIndex = Trace_ReleasePoints[ss].ElementIndex;
						nnode->pt.Scalar = 0.0;
						
						// Link node to lists 
						pnode[ss]->next = nnode;	  
						pnode[ss] = nnode;
						// Reset voxdim0[ss] 
						voxdim0[ss] = 0;
						ti = -N/M;
						Trace_NumTracers++;
					} else {
						if(voxdim0[ss] > TINY)
							ti = Data_LoadedTMin;
						else if(voxdim0[ss] < -TINY) {
							ti = Data_LoadedTMin;
							voxdim0[ss] = 0;
						}
						else
							ti = pnode[ss]->slide.Start_time;
					}
					
					while(1) { 
						//printf("ti = %g, %s", ti);
						
						// Interpolate velocity at time ti 
						uibar = M*ti + N;
						
						// Set target distance to advect point
						d = voxdim - voxdim0[ss];
						if(d < 0) {fprintf(stderr, "Error: d < 0\n"); exit(1);}
						
						// Set up quadratic equation for next release time, tii = (-b + sqrt(b * b - 4 * a * c)) / (2 * a) 
						a = M;
						b = N + uibar - M * ti;
						c = -N * ti - uibar * ti - 2*d;
						
					// 
					//	 Since velocity is linear in time, distance a point travels from time ti to time tii is given by 
					//	 d = 0.5 * (uiibar + uibar) * (tii - ti), where uiibar is the velocity magnitude at time tii. 
					//	 However uiibar = M * tii + N, so we solve d = 0.5 * (M * tii + N + uibar) * (tii - ti), or rewriting, 
					//	 M * tii^s + (N + uibar - M * ti) * tii + (-N * ti - uibar * ti - 2d) = 0 
					//	 
						
						// Check if imaginary roots 
						if(b * b - 4 * a * c < 0) { 
							if(usbar < uebar)  // Flow rate increasing, solution should exist 
								FatalError("Next output time has imaginary roots even with increasing flow");
							else { // Flow rate must be decreasing too fast, break 
								//printf("flow rate decreasing too fast, breaking\n");
								voxdim0[ss] = -1; // Used so next consideration of point ss starts from Data_LoadedTMin 
								break;
							}
						}
						else {
							// Compute next release time 
							tii = (-b + sqrt(b * b - 4 * a * c)) / (2 * a);
							
							//printf("d = %g, tii = %g, ", d, tii); 
							if (tii >= Output_TEnd) {
							
								ptdone[ss] = 1;
								break;
							} else {
							
								if(tii <= ti) 
									FatalError("Next output time not after previous");
								else {
									if(tii >= Trace_ReleaseTMax) { // Release time greater than requested, break 
										//printf("Terminating point %d\n", ss);
										ptdone[ss] = 1;
										break; 
									}
									else if(tii > Data_LoadedTMax) {
										// Save distance point will travel during the remainder of current interval and break 
										voxdim0[ss] = voxdim0[ss] + 0.5 * (uebar + uibar) * (Data_LoadedTMax - ti);
										//printf("tii > Data_LoadedTMax, setting voxdim0[%d] = %.9f, %.9f\n", ss, voxdim0[ss], 0.5 * (uebar + uibar) * (Data_LoadedTMax - ti)); 
										if(voxdim0[ss] < 0) 
											FatalError("voxdim0 < 0");
										break;
									}
									else {
										//printf("creating new release time at %g\n", tii);
										// Create next entry in global linked list Trace_ReleaseList 
										if((nnode = (ReleaseLocation *)malloc(sizeof(ReleaseLocation))) == NULL)
											FatalError("Malloc failed for nnode");
										nnode->slide.Start_time = tii;
										nnode->slide.Stop_time = fmin(Output_TEnd, nnode->slide.Start_time + Trace_IntTLength);
										nnode->slide.Status = UNLAUNCHED;
										for(ii = 0; ii < 3; ii++) 
											nnode->pt.X[ii] = Trace_ReleasePoints[ss].X[ii];
										nnode->pt.LeftDomain = 0;
										nnode->pt.ElementIndex = Trace_ReleasePoints[ss].ElementIndex;
										nnode->pt.Scalar = 0.0;
									
										// Link node to lists 
										pnode[ss]->next = nnode;	  
										pnode[ss] = nnode;
										// Reset voxdim0[ss] 
										voxdim0[ss] = 0;
										// Update ti 
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
		
	
	} // for(df = Data_FirstFrame; df < Data_LastFrame; df++) //

	printf("OK!\n");
	printf("Total number of points to be released is %d\n", Trace_NumTracers);
	fflush(stdout);
	
	// Add terminator to last node of global linked list Trace_ReleaseList 
	for(ss = 0; ss < Trace_NumReleaseLocations; ss++) 
		if(!Trace_ReleasePoints[ss].LeftDomain) 
			pnode[ss]->next = NULL;
	
	// Free local linked list used to connect nodes of Trace_ReleaseList and Trace_ReleasePoints and voxdim0 array 
	free(pnode);
	pnode = NULL;
	free(Trace_ReleasePoints);
	Trace_ReleasePoints = NULL;
	free(voxdim0);
	voxdim0 = NULL;
	free(ptdone);
	ptdone = NULL;
	
	N_Total = Trace_NumTracers;
	
	// Linked List will be converted into a 1-D array called Tracer
	// Allocate memory for Tracer
		
			if((Tracer.x = (double *)malloc(N_Total*sizeof(double))) == NULL) {
				fprintf(stderr, "Malloc failed for Tracer.x \n");
				exit(1);
			}

			if((Tracer.y = (double *)malloc(N_Total*sizeof(double))) == NULL) {
				fprintf(stderr, "Malloc failed for Tracer.y \n");
				exit(1);
			}

		
				if((Tracer.z = (double *)malloc(N_Total*sizeof(double))) == NULL) {
					fprintf(stderr, "Malloc failed for Tracer.z \n");
					exit(1);
				}
		
			if((Tracer.ElementIndex = (int *)malloc(N_Total*sizeof(int))) == NULL) {
				fprintf(stderr, "Malloc failed for Tracer.ElementIndex \n");
				exit(1);
			}
	
			if((Tracer.LeftDomain = (int *)malloc(N_Total*sizeof(int))) == NULL) {
				fprintf(stderr, "Malloc failed for Tracer.leftDomain \n");
				exit(1);
			}

			if((Tracer.Start_time = (double *)malloc(N_Total*sizeof(double))) == NULL) {
				fprintf(stderr, "Malloc failed for Tracer.Start_time \n");
				exit(1);
			}

			if((Tracer.Stop_time = (double *)malloc(N_Total*sizeof(double))) == NULL) {
				fprintf(stderr, "Malloc failed for Tracer.Stop_time \n");
				exit(1);
			}
		
		
			
			
			if((Tracer1.x = (double *)malloc(N_Total*sizeof(double))) == NULL) {
				fprintf(stderr, "Malloc failed for Tracer1.x \n");
				exit(1);
			}

			if((Tracer1.y = (double *)malloc(N_Total*sizeof(double))) == NULL) {
				fprintf(stderr, "Malloc failed for Tracer1.y \n");
				exit(1);
			}

		
				if((Tracer1.z = (double *)malloc(N_Total*sizeof(double))) == NULL) {
					fprintf(stderr, "Malloc failed for Tracer1.z \n");
					exit(1);
				}
		
			if((Tracer1.ElementIndex = (int *)malloc(N_Total*sizeof(int))) == NULL) {
				fprintf(stderr, "Malloc failed for Tracer1.ElementIndex \n");
				exit(1);
			}
	
			if((Tracer1.LeftDomain = (int *)malloc(N_Total*sizeof(int))) == NULL) {
				fprintf(stderr, "Malloc failed for Tracer1.leftDomain \n");
				exit(1);
			}

			if((Tracer1.Start_time = (double *)malloc(N_Total*sizeof(double))) == NULL) {
				fprintf(stderr, "Malloc failed for Tracer1.Start_time \n");
				exit(1);
			}

			if((Tracer1.Stop_time = (double *)malloc(N_Total*sizeof(double))) == NULL) {
				fprintf(stderr, "Malloc failed for Tracer1.Stop_time \n");
				exit(1);
			}
			
			if((index1 = (int *)malloc(N_Total*sizeof(int))) == NULL) {
				fprintf(stderr, "Malloc failed for index \n");
				exit(1);
			}
			
			if((Tracer.Status = (int *)malloc(N_Total*sizeof(int))) == NULL) {
				fprintf(stderr, "Malloc failed for Tracer.Status \n");
				exit(1);
			}
		
		
		int num1 = 0, num2 = 0;
		// Loop over all release locations 
		for(ss = 0; ss < Trace_NumReleaseLocations; ss++) {
			
			num2 += num1;
			num1 = 0;
			
			for (rnode = Trace_ReleaseList[ss]; rnode != NULL; rnode = rnode->next) {
			
				Tracer.x[num2 + num1] = rnode->pt.X[0];
				Tracer.y[num2 + num1] = rnode->pt.X[1];
				Tracer.z[num2 + num1] = rnode->pt.X[2];
				
				Tracer.ElementIndex[num2 + num1] = rnode->pt.ElementIndex;
				Tracer.LeftDomain[num2 + num1] = rnode->pt.LeftDomain;
				
				Tracer.Start_time[num2 + num1] = rnode->slide.Start_time;
				Tracer.Stop_time[num2 + num1] = rnode->slide.Stop_time;
				
				Tracer.Status[num2 + num1] = rnode->slide.Status;
				
				//printf("starttime = %f \t",Tracer.Start_time[num2 + num1]);
				//printf("starttime = %f \n",Tracer.Stop_time[num2 + num1]);
				num1++;		
			}
		}
		// Release Linked list
		ReleaseLocation *n, *p;
		free(Trace_ReleasePoints);
		for(ss = 0; ss < Trace_NumReleaseLocations; ss++) {
			n = Trace_ReleaseList[ss];
			while(n != NULL) {
				p = n->next;
				free(n);
				n = p;
			}
		}
		free(Trace_ReleaseList);
		Trace_ReleaseList = NULL;
		
		
		
	
}


void LoadReleasePoints(double *voxdim) {
	/***
	 Reads in number of release points and release point locations. Determines whether 
	 each release point is outside velocity data domain and, if unstructured velocity data
	 is being used, the tetrahedral element the release point is located in. 
	 ***/
	
	char BinFilePath[LONGSTRING];
	FILE *BinFileID;
	int i, index = -1, found = 0, guess = -1;
	
	// Open binary file for reading 
	sprintf(BinFilePath, "%s%s", Path_Data, Trace_InFile);
	if((BinFileID = fopen(BinFilePath, "rb")) == NULL)
		FatalError("Could not open %s", BinFilePath);
	
	// Read in the number of release locations 
	if(fread(&Trace_NumReleaseLocations, sizeof(int), 1, BinFileID) < 1)
		FatalError("Could not read Trace_NumReleaseLocations from %s", BinFilePath);
	
	printf("Loading %d release locations and determining if any outside of domain...", Trace_NumReleaseLocations);
	fflush(stdout);
	
	// Allocate memory to store material point information for the release points
	if((Trace_ReleasePoints = (LagrangianPoint *)malloc(Trace_NumReleaseLocations * sizeof(LagrangianPoint))) == NULL) 
		FatalError("Malloc failed for Trace_ReleasePoints in function LoadReleasePoints()");
	
	// Read in coordinates 
	for(i = 0; i < Trace_NumReleaseLocations; i++) 
		if(fread(&Trace_ReleasePoints[i].X[0], sizeof(double), 1, BinFileID) < 1 || 
		   fread(&Trace_ReleasePoints[i].X[1], sizeof(double), 1, BinFileID) < 1 || 
		   fread(&Trace_ReleasePoints[i].X[2], sizeof(double), 1, BinFileID) < 1) 
			FatalError("Could not read complete coordinate list from %s", BinFilePath);
	
	// Compute spacing between points (assumed release points are spaced uniformly) 
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
		fscanf(stdin, "%lf", voxdim); 
	} 
	else
		FatalError("Number of release locations cannot be zero");
	
	// Close file
	fclose(BinFileID);
	
	
	for(i = 0; i < Trace_NumReleaseLocations; i++) {
		
		// Check if release point outside bounding box for data 
		Trace_ReleasePoints[i].LeftDomain = TestOutsideDomain_host(Trace_ReleasePoints[i].X);
		
		// If point in bounding box and velocity data is unstrutured, find element that release point is in 
		if(Data_MeshType == UNSTRUCTURED && !Trace_ReleasePoints[i].LeftDomain) { 
			if(!found) {
				// First point located using global search, and result serves as starting point for further (local) searching 
				index = Get_Element_Global_Search(Trace_ReleasePoints[i].X);
				if(index < 0) {
					Trace_ReleasePoints[i].LeftDomain = 1;
					Trace_ReleasePoints[i].ElementIndex = -1;
				}
				else { // Global search succeeded 
					
					// Not using aux mesh
					Trace_ReleasePoints[i].ElementIndex = index;
					Trace_ReleasePoints[i].LeftDomain = 0;
					guess = index;
					found++;
					
				}
			}
			else { // A point has been located
				// Find location of point using local search 
				index = Get_Element_Local_Search(Trace_ReleasePoints[i].X, guess);
				if(index < 0) {
					// If local seeach fails, always double-check with a global search 
					index = Get_Element_Global_Search(Trace_ReleasePoints[i].X);
					if(index < 0) {
						Trace_ReleasePoints[i].LeftDomain = 1;
						Trace_ReleasePoints[i].ElementIndex = -1;
						
					}
					else {
						 // Not using aux mesh 
							
						Trace_ReleasePoints[i].ElementIndex = index;
						Trace_ReleasePoints[i].LeftDomain = 0;
						guess = index;
						found++;
					}
				}
				else { // Local search succeeded 
					 	
					// Not using aux mesh 
					Trace_ReleasePoints[i].ElementIndex = index;
					Trace_ReleasePoints[i].LeftDomain = 0;
					guess = index;
					found++;
					
				} 
			}
		}
		else if(!Trace_ReleasePoints[i].LeftDomain) 
			found++; 
	}
	
	printf("OK!\n");
	fflush(stdout);
	
	//printf("hiii ...\n\n");
	// Write to Output File
	//copyresult2vtk(0, Trace_NumReleaseLocations, Output_TStart);
	//copyresult2bin(0, Trace_NumReleaseLocations, Output_TStart);
	
	// Rewrite release point file with only points inside domain if needed 
	if(found < Trace_NumReleaseLocations) {
		printf("%d located outside of velocity domain, ", Trace_NumReleaseLocations - found);
		printf("Rewriting %s to contain interior points only.\n", Trace_InFile);
		fflush(stdout);
		// Open file to (re)write data 
		if((BinFileID = fopen(BinFilePath, "wb")) == NULL) 
			FatalError("Could not open file %s", BinFilePath);
		// Write out the number of release locations 
		if(fwrite(&found, sizeof(int), 1, BinFileID) < 1) 
			FatalError("Could not write number of points to file %s", BinFilePath);
		// Write out coordinates 
		for(i = 0; i < Trace_NumReleaseLocations; i++) 
			if(!Trace_ReleasePoints[i].LeftDomain) 
				if(fwrite(&Trace_ReleasePoints[i].X[0], sizeof(double), 1, BinFileID) < 1 ||
				   fwrite(&Trace_ReleasePoints[i].X[1], sizeof(double), 1, BinFileID) < 1 ||
				   fwrite(&Trace_ReleasePoints[i].X[2], sizeof(double), 1, BinFileID) < 1) 
					FatalError("Could not write complete coordinate list to file %s", BinFilePath);
		// Close file 
		fclose(BinFileID);   
	}
	
}



// Function used for releasing tracer launches .....(ss equals rrelease number)
void release(int ss) {

	int i;
	
	char Data_OutFilePath[LONGSTRING];
	
	FILE *fileout;
	
	int num;
	
	sprintf(Data_OutFilePath, "%s%s.bin", Path_Output, Temp_OutFilePrefix);
	
	if((fileout = fopen(Data_OutFilePath, "rb")) == NULL) {
		fprintf(stderr, "Could not open file %s \n", Data_OutFilePath);
		fprintf(stderr, "Please make sure Tracer release file is kept at %s and named %s.bin \n", Path_Output, Temp_OutFilePrefix );
		exit(1);
	}
	
//	Point temp;
//	tempalloc(temp);
	
	/* Read in number of tracers from file */ 
	if(fread(&num, sizeof(int), 1, fileout) < 1)
		printf("Could not read N_Frame (Number of tracers) from %s", Data_OutFilePath);
		
	if (ss == 0) {
	
		//Check if num = N_Frame.
		if (num != N_Frame) {
		
			fprintf(stderr, " Error: The tracer information file is not Valid the current settings.\n Total number of tracers in information file is not equal to total number of tracers for current run. \n\n");
			exit(1);
		}
	}
			
	
	
	/* Read Tracer information from file 
	if((int)fread(temp.x, sizeof(double), num, fileout) < num)
		printf("Could not completely read Tracer data from %s", Data_OutFilePath);
		
	for (i = 0 ; i < num; i++) {	
		printf("Tracer.x[] = %f \n",temp.x[i]);	
	}
	
	if((int)fread(temp.y, sizeof(double), num, fileout) < num)
		printf("Could not completely read Tracer data from %s", Data_OutFilePath);
	
	if((int)fread(temp.z, sizeof(double), num, fileout) < num)
		printf("Could not completely read Tracer data from %s", Data_OutFilePath);
		
	if((int)fread(temp.ElementIndex, sizeof(int), num, fileout) < num)
		printf("Could not completely read Tracer data from %s", Data_OutFilePath);		

	if((int)fread(temp.LeftDomain, sizeof(int), num, fileout) < num)
		printf("Could not completely read Tracer data from %s", Data_OutFilePath);
*/	
	for (i = 0 ; i < num; i++) {

		if((int)fread(&Tracer.x[ss*N_Frame +i], sizeof(double), 1, fileout) < 1)
			FatalError("Could not read x nodal coordinates completely from %s", Data_OutFilePath);
	
		if((int)fread(&Tracer.y[ss*N_Frame +i], sizeof(double), 1, fileout) < 1)
			FatalError("Could not read x nodal coordinates completely from %s", Data_OutFilePath);
		
		if((int)fread(&Tracer.z[ss*N_Frame +i], sizeof(double), 1, fileout) < 1)
			FatalError("Could not read x nodal coordinates completely from %s", Data_OutFilePath);
			
		if(fread(&Tracer.ElementIndex[ss*N_Frame +i], sizeof(int), 1, fileout) < 1)
			FatalError("Could not read x nodal coordinates comrace_BinFileArray[0].FilePathpletely from %s", Data_OutFilePath);
			
		if(fread(&Tracer.LeftDomain[ss*N_Frame +i], sizeof(int), 1, fileout) < 1)
			FatalError("Could not read x nodal coordinates completely from %s", Data_OutFilePath);	
		
		//printf("Tracer.ld[] = %d \n",Tracer.LeftDomain[ss*N_Frame +i]);	
		//printf("Tracer.eid[] = %d \n",Tracer.ElementIndex[ss*N_Frame +i]);	
		
	/*
		// Releasing the initial x,y,z and element index of the tracer from our temp1 buffer...
		Tracer.x[ss*N_Frame + i] = temp.x[i];
	printf("Tracer.x[%d] = %f \n",ss*N_Frame +i, Tracer.x[ss*N_Frame +i]);
		Tracer.y[ss*N_Frame + i] = temp.y[i];

		Tracer.z[ss*N_Frame + i] = temp.z[i];

		Tracer.ElementIndex[ss*N_Frame + i] = temp.ElementIndex[i];
		
		Tracer.LeftDomain[ss*N_Frame + i] = temp.LeftDomain[i];
*/		
		// Calculating the start and stop time for this tracer release.
		Tracer.Start_time[ss*N_Frame + i] = Output_TStart + (ss * Trace_LaunchTimeSpacing);

		if( (Tracer.Start_time[ss*N_Frame + i] + Trace_IntTLength) >= Output_TEnd) {
			Tracer.Stop_time[ss*N_Frame + i] = Tracer.Start_time[ss*N_Frame + i] + Output_TEnd;
		} else {
			Tracer.Stop_time[ss*N_Frame + i] = Tracer.Start_time[ss*N_Frame + i] + Trace_IntTLength;
		}


	}

//	tempclean(temp);

	fclose(fileout);
	printf(" Release %d launched successfully from the information file. \n", ss);


}



void initialize_new(void) {

	// Initializing the 0th Tracer launch ......
	N_Launches = 0 ;
	
	int i,j,k;
	//Point	temp;
	

	// initializing Launch time for different tracer launch times.....
	for (i = 0; i < Trace_NumLaunchTimes; i++) {

		Launch_time[i] = Output_TStart + (double)(i * Trace_LaunchTimeSpacing);
	}

	// initializing different output times for output files .......
	for (i = 0; i < Output_TRes; i++) {

		Output_time[i] = Output_TStart + (double)(i* Output_TDelta);
		
	}
	
	
	// Check How many data files will be used for integration only valid if data is periodic
	if (Data_TPeriodic) {
	
		// Get the stop time of last Data file
		double stop_data = Data_TMin + (double)((Data_TRes - 1)* Data_TDelta);
		
		if (Output_time[Output_TRes - 1] > stop_data) {
		
			frames_data = (Output_time[Output_TRes - 1]/Data_TDelta == 0) ? (floor(Output_time[Output_TRes - 1]/Data_TDelta)) : (floor(Output_time[Output_TRes - 1]/Data_TDelta) + 1);
		} else
			frames_data = Data_TRes - 1;
	
	
	}
	
	printf("frames_data = %d \n\n", frames_data);

	// initializing start time and stop time for different data set .......
	//for (i = 0; i < (Data_TRes - 1); i++) {
	
	for (i = 0; i < frames_data; i++) {
	
		DataTime1[i].Start_time = Data_TMin + (double)(i* Data_TDelta);
		DataTime1[i].Stop_time = DataTime1[i].Start_time + Data_TDelta;
		
		//printf("start time for data = %f and end time for data = %f \n",Data_Time1[i].Start_time, Data_Time1[i].Stop_time);
		
	}		

	//if (Output_time[Output_TRes - 1] > DataTime1[Data_TRes - 2].Stop_time) {
	
		
	
	
	//}
	if ( Trace_ReleaseStrategy == 0 ) {

		if (Generate_Frame) {

			// Calculating initial tracer positions. Also storing these positions into our temp buffer so that we can use for the later launches of the tracer grids.
	
			// These loops will initialize only the positions......
			for(i = 0; i < Trace_CartMesh.XRes; i++) {
				for(j = 0; j < Trace_CartMesh.YRes; j++) {
					for(k = 0; k < Trace_CartMesh.ZRes; k++) {

							Tracer.x[getindexwosshost(i, j, k)] = Trace_CartMesh.XMin + (i * Trace_CartMesh.XDelta);


							Tracer.y[getindexwosshost(i, j, k)] = Trace_CartMesh.YMin + (j * Trace_CartMesh.YDelta);

	
							Tracer.z[getindexwosshost(i, j, k)] = Trace_CartMesh.ZMin + (k * Trace_CartMesh.ZDelta);
							//printf("num = %d \t",getindexwosshost(i, j, k) );
							//printf("x = %0.9f  y = %0.9f and z = %0.9f \n",Tracer.x[getindexwosshost(i, j, k)], Tracer.y[getindexwosshost(i, j, k)], Tracer.z[getindexwosshost(i, j, k)] );

							// This value is set to 1 by default. It will be changed to 0 in globalsearch if tracer is found inside the domain.
							Tracer.LeftDomain[getindexwosshost(i, j, k)] = 1;

					}
				}
			}


			// This loop initializes different other properties (such as start time, end time and element index) of our first tracer launch ....
			for (i = 0 ; i < N_Frame; i++) {

				Tracer.Start_time[i] = Output_TStart;
		
				if( (Output_TStart + Trace_IntTLength) >= Output_TEnd) {
					Tracer.Stop_time[i] = Output_TStart  + Output_TEnd;
				} else {
					Tracer.Stop_time[i] = Output_TStart + Trace_IntTLength;
				}

				Tracer.ElementIndex[i] = -1;



			}


			// Get element index of the traccer launch if mesh is unstructured
			if(Data_MeshType == UNSTRUCTURED) {
		
				 // Searching on Host
					search_host();
				
			} else {
				// Cartesian Mesh.
			
				double Point[3];
			
				for(i = 0; i < Trace_CartMesh.XRes; i++) {
					for(j = 0; j < Trace_CartMesh.YRes; j++) {
						for(k = 0; k < Trace_CartMesh.ZRes; k++) {

								Point[0] = Tracer.x[getindexwosshost(i, j, k)];
								Point[1] = Tracer.y[getindexwosshost(i, j, k)];
								Point[2] = Tracer.z[getindexwosshost(i, j, k)];

								// Check if tracer is found inside the domain.
								Tracer.LeftDomain[getindexwosshost(i, j, k)] = TestOutsideDomain_host(Point);
								//if (Tracer.LeftDomain[getindexwosshost(i, j, k)] == 1)
									//printf("Hiiiii \n");
								//printf("x = %0.9f  y = %0.9f \n",Tracer.x[getindexwosshost(i, j, k)], Tracer.y[getindexwosshost(i, j, k)] );
						}
					}
				}
			
			
			
		
			}
			// Store initial frame into output file which can be loaded for different releases...
			temp2file();
		
			printf("Release 0 Launched \n\n");
		} else { // Read in from temp file itself
	
		//	release(0);
			FatalError("This Function of reading in from input file is currently demoted \n\n\n");
		
	
		}
		//for ( i = 0; i < N_Frame; i++)
			//	printf("x = %0.9f  y = %0.9f \n",Tracer.x[i], Tracer.y[i] );

		
	// We will output the first tracer launch to the output file (This will be always done)
	
		Output_Release = 0;	
		copyresult2vtk(0, 1, Output_TStart);
		copyresult2bin(0, 1, Output_TStart);
	
		// Our first tracer frame is launched.
		N_Launches = N_Launches + 1;
	}
}



// function used for launching tracers ....
int Launch_Tracer(double tmin, double tmax, int N_Launches) {

	int i;
	int sum = 0;
	
	// Check for all tracer launches ...
	for (i = 0; i < Trace_NumLaunchTimes; i++) {
		
		// Logic: if they are supposed to be launched during the loaded time (i.e between tmin and tmax) 
		if ( (tmin < Launch_time[i]) && (Launch_time[i] <= tmax)) {
			sum = sum + 1; 
		}
	}

	for (i = 0; i < sum; i++) {

		// Release and set tracers ....
		release(N_Launches);
		
		N_Launches = N_Launches + 1;

	}


	return N_Launches;

}


// This function is no more in use for unstructured grid...
int Output_Tracer(double tmin, double tmax, int Output_Release) {

	int i,j;
	int N[Output_TRes];
	int sum =0;
	int sum1 =0;

	for (i = 0; i < Output_TRes; i++) {
	
		if ( (Output_time[i] > (tmin )) && (Output_time[i] <= (tmax + TINY))) {
			sum = sum + 1; 
			printf("Ans: %0.9f and tmin = %0.9f and tmax = %0.9f \n",Output_time[i] , tmin, tmax);
			for (j = 0; j < Trace_NumLaunchTimes ;j++) {
				if(Launch_time[j] <= (Output_time[i] + TINY))
					sum1++;
			}
			N[(sum -1)] = sum1;
		}
	}

	for (i = 0; i < sum; i++) {
	
		copyresult2vtk(Output_Release + 1, N[i], Output_time[Output_Release + 1]);
		Output_Release++;
		printf("Output release number %d was executed .... \n", Output_Release);
	
	}


	return Output_Release;

}




