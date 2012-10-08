# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <CL/cl.h>

# include "globals.h"
# include "GlobalSearch.h"
# include "velocity.h"
# include "clean.h"



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
			x0 = MeshNodeArray_double.x[MeshElementArray.Node1[i]];
			y0 = MeshNodeArray_double.y[MeshElementArray.Node1[i]];
			z0 = MeshNodeArray_double.z[MeshElementArray.Node1[i]];
			x1 = MeshNodeArray_double.x[MeshElementArray.Node2[i]];
			y1 = MeshNodeArray_double.y[MeshElementArray.Node2[i]];
			z1 = MeshNodeArray_double.z[MeshElementArray.Node2[i]];
			x2 = MeshNodeArray_double.x[MeshElementArray.Node3[i]];
			y2 = MeshNodeArray_double.y[MeshElementArray.Node3[i]];
			z2 = MeshNodeArray_double.z[MeshElementArray.Node3[i]];
			x3 = MeshNodeArray_double.x[MeshElementArray.Node4[i]];
			y3 = MeshNodeArray_double.y[MeshElementArray.Node4[i]];
			z3 = MeshNodeArray_double.z[MeshElementArray.Node4[i]];
			
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

/*
int TestOutsideDomain(double point[3]) {
	  
//	 inside data bounding box --> return 0
//	 outside data bounding box --> return 1
	
	
	if(point[0] < (Data_MeshBounds.XMin - TINY) || point[0] > (Data_MeshBounds.XMax + TINY)) 
		return(1);
	else if(point[1] < (Data_MeshBounds.YMin - TINY) || point[1] > (Data_MeshBounds.YMax + TINY)) 
		return(1);
	else if(point[2] < (Data_MeshBounds.ZMin - TINY) || point[2] > (Data_MeshBounds.ZMax + TINY)) 
		return(1);
	else
		return(0);
	
}
*/

void search_host(void) {

	
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
		
		//printf("Xpoint[0] = %f \n",Xpoint[1] );
		
		// See if coordinates outside of data bounding box
		Tracer.LeftDomain[ii] = TestOutsideDomain_host(Xpoint);
/*		
		if (!Tracer.LeftDomain[ii]) {
		
			seed = Get_Element_Global_Search(Xpoint);
		
			Tracer.ElementIndex[ii] = seed;
			
			if (seed > 0) {
			
				Tracer.LeftDomain[ii] = 0;
				found++;
				foundGS++;
				
			} else {
				Tracer.LeftDomain[ii] = 1;
				
			}
		}
	} */	
	
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
		
	//	printf("eid = %d \n", Tracer.ElementIndex[ii]);
	
	}
	
	if(LocalSearchChecking) 
		printf("  %d of %d located in domain (%d caught by global search checking)\n", found, Trace_CartMesh.XRes * Trace_CartMesh.YRes * Trace_CartMesh.ZRes, foundGS);
	else
		printf("  %d of %d located in domain\n", found, Trace_CartMesh.XRes * Trace_CartMesh.YRes * Trace_CartMesh.ZRes);  
	
	fflush(stdout);
	

	
}








