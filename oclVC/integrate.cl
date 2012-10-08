#if USE_DOUBLE
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#endif

# define TINY 1.0e-10

# define FABS(b) (b >= 0.0 ? b : -b)
#define FMIN(a, b) ((a) < (b) ? (a) : (b))
#define fmax(a, b) ((a) > (b) ? (a) : (b))


__kernel void check_int(__global double *Start_time_dev, __global double *Stop_time_dev, __global int *integrate,__global const double *t1, __global const double *t2, const int ss) {

	// Get thread index
	int tid = get_global_id(0);
	
	if( tid >= ss) { // Redundant thread ...
		return; 
	} else {
	
		double tstart,tstop;
//		double tstart1,tstop1;
		
		double t111, t221;
//		double t112, t222;
		tstart = Start_time_dev[tid];
		tstop = Stop_time_dev[tid];
		t111 = t1[tid];
		t221 = t2[tid];
		
		if ( (tstop < t111) || (tstart > t221) ) {
			integrate[tid] = 0;
		} else {
			integrate[tid] = 1;
		}
		
	
		return;
	}

}

__kernel void initialize_timestep3D( __global double *x_dev, __global double *y_dev, __global double *z_dev, __global double *posx, __global double *posy, __global double *posz, __global double *xn0, __global double *xn1, __global double *xn2, __global int *eid, __global int *ElementIndex_dev, const int ss) {

	int tid = get_global_id(0);
	
	if( tid >= ss) { // Redundant thread ...
		return; 
	} else {

		double a0,b0,c0;
		int d0;		
		a0 = x_dev[tid];
		
		b0 = y_dev[tid];
		
	
		c0 = z_dev[tid];
		
		d0 = ElementIndex_dev[tid];		
		
		posx[tid] = a0;
	
		
		posy[tid] = b0;
	
		
		posz[tid] = c0;
		
		
		eid[tid] = d0;
	
		
		xn0[tid] = 0.00;
		xn1[tid] = 0.00;
		xn2[tid] = 0.00;
		
	

		return;
	}
}



void GetVel_unstruct3D( double tloc, double *dXdt, __global const double *vel_u0, __global const double *vel_u1, __global const double *vel_v0, __global const double *vel_v1, __global const double *vel_w0, __global const double *vel_w1, __global const int *Mesh_Element_Node1, __global const int *Mesh_Element_Node2,__global const int *Mesh_Element_Node3, __global const int *Mesh_Element_Node4,  __global const int *Mesh_Element_Neighborindex1, __global const int *Mesh_Element_Neighborindex2, __global const int *Mesh_Element_Neighborindex3, __global const int *Mesh_Element_Neighborindex4 , const int eid, const double r, const double s, const double t ) {

	if (eid == -1) {  // Point out of Domain ....

		dXdt[0] = 0.0;
		dXdt[1] = 0.0;
		dXdt[2] = 0.0;
		
		return;

	} else {
	
		double U1t, U2t, U3t, U4t, V1t, V2t, V3t, V4t, W1t, W2t, W3t, W4t;
		
		// Interpolate velocity at nodes in time 
		U1t = (1 - tloc) * vel_u0[Mesh_Element_Node1[eid]] + tloc * vel_u1[Mesh_Element_Node1[eid]];
		V1t = (1 - tloc) * vel_v0[Mesh_Element_Node1[eid]] + tloc * vel_v1[Mesh_Element_Node1[eid]];
		W1t = (1 - tloc) * vel_w0[Mesh_Element_Node1[eid]] + tloc * vel_w1[Mesh_Element_Node1[eid]];

		U2t = (1 - tloc) * vel_u0[Mesh_Element_Node2[eid]] + tloc * vel_u1[Mesh_Element_Node2[eid]];
		V2t = (1 - tloc) * vel_v0[Mesh_Element_Node2[eid]] + tloc * vel_v1[Mesh_Element_Node2[eid]];
		W2t = (1 - tloc) * vel_w0[Mesh_Element_Node2[eid]] + tloc * vel_w1[Mesh_Element_Node2[eid]];
				
		U3t = (1 - tloc) * vel_u0[Mesh_Element_Node3[eid]] + tloc * vel_u1[Mesh_Element_Node3[eid]];
		V3t = (1 - tloc) * vel_v0[Mesh_Element_Node3[eid]] + tloc * vel_v1[Mesh_Element_Node3[eid]];
		W3t = (1 - tloc) * vel_w0[Mesh_Element_Node3[eid]] + tloc * vel_w1[Mesh_Element_Node3[eid]];
		
		U4t = (1 - tloc) * vel_u0[Mesh_Element_Node4[eid]] + tloc * vel_u1[Mesh_Element_Node4[eid]];
		V4t = (1 - tloc) * vel_v0[Mesh_Element_Node4[eid]] + tloc * vel_v1[Mesh_Element_Node4[eid]];
		W4t = (1 - tloc) * vel_w0[Mesh_Element_Node4[eid]] + tloc * vel_w1[Mesh_Element_Node4[eid]];
	
		// Get interpolate velocity at point in space using time-interpolated velocity at nodes
//		dXdt[0] = U1t + (U2t - U1t) * fabs(r) + (U3t - U1t) * fabs(s) + (U4t - U1t) * fabs(t);
//		dXdt[1] = V1t + (V2t - V1t) * fabs(r) + (V3t - V1t) * fabs(s) + (V4t - V1t) * fabs(t);
//		dXdt[2] = W1t + (W2t - W1t) * fabs(r) + (W3t - W1t) * fabs(s) + (W4t - W1t) * fabs(t);
		
		dXdt[0] = U1t + (U2t - U1t) * FABS(r) + (U3t - U1t) * FABS(s) + (U4t - U1t) * FABS(t);
		dXdt[1] = V1t + (V2t - V1t) * FABS(r) + (V3t - V1t) * FABS(s) + (V4t - V1t) * FABS(t);
		dXdt[2] = W1t + (W2t - W1t) * FABS(r) + (W3t - W1t) * FABS(s) + (W4t - W1t) * FABS(t);
		//cuPrintf("U1 = %0.9f \n", U1t);
		//cuPrintf("vel.u0 = %0.9f and vel.u1 = %0.9f and tloc = %0.9f\n", vel.u0[MeshElementArray_device.Node1[eid]], vel.u1[MeshElementArray_device.Node1[eid]], tloc);		
		
		

		return;
	}
}



__kernel void LocalSearch3D(__global double *posx, __global double *posy, __global double *posz, __global const int *Mesh_Element_Node1, __global const int *Mesh_Element_Node2 , __global const int *Mesh_Element_Node3, __global const int *Mesh_Element_Node4,  __global const int *Mesh_Element_Neighborindex1, __global const int *Mesh_Element_Neighborindex2, __global const int *Mesh_Element_Neighborindex3, __global const int *Mesh_Element_Neighborindex4 , __global const double *Mesh_Node_x, __global const double *Mesh_Node_y, __global const double *Mesh_Node_z  , __global int *eid, const int ss, __global double *r, __global double *s, __global double *t, __global int *integrate) {

	const int tid = get_global_id(0);
	
	if( tid >= ss) { // Redundant thread ...
		return; 
	} else {
	
		if (integrate[tid] == 1) { // Only Do Local Search If a point has to be integrated over the interval.
		
			double x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3;
			double d,V;
			int ID = eid[tid];
			if(ID == -1) { // Point already left.
		
				return;
		
			} else {
		
			
				while(1) {
				
					// Physical coordinates of nodes of the test element
					x0 = Mesh_Node_x[Mesh_Element_Node1[ID]];
					y0 = Mesh_Node_y[Mesh_Element_Node1[ID]];
					z0 = Mesh_Node_z[Mesh_Element_Node1[ID]];

					x1 = Mesh_Node_x[Mesh_Element_Node2[ID]];
					y1 = Mesh_Node_y[Mesh_Element_Node2[ID]];
					z1 = Mesh_Node_z[Mesh_Element_Node2[ID]];

					x2 = Mesh_Node_x[Mesh_Element_Node3[ID]];
					y2 = Mesh_Node_y[Mesh_Element_Node3[ID]];
					z2 = Mesh_Node_z[Mesh_Element_Node3[ID]];
				
					x3 = Mesh_Node_x[Mesh_Element_Node4[ID]];
					y3 = Mesh_Node_y[Mesh_Element_Node4[ID]];
					z3 = Mesh_Node_z[Mesh_Element_Node4[ID]];
			
					// Determinant of mapping from natural to physical coordinates of test element
					V = (x1 - x0) * ((y2 - y0) * (z3 - z0) - (z2 - z0) * (y3 - y0)) + (x2 - x0) * ((y0 - y1) * (z3 - z0) - (z0 - z1) * (y3 - y0)) + (x3 - x0) * ((y1 - y0) * (z2 - z0) - (z1 - z0) * (y2 - y0));
			
			
					// Natural coordinates of point to be interpolated
					r[tid] = ( (((z3 - z0) * (y2 - y3) - (z2 - z3) * (y3 - y0)) * (posx[tid] - x0)) + ( ((x3 - x0) * (z2 - z3) - (x2 - x3) * (z3 - z0)) * (posy[tid] - y0)) + ( ((y3 - y0) * (x2 - x3) - (y2 - y3) * (x3 - x0)) *(posz[tid] - z0) ))/ V;
			
					s[tid] = ( (((z3 - z0) * (y0 - y1) - (z0 - z1) * (y3 - y0)) * (posx[tid] - x0)) + (((x3 - x0) * (z0 - z1) - (x0 - x1) * (z3 - z0)) * (posy[tid] - y0)) + (((y3 - y0) * (x0 - x1) - (y0 - y1) * (x3 - x0)) * (posz[tid] - z0)) ) / V;
			
			
					t[tid] = ( (((z1 - z2) * (y0 - y1) - (z0 - z1) * (y1 - y2)) * (posx[tid] - x0)) + (((x1 - x2) * (z0 - z1) - (x0 - x1) * (z1 - z2)) * (posy[tid] - y0)) + (((y1 - y2) * (x0 - x1) - (y0 - y1) * (x1 - x2)) * (posz[tid] - z0)) ) / V;
			
			
					d = FMIN(r[tid], FMIN(s[tid], FMIN(t[tid], 1 - r[tid] - s[tid] - t[tid])));

					
				
					if((d + TINY) >= 0) { // Point inside the element .....
						eid[tid] = ID;
						return;

					} else { // Point is outside ....

						if(FABS(r[tid] - d) < TINY) 
							ID = Mesh_Element_Neighborindex1[ID];
						else if(FABS(s[tid] - d) < TINY) 
							ID = Mesh_Element_Neighborindex2[ID];
						else if(FABS(t[tid] - d) < TINY) 
							ID = Mesh_Element_Neighborindex3[ID];
						else if(FABS(d - (1 - r[tid] - s[tid] - t[tid])) < TINY) 
							ID = Mesh_Element_Neighborindex4[ID];
				
						if (ID == -1) {
					
							eid[tid] = ID;
							return;
						
						}
						
					}

				}
			
			
			}
		
			return;
		}
		else
			return;
	
	}
}



__kernel void compute_points_Unstructure3D_1(__global double *posx, __global double *posy, __global double *posz, __global const double *tloc, __global double *xn0, __global double *xn1, __global double *xn2, __global double *x, __global double *y, __global double *z, __global int *eid, __global double *r, __global double *s, __global double *t, __global const double *Vel_U0, __global const double *Vel_U1, __global const double *Vel_V0, __global const double *Vel_V1, __global const double *Vel_W0, __global const double *Vel_W1, __global const int *Mesh_Element_Node1, __global const int *Mesh_Element_Node2,__global const int *Mesh_Element_Node3, __global const int *Mesh_Element_Node4,  __global const int *Mesh_Element_Neighborindex1, __global const int *Mesh_Element_Neighborindex2,__global const int *Mesh_Element_Neighborindex3, __global const int *Mesh_Element_Neighborindex4 , __global const double *h, __global int *integrate, const int stage, const int ss) {

	int tid = get_global_id(0);
	if( tid >= ss) { // Redundant thread ...
		return; 
	} else {
		
		if (integrate[tid] == 1) {
			double k[3];
		
			// Get velocity at current position given by posx, and posy.
			GetVel_unstruct3D(tloc[tid], k, Vel_U0, Vel_U1, Vel_V0, Vel_V1, Vel_W0, Vel_W1, Mesh_Element_Node1, Mesh_Element_Node2, Mesh_Element_Node3, Mesh_Element_Node4,  Mesh_Element_Neighborindex1, Mesh_Element_Neighborindex2, Mesh_Element_Neighborindex3, Mesh_Element_Neighborindex4 , eid[tid], r[tid], s[tid], t[tid] );
			
			
			xn0[tid] = xn0[tid] + (stage == 0 ? 0.16666666666667*h[tid] : 0.33333333333334*h[tid])*k[0];
			xn1[tid] = xn1[tid] + (stage == 0 ? 0.16666666666667*h[tid] : 0.33333333333334*h[tid])*k[1];
			xn2[tid] = xn2[tid] + (stage == 0 ? 0.16666666666667*h[tid] : 0.33333333333334*h[tid])*k[2];
			
			
		
			posx[tid] = x[tid] + (stage == 2 ? h[tid] : 0.5*h[tid])*k[0];
			posy[tid] = y[tid] + (stage == 2 ? h[tid] : 0.5*h[tid])*k[1];
			posz[tid] = z[tid] + (stage == 2 ? h[tid] : 0.5*h[tid])*k[2];
		
			
		
		}
	
	}

	return;
}

__kernel void compute_points_Unstructure3D_2(__global double *posx, __global double *posy, __global double *posz, __global const double *tloc, __global double *xn0, __global double *xn1, __global double *xn2, __global double *x, __global double *y, __global double *z, __global int *eid, __global double *r, __global double *s, __global double *t, __global const double *Vel_U0, __global const double *Vel_U1, __global const double *Vel_V0, __global const double *Vel_V1, __global const double *Vel_W0, __global const double *Vel_W1, __global const int *Mesh_Element_Node1, __global const int *Mesh_Element_Node2,__global const int *Mesh_Element_Node3, __global const int *Mesh_Element_Node4,  __global const int *Mesh_Element_Neighborindex1, __global const int *Mesh_Element_Neighborindex2,__global const int *Mesh_Element_Neighborindex3, __global const int *Mesh_Element_Neighborindex4 , __global const double *h, __global int *integrate, const int ss) {

	int tid = get_global_id(0);;

	if( tid >= ss) { // Redundant thread ...
		return; 
	} else {
	
		if (integrate[tid] == 1) {
			double k[3];
		
			// Get velocity at current position given by posx, and posy.
			GetVel_unstruct3D(tloc[tid], k, Vel_U0, Vel_U1, Vel_V0, Vel_V1, Vel_W0, Vel_W1, Mesh_Element_Node1, Mesh_Element_Node2, Mesh_Element_Node3, Mesh_Element_Node4, Mesh_Element_Neighborindex1, Mesh_Element_Neighborindex2, Mesh_Element_Neighborindex3, Mesh_Element_Neighborindex4 , eid[tid], r[tid], s[tid], t[tid] );
			
			
			xn0[tid] = xn0[tid] + 0.166666666666667*h[tid]*k[0];
			xn1[tid] = xn1[tid] + 0.166666666666667*h[tid]*k[1];
			xn2[tid] = xn2[tid] + 0.166666666666667*h[tid]*k[2];
			
			
		
			x[tid] = x[tid] + xn0[tid];
			y[tid] = y[tid] + xn1[tid];
			z[tid] = z[tid] + xn2[tid];
		
			
		
		}

	}

	return;
}


