#if USE_DOUBLE
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#endif
void GetVel_unstruct3D( double tloc, double *dXdt, __global const double *vel_u0, __global const double *vel_u1, __global const double *vel_v0, __global const double *vel_v1, __global const double *vel_w0, __global const double *vel_w1, __global const int *Mesh_Element_Node1, __global const int *Mesh_Element_Node2,__global const int *Mesh_Element_Node3, __global const int *Mesh_Element_Node4,  __global const int *Mesh_Element_Neighborindex1, __global const int *Mesh_Element_Neighborindex2, __global const int *Mesh_Element_Neighborindex3, __global const int *Mesh_Element_Neighborindex4 ,  int eid, double r,  double s, double t ) {

	if (eid == -1) {  // Point out of Domain ....

		dXdt[0] = 0.0;
		dXdt[1] = 0.0;
		dXdt[2] = 0.0;
		
		return;

	} else {
	
		double U1t, U2t, U3t, U4t, V1t, V2t, V3t, V4t, W1t, W2t, W3t, W4t;
		
		// Interpolate velocity at nodes in time 
		U1t = (1 - tloc) * vel.u0[Mesh_Element_Node1[eid]] + tloc * vel.u1[Mesh_Element_Node1[eid]];
		V1t = (1 - tloc) * vel.v0[Mesh_Element_Node1[eid]] + tloc * vel.v1[Mesh_Element_Node1[eid]];
		W1t = (1 - tloc) * vel.w0[Mesh_Element_Node1[eid]] + tloc * vel.w1[Mesh_Element_Node1[eid]];

		U2t = (1 - tloc) * vel.u0[Mesh_Element_Node2[eid]] + tloc * vel.u1[Mesh_Element_Node2[eid]];
		V2t = (1 - tloc) * vel.v0[Mesh_Element_Node2[eid]] + tloc * vel.v1[Mesh_Element_Node2[eid]];
		W2t = (1 - tloc) * vel.w0[Mesh_Element_Node2[eid]] + tloc * vel.w1[Mesh_Element_Node2[eid]];
				
		U3t = (1 - tloc) * vel.u0[Mesh_Element_Node3[eid]] + tloc * vel.u1[Mesh_Element_Node3[eid]];
		V3t = (1 - tloc) * vel.v0[Mesh_Element_Node3[eid]] + tloc * vel.v1[Mesh_Element_Node3[eid]];
		W3t = (1 - tloc) * vel.w0[Mesh_Element_Node3[eid]] + tloc * vel.w1[Mesh_Element_Node3[eid]];
		
		U4t = (1 - tloc) * vel.u0[Mesh_Element_Node4[eid]] + tloc * vel.u1[Mesh_Element_Node4[eid]];
		V4t = (1 - tloc) * vel.v0[Mesh_Element_Node4[eid]] + tloc * vel.v1[Mesh_Element_Node4[eid]];
		W4t = (1 - tloc) * vel.w0[Mesh_Element_Node4[eid]] + tloc * vel.w1[Mesh_Element_Node4[eid]];
	
		// Get interpolate velocity at point in space using time-interpolated velocity at nodes
		dXdt[0] = U1t + (U2t - U1t) * fabs(r) + (U3t - U1t) * fabs(s) + (U4t - U1t) * fabs(t);
		dXdt[1] = V1t + (V2t - V1t) * fabs(r) + (V3t - V1t) * fabs(s) + (V4t - V1t) * fabs(t);
		dXdt[2] = W1t + (W2t - W1t) * fabs(r) + (W3t - W1t) * fabs(s) + (W4t - W1t) * fabs(t);
	

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
		//	GetVel_unstruct3D(tloc[tid], k, Vel_U0, Vel_U1, Vel_V0, Vel_V1, Vel_W0, Vel_W1, Mesh_Element_Node1, Mesh_Element_Node2, Mesh_Element_Node3, Mesh_Element_Node4,  Mesh_Element_Neighborindex1, Mesh_Element_Neighborindex2, Mesh_Element_Neighborindex3, Mesh_Element_Neighborindex4 , eid[tid], r[tid], s[tid], t[tid] );
			
			
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



