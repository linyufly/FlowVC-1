
# include "unstructured.h"

// This file defines all the kernels used for Cartesian (2D and 3D) Integration. The entire integration algorithm is divided into smaller kernels. Global memory is used to store intermediate information by the integrator. posx, posy and posz keeps the updated x, y, and z position of a particle. xn0, xn1 and xn2 is used as a counter to collect sum of of k1, k2, k3 and k4. This approach reduces the usage of registers. Currently only functions related to velocity interpolation uses more registers. Instruction Level parallelism is used to hide memory latency.


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 					RK4 algorithm																						//
//																														//
// 	k1 = f( t(n), y(n) )																								//
//																														//
//	k2 = f( t(n) + h/2 , y(n) + (h/2 * k1) )																			//
//																														//
//	k3 = f( t(n) + h/2 , y(n) + (h/2 * k2) )																			//
//																														//
//	k4 = f( t(n) + h , y(n) + (h * k3) )																				//
//																														//
//	y(n+1) = y(n) + h/6 * (k1 + 2*k2 + 2*k3 + k4)																		//
//																														//
//																														//
//	xn array stores the sum of (h/6 * (k1 + 2*k2 + 2*k3 + k4))															//
//																														//
//	pos array keeps the intermediate position (y(n), y(n) + (h/2 * k1), y(n) + (h/2 * k2) & y(n) + (h * k3) )			//
//																														//
//	y(n+1) is updated directly from the kernel that computes value of k4												//
//																														//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// 2D Unstructured Kernels

// This function initializes the variables stored in global memory (posx, posy ; xn0, xn1). 
// Approach used for instruction level parallelism: Global memory will be loaded into registers and then copied back to global memory.
__global__ void initialize_timestep2D(Point Tracer_dev, double *posx, double *posy, double *xn0, double *xn1, int ss, int offset) {

	int tid;
	int index;
	// Get thread index
	tid=(2*blockIdx.x*blockDim.x)+threadIdx.x; 
	
	if( tid >= ss) { // Redundant thread ...
		return; 
	} else {
		
		index = tid + offset;
	
		// Load Global memory into registers
		double a0,b0;
		
		double a1,b1;
		
		a0 = Tracer_dev.x[index];
		
		b0 = Tracer_dev.y[index];		
		
		if ((tid + blockDim.x) < ss) {
			
			a1 = Tracer_dev.x[index + blockDim.x];
			b1 = Tracer_dev.y[index + blockDim.x];
				
		}
		
		
		// Copy registers into global memory.
		posx[index] = a0;
		posy[index] = b0;
		xn0[index] = 0.00;
		xn1[index] = 0.00;
		
		
		if ((tid + blockDim.x) < ss) {
		
			posx[index + blockDim.x] = a1;
			posy[index + blockDim.x] = b1;
			xn0[index + blockDim.x] = 0.00;
			xn1[index + blockDim.x] = 0.00;
		
		}
		
		return;
	}
}


// This kernel function computes the velocity of a particle and updates its position based on time-step in posx, posy and posz. This function is utilized to obtain the value of k1, k2 & k3 whose sum is stored in xn0, xn1 and xn2. stage input specifies the constant computing (k1, k2 or k3).

__global__ void compute_points_Unstructure2D_1(double *posx, double *posy, double *tloc, double *xn0, double *xn1, double *x, double *y, int *eid, double *r, double *s, VelData_double vel, Element MeshElementArray_device, double *h, int *integrate, int stage, int ss) {

	int tid;
	
	// Get thread index
	tid=(blockIdx.x*blockDim.x)+threadIdx.x; 
	if( tid >= ss) { // Redundant thread ...
		return; 
	} else {
		
		// Integrate only if this thread need to be integrated.
		if (integrate[tid] == 1) {
			double k[2];
		
			// Get velocity at current position given by posx, and posy.
			GetVel_unstruct2D(tloc[tid], k, vel, MeshElementArray_device, eid[tid], r[tid], s[tid] );
			
			double mul;
			
			// Multiplier used for xn0, xn1 & xn2. See RK4 algorithm for more details.
			mul = (stage == 0) ? (0.1666666666666666666667*h[tid]) : (0.3333333333333333333334*h[tid]);
			
			// Update the xn buffer with the value of k1, k2 or k3
			xn0[tid] = xn0[tid] + mul*k[0];
			xn1[tid] = xn1[tid] + mul*k[1];

			// Multiplier used for updated position of a particle. See RK4 algorithm for more details.
			mul = (stage == 2) ? (h[tid]) : (0.50*h[tid]);

			// Update the intermediate position of a particle.		
			posx[tid] = x[tid] + mul*k[0];
			posy[tid] = y[tid] + mul*k[1];
		
		}
	}

	return;
}

// This function updates the position into the main tracer array. Essentially, this kernel computes the value of k4 and updates the main struct by values xn0 & xn1 respectively.
__global__ void compute_points_Unstructure2D_2(double *posx, double *posy, double *tloc, double *xn0, double *xn1, double *x, double *y, int *eid, double *r, double *s, VelData_double vel, Element MeshElementArray_device, double *h, int *integrate, int ss) {

	int tid;
	
	// Get thread index
	tid=(blockIdx.x*blockDim.x)+threadIdx.x; 
	if( tid >= ss) { // Redundant thread ...
		return; 
	} else {
	
		// Integrate only if this thread need to be integrated.
		if (integrate[tid] == 1) {
			double k[2];
		
			// Get velocity at current position given by posx, and posy.
			GetVel_unstruct2D(tloc[tid], k, vel, MeshElementArray_device, eid[tid], r[tid], s[tid]);
			
			// Add value of k4 into xn array
			xn0[tid] = xn0[tid] + 0.1666666666666666666666666667*h[tid]*k[0];
			xn1[tid] = xn1[tid] + 0.1666666666666666666666666667*h[tid]*k[1];			
		
			// Update the position of a particle in the main struct
			x[tid] = x[tid] + xn0[tid];
			y[tid] = y[tid] + xn1[tid];
		
		}
	}

	return;
}




// 3D Unstructured Kernels
__global__ void initialize_timestep3D(Point Tracer_dev, double *posx, double *posy, double *posz, double *xn0, double *xn1, double *xn2, int ss, int offset) {

	int tid;
	int index;
	// Get thread index
	tid=(2*blockIdx.x*blockDim.x)+threadIdx.x; 
	
	if( tid >= ss) { // Redundant thread ...
		return; 
	} else {

		index = tid + offset;
		
		double a0,b0,c0;
		
		double a1,b1,c1;
	
		// Load Global memory into registers
		a0 = Tracer_dev.x[index];
		
		b0 = Tracer_dev.y[index];
	
		c0 = Tracer_dev.z[index];		
		
		if ((tid + blockDim.x) < ss) {
			
			a1 = Tracer_dev.x[index + blockDim.x];
			b1 = Tracer_dev.y[index + blockDim.x];
			c1 = Tracer_dev.z[index + blockDim.x];
		}
		
		// Copy registers into global memory.
		posx[index] = a0;
		
		posy[index] = b0;
			
		posz[index] = c0;
			
		xn0[index] = 0.00;
		xn1[index] = 0.00;
		xn2[index] = 0.00;
		
		if ((tid + blockDim.x) < ss) {
		
			posx[index + blockDim.x] = a1;
			posy[index + blockDim.x] = b1;
			posz[index + blockDim.x] = c1;
			xn0[index + blockDim.x] = 0.00;
			xn1[index + blockDim.x] = 0.00;
			xn2[index + blockDim.x] = 0.00;
		}
		
		return;
	}
}


__global__ void compute_points_Unstructure3D_1(double *posx, double *posy, double *posz, double *tloc, double *xn0, double *xn1, double *xn2, double *x, double *y, double *z, int *eid, double *r, double *s, double *t, VelData_double vel, Element MeshElementArray_device, double *h, int *integrate, int stage, int ss) {

	int tid;
	
	// Get thread index
	tid=(blockIdx.x*blockDim.x)+threadIdx.x; 
	if( tid >= ss) { // Redundant thread ...
		return; 
	} else {
		// Integrate only if this thread need to be integrated.
		if (integrate[tid] == 1) {
			double k[3];
		
			// Get velocity at current position given by posx, and posy.
			GetVel_unstruct3D(tloc[tid], k, vel, MeshElementArray_device, eid[tid], r[tid], s[tid], t[tid] );
			
			double mul;
			
			// Multiplier used for xn0, xn1 & xn2. See RK4 algorithm for more details.
			mul = (stage == 0) ? (0.1666666666666666666667*h[tid]) : (0.3333333333333333333334*h[tid]);
			
			// Update the xn buffer with the value of k1, k2 or k3
			xn0[tid] = xn0[tid] + mul*k[0];
			xn1[tid] = xn1[tid] + mul*k[1];
			xn2[tid] = xn2[tid] + mul*k[2];

			// Multiplier used for updated position of a particle. See RK4 algorithm for more details.
			mul = (stage == 2) ? (h[tid]) : (0.50*h[tid]);
		
			// Update the intermediate position of particle handled by this thread.
			posx[tid] = x[tid] + mul*k[0];
			posy[tid] = y[tid] + mul*k[1];
			posz[tid] = z[tid] + mul*k[2];
		
			
		
		}
	}

	return;
}



__global__ void compute_points_Unstructure3D_2(double *posx, double *posy, double *posz, double *tloc, double *xn0, double *xn1, double *xn2, double *x, double *y, double *z, int *eid, double *r, double *s, double *t, VelData_double vel, Element MeshElementArray_device, double *h, int *integrate, int ss) {

	int tid;
	
	// Get thread index
	tid=(blockIdx.x*blockDim.x)+threadIdx.x; 
	if( tid >= ss) { // Redundant thread ...
		return; 
	} else {
		
		// Integrate only if this thread need to be integrated.
		if (integrate[tid] == 1) {
			double k[3];
		
			// Get velocity at current position given by posx, and posy.
			GetVel_unstruct3D(tloc[tid], k, vel, MeshElementArray_device, eid[tid], r[tid], s[tid], t[tid] );
			
			// Add value of k4 into xn array
			xn0[tid] = xn0[tid] + 0.1666666666666666666666666667*h[tid]*k[0];
			xn1[tid] = xn1[tid] + 0.1666666666666666666666666667*h[tid]*k[1];
			xn2[tid] = xn2[tid] + 0.1666666666666666666666666667*h[tid]*k[2];
			
			// Update the position of a particle in the main struct
			x[tid] = x[tid] + xn0[tid];
			y[tid] = y[tid] + xn1[tid];
			z[tid] = z[tid] + xn2[tid];
		
			
		
		}
	}

	return;
}


// Kernels used by both 2D and 3D meshes


// This function determines if a thread needs to be integrated from time t1 to t2 and sets the flag array integrate to 0 or 1 accordingly.
// 0 - Do not integrate
// 1 - Integrate
// Algorithm: Dont integrate if start time for a particle is greater than t2 or stop time of a particle is smaller than t1. Integrate in all other cases.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__global__ void check_int_new(Point Tracer_dev, int *integrate, double *t1, double *t2, int ss, int offset) {

	int tid;
	
	// Get thread index
	tid=(2*blockIdx.x*blockDim.x)+threadIdx.x + offset; 
	
	if( tid >= ss) { // Redundant thread ...
		return; 
	} else {
	
		// Load Global memory into registers
		double tstart,tstop;
		double tstart1,tstop1;
		
		double t111, t221;
		double t112, t222;
		tstart = Tracer_dev.Start_time[tid];
		tstop = Tracer_dev.Stop_time[tid];
		t111 = t1[tid];
		t221 = t2[tid];
		if ((tid + blockDim.x) < ss) {
			tstart1 = Tracer_dev.Start_time[tid + blockDim.x];
			tstop1 = Tracer_dev.Stop_time[tid + blockDim.x];
			t112 = t1[tid + blockDim.x];
			t222 = t2[tid + blockDim.x];
		}
		
		// Determine if this thread needs to be integrated during time interval from t1 to t2.
		if ( (tstop < t111) || (tstart > t221) ) {
			integrate[tid] = 0;
		} else {
			integrate[tid] = 1;
		}
		
		if ((tid + blockDim.x) < ss) {
			
			if ( (tstop1 < t112) || (tstart1 > t222) ) {
				integrate[tid+ blockDim.x] = 0;
			} else {
				integrate[tid+ blockDim.x] = 1;
			}
		}
		
		return;
	}

}

