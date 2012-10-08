
# include "cartesian.h"

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


// This function initializes the variables stored in global memory (posx & posy ; xn0 & xn1).
// Approach used for instruction level parallelism: Global memory will be loaded into registers and then copied back to global memory.  
__global__ void initialize_Cart2D(Point Tracer_dev, double *posx, double *posy, double *xn0, double *xn1, int ss, int offset) {


	int tid;
	
	// Get thread index
	tid=(4*blockIdx.x*blockDim.x)+threadIdx.x; 
	
	if( tid >= ss) { // Redundant thread ...
		return; 
	} else {
		int index = tid + offset;
		
		// Load Global memory into registers
		double a0,b0;
		double a1,b1;
		double a2,b2;
		double a3,b3;
		
		a0 = Tracer_dev.x[index];
		b0 = Tracer_dev.y[index];
		
	
		if ((tid + blockDim.x) < ss) {
		
			a1 = Tracer_dev.x[index + blockDim.x];
			b1 = Tracer_dev.y[index + blockDim.x];	
		}
		
		if ((tid + (2*blockDim.x)) < ss) {
		
			a2 = Tracer_dev.x[index + (2*blockDim.x)];
			b2 = Tracer_dev.y[index + (2*blockDim.x)];	
		}
		
		if ((tid + (3*blockDim.x)) < ss) {
		
			a3 = Tracer_dev.x[index + (3*blockDim.x)];
			b3 = Tracer_dev.y[index + (3*blockDim.x)];	
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
		
		if ((tid + (2*blockDim.x)) < ss) {
		
			posx[index + (2*blockDim.x)] = a2;
			posy[index + (2*blockDim.x)] = b2;
			xn0[index + (2*blockDim.x)] = 0.00;
			xn1[index + (2*blockDim.x)] = 0.00;		
		}
		
		if ((tid + (3*blockDim.x)) < ss) {
		
			posx[index + (3*blockDim.x)] = a3;
			posy[index + (3*blockDim.x)] = b3;
			xn0[index + (3*blockDim.x)] = 0.00;
			xn1[index + (3*blockDim.x)] = 0.00;		
		}
		
		
		return;
	}
}


// This function initializes the variables stored in global memory (posx, posy & posz ; xn0, xn1 & xn2). 
// Approach used for instruction level parallelism: Global memory will be loaded into registers and then copied back to global memory. 
__global__ void initialize_Cart3D(Point Tracer_dev, double *posx, double *posy, double *posz, double *xn0, double *xn1, double *xn2, int ss, int offset) {


	int tid;
	int index;
	// Get thread index
	tid=(4*blockIdx.x*blockDim.x)+threadIdx.x; 
	
	if( tid >= ss) { // Redundant thread ...
		return; 
	} else {

		index = tid + offset;

		// Load Global memory into registers
		double a0,b0,c0;
		double a1,b1,c1;
		double a2,b2,c2;
		double a3,b3,c3;
		
		
		a0 = Tracer_dev.x[index];
		b0 = Tracer_dev.y[index];
		c0 = Tracer_dev.z[index];
		
	
		if ((tid + blockDim.x) < ss) {
		
			a1 = Tracer_dev.x[index + blockDim.x];
			b1 = Tracer_dev.y[index + blockDim.x];
			c1 = Tracer_dev.z[index + blockDim.x];
		}
		
		if ((tid + (2*blockDim.x)) < ss) {
		
			a2 = Tracer_dev.x[index + (2*blockDim.x)];
			b2 = Tracer_dev.y[index + (2*blockDim.x)];
			c2 = Tracer_dev.z[index + (2*blockDim.x)];	
		}
		
		if ((tid + (3*blockDim.x)) < ss) {
		
			a3 = Tracer_dev.x[index + (3*blockDim.x)];
			b3 = Tracer_dev.y[index + (3*blockDim.x)];
			c3 = Tracer_dev.z[index + (3*blockDim.x)];	
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
		
		if ((tid + (2*blockDim.x)) < ss) {
		
			posx[index + (2*blockDim.x)] = a2;
			posy[index + (2*blockDim.x)] = b2;
			posz[index + (2*blockDim.x)] = c2;
			xn0[index + (2*blockDim.x)] = 0.00;
			xn1[index + (2*blockDim.x)] = 0.00;
			xn2[index + (2*blockDim.x)] = 0.00;		
		}
		
		if ((tid + (3*blockDim.x)) < ss) {
		
			posx[index + (3*blockDim.x)] = a3;
			posy[index + (3*blockDim.x)] = b3;
			posz[index + (3*blockDim.x)] = c3;
			xn0[index + (3*blockDim.x)] = 0.00;
			xn1[index + (3*blockDim.x)] = 0.00;
			xn2[index + (3*blockDim.x)] = 0.00;		
		}
		
		
		return;
	}
}

// This kernel function computes the velocity of a particle and updates its position based on time-step in posx, posy and posz. This function is utilized to obtain the value of k1, k2 & k3 whose sum is stored in xn0, xn1 and xn2. stage input specifies the constant computing (k1, k2 or k3).

// 3D function
__global__ void compute_points_Cartesian3D_1(double *posx, double *posy, double *posz, double *tloc, double *xn0, double *xn1, double *xn2, double *x, double *y, double *z, int *ld, VelData_double vel, double *h, int *integrate, int stage, int ss) {

	int tid;
	
	// Get thread index
	tid=(blockIdx.x*blockDim.x)+threadIdx.x; 
	if( tid >= ss) { // Redundant thread ...
		return; 
	} else {
		
		// Integrate only if this thread need to be integrated.
		if (integrate[tid] == 1) {
			double k[3];
		
			// Get velocity at current position given by posx, posy and posz.
			GetVel_cartesian3D(posx[tid], posy[tid], posz[tid], tloc[tid], k, vel, ld[tid] );
			
			double mul;
			
			// Multiplier used for xn0, xn1 & xn2. See RK4 algorithm for more details.
			mul = (stage == 0) ? (0.1666666666666666666667*h[tid]) : (0.3333333333333333333334*h[tid]);
			
			xn0[tid] = xn0[tid] + mul*k[0];
			xn1[tid] = xn1[tid] + mul*k[1];
			xn2[tid] = xn2[tid] + mul*k[2];
			
			// Multiplier used for updated position of a particle. See RK4 algorithm for more details.
			mul = (stage == 2) ? (h[tid]) : (0.50*h[tid]);
		
			posx[tid] = x[tid] + mul*k[0];
			posy[tid] = y[tid] + mul*k[1];
			posz[tid] = z[tid] + mul*k[2];
		
			// Check if Point Left Domain
			ld[tid] = TestOutDomain3D_dev(posx[tid], posy[tid], posz[tid]);
		
		}
	}

	return;
}

// 2D function
__global__ void compute_points_Cartesian2D_1(double *posx, double *posy, double *tloc, double *xn0, double *xn1, double *x, double *y, int *ld, VelData_double vel, double *h, int *integrate, int stage, int ss) {

	int tid;
	
	// Get thread index
	tid=(blockIdx.x*blockDim.x)+threadIdx.x; 
	if( tid >= ss) { // Redundant thread ...
		return; 
	} else {
		// Integrate only if this thread need to be integrated.
		if (integrate[tid] == 1) {
			double k[2];
			
			// Get velocity at current position given by posx, posy and posz.
			GetVel_cartesian2D(posx[tid], posy[tid], tloc[tid], k, vel, ld[tid] );
			
			double mul;
			
			// Multiplier used for xn0 & xn1. See RK4 algorithm for more details.
			mul = (stage == 0) ? (0.1666666666666666666667*h[tid]) : (0.3333333333333333333334*h[tid]);
			
			xn0[tid] = xn0[tid] + mul*k[0];
			xn1[tid] = xn1[tid] + mul*k[1];
			
			// Multiplier used for updated position of a particle. See RK4 algorithm for more details.
			mul = (stage == 2) ? (h[tid]) : (0.50*h[tid]);
		
			posx[tid] = x[tid] + mul*k[0];
			posy[tid] = y[tid] + mul*k[1];		
			
			// Check if Point Left Domain
			ld[tid] = TestOutDomain_dev(posx[tid], posy[tid]);
		}
	}

	return;
}

///////////////////////////////////////////////////////////////////////////////////

// This function updates the position into the main tracer array. Essentially, this kernel computes the value of k4 and updates the main struct by values xn0, xn1 & xn2 respectively.
// 3D function
__global__ void compute_points_Cartesian3D_2(double *posx, double *posy, double *posz, double *tloc, double *xn0, double *xn1, double *xn2, double *x, double *y, double *z, int *ld, VelData_double vel, double *h, int *integrate, int ss) {

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
			GetVel_cartesian3D(posx[tid], posy[tid], posz[tid], tloc[tid], k, vel, ld[tid] );
			
			// Add value of k4 into xn array
			xn0[tid] = xn0[tid] + 0.1666666666666666666666666667*h[tid]*k[0];
			xn1[tid] = xn1[tid] + 0.1666666666666666666666666667*h[tid]*k[1];
			xn2[tid] = xn2[tid] + 0.1666666666666666666666666667*h[tid]*k[2];
	
			// Update the position of a particle in the main struct
			x[tid] = x[tid] + xn0[tid];
			y[tid] = y[tid] + xn1[tid];
			z[tid] = z[tid] + xn2[tid];
			
			// Check if Point Left Domain
			ld[tid] = TestOutDomain3D_dev(x[tid], y[tid], z[tid]);
			
		}
	}

	return;
}


// 2D function
__global__ void compute_points_Cartesian2D_2(double *posx, double *posy, double *tloc, double *xn0, double *xn1, double *x, double *y, int *ld, VelData_double vel, double *h, int *integrate, int ss) {

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
			GetVel_cartesian2D(posx[tid], posy[tid], tloc[tid], k, vel, ld[tid] );
			
			// Add value of k4 into xn array
			xn0[tid] = xn0[tid] + 0.1666666666666666666666666667*h[tid]*k[0];
			xn1[tid] = xn1[tid] + 0.1666666666666666666666666667*h[tid]*k[1];			
		
			// Update the position of a particle in the main struct
			x[tid] = x[tid] + xn0[tid];
			y[tid] = y[tid] + xn1[tid];
			
			// Check if Point Left Domain
			ld[tid] = TestOutDomain_dev(x[tid], y[tid]);
		
		}

	}

	return;
}




