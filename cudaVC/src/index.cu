

# include "index.h"

inline int getindexwosshost(int i, int j, int k) {

	int out;
	out= i + (j * (Trace_CartMesh.XRes)) + (k * (Trace_CartMesh.XRes * Trace_CartMesh.YRes));
	return(out);
}

__device__ inline int getindexwoss(int i, int j, int k) {

	int out;
	out= i + (j * (Trace_CartMesh.XRes)) + (k * (Trace_CartMesh.XRes * Trace_CartMesh.YRes));
	
	return(out);
}

inline int getindexhost(int i, int j, int k) {

	int out;
	out= i + (j * Vel_CartMesh.XRes) + (k * (Vel_CartMesh.XRes * Vel_CartMesh.YRes));
	return(out);
}

__device__ inline int getindex(int i, int j, int k) {

	
	int out;
	out= i + (j * RES[0]) + (k * (RES[0] * RES[1]));
	return(out);


}
