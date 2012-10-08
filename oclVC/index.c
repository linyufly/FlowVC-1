# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <CL/cl.h>
# include "globals.h"

# include "settings.h"




# include "structs.h"
# include "index.h"

inline int getindexwosshost(int i, int j, int k) {

	int out;
	out= i + (j * (Trace_CartMesh.XRes)) + (k * (Trace_CartMesh.XRes * Trace_CartMesh.YRes));
	return(out);
}



inline int getindexhost(int i, int j, int k) {

	int out;
	out= i + (j * Vel_CartMesh.XRes) + (k * (Vel_CartMesh.XRes * Vel_CartMesh.YRes));
	return(out);
}


