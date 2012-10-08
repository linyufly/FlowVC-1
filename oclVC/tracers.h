#ifndef INC_TRACERS_H
#define INC_TRACERS_H

#include "structs.h"


extern void GenerateStaggeredRelease(void);
extern void LoadReleasePoints(double *voxdim);

extern void release(int ss);
//extern void initialize(void);
extern void initialize_new(void);

extern int Launch_Tracer(double tmin, double tmax, int N_Launches);


extern int Output_Tracer(double tmin, double tmax, int Output_Release);





#endif
