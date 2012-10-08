#ifndef INC_GLOBALSEARCH_H
#define INC_GLOBALSEARCH_H

# include "structs.h"

# include "settings.h"
 //# include "globals.cu"



extern __global__ void GlobalSearch2D(Element MeshElementArray_device, Node MeshNodeArray_device, int Num_Frame, int *outputid_dev, int Vel_MeshNumElements, int ss, int N_Frame);

extern __global__ void GlobalSearch3D(Element MeshElementArray_device, Node MeshNodeArray_device, int Num_Frame, int *outputid_dev, int Vel_MeshNumElements, int ss, int N_Frame);

extern void globalsearch_device(void);

extern int Get_Element_Global_Search(const double *X);

extern int Get_Element_Local_Search(const double *X, int guess);

extern void search_host(void);
#endif
