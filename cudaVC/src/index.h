#ifndef INC_INDEX_H
#define INC_INDEX_H

# include "structs.h"



extern inline int getindexwosshost(int i, int j, int k);
extern __device__ inline int getindexwoss(int i, int j, int k);

extern inline int getindexhost(int i, int j, int k);
extern __device__ inline int getindex(int i, int j, int k);



#endif
