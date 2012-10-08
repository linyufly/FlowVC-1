#ifndef INC_MEMCOPY_H
#define INC_MEMCOPY_H

# include "structs.h"

// Memcopy common variables which will remain constant throughout ....
extern void host2device_const(void);

// Memcopy functions for global search kernel function
extern void host2device_gs(int ss);
extern void device2host_gs(int ss);
extern void host2const_mem_gs(void);

// Memcopy functions for main kernel function ....
extern void host2device(int ss);
extern void device2host(int ss);


#endif

