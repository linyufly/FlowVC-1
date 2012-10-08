#ifndef INC_GLOBALSEARCH_H
#define INC_GLOBALSEARCH_H

# include "structs.h"

# include "settings.h"
 //# include "globals.cu"


extern int Get_Element_Global_Search(const double *X);

extern int Get_Element_Local_Search(const double *X, int guess);

//extern int TestOutsideDomain(double point[3]);

extern void search_host(void);
#endif
