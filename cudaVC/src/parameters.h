#ifndef INC_PARAMETERS_H
#define INC_PARAMETERS_H

# include "settings.h"
#include <string.h>


extern void ReadInNextValue(FILE *Parameters_InFileID, void *pt, char type);
extern void ReadInParameters(int argc, const char *argv[]);
extern void CheckParameters(void);
extern void SetDerivedParameters(void);

#endif
