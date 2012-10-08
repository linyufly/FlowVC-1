
# include <stdio.h>
# include <math.h>
# include <stdlib.h>
# include <string.h>
#include <stdarg.h>

//# inlude "globals.h"
//# include "parameters.h"

# define LONGSTRING 200
double 					X1[3];
double					X2[3];

int Suffix_TMin, Suffix_TMax, Suffix_TDelta, NumFiles, Data_MeshType;

// This variables stores the number of points loaded onto kernel function for global search purpose.
int 					Number_Points_loaded;

double 					*norm;

// Globally defined file names 

char					Path_Data1[200];

char 					Trace_OutFilePrefix1[200];
char					Path_Output[200];

char					Path_Data2[200];

char 					Trace_OutFilePrefix2[200];




double					Data_Loaded_TMin;
double 					Data_Loaded_TMax;


// Globals related to how much output to generate 

double					Output_TStart;
double					Output_TEnd; 
double					Output_TDelta;
int						Output_TRes;

void FatalError(char *text, ...);
void ReadIn(int ss);
void ReadInNextValue(FILE *Parameters_InFileID, void *pt, char type);
void ReadInParameters(int argc, const char *argv[]);
void CheckParameters(void);
void SetDerivedParameters(void);


int main(int argc, const char *argv[]) {

	// Read in the data file, check the settings and calculated some of the derived parameters....	
	ReadInParameters(argc,argv);
	CheckParameters(); 
	SetDerivedParameters();
	
	// Initialization complete.
	
	int i;
	
	// Allocating memory for norm data.
	if((norm = (double *)malloc(NumFiles*sizeof(double))) == NULL) {
		fprintf(stderr, "Malloc failed for norm buffer \n");
		exit(1);
	}
	
	char OutFile[200];
	FILE *fileout;
	double maxnorm1 = 0.00;
	
	// Open output file
	sprintf(OutFile, "%snorm_data.vtk", Path_Output);
	
	if((fileout = fopen(OutFile, "w")) == NULL) {
		fprintf(stderr, "Could not open file %s \n", OutFile);
		exit(1);
	}
	
	// First file is the initial position which will be the same
	for (i = 1; i < NumFiles; i++) {
	
		// Read In data from the files.
		ReadIn(i);
		
		fprintf(fileout, " %0.9f \n", norm[i]);
		
		printf("Norm for file %d is %0.13f \n", i, norm[i]);
		
		if (norm[i] > maxnorm1)
			maxnorm1 = norm[i];
		
	}
	
	printf("\nAll Norm Data Printed to %s \n\n", OutFile);
	printf("Maximum Norm was found to be %0.13f \n\n", maxnorm1);
	fclose(fileout);
	free(norm);
	


	return 0;
}



void ReadIn(int ss) {

	char InFile1[200],InFile2[200] ;
	int j;
	FILE *filein1, *filein2;
	int num_tracers;
	
	double *norm_temp;
	
	// Open flowVC output file
	sprintf(InFile1, "%s%s.%d.bin",  Path_Data1, Trace_OutFilePrefix1, Suffix_TMin + ss * Suffix_TDelta);
	if((filein1 = fopen(InFile1, "rb")) == NULL) {
		fprintf(stderr, "Could not open file %s \n", InFile1);
		exit(1);
	}
	
	// Open CudaVC output filw
	sprintf(InFile2, "%s%s_out.%d.bin",  Path_Data2, Trace_OutFilePrefix2, ss);
	if((filein2 = fopen(InFile2, "rb")) == NULL) {
		fprintf(stderr, "Could not open file %s \n", InFile2);
		exit(1);
	}
	
	double time1, time2;
	
	// Read time stamp and number of points
    if(fread(&time1, sizeof(double), 1, filein1) < 1)
    	FatalError("Could not read time stamp from file %s", InFile1);
	
	// Read time stamp and number of points
    if(fread(&time2, sizeof(double), 1, filein2) < 1)
    	FatalError("Could not read time stamp from file %s", InFile2);
      
	if (time1 != time2)
		FatalError("Time stamps not equal \n\n");
      
	// Read data from cudaVC output file to determine how many points there are
	if(fread(&num_tracers, sizeof(int), 1, filein2) < 1)
		FatalError("Could not read number of points from file %s", InFile2);
	
	
	// Allocating memory for norm data.
	if((norm_temp = (double *)malloc(num_tracers*sizeof(double))) == NULL) {
		fprintf(stderr, "Malloc failed for norm_temp buffer \n");
		exit(1);
	}
	
	
	printf("Calculating Error for %d points ... \n", num_tracers);
	
	// Get maximum norm out of all points.
	double maxnorm = 0.000;
	
	for(j = 0; j < num_tracers; j++) {
    
    	// Read a point from flowVC
    
    	if(fread(X1, sizeof(double), 3, filein1) < 3)  {
			printf("Could not read location of tracer %d of %d from %s", j + 1, num_tracers, InFile1);
			exit(0);
		}
		
		// Read a point from cudaVC
		
		if(fread(X2, sizeof(double), 3, filein2) < 3)  {
			printf("Could not read location of tracer %d of %d from %s", j + 1, num_tracers, InFile2);
			exit(0);
		}
		//printf("x1 = %f, y1 = %f and z1 = %f \n",X1[0], X1[1], X1[2]);
		//printf("x2 = %f,norm_temp[j] y2 = %f and z2 = %f \n",X2[0], X2[1], X2[2]);
	
		
		norm_temp[j] = sqrt( ((X1[0] - X2[0]) * (X1[0] - X2[0])) +  ( (X1[1] - X2[1]) * (X1[1] - X2[1]) ) + ( (X1[2] - X2[2]) * (X1[2] - X2[2]) )  );
		if (norm_temp[j] > maxnorm)
			maxnorm = norm_temp[j];
		//printf("Norm for Point %d is = %f \n", j, norm_temp[j]);
		
	}
	
	fclose(filein1);
	fclose(filein2);
	
	//printf("\nMaxnorm = %f \n\n", maxnorm);
	norm[ss] = maxnorm;
	free(norm_temp);
}




void ReadInParameters(int argc, const char *argv[]) {

	FILE *Parameters_InFileID;
	
	if(argc == 2) {
		if((Parameters_InFileID = fopen(argv[1], "r")) == NULL) 
			FatalError("Could not open input file %s", argv[1]);
		
		fprintf(stderr, "Reading in parameters... \n\n");

		ReadInNextValue(Parameters_InFileID, Path_Data1, 's');
		ReadInNextValue(Parameters_InFileID, Path_Data2, 's');
		
		ReadInNextValue(Parameters_InFileID, Path_Output, 's');

		ReadInNextValue(Parameters_InFileID, &Data_MeshType, 'd'); 
		
		ReadInNextValue(Parameters_InFileID, &Suffix_TMin, 'd');    
		ReadInNextValue(Parameters_InFileID, &Suffix_TMax, 'd');    
		ReadInNextValue(Parameters_InFileID, &Suffix_TDelta, 'd');  

		
		ReadInNextValue(Parameters_InFileID, Trace_OutFilePrefix1, 's');
		ReadInNextValue(Parameters_InFileID, Trace_OutFilePrefix2, 's'); 
		
		

		fclose(Parameters_InFileID);

		fprintf(stderr,"Path_Data1\t\t=\t%10s\n",Path_Data1);
		fprintf(stderr,"Path_Data2\t\t=\t%10s\n",Path_Data2);
		
		fprintf(stderr,"Path_Output1\t\t=\t%10s\n\n",Path_Output);
		
		
		
		fprintf(stderr,"Data_MeshType\t\t=\t%10d\n",Data_MeshType);       

		fprintf(stderr,"Suffix_TMin\t\t=\t%10d\n",Suffix_TMin);
		fprintf(stderr,"Suffix_TMax\t\t=\t%10d\n",Suffix_TMax);
		fprintf(stderr,"Suffix_TDelta\t\t=\t%10d\n\n",Suffix_TDelta);

		fprintf(stderr,"Trace_OutFilePrefix1\t=\t%10s\n",Trace_OutFilePrefix1);   
		fprintf(stderr,"Trace_OutFilePrefix2\t=\t%10s\n",Trace_OutFilePrefix2);


	}
	else {
		fprintf(stderr, "Please use input file .... \n\n Usage: %s inputfile ...\n", argv[0]);
		exit(1);
	}

}

void ReadInNextValue(FILE *Parameters_InFileID, void *pt, char type) {
	
	int c;
	char buf[LONGSTRING];
	
	while(1) {
		c = fgetc(Parameters_InFileID);
		
		if(c == ' ' || c == '\n' || c == '\r') {
			continue;
		}
		
		if(c == '#') {
			/* Comment line, read until end of line */
			while(c != '\n') {
				c = fgetc(Parameters_InFileID);
			}
		}
		else {
			/* Line should contain variable information of form VARIABLE = VALUE*/
			ungetc(c, Parameters_InFileID);
			fscanf(Parameters_InFileID, "%s %*s", buf);
			if(type == 's')
				fscanf(Parameters_InFileID,"%s\n", (char *) pt);
			else if(type == 'd')
				fscanf(Parameters_InFileID,"%d\n", (int *) pt);
			else if(type == 'f')
				fscanf(Parameters_InFileID,"%lf\n", (double *) pt);
			else 
				printf("Unsupported data type passed to ReadInNextValue().\n");
			break;
		}
	}
}

void CheckParameters(void) {
	
	printf("Checking parameters...\n");
	
	
	
	printf("OK!\n\n");
	fflush(stdout);
	
}




void SetDerivedParameters(void) {
	
	
	//double Data_TMax_Req ;

	printf("Setting derived parameters...\n");
	
	if(strcmp(Path_Output, "pwd")) {
		if(Path_Output[strlen(Path_Output) - 1] != '/') 
			sprintf(Path_Output, "%s/", Path_Output);
	}
	else 
		sprintf(Path_Output, "./");
	
	if(strcmp(Path_Data1, "pwd")) {
		if(Path_Data1[strlen(Path_Data1) - 1] != '/' ) 
			sprintf(Path_Data1, "%s/", Path_Data1);
	}
	else 
		sprintf(Path_Data1, "./");
		
	if(strcmp(Path_Data2, "pwd")) {
		if(Path_Data2[strlen(Path_Data2) - 1] != '/' ) 
			sprintf(Path_Data2, "%s/", Path_Data2);
	}
	else 
		sprintf(Path_Data2, "./");
		
	if(Suffix_TMax != Suffix_TMin) {
		if(Suffix_TDelta < 1)
		  FatalError("Suffix_TDelta < 1");
		NumFiles = (Suffix_TMax - Suffix_TMin) / Suffix_TDelta + 1;
	} else
		NumFiles = 1;
	
	printf("Checking Error from %d files \n", NumFiles);
	
	printf("OK!  \n\n");
	
}

void FatalError(char *text, ...) {
	va_list ap;
	char *p, *sval;
	int ival;
	double dval;
	
	fprintf(stderr, "\nERROR:\n");
	va_start(ap, text);
	for(p = text; *p; p++) {
		if(*p != '%') {
			putc(*p, stderr);
			continue;
		}
		switch (*++p) {
			case 'd':
				ival = va_arg(ap, int);
				fprintf(stderr, "%d", ival);
				break;
			case 'f':
				dval = va_arg(ap, double);
				fprintf(stderr, "%f", dval);
				break;
			case 's':
				for(sval = va_arg(ap, char *); *sval; sval++)
					putc(*sval, stderr);
				break;
			default:
				putc(*p, stderr);
				break;
		}
	}
	va_end(ap);
	
	fprintf(stderr, "\n");
	fflush(stderr);
	
	exit(1);
}


