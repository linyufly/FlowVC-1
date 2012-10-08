# include <stdio.h>
# include <stdlib.h>

# include "settings.h"
# include "structs.h"
# include "globals.h"
# include "fileoutput.h"

void temp2file (void) {

	char Data_OutFilePath[LONGSTRING];
	
	FILE *fileout;
	
	sprintf(Data_OutFilePath, "%s%s.bin", Path_Output,Temp_OutFilePrefix);
	
	if((fileout = fopen(Data_OutFilePath, "wb")) == NULL) {
		fprintf(stderr, "Could not open file %s \n", Data_OutFilePath);
		exit(1);
	}
	
	int i;
	
	if(fwrite(&N_Frame, sizeof(int), 1, fileout) < 1)
		printf("Could not write Trace_NumTracers to %s", Data_OutFilePath);
	for (i = 0; i < N_Frame; i++) {
	// Write Tracer information to file
		if((int)fwrite(&Tracer.x[i], sizeof(double), 1, fileout) < 1)
			printf("Could not completely write TraceIC_Array data to %s", Data_OutFilePath);

		if((int)fwrite(&Tracer.y[i], sizeof(double), 1, fileout) < 1)
			printf("Could not completely write TraceIC_Array data to %s", Data_OutFilePath);	

		if((int)fwrite(&Tracer.z[i], sizeof(double), 1, fileout) < 1)
			printf("Could not completely write TraceIC_Array data to %s", Data_OutFilePath);

		if(fwrite(&Tracer.ElementIndex[i], sizeof(int), 1, fileout) < 1)
			printf("Could not completely write TraceIC_Array data to %s", Data_OutFilePath);

		if(fwrite(&Tracer.LeftDomain[i], sizeof(int), 1, fileout) < 1)
			printf("Could not completely write TraceIC_Array data to %s", Data_OutFilePath);
	}		
	/*
	
		fprintf(fileout, " %0.9f %0.9f %0.9f %d %d \n", temp.x[i], temp.y[i], temp.z[i], temp.ElementIndex[i], temp.LeftDomain[i]);	
	
	}
	*/
	
	

	fclose(fileout);

}

/*
void copyresult2file(int ss) {

	char Data_OutFilePath[LONGSTRING];
	
	FILE *fileout;


	
	
	sprintf(Data_OutFilePath, "%soutput.%d.data",Output_Data, ss);
	if((fileout = fopen(Data_OutFilePath, "w")) == NULL) {
		fprintf(stderr, "Could not open file %s \n", Data_OutFilePath);
		exit(1);
	}
	
	fprintf(fileout,"Points are Integrated from %f to %f \n\n\n",output[ss*N].t.Start_time,output[ss*N].t.Stop_time);
	fprintf(fileout,"x(old)\t\t  x(new)\t\t\ty(old)\t\ty(new) \n\n\n");
	double xold[2];	
	int k;
	for (k=0;k<N;k++)
	 {   
		xold[0] = MIN_X+((double)XDelta*(int)(k%(XRes)));
  		xold[1] = MIN_Y+((double)YDelta*(int)(k/(XRes))); 		
		fprintf(fileout, "%f \t %f\t\t %f\t%f \n",xold[0],output[k + ss*N].X[0],xold[1],output[k + ss*N].X[1]);
	}
	
	 fclose(fileout); 

} */

void copyresult2bin(int ss, int N, double t) {

	char Data_OutFilePath[LONGSTRING];
	
	FILE *fileout;
	int sum = 0;
	double X[3];
	int num_tracers;
	
	if ( Trace_ReleaseStrategy == 0 ) 	
		num_tracers = N*N_Frame;
	else if ( Trace_ReleaseStrategy == 1 )
		num_tracers = N; 
		
	// Check for total number of points to output.
	for (int k = 0; k < num_tracers; k++) {	
	
		if (Data_MeshType == UNSTRUCTURED) {
			// Check element index for unstructured
			if (Tracer.ElementIndex[k] != -1)
				sum++;
		} else {
				// For catresian check if a point is in the domain.
			if (Tracer.LeftDomain[k] != 1)
				sum++;
							
		}
	
	}
	if (sum > 0) {
		sprintf(Data_OutFilePath, "%s%s_out.%d.bin", Path_Output, Trace_OutFilePrefix , ss);
	
		if((fileout = fopen(Data_OutFilePath, "wb")) == NULL) {
			fprintf(stderr, "Could not open file %s \n", Data_OutFilePath);
			exit(1);
		}
		
		if(fwrite(&t, sizeof(double), 1, fileout) < 1)
			FatalError("Could not write timestamp to %s ", Data_OutFilePath );
		
		if(fwrite(&sum, sizeof(int), 1, fileout) < 1) 
			FatalError("Could not write total number of points to %s ", Data_OutFilePath );
		
		
			
		
		for (int k = 0; k < num_tracers; k++) {	
		
			X[0] = Tracer.x[k];
			X[1] = Tracer.y[k];
			X[2] = Tracer.z[k];
		
			if (Data_MeshType == UNSTRUCTURED) {
	
				if(Tracer.ElementIndex[k] >= 0) {
					if(fwrite(X, sizeof(double), 3, fileout) < 3) 
						FatalError("Could not completely write co-ordinates of %d point data to %s", k, Data_OutFilePath);
				}
		
			} else { // Cartesian Grid
		
				if (!Tracer.LeftDomain[k]) {
				
					if(fwrite(X, sizeof(double), 3, fileout) < 3) 
						FatalError("Could not completely write co-ordinates of %d point data to %s", k, Data_OutFilePath);
				
				}
			}
		}

		fclose(fileout);
	}
}

void copyresult2vtk_staggered(int ss, int N, double t, double tmin, double tmax) {

	char Data_OutFilePath[LONGSTRING];
	
	FILE *fileout;
	int sum = 0;
	
	for (int k = 0; k < N; k++) {	
	
		if (Data_MeshType == UNSTRUCTURED) {
			// Check element index for unstructured
			if( (Tracer.Start_time[k] > tmax) || (Tracer.Stop_time[k] < tmin ) ) {
				continue;
			} else {
				if (Tracer.ElementIndex[k] != -1)
					sum++;
			}
		} else {
				// For catresian check if a point is in the domain.
			if((Tracer.Start_time[k] > tmax) || (Tracer.Stop_time[k] < tmin ))	{
				continue;
			} else {
				if (Tracer.LeftDomain[k] != 1)
					sum++;
			}		
		}

	}
	N_Present = sum;
	
//	printf("num = %d particles outputted \n\n", sum);
	
	

		sprintf(Data_OutFilePath, "%s%s_out.%d.vtk", Path_Output, Trace_OutFilePrefix , ss);
		if((fileout = fopen(Data_OutFilePath, "w")) == NULL) {
			fprintf(stderr, "Could not open file %s \n", Data_OutFilePath);
			exit(1);
		}

			fprintf(fileout, "# vtk DataFile Version 3.0\n");
			fprintf(fileout, "t = %.9f\n", t);
			fprintf(fileout, "ASCII\n");
			fprintf(fileout, "DATASET POLYDATA\n");
			fprintf(fileout, "\nPOINTS %d double\n", N_Present);
			
	if (sum > 0) {

		for (int k = 0; k < N; k++) {	
	
			if (Data_MeshType == UNSTRUCTURED) {
	
				if(Tracer.ElementIndex[k] >= 0) {
					
					if((Tracer.Start_time[k] > tmax) || (Tracer.Stop_time[k] < tmin )) {
						//continue;
					} else {
						fprintf(fileout, " %0.9f %0.9f %0.9f \n", Tracer.x[k], Tracer.y[k], Tracer.z[k]);
					}
				}
		
			} else { // Cartesian Grid
		
				if (!Tracer.LeftDomain[k]) {
					//printf("%0.9f %0.9f %0.9f \n", Tracer.x[k], Tracer.y[k], Tracer.z[k]);
					if((Tracer.Start_time[k] > tmax) || (Tracer.Stop_time[k] < tmin )) {
						//continue;
					} else {
						fprintf(fileout, " %0.9f %0.9f %0.9f \n", Tracer.x[k], Tracer.y[k], Tracer.z[k]);	
					}
				}
			}
	
		}

		fclose(fileout);

	}


}

void copyresult2vtk(int ss, int N, double t) {

	char Data_OutFilePath[LONGSTRING];
	
	FILE *fileout;
	int sum = 0;
	//for (int k = 0; k < (N*N_Frame); k++) {	
	
		//printf("x = %0.9f  y = %0.9f \n",Tracer.x[k], Tracer.y[k] );
	//}
	int num_tracers;
	
	if ( Trace_ReleaseStrategy == 0 ) 	
		num_tracers = N*N_Frame;
	else if ( Trace_ReleaseStrategy == 1 ) {
		num_tracers = N; 
	}
	
	for (int k = 0; k < num_tracers; k++) {	
	
		if (Data_MeshType == UNSTRUCTURED) {
			// Check element index for unstructured
			if (Tracer.ElementIndex[k] != -1)
				sum++;
		} else {
				// For catresian check if a point is in the domain.
			if (Tracer.LeftDomain[k] != 1)
				sum++;
							
		}

	}
	N_Present = sum;
//printf("num = %d \n\n", N_Present);

	if (sum > 0) {

		sprintf(Data_OutFilePath, "%s%s_out.%d.vtk", Path_Output, Trace_OutFilePrefix , ss);
		if((fileout = fopen(Data_OutFilePath, "w")) == NULL) {
			fprintf(stderr, "Could not open file %s \n", Data_OutFilePath);
			exit(1);
		}

			fprintf(fileout, "# vtk DataFile Version 3.0\n");
			fprintf(fileout, "t = %.9f\n", t);
			fprintf(fileout, "ASCII\n");
			fprintf(fileout, "DATASET POLYDATA\n");
			fprintf(fileout, "\nPOINTS %d double\n", N_Present);
			


		for (int k = 0; k < (num_tracers); k++) {	
	
			if (Data_MeshType == UNSTRUCTURED) {
	
				if(Tracer.ElementIndex[k] >= 0)
					fprintf(fileout, " %0.9f %0.9f %0.9f \n", Tracer.x[k], Tracer.y[k], Tracer.z[k]);
		
		
			} else { // Cartesian Grid
		
				if (!Tracer.LeftDomain[k]) {
					//printf("%0.9f %0.9f %0.9f \n", Tracer.x[k], Tracer.y[k], Tracer.z[k]);
					fprintf(fileout, " %0.9f %0.9f %0.9f \n", Tracer.x[k], Tracer.y[k], Tracer.z[k]);	
				}
			}
	

		//	if((Tracer.ElementIndex[k] >= 0)  || (!Tracer.LeftDomain)) {
	//			fprintf(fileout, " %0.9f %0.9f %0.9f \n", Tracer.x[k], Tracer.y[k], Tracer.z[k]);	
		//	}
		}

		fclose(fileout);

	}

}
