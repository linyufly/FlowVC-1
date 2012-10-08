/*
 *  main.c
 *  flowVC
 *
 *  Created by Shawn Shadden.
 *  Copyright 2010 Illinois Institute of Technology. All rights reserved.
 *
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <omp.h>
#include "exposuretime.h"
#include "ftle.h"
#include "globals.h"
#include "integration.h"
#include "io.h"
#include "macros.h"
#include "memory.h"
#include "mesh.h"
#include "parameters.h"
#include "residencetime.h"
#include "structs.h"
#include "strainrate.h"
#include "tracers.h"
#include "velocity.h"
#include "velout.h"



/* main function */
int main (int argc, const char * argv[]) {
    double runtime, nextoutputtime;
	int    df, ss, ee;
	char   BinFile[LONGSTRING];
	time_t SimulationStartTime, InitializationTime, SimulationStopTime;
	FILE   *Trace_BinFileID;
	FileData *Trace_BinFileArray;
	//int ii,num_integrated;
	//double t1,t2;
	//LagrangianPoint TempPoint;  
	//ReleaseLocation *rnode;
		
	/* Get simulation start time so we can time how long simulation runs */
	time(&SimulationStartTime);
	
	/* Input paramenters */
	ReadInParameters(argc, argv);
	CheckParameters(); 
	SetDerivedParameters();
	
	/* Load mesh */
	LoadMeshData();    
	
	/* If performing a staggered tracer release, compute release times at each release location based on velocity data */
	if(Trace_Compute && Trace_ReleaseStrategy == STAGGERED) 
		GenerateStaggeredRelease(); 
	
	/* Load velocity for first frame needed to start integration */ 
	if(Data_MeshType == CARTESIAN)  
		LoadCartVelDataFrame(Data_FirstFrame);
	else if(Data_MeshType == UNSTRUCTURED) {
		if(Int_NormalFlow) 
			LoadMeshSurfaceNormals(); /* Mesh normals needed to impose inflow along mesh walls */
		LoadUnstructVelDataFrame(Data_FirstFrame);
	}
	
	/* If interpolating velocity fields onto output mesh, initialize output mesh */
	if(VelOut_Compute) {
		if(VelOut_GenerateMesh) 
			GenerateVelOutMesh(); 
		else
			ReadInVelOutMesh();
	}
	if(FTLE_Compute) {
		InitializeFTLEArray();		// 	***************	OPENMP IMPLEMENTATION Possible************
		Trace_BinFileID = NULL;
		Trace_BinFileArray = NULL;
	}
	else if(Trace_Compute) {
		/* If computing exposure time data, allocate memory and initialize exposure time array */
		if(Trace_CETCompute) {
			if((Trace_CETArray = (CETE *)malloc(CET_MeshNumElements * sizeof(CETE))) == NULL)
				FatalError("Malloc failed for Trace_CETArray");
			for(ee = 0; ee < CET_MeshNumElements; ee++) {
				Trace_CETArray[ee].CETsum = 0; /* Cummulative exposure time for element ee */
				Trace_CETArray[ee].Encounters = 0; /* Numer of encounters for element ee */
			}
		} 
		if(Trace_ReleaseStrategy == STAGGERED) {
			Trace_BinFileID = NULL;
			/* Generate an output file for each output time */
			if((Trace_BinFileArray = (FileData *)malloc(Output_TRes * sizeof(FileData))) == NULL) 
				FatalError("Malloc failed for Trace_BinFileArray"); 
			for(ss = 0; ss < Output_TRes; ss++) {
				/* Open file */
				sprintf(Trace_BinFileArray[ss].FilePath, "%s%s.%d.bin", Path_Output, Trace_OutFilePrefix, ss);
				if((Trace_BinFileArray[ss].FileID = fopen(Trace_BinFileArray[ss].FilePath, "wb")) == NULL)
					FatalError("Could not open file %s", Trace_BinFileArray[ss].FilePath);
				/* Write output time as first entry in file */
				nextoutputtime = Output_TStart + ss * Output_TDelta;
				if(fwrite(&nextoutputtime, sizeof(double), 1, Trace_BinFileArray[ss].FileID) < 1) 
					FatalError("Could not write output time to %s", Trace_BinFileArray[ss].FilePath);
			} 
		}
		else {
			Trace_BinFileArray = NULL;
			
			/* Set up mesh of tracers to be integrated */
			if(!Trace_GenerateMesh)
				ReadInTraceMesh(); // 	***************	OPENMP IMPLEMENTATION Possible ************
			else 
				GenerateTracerMesh(); //	***************	OPENMP IMPLEMENTATION Possible ************
			
			if (TRACER_LOWERLOOP) {
			
				/* Open binary file for tracer output */
				if((Trace_BinFileID = fopen(TraceOUT_BinFilePath, "wb")) == NULL) 
					FatalError("Could not open file %s", TraceOUT_BinFilePath);
			
			}
				
			/* Initialize output counters */
			outputframe = 0;
			nextoutputtime = Output_TStart;
			for(ss = 0; ss < Trace_NumLaunchTimes; ss++)  
				CreateNewTraceLaunch(ss);
		}
		if(Int_NormalFlow && Particle_Radius > TINY)
			LoadBoundaryElementFlags();
			
		/* Allocate memory for array to store number of tracers output each output time frame for each release */
		if((Trace_NumOutput_OMP = (int **)calloc(Output_TRes, sizeof(int*))) == NULL)
			FatalError("Calloc failed for Trace_NumOutput_OMP");

		/* Allocate memory for array to store number of tracers output each output time frame (This number is used for stagerred release only) */
		if((Trace_NumOutput = (int *)calloc(Output_TRes, sizeof(int))) == NULL)
			FatalError("Calloc failed for Trace_NumOutput");	
			
		for (ss = 0; ss < Output_TRes; ss++) {
				
			if((Trace_NumOutput_OMP[ss] = (int *)calloc(Trace_NumLaunchTimes, sizeof(int))) == NULL)
				FatalError("Calloc failed for Trace_NumOutput_OMP1[%d]", ss);
		}
		
		/* If computing activation potential (integrated strain rate), then load first strain rate field */ 
		if(Trace_APCompute) {
			if(Data_MeshType == CARTESIAN)  
				LoadCartStrainRateDataFrame(Data_FirstFrame);
			else
				LoadUnstructStrainRateDataFrame(Data_FirstFrame);
		}  
	}
	else {
		Trace_BinFileID = NULL;
		Trace_BinFileArray = NULL;
	}
	
	/* Output initialization time */
	time(&InitializationTime);
	runtime = difftime(InitializationTime, SimulationStartTime);
	printf("\nInitialization time: %d hr %d min %d s\n", (int)(runtime/3600.0), 
		   (int)(fmod(runtime, 3600.0)/60.0), (int)(fmod(fmod(runtime, 3600.0), 60.0)));  
	fflush(stdout);
	
	
	/*** MAIN LOOP ***/
	for(df = Data_FirstFrame; df != Data_LastFrame; df = df + Int_TimeDirection) { 
	
	
		/* Load next velocity data frame */
		if(Data_MeshType == CARTESIAN) 
			LoadCartVelDataFrame(df + Int_TimeDirection);
		else if(Data_MeshType == UNSTRUCTURED) 
			LoadUnstructVelDataFrame(df + Int_TimeDirection);
		
		/* If computing activation potential, load strain rate data */
		if(Trace_Compute && Trace_APCompute) {
			if(Data_MeshType == CARTESIAN) 
				LoadCartStrainRateDataFrame(df + Int_TimeDirection);
			else if(Data_MeshType == UNSTRUCTURED) 
				LoadUnstructStrainRateDataFrame(df + Int_TimeDirection);
		}
		
		/* Determine time interval of loaded data */
		if(Int_TimeDirection > 0) {
			Data_LoadedTMin = Data_TMin + df * Data_TDelta;
			Data_LoadedTMax = Data_LoadedTMin + Data_TDelta;
		}
		else {
			Data_LoadedTMin = Data_TMin + (df - 1)* Data_TDelta;
			Data_LoadedTMax = Data_LoadedTMin + Data_TDelta;
		}    
		printf("Data loaded in memory spans t = %g to %g\n", Data_LoadedTMin, Data_LoadedTMax);
		
		fflush(stdout);
		/**************************************** Perform Computations *******************************************/
		/* If requested, interpolate velocity data to output mesh */
		if(VelOut_Compute) 
			for(ss = 0; ss < Output_TRes; ss++)  
			/* Check if current output time for ss is in time interval of loaded data and that output hasn't already been completed */
				if((VelOut_OutputTime[ss] - Data_LoadedTMin > -TINY) && (VelOut_OutputTime[ss] - Data_LoadedTMax < TINY) && !VelOut_Complete[ss]) {
					OutputVelOut(VelOut_OutputTime[ss], ss);
					/* OutputStrainOut(VelOut_OutputTime[ss], ss); */
					VelOut_Complete[ss] = 1;
				}
		
        /* If requested, compute FTLE */
         if(FTLE_Compute) {
         	
         	ComputeFTLE(df);         
         
		} else if(Trace_Compute) {
				
			if(Trace_ReleaseStrategy == STAGGERED) {
			
				Tracer_StaggeredRelease(Trace_BinFileID, Trace_BinFileArray);
		
			} else { // Normal release .....
			
				nextoutputtime = Compute_Tracer(&outputframe, nextoutputtime, Trace_BinFileID);
			}
		}
		
	}	
	
	
	/* Do any additional clean up */  
	if(FTLE_Compute) {
		if(FTLE_ComputeVariation) {
			/* Delete temp FTLE files */
			sprintf(BinFile, "%s%s.0.bin", Path_Data, FTLE_OutFilePrefix);
			if(remove(BinFile))
				fprintf(stderr, "Warning: Could not delete file %s\n", BinFile);
		}
	}
	if(Trace_Compute) {
		/* Output RT results */
		if(Trace_RTCompute) 
			OutputRT();
		else if(Trace_CETCompute) 
			OutputCET();
		
		/* Output AP field */
		if(Trace_APCompute)
			OutputAP();
		
		/* Output tracer positions */
		if(Trace_ReleaseStrategy == STAGGERED) {
			for(ss = 0; ss < Trace_NumLaunchTimes; ss++)  
				if(fclose(Trace_BinFileArray[ss].FileID)) 
					fprintf(stderr, "Warning: could not close file %s\n", Trace_BinFileArray[ss].FilePath);
			free(Trace_BinFileArray);
		}
		else {
			/* Remove temporary bin files */
			for(ss = 0; ss < Trace_NumLaunchTimes; ss++) {
				sprintf(BinFile, "%s.%d.bin", TraceN_BinFilePathPrefix, ss);
				if(remove(BinFile))
					fprintf(stderr, "Warning: Could not delete file %s\n", BinFile);
			}
			
			if (TRACER_LOWERLOOP) {
				
				fclose(Trace_BinFileID);
				OutputTracers();
				
			} else {
			/* Convert tracer position data from temp files to output files */
			OutputTracers_OMP();
			
			}
		}
	}
	
	printf("\nCleaning memory\n");
	fflush(stdout);
	
	FreeVelFieldData();
	
	if(FTLE_Compute) {
		FreeFTLEData();
	}
	if(Trace_Compute) {
		FreeTracerData();
		if(Trace_APCompute)
			FreeStrainRateData();
	}
	if(VelOut_Compute) {
		FreeVelOutData();
	}
	
	printf("\nSimulation complete!\n");
	
	time(&SimulationStopTime);
	runtime = difftime(SimulationStopTime, SimulationStartTime);
	printf("Total run time: %d:%d:%d\n", (int)(runtime/3600), (int)(fmod(runtime, 3600)/60), (int)(fmod(fmod(runtime, 3600), 60)));  
	printf("Processor time: %d:%d:%d\n", (int)(clock()/CLOCKS_PER_SEC/3600), (int)(fmod(clock()/CLOCKS_PER_SEC, 3600)/60), 
		   (int)(fmod(fmod(clock()/CLOCKS_PER_SEC, 3600), 60)));  
	
	
	return(0);
}
