

# include "parameters.h"
# include "settings.h"
# include "memcopy.h"
# include "index.cu"


void ReadInParameters(int argc, const char *argv[]) {

	FILE *Parameters_InFileID;
	
	if(argc == 2) {
		if((Parameters_InFileID = fopen(argv[1], "r")) == NULL) 
			FatalError("Could not open input file %s", argv[1]);
		
		fprintf(stderr, "Reading in parameters... \n\n");

		ReadInNextValue(Parameters_InFileID, Path_Data, 's');
		ReadInNextValue(Parameters_InFileID, Path_Output, 's');
		ReadInNextValue(Parameters_InFileID, &Dimensions, 'd');
		ReadInNextValue(Parameters_InFileID, &Data_MeshType, 'd'); 
		ReadInNextValue(Parameters_InFileID, Data_InFilePrefix, 's');
		ReadInNextValue(Parameters_InFileID, &Data_SuffixTMin, 'd');    
		ReadInNextValue(Parameters_InFileID, &Data_SuffixTDelta, 'd');
		ReadInNextValue(Parameters_InFileID, &Data_TRes, 'd');    
		ReadInNextValue(Parameters_InFileID, &Data_TDelta, 'f');    
		ReadInNextValue(Parameters_InFileID, &Data_TMin, 'f');  
		ReadInNextValue(Parameters_InFileID, &Data_TPeriodic, 'd'); 
		ReadInNextValue(Parameters_InFileID, &Data_MeshBounds.XMin, 'f');    
		ReadInNextValue(Parameters_InFileID, &Data_MeshBounds.XMax, 'f');    
		ReadInNextValue(Parameters_InFileID, &Data_MeshBounds.YMin, 'f');    
		ReadInNextValue(Parameters_InFileID, &Data_MeshBounds.YMax, 'f');    
		ReadInNextValue(Parameters_InFileID, &Data_MeshBounds.ZMin, 'f');    
		ReadInNextValue(Parameters_InFileID, &Data_MeshBounds.ZMax, 'f');
		
		ReadInNextValue(Parameters_InFileID, &Output_TStart, 'f');    
		ReadInNextValue(Parameters_InFileID, &Output_TRes, 'd');    
		ReadInNextValue(Parameters_InFileID, &Output_TDelta, 'f');  

		ReadInNextValue(Parameters_InFileID, &Int_Type, 'd');
		ReadInNextValue(Parameters_InFileID, &Int_TimeStep, 'f');
		ReadInNextValue(Parameters_InFileID, &Int_TimeDirection, 'd');   
		
		ReadInNextValue(Parameters_InFileID, &Int_Extrapolate, 'd');
		ReadInNextValue(Parameters_InFileID, &Trace_Compute, 'd');
		
		ReadInNextValue(Parameters_InFileID, &Trace_ReleaseStrategy, 'd');
		ReadInNextValue(Parameters_InFileID, &Trace_ReleaseTMax, 'f');
		
		  
		
		ReadInNextValue(Parameters_InFileID, &Trace_GenerateMesh, 'd');   
		
		ReadInNextValue(Parameters_InFileID, Trace_InFile, 's');
		ReadInNextValue(Parameters_InFileID, &Trace_MultipleInFiles, 'd');
		ReadInNextValue(Parameters_InFileID, &Trace_InFileFormat, 'd');		
		
		ReadInNextValue(Parameters_InFileID, Trace_OutFilePrefix, 's');
		ReadInNextValue(Parameters_InFileID, &Trace_NumLaunchTimes, 'd');
		ReadInNextValue(Parameters_InFileID, &Trace_LaunchTimeSpacing, 'f');    
		ReadInNextValue(Parameters_InFileID, &Trace_IntTLength, 'f'); 

		ReadInNextValue(Parameters_InFileID, &Trace_AlwaysOutput, 'd');

		ReadInNextValue(Parameters_InFileID, &Trace_CartMesh.XMin, 'f');    
		ReadInNextValue(Parameters_InFileID, &Trace_CartMesh.XMax, 'f');   
		ReadInNextValue(Parameters_InFileID, &Trace_CartMesh.YMin, 'f');    
		ReadInNextValue(Parameters_InFileID, &Trace_CartMesh.YMax, 'f');    
		ReadInNextValue(Parameters_InFileID, &Trace_CartMesh.ZMin, 'f');    
		ReadInNextValue(Parameters_InFileID, &Trace_CartMesh.ZMax, 'f');    
		ReadInNextValue(Parameters_InFileID, &Trace_CartMesh.XRes, 'd');    
		ReadInNextValue(Parameters_InFileID, &Trace_CartMesh.YRes, 'd');    
		ReadInNextValue(Parameters_InFileID, &Trace_CartMesh.ZRes, 'd');  


		ReadInNextValue(Parameters_InFileID, &Memory_Usage_GS, 'd');
		ReadInNextValue(Parameters_InFileID, &Memory_Usage_Tracer, 'd');
		ReadInNextValue(Parameters_InFileID, &Device_ID, 'd');
		ReadInNextValue(Parameters_InFileID, &Generate_Frame, 'd');
		ReadInNextValue(Parameters_InFileID, Temp_OutFilePrefix, 's'); 
		ReadInNextValue(Parameters_InFileID, &Keep_Tempfile, 'd');
		ReadInNextValue(Parameters_InFileID, &Searching, 'd');
		ReadInNextValue(Parameters_InFileID, &LocalSearchChecking, 'd');

		fclose(Parameters_InFileID);

		fprintf(stderr,"Path_Data\t\t=\t%10s\n",Path_Data);
		fprintf(stderr,"Path_Output\t\t=\t%10s\n\n",Path_Output);
		
		fprintf(stderr,"Dimensions\t\t=\t%10d\n", Dimensions);
		fprintf(stderr,"Data_MeshType\t\t=\t%10d\n",Data_MeshType);       
		fprintf(stderr,"Data_InFilePrefix\t=\t%10s\n",Data_InFilePrefix);

		fprintf(stderr,"Data_SuffixTMin\t\t=\t%10d\n",Data_SuffixTMin);
		fprintf(stderr,"Data_SuffixTDelta\t=\t%10d\n",Data_SuffixTDelta);

		fprintf(stderr,"Data_TRes\t\t=\t%10d\n",Data_TRes);
		fprintf(stderr,"Data_TDelta\t\t=\t%10g\n",Data_TDelta);
		fprintf(stderr,"Data_TMin\t\t=\t%10g\n",Data_TMin);
		fprintf(stderr,"Data_TPeriodic\t\t=\t%10d\n",Data_TPeriodic);

		fprintf(stderr,"Data_MeshBounds.XMin\t=\t%10g\n",Data_MeshBounds.XMin);
		fprintf(stderr,"Data_MeshBounds.XMax\t=\t%10g\n",Data_MeshBounds.XMax);
		fprintf(stderr,"Data_MeshBounds.YMin\t=\t%10g\n",Data_MeshBounds.YMin);
		fprintf(stderr,"Data_MeshBounds.YMax\t=\t%10g\n",Data_MeshBounds.YMax);
		fprintf(stderr,"Data_MeshBounds.ZMin\t=\t%10g\n",Data_MeshBounds.ZMin);
		fprintf(stderr,"Data_MeshBounds.ZMax\t=\t%10g\n\n",Data_MeshBounds.ZMax);

		fprintf(stderr,"Output_TStart\t\t=\t%10g\n",Output_TStart);
		fprintf(stderr,"Output_TRes\t\t=\t%10d\n",Output_TRes);
		fprintf(stderr,"Output_TDelta\t\t=\t%10g\n\n",Output_TDelta);

		fprintf(stderr,"Int_Type\t\t=\t%10d\n",Int_Type);
		fprintf(stderr,"Int_TimeStep\t\t=\t%10g\n",Int_TimeStep);

		fprintf(stderr,"Int_TimeDirection\t=\t%10d\n",Int_TimeDirection);
		fprintf(stderr,"Int_Extrapolate\t\t=\t%10d\n\n", Int_Extrapolate);

		fprintf(stderr,"Trace_Compute\t\t=\t%10d\n",Trace_Compute);
		
		fprintf(stderr,"Trace_ReleaseStrategy\t=\t%10d\n",Trace_ReleaseStrategy);
		fprintf(stderr,"Trace_ReleaseTMax\t=\t%10g\n",Trace_ReleaseTMax);
		fprintf(stderr,"Trace_GenerateMesh\t=\t%10d\n",Trace_GenerateMesh); 
		
		fprintf(stderr,"Trace_InFile\t\t=\t%10s\n",Trace_InFile);
		fprintf(stderr,"Trace_MultipleInFiles\t=\t%10d\n",Trace_MultipleInFiles);
		fprintf(stderr,"Trace_InFileFormat\t=\t%10d\n",Trace_InFileFormat);
		
		fprintf(stderr,"Trace_OutFilePrefix\t=\t%10s\n",Trace_OutFilePrefix);   
		
		
		fprintf(stderr,"Trace_NumLaunchTimes\t=\t%10d\n",Trace_NumLaunchTimes);
		fprintf(stderr,"Trace_LaunchTimeSpacing\t=\t%10g\n",Trace_LaunchTimeSpacing);
		fprintf(stderr,"Trace_IntTLength\t=\t%10g\n",Trace_IntTLength);
		fprintf(stderr,"Trace_AlwaysOutput\t=\t%10d\n",Trace_AlwaysOutput);

		fprintf(stderr,"Trace_CartMesh.XMin\t=\t%10g\n",Trace_CartMesh.XMin);
		fprintf(stderr,"Trace_CartMesh.XMax\t=\t%10g\n",Trace_CartMesh.XMax);
		fprintf(stderr,"Trace_CartMesh.YMin\t=\t%10g\n",Trace_CartMesh.YMin);
		fprintf(stderr,"Trace_CartMesh.YMax\t=\t%10g\n",Trace_CartMesh.YMax);
		fprintf(stderr,"Trace_CartMesh.ZMin\t=\t%10g\n",Trace_CartMesh.ZMin);
		fprintf(stderr,"Trace_CartMesh.ZMax\t=\t%10g\n",Trace_CartMesh.ZMax);
		fprintf(stderr,"Trace_CartMesh.XRes\t=\t%10d\n",Trace_CartMesh.XRes);
		fprintf(stderr,"Trace_CartMesh.YRes\t=\t%10d\n",Trace_CartMesh.YRes);
		fprintf(stderr,"Trace_CartMesh.ZRes\t=\t%10d\n\n",Trace_CartMesh.ZRes);
		
		if(Memory_Usage_GS == Use_Pinned_Memory)
			fprintf(stderr,"Memory Info: \n  You have opted to use Page-locked memory for global search.\n  Usage of Pinned memory can affect the performance of the system. \n\n");
		
		if(Memory_Usage_Tracer == Use_Pinned_Memory)
			fprintf(stderr,"Memory Info: \n  You have opted to use Page-locked memory for global search.\n  Usage of Pinned memory can affect the performance of the system. \n\n");
			
		fprintf(stderr,"The code will run on Device ID = %10d\n\n",Device_ID);
		
		fprintf(stderr,"Generate_Frame \t\t=\t%10d\n\n",Generate_Frame);
		fprintf(stderr,"Temp_OutFilePrefix\t=\t%10s\n",Temp_OutFilePrefix); 
		
		fprintf(stderr,"Keep_Tempfile\t\t=\t%10d\n\n",Keep_Tempfile); 
		fprintf(stderr,"Searching\t\t=\t%10d\n\n",Searching);
		fprintf(stderr,"LocalSearchChecking\t=\t%10d\n\n",LocalSearchChecking);
		

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
				FatalError("Unsupported data type passed to ReadInNextValue().\n");
			break;
		}
	}
}

void CheckParameters(void) {
	
	printf("Checking parameters...\n");
	
	if(Dimensions != 2 && Dimensions != 3) 
		FatalError("Dimensions must be equal to 2 or 3");
	if(Data_MeshType != CARTESIAN && Data_MeshType != UNSTRUCTURED)
		FatalError("Unsupported Data_MeshType");
	if(Data_TRes < 2)
		FatalError("Data_TRes must be 2 or greater.");
	if(Data_TDelta <= 0)
		FatalError("Data_TDelta <= 0");
	if(Data_TPeriodic != 0 && Data_TPeriodic != 1)
		FatalError("Data_TPeriodic should be 0 or 1");
	else if(Data_TPeriodic)
		printf("  Data specified as periodic: Ensure first and last data files correspond to same point in cycle.\n");
	if(Data_MeshBounds.XMin >= Data_MeshBounds.XMax)
		FatalError("Data_MeshBounds.XMin >= Data_MeshBounds.XMax");
	if(Data_MeshBounds.YMin >= Data_MeshBounds.YMax)
		FatalError("Data_MeshBounds.YMin >= Data_MeshBounds.YMax");
	if(Data_MeshBounds.ZMin > Data_MeshBounds.ZMax)
		FatalError("Data_MeshBounds.ZMin > Data_MeshBounds.ZMax");
	if(Output_TRes < 1)
		FatalError("Output_TRes < 1");
	if(Output_TDelta < TINY && Output_TRes > 1)
		FatalError("Output_TDelta must be positive");
	if(Int_Type != 0 && Int_Type != 1 && Int_Type != 2)
		fprintf(stderr, "Warning: Unrecognized value for Int_Type, setting to 0\n");
	if(Int_TimeStep < TINY && Int_Type != 2)
		FatalError("Int_TimeStep must be positive");
	if(Int_TimeDirection != 1 && Int_TimeDirection != -1)
		FatalError("Unsupported Int_TimeDirection value, must be -1 or 1");
	
	if(Trace_Compute) {
		if(Trace_GenerateMesh) {
			if(Trace_CartMesh.XMin > Trace_CartMesh.XMax)
				FatalError("Trace_CartMesh.XMin > Trace_CartMesh.XMax");
			if(Trace_CartMesh.YMin > Trace_CartMesh.YMax)
				FatalError("Trace_CartMesh.YMin > Trace_CartMesh.YMax");
			if(Trace_CartMesh.ZMin > Trace_CartMesh.ZMax)
				FatalError("Trace_CartMesh.ZMin > Trace_CartMesh.ZMax");
			if(Trace_CartMesh.XRes < 1)
				FatalError("Trace_CartMesh.XRes < 1");
			if(Trace_CartMesh.YRes < 1)
				FatalError("Trace_CartMesh.YRes < 1");
			if(Trace_CartMesh.ZRes < 1)
				FatalError("Trace_CartMesh.ZRes < 1");
		}
		
		if(Trace_NumLaunchTimes < 1)
			FatalError("Trace_NumLaunchTimes < 1");
		if(Trace_LaunchTimeSpacing < TINY && Trace_NumLaunchTimes > 1)
			FatalError("Trace_LaunchTimeSpacing must be positive");
		if(Trace_IntTLength < 0)
			FatalError("Trace_IntTLength < 0");
	
	}   
	
	printf("OK!\n\n");
	fflush(stdout);
	
}




void SetDerivedParameters(void) {


	printf("Setting derived parameters...\n");
	
	Data_TMax = Data_TMin + (Data_TRes - 1) * Data_TDelta;

	if((Output_TStart < Data_TMin || Output_TStart > Data_TMax) && !Data_TPeriodic)
		FatalError("Output_TStart outside data range");
	
	Output_TEnd = Output_TStart + (Output_TRes - 1) * Int_TimeDirection * Output_TDelta;

	if((Output_TEnd > Data_TMax || Output_TEnd < Data_TMin) && !Data_TPeriodic)
		FatalError("Output_TEnd outside data range");
	else
		printf("  Output_TEnd\t\t=\t%10g\n", Output_TEnd);
	
	Data_FirstFrame = (Int_TimeDirection > 0) ? floor((Output_TStart - Data_TMin) / Data_TDelta) : ceil((Output_TStart - Data_TMin) / Data_TDelta);
	Data_LastFrame = (Int_TimeDirection > 0) ? ceil((Output_TEnd - Data_TMin) / Data_TDelta) : floor((Output_TEnd - Data_TMin) / Data_TDelta);

	
	// Unlike flowVC computation is partitioned according to output frames ...
	frames = Output_TRes - 1;

	printf("  Number of frames\t=\t%10d \n",frames);
	
	printf("  Data_FirstFrame\t=\t%10d \n", Data_FirstFrame);
	printf("  Data_LastFrame\t=\t%10d \n", Data_LastFrame);

	if(Trace_Compute) {
	
		if(Trace_GenerateMesh) {

			if(Trace_CartMesh.XRes == 1)
				Trace_CartMesh.XDelta = 0.0;
			else
				Trace_CartMesh.XDelta = (Trace_CartMesh.XMax - Trace_CartMesh.XMin) / ((double)Trace_CartMesh.XRes - 1.0);
			if(Trace_CartMesh.YRes == 1)
				Trace_CartMesh.YDelta = 0.0;
			else 
				Trace_CartMesh.YDelta = (Trace_CartMesh.YMax - Trace_CartMesh.YMin) / ((double)Trace_CartMesh.YRes - 1.0);
			if(Trace_CartMesh.ZRes == 1)
				Trace_CartMesh.ZDelta = 0.0;
			else {
				if(Dimensions == 3)
					Trace_CartMesh.ZDelta = (Trace_CartMesh.ZMax - Trace_CartMesh.ZMin) / ((double)Trace_CartMesh.ZRes - 1.0);
				else {
					Trace_CartMesh.ZMin = 0;
					Trace_CartMesh.ZMax = 0;
					Trace_CartMesh.ZDelta = 0;
					fprintf(stderr, "\nWarning: Trace_CartMesh.ZRes being set to 1\n");
					Trace_CartMesh.ZRes = 1;
				}
			} 
		}
	
		// Each tracer launch will have same number of particles equal to this number. There will be Trace_NumLaunchTimes number of launches.
		N_Frame = Trace_CartMesh.XRes * Trace_CartMesh.YRes * Trace_CartMesh.ZRes ;

		// Check how many out of Trace_NumLaunchTimes will be actually released.
		
		int i,num = 0;
		double time1, timeend;
		// Time at the last frame.
		timeend =  Output_TStart + (Output_TRes - 1) *Output_TDelta;
		//printf("timeend = %f \n\n",timeend);
		for(i = 0; i < Trace_NumLaunchTimes; i++) {
		
			time1 = Output_TStart + (i*Trace_LaunchTimeSpacing);
			
			if (time1 <= timeend)
				num++;
		
		}
		if (num < Trace_NumLaunchTimes) {
			printf("\n  Out of %d launches, only %d will be actually launched.\n\n",Trace_NumLaunchTimes,num);
		}
		
		N_Total = N_Frame * num;	

		printf("  Number of points in starting frame = %d \n\n  Number of points in ending frame = %d \n",N_Frame,N_Total);
		
	}
	
	if(strcmp(Path_Output, "pwd")) {
		if(Path_Output[strlen(Path_Output) - 1] != '/') 
			sprintf(Path_Output, "%s/", Path_Output);
	}
	else 
		sprintf(Path_Output, "./");
	
	if(strcmp(Path_Data, "pwd")) {
		if(Path_Data[strlen(Path_Data) - 1] != '/' ) 
			sprintf(Path_Data, "%s/", Path_Data);
	}
	else 
		sprintf(Path_Data, "./");
		
		
	// Set our target cuda device number on which we will be running our code......
	err = cudaSetDevice(Device_ID);
	if(err != cudaSuccess) {
		fprintf(stderr, "Something went terribly wrong in starting the device. Please make sure you are running on cuda-enabled device...\n");
	    	printf("CUDA Error: %s \n\n", cudaGetErrorString(err));
	    	exit(1);
	}
	
	
	
	
	printf("OK!  \n\n");
	
}


