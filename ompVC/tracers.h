/*
 *  tracers.h
 *  flowVC
 *
 *  Created by Shawn Shadden.
 *  Copyright 2010 Flow Physics Group. All rights reserved.
 *
 */



#ifndef INC_TRACERS_H
#define INC_TRACERS_H

void Tracer_StaggeredRelease(FILE *Trace_BinFileID, FileData *Trace_BinFileArray);
//void Tracer_NormalRelease(FILE *Trace_BinFileID, FileData *Trace_BinFileArray, double nextoutputtime);
double Compute_Tracer (int *outputframe, double nextoutputtime, FILE *Trace_BinFileID);

void GenerateStaggeredRelease(void);
void LoadReleasePoints(double *voxdim);
void ReadInTraceMesh(void);
void GenerateTracerMesh(void);
void CreateNewTraceLaunch(int ss);
LagrangianPoint* ReadInTraceLaunch(int ss, LagrangianPoint *Trace_MeshPt);
void WriteOutTraceLaunch(int ss, LagrangianPoint *Trace_MeshPt);
void OutputTracers(void);
void FreeTracerData(void);

void OutputTracers_OMP(void);

#endif
