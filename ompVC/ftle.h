/*
 *  ftle.h
 *  flowVC
 *
 *  Created by Shawn Shadden.
 *  Copyright 2010 Flow Physics Group. All rights reserved.
 *
 */


#ifndef INC_FTLE_H
#define INC_FTLE_H

#include "structs.h"


void ComputeFTLE(int df);

void InitializeFTLEArray(void);
void ReadInFTLELaunch(int ss, FTLEPoint ***FTLE_MeshPt);
void WriteOutFTLELaunch(int ss, FTLEPoint ***FTLE_MeshPt);
void GetFTLEForPointEarly(int ss, int i, int j, int k, double t1, double te, FTLEPoint ***FTLE_MeshPt);
void GetFTLEForPoint(int i,int j, int k, double IntTime, FTLEPoint ***FTLE_MeshPt);
LagrangianPoint Advect_FTLEPoint(int i, int j, int k, double t1, double t2, FTLEPoint ***FTLE_MeshPt);
void UpdateFTLELocations(FTLEPoint ***FTLE_MeshPt, double ****FTLE_NewArray);
void OutputFTLE(int ss, double IntT, FTLEPoint ***FTLE_MeshPt);
void FreeFTLEData(void);

#endif
