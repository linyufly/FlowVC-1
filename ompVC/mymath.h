/*
 *  mymath.h
 *  flowVC
 *
 *  Created by Shawn Shadden.
 *  Copyright 2010 Flow Physics Group. All rights reserved.
 *
 */


#ifndef INC_MYMATH_H
#define INC_MYMATH_H

double 	vdot(double *v1, double *v2, int dim);
void 	cross(const double *u, const double *v, double *result);
double 	dist(const double *v1, const double *v2, const int dim);
double	pythag(double a, double b);
double 	GetMaxEigenvalue(double a[][3]);

#endif

