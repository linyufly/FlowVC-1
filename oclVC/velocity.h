#ifndef INC_VELOCITY_H
#define INC_VELOCITY_H



extern void readveldata(int ss);

extern int TestOutsideDomain_host(double point[3]);

extern int TestOutsideCartVelDomain(double X[3]);

extern void GetVelocity_Cartesian(double tq, LagrangianPoint *pt, double *dXdt);

extern void GetVelocity_Unstructured(const double tq, LagrangianPoint *pt, double *dXdt);




#endif
