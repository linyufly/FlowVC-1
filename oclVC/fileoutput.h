#ifndef INC_FILEOUTPUT_H
#define INC_FILEOUTPUT_H


extern void temp2file (void);
//extern void copyresult2file(int ss);
extern void copyresult2vtk(int ss, int N, double t);
extern void copyresult2vtk_staggered(int ss, int N, double t, double tmin, double tmax);
extern void copyresult2bin(int ss, int N, double t);
#endif
