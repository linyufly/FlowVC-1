


// Cuda Kernel Block size
#define POINTS_BLOCKSIZE_GS 512


#define POINTS_BLOCKSIZE_MAIN 512
#define FTLE_BLOCKSIZE 100


// Intervals for multistep methods.
#define INTERVAL 2 


#define TINY 1.0e-8
#define TINY_GS 1.0e-6
#define LONGSTRING 200
#define SHORTSTRING 100

#define UNLAUNCHED 0
#define LAUNCHED 1
#define COMPLETE 2

#define CARTESIAN 0
#define UNSTRUCTURED 1
#define Use_Pinned_Memory 1

#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define fmin(a, b) ((a) < (b) ? (a) : (b))
#define fmax(a, b) ((a) > (b) ? (a) : (b))

#define MAX_CONSTANT_MEMORY_GS 65536
#define POINT1_MEMORY_2D 8
#define POINT1_MEMORY_3D 12

#define CONSTANT_MEMORY 5312 // This number is chosen carefully. Do not change this unless hardware changes. Currently all GPUs support 64kb of constant memory. We load tracers position (x, y and z co-ordinate) into the constant memory. 

#define POINTS_MAIN 10000  // This number is currently 50K. change it according to our graphics card specs. 



#define CONSTANT_MEMORY_MAIN 4000




