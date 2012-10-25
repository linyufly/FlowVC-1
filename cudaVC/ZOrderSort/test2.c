#include <stdlib.h> 
#include <stdio.h>
#include <math.h>

int main (int argc, const char * argv[]) {

  FILE *fin;
  char fname[100];
  int i, numOfCells = 0;
  int *adjacencies;

  sprintf(fname, "%s_adjacency.bin", argv[1]);
  if((fin = fopen(fname, "rb")) == NULL) { 
    fprintf(stderr,"Could not open file %s",fname);
    exit(1);
  }

  if (fread(&numOfCells, sizeof(int), 1, fin) != 1){fprintf(stderr, "Failed to read numOfCells"); exit(1);}
  printf("numOfCells = %d\n", numOfCells);
  if((adjacencies = (int *)calloc(4 * numOfCells, sizeof(int))) == NULL) {fprintf(stderr,"Calloc failed"); exit(1);}
  if (fread(adjacencies, sizeof(int), numOfCells * 4, fin) != numOfCells * 4){ printf("Failed to read adjacencies"); exit(1);}

  /// DEBUG ///
  for (i = 0; i < numOfCells * 4; i++)
    if (adjacencies[i] >= numOfCells)
      printf("%d(%d) ", i, adjacencies[i]);
  printf("\n");

  fclose(fin);

}
