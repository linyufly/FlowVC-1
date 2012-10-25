/*******************************************
Program		:	Z-Order Sort
Author		:	Mingcheng Chen
Last Update	:	October 24th, 2012
*******************************************/

#include <cstdio>
#include <cstdlib>
#include <algorithm>

char *coordinatesFile, *connectivitiesFile, *adjacenciesFile;
char *velocitiesPrefix;
int velStart, velEnd, velStepSize;

char fileName[1000];

double *coordinates, *velocities;
double *xValues, *yValues, *zValues;
int *intCoordinates;
int *connectivities, *adjacencies;
int *sortedPointIDs, *newPointIDs;
int *sortedCellIDs, *newCellIDs;
int numOfPoints, numOfCells, lx, ly, lz;

void Error(const char *message) {
	printf("Error: %s\n", message);
	exit(0);
}

void Arrange(double *values, int &length, int originalLength) {
	std::sort(values, values + originalLength);
	length = std::unique(values, values + originalLength) - values;
}

void Init() {
	// Read coordinates
	FILE *fin = fopen(coordinatesFile, "rb");

	if (fread(&numOfPoints, sizeof(int), 1, fin) != 1) Error("Fail to read numOfPoints");

	printf("numOfPoints = %d\n", numOfPoints);
	
	coordinates = new double [numOfPoints * 3];

	if (fread(coordinates, sizeof(double), numOfPoints * 3, fin) != numOfPoints * 3) Error("Fail to read coordinates");

	fclose(fin);

	// Read connectivities
	fin = fopen(connectivitiesFile, "rb");

	if (fread(&numOfCells, sizeof(int), 1, fin) != 1) Error("Fail to read numOfCells");

	printf("numOfCells = %d\n", numOfCells);

	connectivities = new int [numOfCells * 4];

	if (fread(connectivities, sizeof(int), numOfCells * 4, fin) != numOfCells * 4) Error("Fail to read connectivities");

	fclose(fin);

	// Read adjacencies
	fin = fopen(adjacenciesFile, "rb");

	int _numOfCells;

	if (fread(&_numOfCells, sizeof(int), 1, fin) != 1) Error("Fail to read _numOfCells");

	printf("_numOfCells = %d\n", _numOfCells);

	if (_numOfCells != numOfCells) Error("numOfCells is not consistent.");

	adjacencies = new int [numOfCells * 4];

	if (fread(adjacencies, sizeof(int), numOfCells * 4, fin) != numOfCells * 4) Error("Fail to read adjacencies");

	/// DEBUG ///
	//for (int i = 0; i < numOfCells * 4; i++)
	//	if (adjacencies[i] < 0 || adjacencies[i] >= numOfCells)
	//		printf("%d(%d) ", i, adjacencies[i]);
	//printf("\n");

	fclose(fin);
}

bool CompareByZOrder(int id1, int id2) {
	int a[3] = {intCoordinates[id1 * 3], intCoordinates[id1 * 3 + 1], intCoordinates[id1 * 3 + 2]};
	int b[3] = {intCoordinates[id2 * 3], intCoordinates[id2 * 3 + 1], intCoordinates[id2 * 3 + 2]};

	int j = 0, x = 0;

	for (int k = 0; k < 3; k++) {
		int y = a[k] ^ b[k];
		if (x < y && x < (x ^ y)) {
			j = k;
			x = y;
		}
	}

	return a[j] < b[j];
}

void Work() {
	// Build xValues, yValues and zValues
	xValues = new double [numOfPoints];
	yValues = new double [numOfPoints];
	zValues = new double [numOfPoints];

	for (int i = 0; i < numOfPoints; i++) {
		xValues[i] = coordinates[i * 3];
		yValues[i] = coordinates[i * 3 + 1];
		zValues[i] = coordinates[i * 3 + 2];
	}

	Arrange(xValues, lx, numOfPoints);
	Arrange(yValues, ly, numOfPoints);
	Arrange(zValues, lz, numOfPoints);

	// Build int coordinates
	intCoordinates = new int [numOfPoints * 3];
	for (int i = 0; i < numOfPoints; i++) {
		intCoordinates[i * 3] = std::lower_bound(xValues, xValues + lx, coordinates[i * 3]) - xValues;
		intCoordinates[i * 3 + 1] = std::lower_bound(yValues, yValues + ly, coordinates[i * 3 + 1]) - yValues;
		intCoordinates[i * 3 + 2] = std::lower_bound(zValues, zValues + lz, coordinates[i * 3 + 2]) - zValues;	
	}

	// Build sorted point ID list
	sortedPointIDs = new int [numOfPoints];

	for (int i = 0; i < numOfPoints; i++)
		sortedPointIDs[i] = i;

	std::sort(sortedPointIDs, sortedPointIDs + numOfPoints, CompareByZOrder);

	newPointIDs = new int [numOfPoints];

	for (int i = 0; i < numOfPoints; i++)
		newPointIDs[sortedPointIDs[i]] = i;

	// Output new coordinates
	sprintf(fileName, "%s.new", coordinatesFile);

	FILE *fout = fopen(fileName, "wb");

	if (fwrite(&numOfPoints, sizeof(int), 1, fout) != 1) Error("Fail to write numOfPoints");

	for (int i = 0; i < numOfPoints; i++) {
		int id = sortedPointIDs[i];

		if (fwrite(coordinates + id * 3, sizeof(double), 3, fout) != 3) Error("Fail to write coordinates");
	}

	fclose(fout);

	// Clean xVAlues, yValues, zValues and intCoordinates
	delete [] xValues;
	delete [] yValues;
	delete [] zValues;
	delete [] intCoordinates;

	// Build xValues, yValues and zValues
	xValues = new double [numOfCells];
	yValues = new double [numOfCells];
	zValues = new double [numOfCells];

	for (int i = 0; i < numOfCells; i++) {
		xValues[i] = 0;
		yValues[i] = 0;
		zValues[i] = 0;

		for (int j = 0; j < 4; j++) {
			int id = connectivities[i * 4 + j];
			xValues[i] += coordinates[id * 3];
			yValues[i] += coordinates[id * 3 + 1];
			zValues[i] += coordinates[id * 3 + 2];
		}

		xValues[i] /= 4;
		yValues[i] /= 4;
		zValues[i] /= 4;
	}

	Arrange(xValues, lx, numOfCells);
	Arrange(yValues, ly, numOfCells);
	Arrange(zValues, lz, numOfCells);

	// Build int coordinates
	intCoordinates = new int [numOfCells * 3];
	for (int i = 0; i < numOfCells; i++) {
		double x = 0, y = 0, z = 0;

		for (int j = 0; j < 4; j++) {
			int id = connectivities[i * 4 + j];
			x += coordinates[id * 3];
			y += coordinates[id * 3 + 1];
			z += coordinates[id * 3 + 2];
		}

		x /= 4;
		y /= 4;
		z /= 4;

		intCoordinates[i * 3] = std::lower_bound(xValues, xValues + lx, x) - xValues;
		intCoordinates[i * 3 + 1] = std::lower_bound(yValues, yValues + ly, y) - yValues;
		intCoordinates[i * 3 + 2] = std::lower_bound(zValues, zValues + lz, z) - zValues;	
	}

	// Build sorted cell ID list
	sortedCellIDs = new int [numOfCells];
	
	for (int i = 0; i < numOfCells; i++)
		sortedCellIDs[i] = i;

	std::sort(sortedCellIDs, sortedCellIDs + numOfCells, CompareByZOrder);

	newCellIDs = new int [numOfCells];

	for (int i = 0; i < numOfCells; i++)
		newCellIDs[sortedCellIDs[i]] = i;

	// Output new connectivities
	sprintf(fileName, "%s.new", connectivitiesFile);

	fout = fopen(fileName, "wb");

	fwrite(&numOfCells, sizeof(int), 1, fout);

	for (int i = 0; i < numOfCells; i++) {
		int id = sortedCellIDs[i];

		for (int j = 0; j < 4; j++) {
			int pointID = newPointIDs[connectivities[id * 4 + j]];

			fwrite(&pointID, sizeof(int), 1, fout);
		}
	}

	fclose(fout);

	printf("Finish new connectivities\n");

	// Output new adjacencies
	sprintf(fileName, "%s.new", adjacenciesFile);

	fout = fopen(fileName, "wb");

	fwrite(&numOfCells, sizeof(int), 1, fout);

	for (int i = 0; i < numOfCells; i++) {
		int id = sortedCellIDs[i];

		for (int j = 0; j < 4; j++) {
			int adj = adjacencies[id * 4 + j];

			//if (adj < 0 || adj >= numOfCells) printf("adj = %d\n", adj);

			int neighborID = adj == -1 ? -1 : newCellIDs[adj];

			fwrite(&neighborID, sizeof(int), 1, fout);
		}
	}

	fclose(fout);

	printf("Finish new adjacencies\n");

	// Output new velocities
	velocities = new double [numOfPoints * 3];

	for (int frameID = velStart; frameID <= velEnd; frameID += velStepSize) {
		sprintf(fileName, "%s.%d.bin", velocitiesPrefix, frameID);

		FILE *fin = fopen(fileName, "rb");

		double timeStamp;

		fread(&timeStamp, sizeof(double), 1, fin);

		if (fread(velocities, sizeof(double), numOfPoints * 3, fin) != numOfPoints * 3) Error("Fail to read velocities");

		fclose(fin);

		sprintf(fileName, "%s.%d.bin.new", velocitiesPrefix, frameID);

		fout = fopen(fileName, "wb");

		fwrite(&timeStamp, sizeof(double), 1, fin);

		for (int i = 0; i < numOfPoints; i++) {
			int id = sortedPointIDs[i];

			fwrite(velocities + id * 3, sizeof(double), 3, fout);	
		}

		fclose(fout);
	}

	printf("Finish new velocities\n");
}

int main(int argc, char **argv) {
	coordinatesFile = argv[1];
	connectivitiesFile = argv[2];
	adjacenciesFile = argv[3];

	velocitiesPrefix = argv[4];
	sscanf(argv[5], "%d", &velStart);
	sscanf(argv[6], "%d", &velEnd);
	sscanf(argv[7], "%d", &velStepSize);

	Init();
	Work();

	return 0;
}
