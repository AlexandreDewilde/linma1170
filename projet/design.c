
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmshc.h>
#include <errno.h>

void designTuningFork(double const r1, double const r2, double const e, double const l, double const meshSizeFactor, char const* filename) {
	/**
	 * r1 = inner radius (half-distance between prongs)
	 * r2 = outer radius (half-width of fork)
	 * e  = length of handle
	 * l  = length of prongs
	 * meshSizeFactor = meshSize / width of prongs
	 * if `filename` is not NULL, save to file
	*/
	int ierr;
	gmshClear(&ierr);

	double const h = r2 - r1; // width of prongs
	double const meshSize = h * meshSizeFactor;

	// Add points
	double x = 0;
	double y = 0;
	double z = 0;

	gmshModelOccAddPoint(x,y,z,meshSize,1,&ierr);
	x += h;
	gmshModelOccAddPoint(x,y,z,meshSize,2,&ierr);
	y += e;
	gmshModelOccAddPoint(x,y,z,meshSize,3,&ierr);
	x += r2;
	y += r2;
	gmshModelOccAddPoint(x,y,z,meshSize,4,&ierr);
	y += l;
	gmshModelOccAddPoint(x,y,z,meshSize,5,&ierr);
	x -= h;
	gmshModelOccAddPoint(x,y,z,meshSize,6,&ierr);
	y -= l;
	gmshModelOccAddPoint(x,y,z,meshSize,7,&ierr);
	x -= r1;
	y -= r1;
	gmshModelOccAddPoint(x,y,z,meshSize,8,&ierr);
	x -= h;
	gmshModelOccAddPoint(x,y,z,meshSize,9,&ierr);
	x -= r1;
	y += r1;
	gmshModelOccAddPoint(x,y,z,meshSize,10,&ierr);
	y += l;
	gmshModelOccAddPoint(x,y,z,meshSize,11,&ierr);
	x -= h;
	gmshModelOccAddPoint(x,y,z,meshSize,12,&ierr);
	y -= l;
	gmshModelOccAddPoint(x,y,z,meshSize,13,&ierr);
	x += r2;
	y -= r2;
	gmshModelOccAddPoint(x,y,z,meshSize,14,&ierr);
	y += (h+r1);
	gmshModelOccAddPoint(x,y,z,meshSize,15,&ierr);
	x += h;
	gmshModelOccAddPoint(x,y,z,meshSize,16,&ierr);
	
	// Add curves
	gmshModelOccAddLine(1,2,1,&ierr);
	gmshModelOccAddLine(2,3,2,&ierr);
	gmshModelOccAddCircleArc(3,16,4,3,&ierr);
	gmshModelOccAddLine(4,5,4,&ierr);
	gmshModelOccAddLine(5,6,5,&ierr);
	gmshModelOccAddLine(6,7,6,&ierr);
	gmshModelOccAddCircleArc(7,16,8,7,&ierr);
	gmshModelOccAddLine(8,9,8,&ierr);
	gmshModelOccAddCircleArc(9,15,10,9,&ierr);
	gmshModelOccAddLine(10,11,10,&ierr);
	gmshModelOccAddLine(11,12,11,&ierr);
	gmshModelOccAddLine(12,13,12,&ierr);
	gmshModelOccAddCircleArc(13,15,14,13,&ierr);
	gmshModelOccAddLine(14,1,14,&ierr);

	// Add wire (closed curve)
	int curveTags[14] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14};
	gmshModelOccAddWire(curveTags, 14, 1, 1, &ierr);

	// Add surface
	int wireTags[1] = {1};
	gmshModelOccAddPlaneSurface(wireTags, 1, 100, &ierr);

	// Sync
	gmshModelOccSynchronize(&ierr);

	// Create physical group for surface
	int surfaceTags[1] = {100};
	gmshModelAddPhysicalGroup(2, surfaceTags, 1, -1, "bulk", &ierr);

	// Create physical group for clamped curves
	int clampedCurveTags[3] = {1, 2, 14};
	gmshModelAddPhysicalGroup(1, clampedCurveTags, 3, -1, "clamped", &ierr);

	gmshModelMeshGenerate(2, &ierr);

	if(filename != NULL)
		gmshWrite(filename, &ierr);
}

void designTuningForkSymmetric(double const r1, double const r2, double const e, double const l, double const meshSizeFactor, char const* filename) {
	int ierr;

	gmshClear(&ierr);

	double const h = r2 - r1; // width of prongs
	double const meshSize = h * meshSizeFactor;

	// Add points
	double x = 0;
	double y = 0;
	double z = 0;

	gmshModelOccAddPoint(x, y , z, meshSize, 1, &ierr);
	x += h / 2;
	gmshModelOccAddPoint(x, y, z, meshSize, 2, &ierr);
	y += e;
	gmshModelOccAddPoint(x, y, z, meshSize, 3, &ierr);
	x += r2;
	y += r2;
	gmshModelOccAddPoint(x, y, z, meshSize, 4, &ierr);
	y += l;
	gmshModelOccAddPoint(x, y, z, meshSize, 5, &ierr);
	x -= h;
	gmshModelOccAddPoint(x, y, z, meshSize, 6, &ierr);
	y -= l;
	gmshModelOccAddPoint(x, y, z, meshSize, 7, &ierr);
	x -= r1;
	y -= r1;
	gmshModelOccAddPoint(x, y, z, meshSize, 8, &ierr);
	x = 0;
	gmshModelOccAddPoint(x, y, z, meshSize, 9, &ierr);
	x += h / 2;
	y += r1;
	gmshModelOccAddPoint(x, y, z, meshSize, 11, &ierr);

	gmshModelOccAddLine(1, 2, 1, &ierr);
	gmshModelOccAddLine(2, 3, 2, &ierr);
	gmshModelOccAddCircleArc(3, 11, 4, 3, &ierr);
	gmshModelOccAddLine(4, 5, 4, &ierr);
	gmshModelOccAddLine(5, 6, 5, &ierr);
	gmshModelOccAddLine(6, 7, 6, &ierr);
	gmshModelOccAddCircleArc(7, 11, 8, 7, &ierr);
	gmshModelOccAddLine(8, 9, 8, &ierr);
	gmshModelOccAddLine(9, 1, 9, &ierr);

	int curveTags[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
	gmshModelOccAddWire(curveTags, 9, 1, 1, &ierr);
	int wireTags[1] = {1};
	gmshModelOccAddPlaneSurface(wireTags, 1, 100, &ierr);
	gmshModelOccSynchronize(&ierr);

	// Create physical group for surface
	int surfaceTags[1] = {100};
	gmshModelAddPhysicalGroup(2, surfaceTags, 1, -1, "bulk", &ierr);

	int clampedCurveTags[3] = {1, 2, 9};
	gmshModelAddPhysicalGroup(1, clampedCurveTags, 3, -1, "clamped", &ierr);

	gmshModelMeshGenerate(2, &ierr);
}

static int designRec(size_t current, size_t idx, double x, double y, double const e, double const* b, double const* d, double const* h, double const* l, size_t const n, double const meshSize) {
	int ierr;
	if (current == n - 1) {
		gmshModelOccAddPoint(x, y, 0, meshSize, idx++, &ierr);
		x += b[current] + d[current];
		gmshModelOccAddPoint(x, y, 0, meshSize, idx++, &ierr);
		y += l[current] + h[current];
		gmshModelOccAddPoint(x, y, 0, meshSize, idx++, &ierr);
		x -= d[current];
		gmshModelOccAddPoint(x, y, 0, meshSize, idx++, &ierr);
		y -= l[current];
		gmshModelOccAddPoint(x, y, 0, meshSize, idx++, &ierr);
		x -= 2 * b[current] + e;
		gmshModelOccAddPoint(x, y, 0, meshSize, idx++, &ierr);
		y += l[current];
		gmshModelOccAddPoint(x, y, 0, meshSize, idx++, &ierr);
		x -= d[current];
		gmshModelOccAddPoint(x, y, 0, meshSize, idx++, &ierr);
		y -= l[current] + h[current];
		gmshModelOccAddPoint(x, y, 0, meshSize, idx++, &ierr);
		return idx;
	}
	gmshModelOccAddPoint(x, y, 0, meshSize, idx++, &ierr);
	x += b[current] + d[current];
	gmshModelOccAddPoint(x, y, 0, meshSize, idx++, &ierr);
	y += l[current] + h[current];
	idx = designRec(current + 1, idx, x, y, d[current], b, d, h, l, n, meshSize);
	x -= d[current];
	gmshModelOccAddPoint(x, y, 0, meshSize, idx++, &ierr);
	y -= l[current];
	gmshModelOccAddPoint(x, y, 0, meshSize, idx++, &ierr);
	x -= 2 * b[current] + e;
	gmshModelOccAddPoint(x, y, 0, meshSize, idx++, &ierr);
	y += l[current];
	idx = designRec(current + 1, idx, x, y, d[current], b, d, h, l, n, meshSize);
	x -= d[current];
	gmshModelOccAddPoint(x, y, 0, meshSize, idx++, &ierr);
	y -= l[current] + h[current];
	gmshModelOccAddPoint(x, y, 0, meshSize, idx++, &ierr);
	return idx;
}

void designTuningForkNLayer(double const lh, double const wh, double const* d, double const* dec, double const* h, double const* l, size_t n, double meshSizeFactor, char const* filename) {
	int ierr;
	gmshClear(&ierr);

	double b[n];
	b[n - 1] = 0.;
	for (int i = n - 1; i > 0; i--)
		b[i - 1] = d[i] + dec[i] + (n - i) * b[i];

	double meshSize = d[0] * meshSizeFactor;

	size_t idx = 0;
	double x = 0.;
	double y = 0.;
	gmshModelOccAddPoint(x, y, 0, meshSize, idx++, &ierr);
	x += wh;
	gmshModelOccAddPoint(x, y, 0, meshSize, idx++, &ierr);
	y += lh;
	// gmshModelOccAddPoint(x, y, 0, meshSize, idx++, &ierr);
	idx = designRec(0, idx, x, y, wh, b, d, h, l, n, meshSize);
	x -= wh;
	gmshModelOccAddPoint(x, y, 0, meshSize, idx++, &ierr);
	for (size_t i = 0; i < idx; i++)
		gmshModelOccAddLine(i, (i + 1) % idx, i, &ierr);

	// Add wire (closed curve)
	int curveTags[idx];
	for (size_t i = 0; i < idx; i++)
		curveTags[i] = i;
	gmshModelOccAddWire(curveTags, idx, 1, 1, &ierr);

	// Add surface
	int wireTags[1] = {1};
	gmshModelOccAddPlaneSurface(wireTags, 1, 100, &ierr);
	
	gmshModelOccSynchronize(&ierr);

	// Create physical group for surface
	int surfaceTags[1] = {100};
	gmshModelAddPhysicalGroup(2, surfaceTags, 1, -1, "bulk", &ierr);

	int clampedCurveTags[3] = {0, 1, idx - 1};
	gmshModelAddPhysicalGroup(1, clampedCurveTags, 3, -1, "clamped", &ierr);

	gmshModelMeshGenerate(2, &ierr);

	if (filename != NULL)
		gmshWrite(filename, &ierr);
}

void designTuningForkSymmetricNLayer(double const lh, double const wh, double const* d, double const* dec, double const* h, double const* l, size_t n, double const meshSizeFactor, char const* filename) {
	int ierr;
	gmshClear(&ierr);

	double b[n];
	b[n - 1] = 0.;
	for (int i = n - 1; i > 0; i--) 
		b[i - 1] = d[i] + dec[i] + (n - i) * b[i];

	double meshSize = d[0] * meshSizeFactor;
	size_t idx = 0;
	double x = 0.;
	double y = 0.;
	gmshModelOccAddPoint(x, y, 0, meshSize, idx++, &ierr);
	x += wh / 2;
	gmshModelOccAddPoint(x, y, 0, meshSize, idx++, &ierr);
	y += lh;
	gmshModelOccAddPoint(x, y, 0, meshSize, idx++, &ierr);
	x += d[0] + b[0];
	gmshModelOccAddPoint(x, y, 0, meshSize, idx++, &ierr);
	y += l[0] + h[0];
	// gmshModelOccAddPoint(x, y, 0, meshSize, idx++, &ierr);
	idx = designRec(1, idx, x, y, d[0], b, d, h, l, n, meshSize);
	x -= d[0];
	gmshModelOccAddPoint(x, y, 0, meshSize, idx++, &ierr);
	y -= l[0];
	gmshModelOccAddPoint(x, y, 0, meshSize, idx++, &ierr);
	x -= b[0] + wh / 2;
	gmshModelOccAddPoint(x, y, 0, meshSize, idx++, &ierr);
	for (size_t i = 0; i < idx; i++)
		gmshModelOccAddLine(i, (i + 1) % idx, i, &ierr);

	// Add wire (closed curve)
	int curveTags[idx];
	for (size_t i = 0; i < idx; i++)
		curveTags[i] = i;
	gmshModelOccAddWire(curveTags, idx, 1, 1, &ierr);

	// Add surface
	int wireTags[1] = {1};
	gmshModelOccAddPlaneSurface(wireTags, 1, 100, &ierr);
	
	gmshModelOccSynchronize(&ierr);

	// Create physical group for surface
	int surfaceTags[1] = {100};
	gmshModelAddPhysicalGroup(2, surfaceTags, 1, -1, "bulk", &ierr);

	int clampedCurveTags[3] = {0, 1, idx - 1};
	gmshModelAddPhysicalGroup(1, clampedCurveTags, 3, -1, "clamped", &ierr);

	gmshModelMeshGenerate(2, &ierr);

	if (filename != NULL)
		gmshWrite(filename, &ierr);
}