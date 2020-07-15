#ifndef MARCHINGCUBES_H_
#define MARCHINGCUBES_H_

#include "mp_vector.h"//vector operations
#include "mp_cubic.h"
#include "mc_table.h"//tables used by Marching Cubes (edgeTable and triTable)

#include<vector> 
#include<algorithm>

//used to save triangles - 3 vertices and a normal vector
struct TRIANGLE {
    mpVector p[3];
    mpVector norm;
    bool operator<(const TRIANGLE& t) const;
    bool operator==(const TRIANGLE& t) const;
};

//does Linear Interpolation between points p1 and p2 (they already contain their computed values)
mpVector LinearInterp(mp4Vector p1, mp4Vector p2, float value);

////////////////////////////////////////////////////////////////////////////////////////
//POINTERS TO FUNCTIONS////
//pointer to function which computes the value at a point////
typedef float (*FORMULA)(mpVector);////
////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////
///// the MARCHING CUBES algorithm /////////
////////////////////////////////////////////////////////////////////////////////////////

//Version 1
//Runs Marching Cubes on 3D array of points which already defined coordinates and values.
// ncellsX, ncellsY, ncellsZ):  number of cells to subdivide on each axis
// minValue: minimum to which value at every vertex is compared to
// points: array of length (ncellsX+1)(ncellsY+1)(ncellsZ+1) of mp4Vector points containing coordinates and values
//returns pointer to triangle array and the number of triangles in numTriangles
//note: array of points is first taken on z axis, then y and then x. So for example, if u iterate through it in a
//       for loop, have indexes i, j, k for x, y, z respectively, then to get the point you will have to make the
// following index: i*(ncellsY+1)*(ncellsZ+1) + j*(ncellsZ+1) + k .
//Also, the array starts at the minimum on all axes. (point at indices (i,j,k) is the minimum)
TRIANGLE* MarchingCubesCross(int ncellsX, int ncellsY, int ncellsZ,
			     float minValue, mp4Vector * points ,int &numTriangles);

//Version 2
//First computes the 3D array of coordintes and values and then runs the above function on the data
//takes dimensions (minx,maxx,miny,...) and the number of cells (ncellsX,...) to subdivide on each axis
// minValue: minimum to which value at every vertex is compared to
// function of type float (mpVector p) formula, which computes value of p at its coordinates
// function of type mpVector (mp4Vector p1, mp4Vector p2) intersection, which determines the
//  point of intersection of the surface and the edge between points p1 and p2
// saves number of triangles in numTriangles and the pointer to them is returned
// (note: mins' and maxs' are included in the algorithm)
TRIANGLE* MarchingCubesCross(float mcMinX, float mcMaxX, float mcMinY, float mcMaxY, float mcMinZ, float mcMaxZ,
			     int ncellsX, int ncellsY, int ncellsZ, float minValue,
			     FORMULA formula, int &numTriangles);


//Version 3
std::vector<TRIANGLE> MarchingCubesCross(std::vector<mpCubic> cubics, float minValue);

#endif