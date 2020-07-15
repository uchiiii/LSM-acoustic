#include "marching_cubes_cross.h"

bool TRIANGLE::operator<(const TRIANGLE& t) const {
    if(p[0] < t.p[0]) return true;
    else if(p[0] > t.p[0]) return false;
    else {
        if(p[1] < t.p[1]) return true;
        else if(p[1] > t.p[1]) return false;
        else {
            if(p[2] < t.p[2]) return true;
            else return false;
        }
    }   
}

bool TRIANGLE::operator<(const TRIANGLE& t) const {
    if(p[0] == t.p[0] and p[1] == t.p[1] and p[2] == t.p[2]) return true;
    else return false;
}

//Linear Interpolation function
mpVector LinearInterp(mp4Vector p1, mp4Vector p2, float value)
{
  mpVector p;
  if(p1.val != p2.val)
    p = (mpVector)p1 + ((mpVector)p2 - (mpVector)p1)/(p2.val - p1.val)*(value - p1.val);
  else
    p = (mpVector)p1;
  return p;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
//MARCHING CUBES//
//////////////////////////////////////////////////////////////////////////////////////////////////////

//  VERSION  1  //
TRIANGLE* MarchingCubesCross(int ncellsX, int ncellsY, int ncellsZ,
			     float minValue, mp4Vector * points, int &numTriangles)
{
  //this should be enough space, if not change 3 to 4
  TRIANGLE * triangles = new TRIANGLE[3*ncellsX*ncellsY*ncellsZ];
  //TRIANGLE * triangles = new TRIANGLE[ncellsX*20*ncellsZ]; //  to save space
  numTriangles = int(0);

  int YtimeZ = (ncellsY+1)*(ncellsZ+1);//for little extra speed
  int ni, nj;

  //go through all the points
  for(int i=0; i < ncellsX; i++){//x axis
    ni = i*YtimeZ;
    for(int j=0; j < ncellsY; j++){//y axis
      nj = j*(ncellsZ+1);
      for(int k=0; k < ncellsZ; k++)//z axis
	{
	  //initialize vertices
	  mp4Vector verts[8];
	  int ind = ni + nj + k;
	  /*(step 3)*/ verts[0] = points[ind];
	  verts[1] = points[ind + YtimeZ];
	  verts[2] = points[ind + YtimeZ + 1];
	  verts[3] = points[ind + 1];
	  verts[4] = points[ind + (ncellsZ+1)];
	  verts[5] = points[ind + YtimeZ + (ncellsZ+1)];
	  verts[6] = points[ind + YtimeZ + (ncellsZ+1) + 1];
	  verts[7] = points[ind + (ncellsZ+1) + 1];

	  //get the index
	  int cubeIndex = int(0);
	  for(int n=0; n < 8; n++)
	    /*(step 4)*/if(verts[n].val <= minValue) cubeIndex |= (1 << n);

	  //check if its completely inside or outside
	  /*(step 5)*/ if(!edgeTable[cubeIndex]) continue;

	  //get linearly interpolated vertices on edges and save into the array
	  mpVector intVerts[12];
	  /*(step 6)*/ if(edgeTable[cubeIndex] & 1) intVerts[0] = LinearInterp(verts[0], verts[1], minValue);
	  if(edgeTable[cubeIndex] & 2) intVerts[1] = LinearInterp(verts[1], verts[2], minValue);
	  if(edgeTable[cubeIndex] & 4) intVerts[2] = LinearInterp(verts[2], verts[3], minValue);
	  if(edgeTable[cubeIndex] & 8) intVerts[3] = LinearInterp(verts[3], verts[0], minValue);
	  if(edgeTable[cubeIndex] & 16) intVerts[4] = LinearInterp(verts[4], verts[5], minValue);
	  if(edgeTable[cubeIndex] & 32) intVerts[5] = LinearInterp(verts[5], verts[6], minValue);
	  if(edgeTable[cubeIndex] & 64) intVerts[6] = LinearInterp(verts[6], verts[7], minValue);
	  if(edgeTable[cubeIndex] & 128) intVerts[7] = LinearInterp(verts[7], verts[4], minValue);
	  if(edgeTable[cubeIndex] & 256) intVerts[8] = LinearInterp(verts[0], verts[4], minValue);
	  if(edgeTable[cubeIndex] & 512) intVerts[9] = LinearInterp(verts[1], verts[5], minValue);
	  if(edgeTable[cubeIndex] & 1024) intVerts[10] = LinearInterp(verts[2], verts[6], minValue);
	  if(edgeTable[cubeIndex] & 2048) intVerts[11] = LinearInterp(verts[3], verts[7], minValue);

	  //now build the triangles using triTable
	  for (int n=0; triTable[cubeIndex][n] != -1; n+=3) {
	    /*(step 7)*/ triangles[numTriangles].p[0] = intVerts[triTable[cubeIndex][n+2]];
	    triangles[numTriangles].p[1] = intVerts[triTable[cubeIndex][n+1]];
	    triangles[numTriangles].p[2] = intVerts[triTable[cubeIndex][n]];
	    //Computing normal as cross product of triangle's edges
	    /*(step 8)*/ triangles[numTriangles].norm = ((triangles[numTriangles].p[1] -
							  triangles[numTriangles].p[0]).Cross(triangles[numTriangles].p[2] -
											      triangles[numTriangles].p[0])).Normalize();
	    numTriangles++;
	  }

	}//END OF Z FOR LOOP
    }//END OF Y FOR LOOP
  }//END OF X FOR LOOP

  //free all the wasted space
  TRIANGLE * retTriangles = new TRIANGLE[numTriangles];
  for(int i=0; i < numTriangles; i++)
    retTriangles[i] = triangles[i];
  delete [] triangles;

  return retTriangles;
}


//VERSION  2  //
TRIANGLE* MarchingCubesCross(float mcMinX, float mcMaxX, float mcMinY, float mcMaxY, float mcMinZ, float mcMaxZ,
			     int ncellsX, int ncellsY, int ncellsZ, float minValue,
			     FORMULA formula, int &numTriangles)
{
  //space is already defined and subdivided (mcMinX,...), staring with step 3
  //first initialize the points
  mp4Vector * mcDataPoints = new mp4Vector[(ncellsX+1)*(ncellsY+1)*(ncellsZ+1)];
  mpVector stepSize((mcMaxX-mcMinX)/ncellsX, (mcMaxY-mcMinY)/ncellsY, (mcMaxZ-mcMinZ)/ncellsZ);

  int YtimesZ = (ncellsY+1)*(ncellsZ+1);//for little extra speed
  for(int i=0; i < ncellsX+1; i++) {
    int ni = i*YtimesZ;//for little extra speed
    float vertX = mcMinX + i*stepSize.x;
    for(int j=0; j < ncellsY+1; j++) {
      int nj = j*(ncellsZ+1);//for little extra speed
      float vertY = mcMinY + j*stepSize.y;
      for(int k=0; k < ncellsZ+1; k++) {
	mp4Vector vert(vertX, vertY, mcMinZ + k*stepSize.z, 0);
	vert.val = formula((mpVector)vert);
	/*(step 3)*/ mcDataPoints[ni + nj + k] = vert;
      }
    }
  }
  //then run Marching Cubes (version 1) on the data
  return MarchingCubesCross(ncellsX, ncellsY, ncellsZ, minValue, mcDataPoints, numTriangles);
}

//VERSION  3 //
std::vector<TRIANGLE> MarchingCubesCross(std::vector<mpCubic> cubics, float minValue) {
    std::vector<TRIANGLE> trgs;

    for(auto cubic: cubics) {
        
        int cubeIndex = 0;
        for(int i=0; i<=7; i++) {
            if(cubic.vertices[i].val <= minValue) cubeIndex |= (1<<i);
        }

        if(!edgeTable[cubeIndex]) continue; // this is not required.

        mpVector intVerts[12];

        if(edgeTable[cubeIndex] & 1) intVerts[0] = LinearInterp(cubic.vertices[0], cubic.vertices[1], minValue);
	    if(edgeTable[cubeIndex] & 2) intVerts[1] = LinearInterp(cubic.vertices[1], cubic.vertices[2], minValue);
	    if(edgeTable[cubeIndex] & 4) intVerts[2] = LinearInterp(cubic.vertices[2], cubic.vertices[3], minValue);
	    if(edgeTable[cubeIndex] & 8) intVerts[3] = LinearInterp(cubic.vertices[3], cubic.vertices[0], minValue);
	    if(edgeTable[cubeIndex] & 16) intVerts[4] = LinearInterp(cubic.vertices[4], cubic.vertices[5], minValue);
	    if(edgeTable[cubeIndex] & 32) intVerts[5] = LinearInterp(cubic.vertices[5], cubic.vertices[6], minValue);
	    if(edgeTable[cubeIndex] & 64) intVerts[6] = LinearInterp(cubic.vertices[6], cubic.vertices[7], minValue);
	    if(edgeTable[cubeIndex] & 128) intVerts[7] = LinearInterp(cubic.vertices[7], cubic.vertices[4], minValue);
	    if(edgeTable[cubeIndex] & 256) intVerts[8] = LinearInterp(cubic.vertices[0], cubic.vertices[4], minValue);
	    if(edgeTable[cubeIndex] & 512) intVerts[9] = LinearInterp(cubic.vertices[1], cubic.vertices[5], minValue);
	    if(edgeTable[cubeIndex] & 1024) intVerts[10] = LinearInterp(cubic.vertices[2], cubic.vertices[6], minValue);
	    if(edgeTable[cubeIndex] & 2048) intVerts[11] = LinearInterp(cubic.vertices[3], cubic.vertices[7], minValue);

        for (int n=0; triTable[cubeIndex][n] != -1; n+=3) {
            TRIANGLE trg;
	        trg.p[0] = intVerts[triTable[cubeIndex][n+2]];
	        trg.p[1] = intVerts[triTable[cubeIndex][n+1]];
	        trg.p[2] = intVerts[triTable[cubeIndex][n]];
	        // Computing normal as cross product of triangle's edges
	        trg.norm = ((trg.p[1] - trg.p[0]).Cross(trg.p[2] - trg.p[0])).Normalize();
            trgs.emplace_back(trg);
	    }

        // added triangles with faces here
        
    }

    sort(trgs.begin(), trgs.end());

    // remove duplicates 
    for(auto itr = trgs.begin(); itr != trgs.end();) {
        auto left = itr; itr++;
        while(*itr == *left) itr++; 
        itr = trgs.erase(left, itr);
    }

    return trgs;
}

