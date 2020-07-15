#include "marching_cubes_cross.h"

#include<vector>
#include<cmath>


static const int INSIDE = 1;
static const int BOUNDARY = 2;
static const int OUTSIDE = 3;
static const int EXTERIOR = 4;

/////////    constant definition
const long double MINVALUE = 0.0;
const long double LENGTH = 2.5;
const int NUM = 50;

const long double DIFF = LENGTH/(double) NUM;

const int NUM_ITER = 20;


// The fixed design domain D is divided into finite elements.
long double phi[NUM+1][NUM+1][NUM+1];
// bool place[NUM+1][NUM+1][NUM+1]; // for points

long double distance(long double x0, long double y0, long double z0, long double x1, long double y1, long double z1) {
    return std::sqrt( (x1-x0)*(x1-x0) + (y1-y0)*(y1-y0) + (z1-z0)*(z1-z0) );
}

void put_exterior() {
    for(int i=0;i<=NUM;i++) {
        for(int j=0;j<=NUM;j++) {
            for(int k=0;k<=NUM;k++) {
                if(false) place[i][j][k] = EXTERIOR; 
            }
        }
    }
}

void initial_phi() {

    // center 
    long double cent_x = LENGTH/2;
    long double cent_y = LENGTH/2;
    long double cent_z = LENGTH/2;

    // radius
    long double rad = 0.2;

    // phi inside
    long double inside = 0.1;
    long double outside = -0.1;

    for(int i=0;i<=NUM;i++) {
        for(int j=0;j<=NUM;j++) {
            for(int k=0;k<=NUM;k++) {
                if(distance(cent_x, cent_y, cent_z, DIFF*i, DIFF*j, DIFF*k) < rad) {
                    phi[i][j][k] = inside;
                } else {
                    phi[i][j][k] = outside;
                }
            }
        }
    }
}

std::pair<bool, bool> minus_plus(int i, int j , int k) {
    bool minus=false, plus=false;
    phi[i][j][k] <= MINVALUE? minus=true: plus=true;
    phi[i][j][k+1] <= MINVALUE? minus=true: plus=true;
    phi[i][j+1][k] <= MINVALUE? minus=true: plus=true;
    phi[i][j+1][k+1] <= MINVALUE? minus=true: plus=true;
    phi[i+1][j][k] <= MINVALUE? minus=true: plus=true;
    phi[i+1][j][k+1] <= MINVALUE? minus=true: plus=true;
    phi[i+1][j+1][k] <= MINVALUE? minus=true: plus=true;
    phi[i+1][j+1][k+1] <= MINVALUE? minus=true: plus=true;
    return {minus, plus};
}

std::vector<mpCubic> create_boundary_cubic() {
    std::vector<mpCubic> cbs;
    std::vector<std::tuple<int, int, int>> idx;
    for(int i=0;i<NUM;i++) {
        for(int j=0;j<NUM;j++) {
            for(int k=0;k<NUM;k++) {
                
                auto [minus, plus] = minus_plus(i,j,k);

                // if(place[i][j][k] == EXTERIOR) continue;

                // if(!minus and plus) place[i][j][k] = INSIDE;
                // if(minus and plus) place[i][j][k] = BOUNDARY;
                // if(minus and !plus) place[i][j][k] = OUTSIDE;

                if(minus and plus) {
                    idx.emplace_back(i, j, k);

                    // prepare vertices of cubic
                    std::vector<mp4Vector> vertices;

                    for(int in_i=0; in_i<=1; in_i++) {
                        for(int in_j=0; in_j<=1; in_j++) {
                            for(int in_k=0; in_k<=1; in_k++) {
                                mp4Vector cur = mp4Vector(DIFF*(i+in_i), DIFF*(j+in_j), DIFF*(k+in_k), phi[i+in_i][j+in_j][k+in_k]);
                                vertices.emplace_back(cur);
                            }
                        }
                    }

                    mpCubic cb = mpCubic(vertices, 0);
                    cbs.emplace_back(cb);
                } 
            }
        }
    } 

    // add ext
    for(auto itr=idx.begin();itr!=idx.end();itr++) {
        int ext = 0;
        for(int in_i=-1; in_i<=1; in_i+=2) {
            if(i+in_i < 0 and i+in_i >= NUM) 
        }
        for(int in_j=-1; in_j<=1; in_j+=2) {
                
        }
        for(int in_k=-1; in_k<=1; in_k+=2) {
                    
        }
    }
}

int main() {

    // here define design domain with EXTERIOR
    // put_exterior();

    // An initial distribution of the level set function φ0 is given on nodes of the finite elements.
    initial_phi();

    for(int cur_iter=0; cur_iter<=NUM_ITER;cur_iter++) {
        // create cubic which is on the boundary
        auto cubics = create_boundary_cubic();

        // Boundary elements are generated from the distribution of φ.
        std::vector<TRIANGLE> bnd = MarchingCubesCross(cubics, MINVALUE);

        // The forward boundary value problem in Eqs. (3)–(5) is solved by the FMBEM.
        

        // The objective function is evaluated by the FMBEM

        // The adjoint boundary value problem in Eqs. (18)–(20) is solved by the FMBEM.

        // Topological derivative in Eq. (42) is evaluated on the nodes of the finite elements.

        // Distribution of the level set function φ is obtained as the solution of the boundary value problem in Eqs. (61) and (60) by the FEM

    }
}