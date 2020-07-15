using Meshing
using FileIO # MeshIO should also be installed
using GeometryBasics

####### Patameters 
const LENGTH = 1.0
const DEVIDES = 100

const DIFF = LENGTH/DEVIDES

const OBJ_FUNC_EPS = 1e-5

const CATG = ["u", "lambda"]

const NUM_ITERATION = 5


###### Other Physical Parameters
const DENSITY = 1.1839e+00 #density of the propagation medium
const SOUND_SPEED = 346.18

####### Level-set method Parameters
K = 1e5
tau = 2.5*2.5*1e-4
delta_t = 0.005

####### Frequency info
frqStart = 0.0 #kHz
frqStep = SOUND_SPEED/1000 #kHz
frqNums = 1

####### Evaluation mesh (place of microphones)
num_microphones = 2
m1 = Point(1.25, 1.25, 5.0, 0)
m2 = Point(1.25, 1.25, -2.5, 0)
microphones = [m1, m2]

mic_p = Array{Complex{Float64}, 1}(undef, num_microphones)
mic_v = Array{Complex{Float64}, 1}(undef, num_microphones)

####### Place of Sound Source
num_source = 2
INTENSITY = 70
s1 = Point(5.0, 1.25, 1.25, INTENSITY)
s2 = Point(-2.5, 1.25, 1.24, INTENSITY)
sources = [s1, s2]

####### Create Domain
phi = Array{Float64, 3}(undef, DEVIDES, DEVIDES, DEVIDES)

####### Array
u_p = Array{Complex{Float64}, 3}(undef, DEVIDES, DEVIDES, DEVIDES)
u_v = Array{Complex{Float64}, 3}(undef, DEVIDES, DEVIDES, DEVIDES)
lambda_p = Array{Complex{Float64}, 3}(undef, DEVIDES, DEVIDES, DEVIDES)
lambda_v = Array{Complex{Float64}, 3}(undef, DEVIDES, DEVIDES, DEVIDES)

####### initial distribution of phi_{0} 
const EPS = 1e-10
const center1_x = 0.20
const center1_y = 0.20
const center1_z = 0.20
const center2_x = 0.80
const center2_y = 0.80
const center2_z = 0.80
const radius = 0.1

println(DIFF)

for i = 1:DEVIDES
    for j = 1:DEVIDES
        for k = 1:DEVIDES
            dt1 = sqrt( (i*DIFF - center1_x)*(i*DIFF - center1_x) + (j*DIFF - center1_y) * (j*DIFF - center1_y) + (k*DIFF - center1_z) * (k*DIFF - center1_z))
            dt2 = sqrt( (i*DIFF - center2_x)*(i*DIFF - center2_x) + (j*DIFF - center2_y) * (j*DIFF - center2_y) + (k*DIFF - center2_z) * (k*DIFF - center2_z))
            if  dt1 < radius - EPS || dt2 < radius - EPS
                phi[i,j,k] = -0.1
            else
                phi[i,j,k] = 0.1
            end
            # if (i == 20 && j == 20 && k == 20 ) || (i == 80 && j == 80 && k == 80 )
            #     phi[i,j,k] = 0.2
            # else
            #     phi[i,j,k] = -0.2
            # end
        end
    end
end

# generate directly using GeometryBasics API
# Rect specifies the sampling intervals
# gy_mesh = Mesh(phi, MarchingCubes(iso=0.0), origin=(DIFF,DIFF,DIFF), widths=(LENGTH, LENGTH, LENGTH))
gy_mesh = Mesh(phi, MarchingTetrahedra(iso=0.0, insidepositive=true), origin=(DIFF,DIFF,DIFF), widths=(LENGTH, LENGTH, LENGTH))
# gy_mesh = Mesh(gyroid_shell, Rect(Vec(0,0,0),Vec(pi*4,pi*4,pi*4)),
#                        MarchingCubes())
print(typeof(gy_mesh))
save("gyroid.ply", gy_mesh)

# view with Makie
# import Makie
# using LinearAlgebra
# m = Makie.mesh(gy_mesh, color=[norm(v) for v in coordinates(gy_mesh)])
# save("v2.png", m)