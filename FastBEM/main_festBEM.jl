include("type.jl")
# include("utility.jl")
using Meshing
using FileIO # MeshIO should also be installed
using GeometryBasics
import Makie
using LinearAlgebra

function read_ouput(which, cur_freq)
    if which == "u"     
        open("$(which)/output_result.dat") do file
            for l in enumerate(eachline(file))
                if l[1] <= DEVIDES*DEVIDES*DEVIDES
                    idx = l[1] - 1
                    a = parse.(Float64, split(l[2]))
                    u_p[floor(Int, idx/DEVIDES/DEVIDES)+1, floor(Int, idx/DEVIDES)%DEVIDES+1, idx%DEVIDES+1] = a[2] + a[3]*im
                else
                    idx = l[1] - DEVIDES*DEVIDES*DEVIDES
                    a = parse.(Float64, split(l[2]))
                    mic_p[idx] = a[2] + a[3]*im
                end                     
            end
        end

        open("$") do file
            for l in enumerate(eachline(file))
                if l[1] <= DEVIDES*DEVIDES*DEVIDES
                    idx = l[1] - 1
                    a = parse.(Float64, split(l[2]))
                    u_v[floor(Int, idx/DEVIDES/DEVIDES)+1, floor(Int, idx/DEVIDES)%DEVIDES+1, idx%DEVIDES+1] = a[2] + a[3]*im
                else
                    idx = l[1] - DEVIDES*DEVIDES*DEVIDES
                    a = parse.(Float64, split(l[2]))
                    mic_v[idx] = a[2] + a[3]*im
                    # print(mic_v[idx])
                end 
                    
            end
        end
    elseif which == "lambda"
        open("settings/$(which)/NumCalc/CPU_1_Core_1/be.out/be.$(cur_freq)/pEvalGrid") do file
            for l in enumerate(eachline(file))
                if l[1] <= DEVIDES*DEVIDES*DEVIDES
                    idx = l[1] - 1
                    a = parse.(Float64, split(l[2]))
                    lambda_p[floor(Int, idx/DEVIDES/DEVIDES)+1, floor(Int, idx/DEVIDES)%DEVIDES+1, idx%DEVIDES+1] = a[2] + a[3]*im
                end
            end
        end

        open("settings/$(which)/NumCalc/CPU_1_Core_1/be.out/be.$(cur_freq)/vEvalGrid") do file
            for l in enumerate(eachline(file))
                if l[1] <= DEVIDES*DEVIDES*DEVIDES
                    idx = l[1] - 1
                    a = parse.(Float64, split(l[2]))
                    lambda_v[floor(Int, idx/DEVIDES/DEVIDES)+1, floor(Int, idx/DEVIDES)%DEVIDES+1, idx%DEVIDES+1] = a[2] + a[3]*im
                end
            end
        end
    end

end

function out_exec(which)
    ###### The forward boundary value problem (FMBEM)
end


function write_lambda_dat(points, faces, p_microphone)
    if ispath("lambda/input.dat")
        rm("lambda/input.dat")
    end
    f = open("lambda/input.dat", "w")

    # Line 1
    write(f, "A Sphere Model for Acoustic Scattering Analysis\n")

    # Line 2 (Job Type (Complete/Field Only/ATM/Use ATM); Solver Type (1=FMM/2=ACA/3=CBEM/4=HFBEM); No. Thread)
    write(f, "Complete 1 4\n")

    # Line 3
    write(f, "Full 0 0.d0\n")

    # Line 4 (Nos. of Boundary Elements/Nodes, Nos. of Field Points/Cells, and No. Panels)
    write(f, "$(size(points, 1)) $(size(faces, 1)) $(DEVIDES*DEVIDES*DEVIDES+num_microphones) 0 0\n")

    # Line 5 (No. of plane waves (nplane); User defined sources (1=Yes/0=No))
    write(f, "0 0\n")

    # Line 6 (Complex amplitude and direction vector of the plane wave(s))
    # skip

    # Line 7 (No. of monopole sources, and no. of dipole sources)
    write(f, "$num_microphones 0\n")
    for (i, p) = emumerate(microphones)
        write(f, "($(real(p_microphone[i])), $(-imag(p_microphone[i]))) $(p.x) $(p.y) $(p.z)\n")
    end

    # Line 8 (cspeed, density, ref. pressure, ref. intensity, wavenumber k ratio)
    write(f, "$SOUND_SPEED $DENSITY 1.0e-6 1.0e-12  0\n")

    # Line 9
    write(f, "$frqStart $(frqStart+frqStep*(frqNums-1)) $frqNums 0 0\n")

    # Line 10
    write(f, "0 3 1 0 0 One\n")

    write(f, "\$ Nodes:\n")
    sz = size(points, 1)
    for i = 1:sz
        write(f, "$(i) $(points[i][1]) $(points[i][2]) $(points[i][3])\n")
    end 

    write(f, "\$ Elements and Boundary Conditions:\n")
    sz = size(faces, 1)
    for i = 1:sz
        write(f, "$(i) $(faces[i][1]) $(faces[i][2]) $(faces[i][3]) 1 ( 0.0, 0.0 )\n") #counter-clock wise
    end


    write(f, "\$ Field Points\n")
    idx = 1
    for i = 1:DEVIDES, j = 1:DEVIDES, k = 1:DEVIDES
        write(f, "$(idx) $(DIFF*i) $(DIFF*j) $(DIFF*k)\n")
        idx += 1
    end
    for p = microphones
        write(f, "$(idx) $(p.x) $(p.y) $(p.z)\n")
        idx += 1
    end 
    
    write(f, "\$ Field Cells\n")

    write(f, "\$ End of the File\n")
end

function write_u_dat(points, faces)
    if ispath("u/input.dat")
        rm("u/input.dat")
    end
    f = open("u/input.dat", "w")

    # Line 1
    write(f, "A Sphere Model for Acoustic Scattering Analysis\n")

    # Line 2 (Job Type (Complete/Field Only/ATM/Use ATM); Solver Type (1=FMM/2=ACA/3=CBEM/4=HFBEM); No. Thread)
    write(f, "Complete 1 4\n")

    # Line 3
    write(f, "Full 0 0.d0\n")

    # Line 4 (Nos. of Boundary Elements/Nodes, Nos. of Field Points/Cells, and No. Panels)
    write(f, "$(size(points, 1)) $(size(faces, 1)) $(DEVIDES*DEVIDES*DEVIDES+num_microphones) 0 0\n")

    # Line 5 (No. of plane waves (nplane); User defined sources (1=Yes/0=No))
    write(f, "0 0\n")

    # Line 6 (Complex amplitude and direction vector of the plane wave(s))
    # skip

    # Line 7 (No. of monopole sources, and no. of dipole sources)
    write(f, "$num_source 0\n")
    for p = sources
        write(f, "($(p.val), 0) $(p.x) $(p.y) $(p.z)\n")
    end

    # Line 8 (cspeed, density, ref. pressure, ref. intensity, wavenumber k ratio)
    write(f, "$SOUND_SPEED $DENSITY 1.0e-6 1.0e-12  0\n")

    # Line 9
    write(f, "$frqStart $(frqStart + frqStep*(frqNums-1)) $frqNums 0 0\n")

    # Line 10
    write(f, "0 3 1 0 0 One\n")

    write(f, "\$ Nodes:\n")
    sz = size(points, 1)
    for i = 1:sz
        write(f, "$(i) $(points[i][1]) $(points[i][2]) $(points[i][3])\n")
    end 

    write(f, "\$ Elements and Boundary Conditions:\n")
    sz = size(faces, 1)
    for i = 1:sz
        write(f, "$(i) $(faces[i][1]) $(faces[i][2]) $(faces[i][3]) 1 ( 0.0, 0.0 )\n") #counter-clock wise
    end


    write(f, "\$ Field Points\n")
    idx = 1
    for i = 1:DEVIDES, j = 1:DEVIDES, k = 1:DEVIDES
        write(f, "$(idx) $(DIFF*i) $(DIFF*j) $(DIFF*k)\n")
        idx += 1
    end
    for p = microphones
        write(f, "$(idx) $(p.x) $(p.y) $(p.z)\n")
        idx += 1
    end 
    
    write(f, "\$ Field Cells\n")

    write(f, "\$ End of the File\n")

end


function remove_files()
    for which = CATG
        rm("$(which)/")
    end
    rm("settings/ObjectMeshes/", recursive=true)
end

####### Patameters 
const LENGTH = 1.0
const DEVIDES = 50

const DIFF = LENGTH/DEVIDES

const OBJ_FUNC_EPS = 1e-5

const CATG = ["u", "lambda"]

const NUM_ITERATION = 1 


###### Other Physical Parameters
const DENSITY = 1.1839e+00 #density of the propagation medium
const SOUND_SPEED = 346.18

####### Level-set method Parameters
K = 1e5
tau = 2.5*2.5*1e-4
delta_t = 0.005

####### Frequency info
frqStart = SOUND_SPEED #Hz
frqStep = 100.0 #Hz
frqNums = 1

####### Evaluation mesh (place of microphones)
num_microphones = 2
m1 = Point(LENGTH/2, LENGTH/2, 2*LENGTH, 0)
m2 = Point(LENGTH/2, LENGTH/2, -LENGTH, 0)
microphones = [m1, m2]

mic_p = Array{Complex{Float64}, 1}(undef, num_microphones)
mic_v = Array{Complex{Float64}, 1}(undef, num_microphones)

####### Place of Sound Source
num_source = 2
INTENSITY = 70
s1 = Point(2*LENGTH, LENGTH/2, LENGTH/2, INTENSITY)
s2 = Point(-LENGTH, LENGTH/2, LENGTH/2, INTENSITY)
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
const center_x = 0.50
const center_y = 0.50
const center_z = 0.50
const radius = 0.10

for i = 1:DEVIDES, j = 1:DEVIDES, k = 1:DEVIDES
    dt = sqrt( (i*DIFF - center_x)*(i*DIFF - center_x) + (j*DIFF - center_y) * (j*DIFF - center_y) + (k*DIFF - center_z) * (k*DIFF - center_z))
    # phi[i,j,k] = radius - dt # inside positive
    if dt < radius - EPS
        phi[i,j,k] = -0.05
    else
        phi[i,j,k] = 0.05
    end
end

println("--------------------------------------------------------")
println("Domain        : Cubic [0, $LENGTH] * [0, $LENGTH] * [0, $LENGTH]")
println("Discretization: $DEVIDES per one axis")
println("Frequency     : ", [frqStart+i*frqStep for i = 1:frqNums])
println("--------------------------------------------------------")
print("Initialization DONE\n\n")

write_eval_mesh(DEVIDES*DEVIDES*DEVIDES*15)
print("Evaluation Grid WRITTEN\n\n")

prev_obj_func_val = typemax(Float64)
cur_obj_func_val = typemax(Float64)

for cur_iter = 1:NUM_ITERATION
    print("Current Iteration: $cur_iter\n\n")
    global phi
    global cur_obj_func_val, u_v, u_p, lambda_p, lambda_v
    ###### generate boundary mesh
    
    # for Marching Cubes method
    # gy_mesh = Mesh(phi, MarchingCubes(iso=0.0), origin=(DIFF,DIFF,DIFF), widths=(LENGTH, LENGTH, LENGTH))
    # save("gyroid_$cur_iter.ply", gy_mesh)
    # points,faces = isosurface(phi, MarchingCubes(iso=0.0), origin=(DIFF,DIFF,DIFF), widths=(LENGTH, LENGTH, LENGTH))

    # for Marching Tetrahedra
    gy_mesh = Mesh(phi, MarchingTetrahedra(iso=0.0, insidepositive=true), origin=(DIFF,DIFF,DIFF), widths=(LENGTH, LENGTH, LENGTH))
    save("gyroid_tetra_$cur_iter.ply", gy_mesh)
    points,faces = isosurface(phi, MarchingTetrahedra(iso=0.0, insidepositive=true), origin=(DIFF,DIFF,DIFF), widths=(LENGTH, LENGTH, LENGTH))
    # points,faces = isosurface(phi, MarchingTetrahedra(iso=0.0, eps=1e-4, insidepositive=true), origin=(DIFF,DIFF,DIFF), widths=(LENGTH, LENGTH, LENGTH))

    println("Meshed : (",size(points, 1), " points, ", size(faces,1), " faces)")
    
    write_u_dat(points, faces)

    #Threads.@threads 
    cd(()-> run(`./NumCalc`), "./u/")


    for cur_freq_c = 1:frqNums
        read_mic()
    end

    prev_obj_func_val = cur_obj_func_val
    cur_obj_func_val = 0
    for val = mic_p
        cur_obj_func_val += 1/2 * abs2(val)
    end
    
    if abs(cur_obj_func_val - prev_obj_func_val) < OBJ_FUNC_EPS
        break
    end
    

    write_lambda_dat(points, faces, mic_p)

    tp_deriv = zeros(Float64, (DEVIDES, DEVIDES, DEVIDES))
    ###### The objective function is evaluated by the FMBEM
    for cur_freq_c = 1:frqNums
        K_w = 2*pi*(frqStart + (cur_freq_c-1)*frqStep)*1000/SOUND_SPEED # wave number
        for s = CATG
            read_ouput(s, cur_freq_c)
        end

        ###### Topological derivative in Eq. (42) is evaluated on the nodes of the finite elements.
        for k = 1:DEVIDES, j = 1:DEVIDES, i = 1:DEVIDES
            tp_deriv[i,j,k] +=  real(3/2*lambda_v[i,j,k]*u_v[i,j,k] - K_w*K_w*lambda_p[i,j,k]*u_p[i,j,k])
            # if tp_deriv[i,j,k] != 0.0
            #     print(tp_deriv[i,j,k])
            # end
        end
    end


    next_phi = copy(phi)
    ###### Distribution of the level set function φ is obtained as the solution of the boundary value problem
    for k = 2:DEVIDES-1, j = 2:DEVIDES-1, i = 2:DEVIDES-1
        # ここで複雑な形状の時はupdateしない部分を作ればよい.
        laplacian = phi[i+1,j,k] + phi[i-1,j,k] + phi[i,j+1,k] + phi[i,j-1,k] + phi[i,j,k+1] + phi[i,j,k-1] - 6*phi[i,j,k]
        next_phi[i,j,k] = phi[i,j,k] + delta_t * (-K * tp_deriv[i,j,k] + tau * laplacian)
    end
    phi = next_phi
    ##### remove the directory
    remove_files()
end

gy_mesh = Mesh(phi, MarchingCubes(iso=0.0), origin=(DIFF,DIFF,DIFF), widths=(LENGTH, LENGTH, LENGTH))
m = Makie.mesh(gy_mesh)
save("v2.png", m)