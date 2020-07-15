include("type.jl")
# include("utility.jl")
using Meshing
using FileIO # MeshIO should also be installed
using GeometryBasics
import Makie
using LinearAlgebra

function read_ouput(which, cur_freq)
    # fp_out = open("settings/NumCalc/CPU_1_Core_1/be.out/be.$(idx)/pEvalGrid", "r")
    # fv_out = open("settings/NumCalc/CPU_1_Core_1/be.out/be.$(idx)/vEvalGrid", "r")
    if which == "u"     
        open("settings/$(which)/NumCalc/CPU_1_Core_1/be.out/be.$(cur_freq)/pEvalGrid") do file
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

        open("settings/$(which)/NumCalc/CPU_1_Core_1/be.out/be.$(cur_freq)/vEvalGrid") do file
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

function out_exec(which, num_object_points, num_object_faces)
    write_NC(which, num_object_points, num_object_faces)
    print("NU.inp for $which is written.\n")

    ###### The forward boundary value problem (FMBEM)
    cd(()-> run(`./NumCalc`), "./settings/$(which)/NumCalc/CPU_1_Core_1/")
end

function write_object(points, faces)
    if isdir("settings/ObjectMeshes/")
        rm("settings/ObjectMeshes/", recursive=true)
    end
    makeDirec = `mkdir -p settings/ObjectMeshes/Reference/`
    run(makeDirec)
    # write node info.
    f_node = open("settings/ObjectMeshes/Reference/Nodes.txt", "w")
    sz = size(points, 1)
    write(f_node, "$sz\n")
    for i = 0:sz-1
        write(f_node, "$(i+1) $(points[i+1][1]) $(points[i+1][2]) $(points[i+1][3])\n")
    end
    close(f_node)

    # write faces info.
    f_elem = open("settings/ObjectMeshes/Reference/Elements.txt", "w")
    sz = size(faces, 1)
    write(f_elem, "$sz\n")
    for i = 0:sz-1
        write(f_elem, "$(i+1) $(faces[i+1][1]) $(faces[i+1][2]) $(faces[i+1][3]) 0 0 0\n") #counter-clock wise
    end
    close(f_elem)
    
end

function write_eval_mesh(base)
    if isdir("settings/EvaluationGrids/8_HPlane/")
        rm("settings/EvaluationGrids/8_HPlane/", recursive=true)
    end
    makeDirec = `mkdir -p settings/EvaluationGrids/8_HPlane/`
    run(makeDirec)
    f_node = open("settings/EvaluationGrids/8_HPlane/Nodes.txt", "w")
    f_elem = open("settings/EvaluationGrids/8_HPlane/Elements.txt", "w")
    #Number of all element groups
    write(f_node, "$(DEVIDES*DEVIDES*DEVIDES+num_microphones)\n")
    write(f_elem, "$(DEVIDES*DEVIDES*DEVIDES+num_microphones)\n") #elementも書き出してみる

    idx = base
    for i = 1:DEVIDES, j = 1:DEVIDES, k = 1:DEVIDES
        write(f_node, "$(idx) $(DIFF*i) $(DIFF*j) $(DIFF*k)\n")
        write(f_elem, "$idx $idx $idx $idx 2 0 1\n")
        idx += 1
    end
    for p = microphones
        write(f_node, "$(idx) $(p.x) $(p.y) $(p.z)\n")
        write(f_elem, "$idx $idx $idx $idx 2 0 1\n")
        idx += 1
    end

    close(f_node)
    close(f_elem)
end


function write_NC(which, num_object_points, num_object_faces)
    # makeDirec = `mkdir -p settings/$(which)/NumCalc/CPU_1_Core_1/`
    # run(makeDirec)
    if ispath("settings/$(which)/NumCalc/CPU_1_Core_1/NC.inp")
        rm("settings/$(which)/NumCalc/CPU_1_Core_1/NC.inp")
    end
    f = open("settings/$(which)/NumCalc/CPU_1_Core_1/NC.inp", "w")

    #Identifier and Versionnumber
    #Problemname (OK)
    write(f, "Mesh2HRTF 0.2.0\n##\nHead-Related Transfer Functions\n##\n")

    #Controlparameter I (OK)
    write(f, "## Controlparameter I\n")
    write(f, "0 0 0 0 7 0\n")

    #Controlparameter II (OK)
    write(f, "## Controlparameter II\n")
    write(f, "1 $frqNums $frqStep 0.00e+00 1 0 0\n")

    #Frequency Curve (OK)
    write(f, """## Load Frequency Curve
    0 11
    0.000000e+00 0.000000e+00 0.0
    0.100000e+00 0.010000e+04 0.0
    0.200000e+00 0.020000e+04 0.0
    0.300000e+00 0.030000e+04 0.0
    0.400000e+00 0.040000e+04 0.0
    0.500000e+00 0.050000e+04 0.0
    0.600000e+00 0.060000e+04 0.0
    0.700000e+00 0.070000e+04 0.0
    0.800000e+00 0.080000e+04 0.0
    0.900000e+00 0.090000e+04 0.0
    1.000000e+00 0.100000e+04 0.0
    """)

    #Main Parameters I
    num_points = num_object_points + DEVIDES*DEVIDES*DEVIDES + num_microphones
    num_faces = num_object_faces + DEVIDES*DEVIDES*DEVIDES + num_microphones
    write(f, "## 1. Main Parameters I\n")
    write(f, "2 $num_faces $num_points 0 0 2 1 4 0\n")

    #MainParameters II
    write(f, "## 2. Main Parameters II\n")
    write(f, "0 $(num_source) 0 0.0000e+00 0 0 0\n") # # of point sources is 2 here.

    #Main Parameters III (OK)
    write(f, """## 3. Main Parameters III
    0 0 0 0
    """)

    #Main Parameters IV (OK)
    write(f, """## 4. Main Parameters IV
    $(SOUND_SPEED) $(DENSITY) 1.0 0.0e+00 0.0e+00 0.0e+00 0.0e+00
    ##
    """)

    #NODES (OK)
    write(f, "NODES\n")
    write(f, "../../../ObjectMeshes/Reference/Nodes.txt\n")
    write(f, "../../../EvaluationGrids/8_HPlane/Nodes.txt\n")
    write(f, "##\n")

    #ELEMENTS (OK)
    write(f, "ELEMENTS\n")
    write(f, "../../../ObjectMeshes/Reference/Elements.txt\n")
    write(f, "../../../EvaluationGrids/8_HPlane/Elements.txt\n")
    write(f, "##\n")

    #SYMMETRY (OK)
    write(f, """
    # SYMMETRY
    # 0 0 0
    # 0.0000e+00 0.0000e+00 0.0000e+00
    ##
    """)

    #BOUNDARY
    write(f, "BOUNDARY\n")
    # write nothing which is equivalent to "FROM FirstElement TO LastElement VELO 0.0 -1 0.0 -1"
    write(f, "RETU\n")

    #PLANE WAVES (OK)
    write(f, """
    # PLANE WAVES
    # 0 0.0000e+00 -1.0000e+00 0.0000e+00 1.0000e-6 -1 0.0000e+00 -1
    """)

    #POINT SOURCES
    write(f, "POINT SOURCES\n")
    for (k,p) = enumerate(sources)
        write(f, "$k $(p.x) $(p.y) $(p.z) 0.0 -1 $(p.val) -1\n") #strength should be 1/2*freq*p.val, but p.val in order to calculate diffrent freq simulaneously
    end

    #CURVES
    write(f, """# CURVES
    # Frequency Factor 0.0
    """)

    #POSTPROCESS (OK)
    write(f, "POST PROCESS\n")
    #END (OK)
    write(f, "END\n")

    close(f)
end

function remove_files()
    for which = CATG
        rm("settings/$(which)/NumCalc/CPU_1_Core_1/NC.inp")
        rm("settings/$(which)/NumCalc/CPU_1_Core_1/be.out/", recursive=true)
    end
    rm("settings/ObjectMeshes/", recursive=true)
end

####### Patameters 
const LENGTH = 1.0
const DEVIDES = 50

const DIFF = LENGTH/DEVIDES

const OBJ_FUNC_EPS = 1e-5

const CATG = ["u", "lambda"]

const NUM_ITERATION = 100


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
    
    write_object(points, faces)

    #Threads.@threads 
    for s = CATG #threading 
        out_exec(s, size(points, 1), size(faces, 1))
    end

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

    prev_obj_func_val = cur_obj_func_val
    cur_obj_func_val = 0
    for val = mic_p
        cur_obj_func_val += 1/2 * abs2(val)
    end
    
    # if abs(cur_obj_func_val - prev_obj_func_val) < OBJ_FUNC_EPS
    #     break
    # end

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