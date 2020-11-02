using Meshing
using FileIO # MeshIO should also be installed
using GeometryBasics
import Makie
using MeshIO
using LinearAlgebra
using ProgressBars
using Printf
import GR
using Plots

"""
FastBEMの出力を読み取る関数.
output_result.datのformatに依存するので、まだ完成してない.
"""
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

        open("$(which)/output_result.dat") do file
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
    end

end

"""
FastBEMを実行する関数.
whichは"lambda" or "u".
"""
function out_exec(which)
    println("FastBEM for $which started.")

    cd(()-> run(`FastBEM_Acoustics.exe .`), "./$(which)/")
end


# function write_lambda_dat(points, faces, p_microphone)
#     if ispath("lambda/input.dat")
#         rm("lambda/input.dat")
#     end
#     f = open("lambda/input.dat", "w")

#     # Line 1
#     write(f, "A Sphere Model for Acoustic Scattering Analysis\n")

#     # Line 2 (Job Type (Complete/Field Only/ATM/Use ATM); Solver Type (1=FMM/2=ACA/3=CBEM/4=HFBEM); No. Thread)
#     write(f, "Complete 1 4\n")

#     # Line 3
#     write(f, "Full 0 0.d0\n")

#     # Line 4 (Nos. of Boundary Elements/Nodes, Nos. of Field Points/Cells, and No. Panels)
#     write(f, "$(size(points, 1)) $(size(faces, 1)) $(DEVIDES*DEVIDES*DEVIDES+num_microphones) 0 0\n")

#     # Line 5 (No. of plane waves (nplane); User defined sources (1=Yes/0=No))
#     write(f, "0 0\n")

#     # Line 6 (Complex amplitude and direction vector of the plane wave(s))
#     # skip

#     # Line 7 (No. of monopole sources, and no. of dipole sources)
#     write(f, "$num_microphones 0\n")
#     for (i, p) = emumerate(microphones)
#         write(f, "($(real(p_microphone[i])), $(-imag(p_microphone[i]))) $(p.x) $(p.y) $(p.z)\n")
#     end

#     # Line 8 (cspeed, density, ref. pressure, ref. intensity, wavenumber k ratio)
#     write(f, "$SOUND_SPEED $DENSITY 1.0e-6 1.0e-12  0\n")

#     # Line 9
#     write(f, "$frqStart $(frqStart+frqStep*(frqNums-1)) $frqNums 0 0\n")

#     # Line 10
#     write(f, "0 3 1 0 0 One\n")

#     write(f, "\$ Nodes:\n")
#     sz = size(points, 1)
#     for i = 1:sz
#         write(f, "$(i) $(points[i][1]) $(points[i][2]) $(points[i][3])\n")
#     end 

#     write(f, "\$ Elements and Boundary Conditions:\n")
#     sz = size(faces, 1)
#     for i = 1:sz
#         write(f, "$(i) $(faces[i][1]) $(faces[i][2]) $(faces[i][3]) 1 ( 0.0, 0.0 )\n") #counter-clock wise
#     end


#     write(f, "\$ Field Points\n")
#     idx = 1
#     for i = 1:DEVIDES, j = 1:DEVIDES, k = 1:DEVIDES
#         write(f, "$(idx) $(DIFF*i) $(DIFF*j) $(DIFF*k)\n")
#         idx += 1
#     end
#     for p = microphones
#         write(f, "$(idx) $(p.x) $(p.y) $(p.z)\n")
#         idx += 1
#     end 
    
#     write(f, "\$ Field Cells\n")

#     write(f, "\$ End of the File\n")
# end

function write_u_dat(points, faces)
    if ispath("u/input.dat")
        rm("u/input.dat")
    end
    f = open("u/input.dat", "w")

    # Line 1
    write(f, "A Sphere Model for Acoustic Scattering Analysis\n")

    # Line 2 (Job Type (Complete/Field Only/ATM/Use ATM); Solver Type (1=FMM/2=ACA/3=CBEM/4=HFBEM); No. Thread)
    write(f, "Complete 1 1\n")

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
const DEVIDES = 10

const DIFF = LENGTH/DEVIDES

const OBJ_FUNC_EPS = 1e-5

const CATG = ["u", "lambda"]

const NUM_ITERATION = 1


###### Other Physical Parameters
const DENSITY = 1.1839e+00 #density of the propagation medium
const SOUND_SPEED = 346
const SKIPLINE = 4
const IMPEDANCE = 1e7 + 0im

####### Level-set method Parameters
K = 0.1
tau = 2.5*2.5*1e-4
delta_t = 0.005

####### Frequency info
# frqStart = 0.0 #Hz
# frqEnd = 1000.0
# frqStep = 100.0 #Hz
frqStart = SOUND_SPEED
frqEnd = SOUND_SPEED
frqStep = 1.0
FREQ = frqStart:frqStep:frqEnd

####### Evaluation mesh (place of microphones)
num_microphones = 5
microphones = Array{Float64, 2}(undef, num_microphones, 3) #[x,y,z]
# microphones[1,1:3] = [LENGTH/2, LENGTH/2, 2*LENGTH]
# microphones[2,1:3] = [LENGTH/2, LENGTH/2, -LENGTH]
const MICDIFF = 0.02
microphones[1,1:3] = [LENGTH/2, LENGTH/2, -LENGTH]
microphones[2,1:3] = [LENGTH/2+MICDIFF, LENGTH/2, -LENGTH]
microphones[3,1:3] = [LENGTH/2, LENGTH/2+MICDIFF, -LENGTH]
microphones[4,1:3] = [LENGTH/2-MICDIFF, LENGTH/2, -LENGTH]
microphones[5,1:3] = [LENGTH/2, LENGTH/2-MICDIFF, -LENGTH]

mic_p = Array{Complex{Float64}, 1}(undef, num_microphones)
mic_v = Array{Complex{Float64}, 2}(undef, num_microphones, 3)

####### Place of Sound Source
num_source = 1
INTENSITY = 10.0
sources = Array{Float64, 2}(undef, num_source, 4) #[x,y,z,intensity]
# sources[1,1:4] = [2*LENGTH, LENGTH/2, LENGTH/2, INTENSITY] 
# sources[2,1:4] = [-LENGTH, LENGTH/2, LENGTH/2, INTENSITY] 
sources[1,1:4] = [LENGTH/2, LENGTH/2, 2*LENGTH, INTENSITY]
plot_setting()

####### Create Domain
phi = Array{Float64, 3}(undef, DEVIDES, DEVIDES, DEVIDES)

####### Array
u_p = Array{Complex{Float64}, 3}(undef, DEVIDES, DEVIDES, DEVIDES)
u_v = Array{Complex{Float64}, 4}(undef, DEVIDES, DEVIDES, DEVIDES, 3)
lambda_p = Array{Complex{Float64}, 3}(undef, DEVIDES, DEVIDES, DEVIDES)
lambda_v = Array{Complex{Float64}, 4}(undef, DEVIDES, DEVIDES, DEVIDES, 3)


###### Cost 
costs = Array{Float64, 1}(undef, NUM_ITERATION)

####### initial distribution of phi_{0} 
const EPS = 1e-10
const center_x = LENGTH/2.0
const center_y = LENGTH/2.0
const center_z = LENGTH/2.0
const radius = 0.10

for i = 1:DEVIDES, j = 1:DEVIDES, k = 1:DEVIDES
    dt = sqrt( (i*DIFF - center_x)*(i*DIFF - center_x) + (j*DIFF - center_y) * (j*DIFF - center_y) + (k*DIFF - center_z) * (k*DIFF - center_z))
    phi[i,j,k] = radius - dt # inside positive
    # if dt < radius - EPS
    #     phi[i,j,k] = -0.05
    # else
    #     phi[i,j,k] = 0.05
    # end
end

write_observation_points()
write_command_file_for_u()

printstyled("--------------------------------------------------------\n"; color = :red)
println("Domain        : Cubic [0, $LENGTH] * [0, $LENGTH] * [0, $LENGTH]")
println("Discretization: $DEVIDES per one axis")
println("Frequency     : ", collect(FREQ))
printstyled("--------------------------------------------------------\n"; color = :red)
print("Initialization DONE\n\n")

print("Evaluation Grid WRITTEN\n\n")

prev_obj_func_val = typemax(Float64)
cur_obj_func_val = typemax(Float64)

for cur_iter = ProgressBar(1:NUM_ITERATION)
    global prev_obj_func_val, cur_obj_func_val
    printstyled("--------------------------------------------------------\n"; color = :red)
    print("Current Iteration: $cur_iter\n\n")
    printstyled("--------------------------------------------------------\n"; color = :red)
    global phi
    global cur_obj_func_val, u_v, u_p, lambda_p, lambda_v


    ###### remove files created in the last iteration ########
    remove_files()

    ###### generate boundary mesh
    # for Marching Tetrahedra  ### eps=SOUND_SPEED/frqEnd/6.0
    gy_mesh = Mesh(phi, MarchingTetrahedra(iso=0.0, eps=LENGTH/DEVIDES/4.0, insidepositive=true), origin=(DIFF,DIFF,DIFF), widths=(LENGTH, LENGTH, LENGTH))
    # gy_mesh = Mesh(phi, MarchingTetrahedra(iso=0.0), origin=(DIFF,DIFF,DIFF), widths=(LENGTH, LENGTH, LENGTH))
    xyplot(cur_iter)
    save("./WAON/topology.stl", gy_mesh)
    save("./topology_output/topology$cur_iter.stl", gy_mesh)

    # waon simulation of u
    printstyled("--------------------------------------------------------\n"; color = :red)
    out_exec("u")
    printstyled("--------------------------------------------------------\n"; color = :red)


    prev_obj_func_val = cur_obj_func_val
    cur_obj_func_val = 0
    for (idx, freq) = enumerate(FREQ)

        #Initialization
        for i = 1:num_microphones
            mic_p[i] = 0.0
        end
        read_ouput_only_mic_u(idx)

        for val = mic_p
            cur_obj_func_val += 1/2 * abs2(val)
        end
    end

    costs[cur_iter] = cur_obj_func_val
    save_costs(costs[1:cur_iter], cur_iter)

    # println("currect value of the objective function: $(cur_obj_func_val)")
    # if abs(cur_obj_func_val - prev_obj_func_val) < OBJ_FUNC_EPS
    #     println("finished after $(cur_iter) iterations.")
    #     break
    # end

    # waon simulation of lambda
    write_command_file_for_lambda()
    printstyled("--------------------------------------------------------\n"; color = :red)
    out_exec("lambda")
    printstyled("--------------------------------------------------------\n"; color = :red)

    tp_deriv = zeros(Float64, (DEVIDES, DEVIDES, DEVIDES))
    ###### The objective function is evaluated by the FMBEM
    for (idx, freq) = enumerate(FREQ)
        K_w = 2*pi*freq/SOUND_SPEED # wave number

        read_ouput_u(idx)
        read_ouput_lambda(idx)

        ###### Topological derivative in Eq. (42) is evaluated on the nodes of the finite elements.
        for i = 1:DEVIDES, j = 1:DEVIDES, k = 1:DEVIDES
            sm = complex(0.0)
            for idx = 1:3
                sm += lambda_v[i,j,k,idx] * u_v[i,j,k,idx]
            end
            tp_deriv[i,j,k] +=  real(3.0/2.0*sm - K_w*K_w*lambda_p[i,j,k]*u_p[i,j,k])
            # println(tp_deriv[i,j,k])
            if isnan(tp_deriv[i,j,k])
                println(DEVIDES*i, " ", DEVIDES*j, " ", DEVIDES*k)
                println(num_microphones+(i*DEVIDES*DEVIDES+j*DEVIDES+k)*3+SKIPLINE)
                for idx = 1:3
                    println("lambda_v_$idx", lambda_v[i,j,k,idx], ", u_v_$idx", u_v[i,j,k,idx])
                end
                println("sm = ", sm)
                println("lambda_p = ", lambda_p[i,j,k])
                println("u_p = ", u_p[i,j,k])
            end
        end
    end

    next_phi = copy(phi)
    ###### Distribution of the level set function φ is obtained as the solution of the boundary value problem
    for k = 2:DEVIDES-1, j = 2:DEVIDES-1, i = 2:DEVIDES-1
        # ここで複雑な形状の時はupdateしない部分を作ればよい.
        laplacian = (phi[i+1,j,k] + phi[i-1,j,k] + phi[i,j+1,k] + phi[i,j-1,k] + phi[i,j,k+1] + phi[i,j,k-1] - 6*phi[i,j,k]) / (DIFF*DIFF)
        next_phi[i,j,k] = phi[i,j,k] + delta_t * (-K * tp_deriv[i,j,k] + tau * laplacian)
        # if next_phi[i,j,k] > 100 
        #     println(DEVIDES*i, DEVIDES*j, DEVIDES*k)
        # end
        if isnan(next_phi[i,j,k])
            println("DETECT NaN: $i,$j,$k")
        end
    end
    phi = next_phi

end