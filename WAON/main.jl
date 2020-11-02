# include("type.jl")
# include("utility.jl")
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
gr()  # for plot 
# ioff() # not to plot interactively

function read_ouput_u(freq)
    """
    freq is not actual frequency but index.
    """
    which = "u"
    for i = 1:2 # 実部と虚部
        filename = @sprintf("./WAON/%s/output/oFPM%04d.s0%d", which, freq, i)
        open(filename) do file
            for l in enumerate(eachline(file))
                if l[1] <= SKIPLINE
                    # println(l[2])
                    continue
                elseif l[1] <= SKIPLINE + num_microphones
                    # println(l[2])
                    # idx = l[1] - SKIPLINE
                    # a = parse.(Float64, split(l[2]))
                    # if i==1
                    #     mic_p[idx] = a[1]
                    # else 
                    #     mic_p[idx] += a[1]*im
                    # end
                elseif l[1] <= SKIPLINE + num_microphones + DEVIDES*DEVIDES*DEVIDES
                    idx = l[1] - SKIPLINE - num_microphones - 1
                    a = parse.(Float64, split(l[2]))
                    # println(floor(Int, idx/DEVIDES/DEVIDES)+1, floor(Int, idx/DEVIDES)%DEVIDES+1, idx%DEVIDES+1)
                    if i==1
                        u_p[floor(Int, idx/DEVIDES/DEVIDES)+1, floor(Int, idx/DEVIDES)%DEVIDES+1, idx%DEVIDES+1] = a[1]
                    else
                        u_p[floor(Int, idx/DEVIDES/DEVIDES)+1, floor(Int, idx/DEVIDES)%DEVIDES+1, idx%DEVIDES+1] += a[1]*im
                    end
                else
                    println(l[1], " ", l[2])
                end
            end
        end
    end
    for i = 1:2
        filename = @sprintf("./WAON/%s/output/oFPM%04d.v0%d", which, freq, i)
        open(filename) do file
            for l in enumerate(eachline(file))
                if l[1] <= SKIPLINE
                    continue
                elseif l[1] <= SKIPLINE + num_microphones * 3
                    # idx = l[1] - SKIPLINE
                    # a = parse.(Float64, split(l[2]))
                    # if i==1
                    #     mic_v[idx] = a[1]
                    # else 
                    #     mic_v[idx] += a[1]*im
                    # end
                else
                    idx = div(l[1] - SKIPLINE - num_microphones * 3 - 1, 3)
                    xyz = rem(l[1] - SKIPLINE - num_microphones * 3 - 1, 3)
                    a = parse.(Float64, split(l[2]))
                    if i==1
                        u_v[floor(Int, idx/DEVIDES/DEVIDES)+1, floor(Int, idx/DEVIDES)%DEVIDES+1, idx%DEVIDES+1, xyz+1] = a[1]
                    else
                        u_v[floor(Int, idx/DEVIDES/DEVIDES)+1, floor(Int, idx/DEVIDES)%DEVIDES+1, idx%DEVIDES+1, xyz+1] += a[1]*im
                    end
                end
            end
        end
    end
end


function read_ouput_lambda(freq)
    """
    freq is not actual frequency but index.
    """
    which = "lambda"
    for i = 1:2
        filename = @sprintf("./WAON/%s/output/oFPM%04d.s0%d", which, freq, i)
        open(filename) do file
            for l in enumerate(eachline(file))
                if l[1] <= SKIPLINE
                    continue
                elseif l[1] <= SKIPLINE + num_microphones
                    # idx = l[1] - SKIPLINE
                    # a = parse.(Float64, split(l[2]))
                    # if i==1
                    #     mic_p[idx] = a[1]
                    # else 
                    #     mic_p[idx] += a[1]*im
                    # end
                else
                    idx = l[1] - SKIPLINE - num_microphones - 1
                    a = parse.(Float64, split(l[2]))
                    if i==1
                        lambda_p[floor(Int, idx/DEVIDES/DEVIDES)+1, floor(Int, idx/DEVIDES)%DEVIDES+1, idx%DEVIDES+1] = a[1]
                    else
                        lambda_p[floor(Int, idx/DEVIDES/DEVIDES)+1, floor(Int, idx/DEVIDES)%DEVIDES+1, idx%DEVIDES+1] += a[1]*im
                    end
                end
            end
        end
    end
    for i = 1:2
        filename = @sprintf("./WAON/%s/output/oFPM%04d.v0%d", which, freq, i)
        open(filename) do file
            for l in enumerate(eachline(file))
                if l[1] <= SKIPLINE
                    continue
                elseif l[1] <= SKIPLINE + num_microphones * 3
                    # idx = l[1] - SKIPLINE
                    # a = parse.(Float64, split(l[2]))
                    # if i==1
                    #     mic_v[idx] = a[1]
                    # else 
                    #     mic_v[idx] += a[1]*im
                    # end
                else
                    idx = div(l[1] - SKIPLINE - num_microphones * 3 - 1, 3)
                    xyz = rem(l[1] - SKIPLINE - num_microphones * 3 - 1, 3)
                    a = parse.(Float64, split(l[2]))
                    if i==1
                        lambda_v[floor(Int, idx/DEVIDES/DEVIDES)+1, floor(Int, idx/DEVIDES)%DEVIDES+1, idx%DEVIDES+1, xyz+1] = a[1]
                    else
                        lambda_v[floor(Int, idx/DEVIDES/DEVIDES)+1, floor(Int, idx/DEVIDES)%DEVIDES+1, idx%DEVIDES+1, xyz+1] += a[1]*im
                    end
                end
            end
        end
    end
end

function read_ouput_only_mic_u(freq)
    """
    freq is not actual frequency but index.
    """
    global mic_p
    which = "u"
    for i = 1:2
        filename = @sprintf("./WAON/%s/output/oFPM%04d.s0%d", which, freq, i)
        open(filename) do file
            for l in enumerate(eachline(file))
                if l[1] <= SKIPLINE
                    continue
                elseif l[1] <= SKIPLINE + num_microphones
                    idx = l[1] - SKIPLINE
                    a = parse.(Float64, split(l[2]))
                    if i==1
                        mic_p[idx] = a[1]
                    else 
                        mic_p[idx] += a[1]*im
                    end
                    # println(mic_p[idx])
                else
                    break
                end
            end
        end
    end
end

function out_exec(which)
    println("WAON for $which started.")

    ###### The forward boundary value problem (FMBEM)
    #cd(()-> run(`ls`), "./WAON")

    cd(()-> run(`waonNogui.bat ./$(which)/input.txt`), "./WAON")
end

function write_command_file_for_u()
    """
    uに関してコマンドファイルを吐き出す. 
    !!! 実行時の始めに一度だけ実行する.
    つまり、ここの設定はイテレーションにおいて不変.
    """
    label = "u"
    if ispath("./WAON/$label/input.txt")
        rm("./WAON/$label/input.txt")
    end
    if ispath("./WAON/$label/input.wdb")
        rm("./WAON/$label/input.wdb")
    end

    f = open("./WAON/$label/input.txt", "w")
    
    # chapter 2
    write(f, "NewDatabase File \"./$label/input.wdb\" Return\n")

    # chapter 3
    write(f, "
    Model 
        Formulation DLF
        Return
    ")

    # chapter 4
    write(f, "
    Mesh
        File \"topology.stl\"
        Format \"STL\"
        Unit m
        Return
    ")

    # chpter 5
    write(f, "
    ")

    # chapter 6
    write(f, "
    Material
        Name \"Air\"
        Speed 346.180+i*0.0
        Rho 1.1839+i*0.0
        Return
    ")

    # chapter 7
    write(f, "
    Assign
        Elements All
        PositiveMaterial \"Air\"
        NegativeMaterial \"Air\"
        Return
    ")

    # chapter 8
    write(f, "
    ")

    # chapter 9
    write(f, "
    ")

    # chapter 10
    for i = 1:num_source
        write(f, "
        Source
            Name \"Source$i\"
            Type SphericalWave
            Location $(sources[i,1]), $(sources[i,2]), $(sources[i,3]) 
            Value $(sources[i,4])
            Material \"Air\"
            Return
        ")
    end
    # write(f, "
    #     Source 
    #         Name \"boundary1\"
    #         Type Velocity
    #         Elements All
    #         Both
    #         Value 0.0+i*0.0
    #         Return
    # ")

    # chapter 11
    write(f, "
    Absorbent
        Name \"ab1\"
        Type Impedance
        Elements All
        Value $(real(IMPEDANCE))+i*$(imag(IMPEDANCE))
        positive
        Return
    ")

    # chapter 12
    write(f, "
    Point
        File \"./observation_points.csv\"
        Format \"CSV\"
        Unit m
        Material \"Air\"
        Return
    ")

    # chapter 13
    write(f, "
    Parameter
        Iteration 5000
        Tolerance 1E-6
        EstimatedIteration 1000
        IntegrationPoints
            Main 3.0 4 2
            Post 2.0 3 2
        Return
    ")
    write(f, "
    Save Return
    ")

    # chapter 14
    write(f, "
    Analysis
        AllProcess FMBEM
        Frequency "
    )
    for frq = FREQ
        write(f, "$(frq) ")
    end
    
    write(f, "
        Solver GPBiCG
        Parallel
            SMP 2
            FLP 1
            DMP 1
            Return
        Memory 100000
        Output
            Format EnsightAscii
            File ./$label/output/o
            Return
        Return
    ")
    

    # chapter 15
    write(f, "
    Save Return
    ")

    # chapter 16
    write(f, "
    Result
        Output
        File \"./$label/output/main\"
        Field
        Frequency ")
    for frq = FREQ
        write(f, "$(frq) ")
    end
    write(f, "
        Coordinate
        Return
    ")

    write(f, "
    Exit
    ")

    close(f)
end

function write_command_file_for_lambda()
    """
    lambdaに関してコマンドファイルを吐き出す. 
    ここはイテレーションの間変更しなければならない
    なぜなら、sourceの値がuに依存するためである.
    """
    label = "lambda"
    if ispath("./WAON/$label/input.txt")
        rm("./WAON/$label/input.txt")
    end
    if ispath("./WAON/$label/input.wdb")
        rm("./WAON/$label/input.wdb")
    end

    f = open("./WAON/$label/input.txt", "w")
    
    # chapter 2
    write(f, "NewDatabase File \"./$label/input.wdb\" Return\n")

    # chapter 3
    write(f, "
    Model
        Formulation DLF
        Return
    ")

    # chapter 4
    write(f, "
    Mesh
        File \"topology.stl\"
        Format \"STL\"
        Unit m
        Return
    ")

    # chpter 5
    write(f, "
    ")

    # chapter 6
    write(f, "
    Material
        Name \"Air\"
        Speed 346.180+i*0.0
        Rho 1.1839+i*0.0
        Return
    ")

    # chapter 7
    write(f, "
    Assign
        Elements All
        PositiveMaterial \"Air\"
        NegativeMaterial \"Air\"
        Return
    ")

    # chapter 8
    write(f, "
    ")

    # chapter 9
    write(f, "
    ")

    # chapter 10
    for i = 1:num_microphones
        write(f, "
        Source
            Name \"Source$i\"
            Type SphericalWave
            Location $(microphones[i,1]), $(microphones[i,2]), $(microphones[i,3]) 
            Value $(real(mic_p[i]))+i*$(-imag(mic_p[i]))
            Material \"Air\"
            Return
        ")
    end

    # write(f, "
    #     Source 
    #         Name \"boundary1\"
    #         Type Velocity
    #         Elements All
    #         Both
    #         Value 0.0+i*0.0
    #         Return
    # ")

    # chapter 11
    write(f, "
    Absorbent
        Name \"abab1\"
        Type Impedance
        Elements All
        Value $(real(IMPEDANCE))+i*$(imag(IMPEDANCE))
        positive
        Return
    ")

    # chapter 12
    write(f, "
    Point
        File \"./observation_points.csv\"
        Format \"CSV\"
        Unit m
        Material \"Air\"
        Return
    ")

    # chapter 13
    write(f, "
    Parameter
        Iteration 5000
        Tolerance 1E-6
        EstimatedIteration 1000
        IntegrationPoints
            Main 3.0 4 2
            Post 2.0 3 2
        Return
    ")
    write(f, "
    Save Return
    ")

    # chapter 14
    write(f, "
    Analysis
        AllProcess FMBEM
        Frequency "
    )
    for frq = FREQ
        write(f, "$(frq) ")
    end
    
    write(f, "
        Solver GPBiCG
        Parallel
            SMP 2
            FLP 1
            DMP 1
            Return
        Output
            Format EnsightAscii
            File ./$label/output/o
            Return
        Memory 100000
        Return
    ")
    

    # chapter 15
    write(f, "
    Save Return
    ")

    # chapter 16
    write(f, "
    Result
        Output
        File \"./$label/output/main\"
        Field
        Frequency ")
    for frq = FREQ
        write(f, "$(frq) ")
    end
    write(f, "
        Coordinate
        Return
    ")


    write(f, "
    Exit
    ")

    close(f)
end

function write_observation_points()
    """
    ./WAONフォルダにobservation_points.csvを吐き出す.
    !!! これは実行時のはじめに一度呼び出すだけ

    """
    if ispath("./WAON/observation_points.csv")
        rm("./WAON/observation_points.csv")
    end

    f = open("./WAON/observation_points.csv", "w")

    for i = 1:num_microphones
        write(f, "$(microphones[i,1]), $(microphones[i,2]), $(microphones[i,3])\n")
    end
    for i = 1:DEVIDES, j = 1:DEVIDES, k = 1:DEVIDES
        write(f, "$(i*DIFF), $(j*DIFF), $(k*DIFF)\n")
    end

    close(f)
end

function remove_files()
    for which = CATG
        if ispath("WAON/$(which)/output")
            rm("WAON/$(which)/output", recursive=true)
        end
        if ispath("WAON/$(which)/input.wdb")
            rm("WAON/$(which)/input.wdb")
        end
    end
    if ispath("WAON/lambda/input.txt")
        rm("WAON/lambda/input.txt")
    end
end

function save_costs(csts, itr)
    #### draw cost function
    plot(csts, marker=:circle)
    savefig("cost_plot/costs$itr.png")
end 

function xyplot(cur_iter)
    println("ploting xy plane...")
    for i = 1:DEVIDES
        img = heatmap(phi[:,:,i])
        savefig(img, "./topology_output/topology$(cur_iter)_zindex$i.png")
    end
    println("DONE plotting xy plane!!!")
end

function xyzplot(cur_iter)
    println("ploting xy plane...")
    x_plot = Array{Float32, 1}(undef,DEVIDES*DEVIDES*DEVIDES)
    y_plot = Array{Float32, 1}(undef,DEVIDES*DEVIDES*DEVIDES)
    z_plot = Array{Float32, 1}(undef,DEVIDES*DEVIDES*DEVIDES)
    val_plot = Array{Float32, 1}(undef,DEVIDES*DEVIDES*DEVIDES)
    for i=1:DEVIDES, j=1:DEVIDES, k=1:DEVIDES
        x_plot[(i-1)*DEVIDES*DEVIDES+(j-1)*DEVIDES+k] = i*DIFF
        y_plot[(i-1)*DEVIDES*DEVIDES+(j-1)*DEVIDES+k] = j*DIFF
        z_plot[(i-1)*DEVIDES*DEVIDES+(j-1)*DEVIDES+k] = k*DIFF
        val_plot[(i-1)*DEVIDES*DEVIDES+(j-1)*DEVIDES+k] = phi[i,j,k]
    end
    img = GR.scatter3(x_plot, y_plot, z_plot, val_plot)
    GR.savefig("./topology_output/topology$(cur_iter).png")
    println("DONE plotting xy plane!!!")
    
end

function plot_setting()
    img = scatter(append!(microphones[:,1], sources[:,1]), append!(microphones[:,2], sources[:,2]), append!(microphones[:,3], sources[:,3]), grid=true, gridstyle=:dash)
    savefig(img, "./setting.png")
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


# gy_mesh = Mesh(phi, MarchingCubes(iso=0.0), origin=(DIFF,DIFF,DIFF), widths=(LENGTH, LENGTH, LENGTH))
# m = Makie.mesh(gy_mesh)
# save("v2.png", m)