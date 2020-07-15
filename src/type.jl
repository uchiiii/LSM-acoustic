########## Type definition  ############

struct Point
    x::Float64
    y::Float64
    z::Float64
    val::Float64
end

Cubic = Vector{Point} # size = 8
Cubics = Vector{Cubic}

Triangle = Vector{Point} # size = 3
Triangles = Vector{Triangle}

Face = Vector{Point} # size = 4
   
Domain = Array{Point, 3}