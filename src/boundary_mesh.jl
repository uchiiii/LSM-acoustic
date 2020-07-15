include("type.jl")
include("utility.jl")

function generate_mesh(dom::Domain)
    # pick up boundary Cubics
    cbs = pick_boundary_boxes(dom)

    # for each cubic and for each surface

    # 
end


function pick_boundary_boxes(dom::Domain)
    (X,Y,Z) = size(dom)
    cbs = Cubics(undef, 0)
    inside = Vector{Bool, 3}(undef, X-1, Y-1, Z-1)
    for i = 1:X-1
        for j = 1:Y-1
            for k = 1:Z-1
                
                posit = false
                negat = false
                for l = 0:1
                    for m = 0:1
                        for n = 0:1
                            if dom[i+l, j+m, k+n].val < 0
                                negat = true
                            end
                            if dom[i+l, j+m, k+n].val > 0
                                posit = true
                            end
                        end
                    end
                end
                
                cb = Cubic(undef, 0)
                if posit == true && negat == true
                    # add cubic
                    for l = 0:1
                        for m = 0:1
                            for n = 0:1
                                push!(cb, dom[i+l, j+m, k+n])
                            end
                        end
                    end
                    push!(cbs, cb)
                end
            end
        end
    end

    cbs
end


function added_triangles(cbs::Cubics)
    trs = Triangles(undef, 0)
    mdps = Vector{Point}(undef, 0)
    for cb = cbs
        
        # for vertex 1,2,3,4
        (cur_trs, cur_mdps) = generate_triangle_face(cb[1], cb[2], cb[4], cb[3])
        append!(trs, cur_trs)
        append!(mdps, cur_mdps)

        # for vertex
        if cb[1] <= 0 && cb[2] > 0 && cb[3] > 0 && cb[4] > 0 && cb[5] > 0 && cb[6] > 0 && cb[7] > 0 && cb[8] > 0
            push!(trs, []
    end
end

# points are assumed to be aligned in the clock-wise or counter-clock-wise order.
function generate_triangle_face(p1::Point, p2::Point, p3::Point, p4::Point)
    trs = Triangles(undef, 0)
    mdps = Vector{Point}(undef, 0)

    if p1.val > 0 && p3.val > 0 && p2.val < 0 && p4.val < 0

    end
    if p1.val * p2.val < 0 
        mdp = middle_point(p1,p2)
    end

    if p2.val * p3.val < 0
        mdp = middle_point(p2,p3)
    end
end

function middle_point(p1::Point, p2::Point)
    internal_division_point(p1, p2)
end
