module utilities

using LinearAlgebra, JuLIP

export find_refnodes, scaled_positions

function halton_base_seq(base, num_points, vol)
    seq = []
    for i in 0:(num_points-1)
        x = (20^3 / vol)^(1/3) * 0.2 # heuristic
        f = 1.0 / base
        while i > 0
            x += f * (i % base)
            i ÷= base
            f /= base
        end
        push!(seq, x)
    end
    return seq
end

function halton_sequence(dim, num_points, vol)
    """Places refnode points quasi-evenly spaced in non-special places
    based on volume heuristic"""
    
    sequence = []
    for d in 0:(dim-1)
        base = 2 + d  # Use different prime bases for each dimension
        push!(sequence, halton_base_seq(base, num_points, vol))
    end
    
    sequence = hcat(sequence...)
    return sequence
end

function place_points_in_cube(num_points, vol)
    dim = 3  # Number of dimensions
    points = halton_sequence(dim, num_points, vol)
    return points
end

function normalise_fraction(frac::Float64)
    res = frac % 1.0
    if res < 0.0
        res += 1.0
    end
    return res
end

function scaled_positions(atoms)
    fractional = [normalise_fraction.(inv(atoms.cell') * i) for i in atoms.X]
    return fractional
end

function find_refnodes(ats)
    
    vol = volume(ats)
    num_points = ceil(Int, (vol/(20^3))) + 2
    
    points = place_points_in_cube(num_points, vol)
    spos = scaled_positions(ats)

    refnodes = []
    for v in eachrow(points)
        diff = [i .- v for i in spos]
        refnode = argmin(vec(norm.(diff)))
        push!(refnodes, refnode)
    end
    
    push!(refnodes, 1) # just a convenience padding element, not used
    refnodes = convert(Array{Int}, refnodes)
    @show refnodes

    return refnodes

end

end