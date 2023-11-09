
module julia_rings

using JuLIP, LinearAlgebra
include(joinpath(@__DIR__,"rings.jl"))
include(joinpath(@__DIR__,"utilities.jl"))
using .rings, .utilities

export get_nodlnkd, ring_statistics, read_nodlnkd


function get_nodlnkd(atoms::AbstractAtoms, 
                     cutoff::Dict{Tuple{Int, Int},Float64})

    for key in keys(cutoff)
        if (key[2], key[1]) ∉ keys(cutoff)
            cutoff[(key[2], key[1])] = cutoff[key]
        end
    end

    max_cutoff = maximum(values(cutoff))
    nl = neighbourlist(atoms, max_cutoff)
    
    nodlnkd = [Int[] for m in 1:length(atoms.Z)]
    for (i, j, rr) in pairs(nl)
        if norm(rr) < cutoff[(atoms.Z[i], atoms.Z[j])]
            push!(nodlnkd[i], j)
        end
    end

    return nodlnkd
end

function ring_statistics(atoms::AbstractAtoms, cutoff::Dict{Tuple{Int, Int},Float64};
                            refnodes::Union{String, Vector{Int}}="auto",
                            maxlevel=12, no_supercell=false, verbosity=0, kwargs...)
    
    nodlnkd = get_nodlnkd(atoms, cutoff)
    refnodes = find_refnodes(atoms)
    
    return rings.ring_statistics(length(atoms),
                           nodlnkd;
                           refnodes,
                           maxlvl=maxlevel,
                           verbosity=verbosity,
                           kwargs...)
end

end