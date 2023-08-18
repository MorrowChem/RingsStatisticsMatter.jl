
module rings

    """
    Implementation of algorithm described in:
    
    Yuan, X. & Cormack, A. N. 
    "Efficient Algorithm for Primitive Ring Statistics in Topological Networks." 
    Computational Materials Science, vol. 24, pp. 343-360, 2002.
    """

    using PyCall, JSON, Formatting
    export ring_statistics

    function read_nodlnkd(filename)
        """Helper for reading neighbourlist from a .csv file"""
        file = open(filename, "r")
        nodlnkd = Vector{Vector{Int64}}()
    
        for line in eachline(file)
            values = split(line, ',')
            int_values = [parse(Int64, v) for v in values]
            push!(nodlnkd, int_values)
        end
    
        close(file)
    
        return nodlnkd
    end


    function dijkstra_nonwgt(nodsrc::Int64,
                             lvlreq::Int64,
                             lnks::Vector{Int64},
                             nodlnkd::Vector{Vector{Int64}})
        """
        Perform Dijkstra's algorithm for non-weighted graphs to calculate shortest distances.

        # Arguments
        - `nodsrc::Int64`: Source node index.
        - `lvlreq::Int64`: Required level of distances.
        - `lnks`: Array of integers representing the number of links for each node.
        - `nodlnkd`: Array of arrays representing the linked nodes for each node.
        
        # Returns
        - `lvldist::Vector{Int64}`: Array of shortest distances from the source node.
        
        This function performs Dijkstra's algorithm on a non-weighted graph
        to calculate the shortest distances from a specified source node to all other nodes.
        The algorithm considers the required level of distances (`lvlreq`) 
        and returns an array of shortest distances.
        """
        # initialise
        lvldist = lvlreq + 2 .+ zeros(Int64, length(lnks))
        lvldist[nodsrc] = 0
        quebgn = 0
        quend = 1
        queue = Vector{Int64}(undef, length(lnks))
        queue[1] = nodsrc

        while quebgn < quend
            quebgn += 1
            nodcrt = queue[quebgn]
            lvlprb = lvldist[nodcrt] + 1

            for lnlscrt = 1:lnks[nodcrt]
                nodprb = nodlnkd[nodcrt][lnlscrt]
                if lvldist[nodprb] > lvlprb
                    lvldist[nodprb] = lvlprb

                    if lvlprb < lvlreq
                        quend += 1
                        queue[quend] = nodprb
                    end

                end
            end
        end
        return lvldist
    end


    function strpth_record(nodcrt::Int64,
                           lvlcrt::Int64,
                           lvlprim::Int64, 
                           pths::Int64,
                           strpth::Matrix{Int64},
                           strpthx::Vector{Int64},
                           lvldist::Vector{Int64},
                           lnks::Vector{Int64},
                           nodlnkd::Vector{Vector{Int64}})
        """
        Record shortest paths recursively from a current node to a set of prime nodes.

        # Arguments
        - `nodcrt::Int64`: Current node index.
        - `lvlcrt::Int64`: Dijkstra-level of the current node.
        - `lvlprim::Int64`: Dijkstra-level of the prime nodes.
        - `pths::Int64`: Number of paths recorded.
        - `strpth`: Array storing recorded paths.
        - `strpthx`: Auxiliary array for path recording.
        - `lvldist`: Array of shortest distances for source node.
        - `lnks`: Array of integers representing the number of links for each node.
        - `nodlnkd`: Array of arrays representing the linked nodes for each node.

        # Returns
        - `pths::Int64`: Updated number of paths recorded.
        - `strpth`: Updated array storing recorded paths.
        - `strpthx`: Updated auxiliary array for path recording.

        This function recursively records the shortest paths from a given current node to a set of prime nodes. The prime nodes are determined by the `lvlprim` parameter. The recorded paths are stored in the `strpth` array, and an auxiliary array `strpthx` is used for path recording. The function returns the updated number of recorded paths and arrays.
        """
        for lnkscrt in 1:lnks[nodcrt]
            nodprb = nodlnkd[nodcrt][lnkscrt]
            lvlprb = lvldist[nodprb]

            if lvlprb == 0
                pths += 1
                for lvl in 1:lvlprim
                    strpth[pths, lvl] = strpthx[lvl]
                end

            elseif lvlprb == lvlcrt - 1
                strpthx[lvlprb] = nodprb
                pths, strpth, strpthx = strpth_record(nodprb, lvlprb, lvlprim, 
                                                      pths, strpth, strpthx, 
                                                      lvldist, lnks, nodlnkd)
            end
        end
    return pths, strpth, strpthx
    end


    function prime_ring(iodd::Int64, 
                        lvlprim::Int64, 
                        lnks::Vector{Int64}, 
                        lvldist::Vector{Int64}, 
                        lvlref::Vector{Vector{Int64}},
                        ringstat::Vector{Float64}, 
                        nodlnkd::Vector{Vector{Int64}}, 
                        querng::Vector{Vector{Int64}}, 
                        srtpth1::Matrix{Int64},
                        srtpth2::Matrix{Int64}, 
                        prime1::Int64, 
                        prime2::Int64, 
                        nodsrc::Int64, 
                        ngf::Vector{Int64}, 
                        rings::Vector{Vector{Vector{Int64}}},
                        maxlvl::Int64)
        """
        Decide if a proposed ring is a prime ring

        Arguments:
        - `iodd::Int64`: 1 for odd, 0 for even rings
        - `lvlprim::Int64`: Dijkstra level of the prime-mid-node.
        - `lnks`: array for the order of each node (i.e. number of bonds on each atom)
        - `lvldist`: array of Dijkstra-levels for each node relative to current source
        - `lvlref`: array of Dijkstra-levels for each node relative to reference nodes
        - `ringstat`: Ring statistics array for updating, if prime ring found
        - `nodlnkd`: neighbourlist for each node
        - `querng`: Queue for rings to check in this run of prime_ring
        - `srtpth1`: Shortest-distance path for prime-mid-node no. 1
        - `srtpth2`: Shortest-distance path for prime-mid-node no. 2 (same as srtpth1 if even ring)
        - `prime1::Int64`: Prime-mid-node 1 index
        - `prime2::Int64`: Prime-mid-node 2 index (same as 1 if even)
        - `nodsrc::Int64`: Current source node
        - `ngf`: rejection count for odd and even for updating
        - `rings`: lists of rings of each size for updating
        - `maxlvl::Int64`: Max level for Dijkstra

        Returns:
        - A tuple `(ngf, rings, ringstat)`.
        -- ngf - number of rejected rings (even, odd)
        -- rings - set of size-separated rings found described by nodes
        -- ringstat - counts of each ring, normalised by size of ring (each ring is found n times where n=ring_size)

        """
        rngs = size(querng)[1]
        for irng in 1:rngs
            pth1 = querng[irng][1]
            pth2 = querng[irng][2]
            
            if pth1 > 0
                for lvlmax in lvlprim:lvlprim+iodd
                    for lvlchk in 1:lvlmax-1
                        nodchk=srtpth1[pth1,lvlchk]
                        nodmid=srtpth2[pth2,lvlmax-lvlchk]

                        if nodmid == 0
                            nodmid = prime2
                        end
                        if nodchk == 0
                            nodchk = prime1
                        end

                        
                        lr = length(lvlref)
                        limit = Matrix{Int64}(undef, lr, 2)

                        limit[1:lr, 1] = [lvlref[i][nodmid] > maxlvl ?  10000 : lvlref[i][nodmid]+lvlprim-1 for i in 1:lr]
                        limit[1:lr, 2] = [lvlref[i][nodmid] > maxlvl ? -10000 : lvlref[i][nodmid]-lvlprim+1 for i in 1:lr]

                        limit[lr, 1] = lvldist[nodmid] + lvlprim - 1
                        limit[lr, 2] = lvldist[nodmid] - lvlprim + 1

                        goal_found = false
                        goal_found = pair_search(nodchk, nodmid, limit,
                                                 goal_found, lvlref, lvldist,
                                                 lnks, nodlnkd, maxlvl)

                        if goal_found
                            goal_found = false
                            ngf[iodd+1] += 1

                            for irgx in irng+1:rngs
                                p1x = querng[irgx][1]
                                p2x = querng[irgx][2]
                                if p1x > 0
                                    if iodd==1
                                        if srtpth1[p1x, lvlchk]==nodchk && srtpth2[p2x, lvlmax-lvlchk]==nodmid
                                            querng[irgx][1] = 0
                                        end
                                    else
                                        if srtpth1[p1x, lvlchk]==nodchk && srtpth1[p2x, lvlmax-lvlchk]==nodmid
                                            querng[irgx][1] = 0
                                        end
                                        if srtpth1[p2x, lvlchk]==nodchk && srtpth1[p1x, lvlmax-lvlchk]==nodmid
                                            querng[irgx][1] = 0
                                        end
                                    end
                                end
                            end
                            @goto next_ring
                        end
                    end
                end
                rs = 2*lvlprim+iodd
                ringstat[rs] += 1
                
                # get a list of the ring indices to save
                if iodd==1
                    primes = [prime1, prime2]
                else
                    primes = [prime1]
                end
                r = collect(Iterators.flatten( # collect ring nodes for master list
                    [[nodsrc,],  
                    srtpth1[pth1, 1:Int64((rs-(2+iodd))/2)], 
                    primes, 
                    srtpth2[pth2, 1:Int64((rs-(2+iodd))/2)]]))
                push!(rings[rs], r)
            end 
        @label next_ring
        end
    return ngf, rings, ringstat
    end


    function pair_search(nodcrt::Int64,
                         nodgoal::Int64,
                         limit::Matrix{Int64},
                         goal_found::Bool,
                         lvlref::Vector{Vector{Int64}},
                         lvldist::Vector{Int64},
                         lnks::Vector{Int64},
                         nodlnkd::Vector{Vector{Int64}},
                         maxlvl::Int64)
        """
        Recursively search for shorter paths between nodes, accelerated by referential distance maps
            
            Arguments:
            - `nodcrt::Int64`: The current node being examined.
            - `nodgoal::Int64`: The target node to be found.
            - `limit`: A matrix of limits for level references and distances.
            - `goal_found::Bool`: A boolean indicating whether the goal node has been found with a shorter path
            - `lvlref`: Array of referential distance maps
            - `lvldist`: Distances of each node to current source node
            - `lnks`: Number of bonds at each node
            - `nodlnkd`: neighbour list
            - `maxlvl`: The maximum Dijkstra-level
            
            Returns:
            - `goal_found::Bool`
        """

        lr = length(lvlref)
        if nodcrt == nodgoal
            goal_found = true
            return goal_found
        end

        for lnkscrt in 1:lnks[nodcrt]
            nodprb = nodlnkd[nodcrt][lnkscrt]

            for iref in 1:lr-1
                if lvlref[iref][nodprb] >= maxlvl
                    continue
                end
                if ( lvlref[iref][nodprb] >= limit[iref, 1] || 
                    lvlref[iref][nodprb] <= limit[iref, 2] )
                    @goto next_search
                end
            end
            
            if ( lvldist[nodprb] >= limit[lr, 1] || 
                lvldist[nodprb] <= limit[lr, 2] )
                @goto next_search
            end
            
            lmtx = Matrix{Int64}(undef, lr, 2)
            for iref in 1:lr
                lmtx[iref, 1] = limit[iref, 1] - 1
                lmtx[iref, 2] = limit[iref, 2] + 1
            end

            goal_found = pair_search(nodprb, nodgoal, lmtx,
                                     goal_found, lvlref, lvldist,
                                     lnks, nodlnkd, maxlvl)
            
            if goal_found
                return goal_found
            end
            @label next_search
        end

    return goal_found
    end


    # Program itself - main function
    function ring_statistics_single(numatoms::Int64,
                                    nodlnkd::Vector{Vector{Int64}},
                                    lvlref::Vector{Vector{Int64}},
                                    lnks::Vector{Int64},
                                    maxlvl::Int64=12,
                                    mxpths::Int64=100,
                                    progress::Bool=true,
                                    nods=nothing,
                                    parallel::Bool=false,
                                    verbosity::Int64=0)
        """
        Same as below, except
        - `nods=nothing`: Indices of atoms to be considered by a thread. If not provided, defaults to all atoms. 
        """
        
        if (nods == nothing)
            nods = collect(range(1,numatoms))
        end

        # fixed array initialisations
        ringstat = zeros(Float64, 2*maxlvl) # array for the ring stats
        ngf = [0, 0] # counts rejections
        # list for ring node indices
        rings = [Vector{Vector{Int64}}() for i in range(1,2*maxlvl)]

        if verbosity>0 && progress && Threads.threadid() == 1
             println("Progress: ")
        end

        for (ct, nodsrc) in enumerate(nods)
            if progress && Threads.threadid()==1
                if nodsrc % floor(length(nods)/20) == 0
                    print(" ..$(Int64(floor(ct/length(nods) * 100)))%")
                    if nodsrc == numatoms-1
                        println()
                    end
                end
            end

            # dist
            lvldist = dijkstra_nonwgt(nodsrc, maxlvl, lnks, nodlnkd)
            
            if verbosity>0 && progress && Threads.threadid() == 1
                println("Past lvldist $nodsrc")
            end


            # find all the prime-mid-nodes
            primes = findall(x -> x<(maxlvl-maxlvl%2)/2 && x>=1, lvldist)
            if progress && Threads.threadid()==1 && verbosity>0
                println("Found primes")
            end
            # check each prime-mid-node
            for prime in primes
                

                # even ring case here
                if size(findall(
                    [lvldist[n] == lvldist[prime]-1 for n in nodlnkd[prime]]
                        ))[1] > 1

                    # find shortest paths to source from this prime-mid-node
                    pths = 0
                    strpth = zeros(Int64, mxpths, Int64(ceil(maxlvl/2)))
                    strpthx = zeros(Int64, Int64(ceil(maxlvl/2)))
                    querng = Vector{Int64}[]
                    
                    pths, strpth, strpthx = strpth_record(prime, lvldist[prime], lvldist[prime],
                                                          pths, strpth, strpthx, lvldist,
                                                          lnks, nodlnkd)
                    # check if there is a ring containing this prime-mid-node and source
                    for pth1 in 1:pths
                        for pth2 in pth1:pths
                            inter = intersect(strpth[pth1, :], strpth[pth2, :])
                            if inter in [[0], []]
                                push!(querng, [pth1, pth2])
                            end
                        end
                    end

                    # check if is prime ring
                    ngf, rings, ringstat = prime_ring(0, lvldist[prime], lnks,
                                                      lvldist, lvlref, ringstat,
                                                      nodlnkd, querng, strpth,
                                                      strpth, prime, prime,
                                                      nodsrc, ngf, rings, maxlvl)
                end

                # odd ring case here
                equal_lvls = findall(
                    [lvldist[n] == lvldist[prime] for n in nodlnkd[prime]]
                    )
                if size(equal_lvls)[1] > 0

                    prime1 = prime
                    # find shortest paths to source from this prime-mid-node
                    pths = 0
                    strpth = zeros(Int64, mxpths, Int64(ceil(maxlvl/2)))
                    strpthx = zeros(Int64, Int64(ceil(maxlvl/2)))
                    
                    pths, strpth, strpthx = strpth_record(prime1, lvldist[prime1], lvldist[prime1],
                                                          pths, strpth, strpthx,
                                                          lvldist, lnks, nodlnkd)
                    for prime2 in nodlnkd[prime1][equal_lvls]
                        if prime2 < prime1
                            continue
                        end
                        querng = Vector{Int64}[]
                        
                        pths2 = 0
                        strpth2 = zeros(Int64, mxpths, Int64(ceil(maxlvl/2)))
                        strpthx2 = zeros(Int64, Int64(ceil(maxlvl/2)))
                        pths2, strpth2, strpthx2 = strpth_record(prime2, lvldist[prime2], lvldist[prime2],
                                                                 pths2, strpth2, strpthx2,
                                                                 lvldist, lnks, nodlnkd)

                        # check both sets of shortest paths for rings
                        for pth1 in 1:pths
                            for pth2 in 1:pths2
                                inter = intersect([prime1; strpth[pth1, :]], [prime2; strpth2[pth2, :]])
                                if length(inter) == 0 || (length(inter) == 1 && inter[1] == 0)
                                    push!(querng, [pth1, pth2])
                                end
                            end
                        end 
                        
                        # check if is prime ring
                        ngf, rings, ringstat = prime_ring(1, lvldist[prime], lnks,
                                                          lvldist, lvlref, ringstat,
                                                          nodlnkd, querng, strpth,
                                                          strpth2, prime1, prime2,
                                                          nodsrc, ngf, rings, maxlvl)
                    end
                end
            end
        end
        return ringstat, ngf, rings
    end


    function process_refnodes_chunk(chunk, refnodes::Vector{Int64}, maxlvl::Int64,
                                    lnks::Vector{Int64}, nodlnkd::Vector{Vector{Int64}})
        lvlref_tmp = Vector{Int64}[]
        for i in chunk
            lvlref_i = dijkstra_nonwgt(refnodes[i], maxlvl, lnks, nodlnkd)
            push!(lvlref_tmp, lvlref_i)
        end
        return lvlref_tmp
    end

    function ring_statistics(numatoms::Int64,
                             nodlnkd::Vector{Vector{Int64}},
                             refnodes::Vector{Int64};

                             maxlvl::Int64=12,
                             mxpths::Int64=100,
                             progress::Bool=true,
                             rings_out::Bool=true,
                             outfile::String="rings_out.json",
                             verbosity::Int64=0)
        """
        Compute ring statistics for a set of atoms.
        
        # Arguments
        - `numatoms::Int`: Total number of atoms in the system.
        - `nodlnkd::Vector{Vector{Int64}}`: A list of vectors, where `nodlnkd[i]` contains the indices of atoms bonded to atom `i`.
        - `refnodes::Vector{Int64}`: Indices of reference atoms used for distance calculations.
        - `maxlvl::Int=12`: Maximum depth-level from source and reference nodes for Dijkstra
        - `mxpths::Int64=100`: Maximum number of shortest paths explored for each source - 100 is plenty for 100k of amorphous Si
        - `progress::Bool=true`: Whether to display progress information.
        - `parallel::Bool=false`: Whether to use parallel execution.
        
        # Returns
        - `ringstat::Vector{Float64}`: Array for storing the ring statistics.
        - `ngf::Vector{Int}`: Counts of rejections for even and odd rings respectively
        - `rings::Vector{Vector{Vector{Int64}}}`: List of ring node indices, separated by ring size. Should be n duplicates per ring (discovered by considering each source)
        
        This function calculates ring statistics for a single set of atoms using Dijkstra's algorithm for distance calculations. It considers even and odd ring cases and returns statistics, rejection counts, and information about detected rings.

        Implementation of algorithm described in:
        
        Yuan, X. & Cormack, A. N. 
        "Efficient Algorithm for Primitive Ring Statistics in Topological Networks." 
        Computational Materials Science, vol. 24, pp. 343-360, 2002.
        
        """

        println("Running ring_statistics with $(Threads.nthreads()) threads, verbosity level $verbosity")

        nods = collect(range(1,numatoms))
        lnks = Int64[length(nodlnkd[m]) for m in 1:numatoms] # no. bonds at each atom
        
        
        # referential distance maps using refnodes
        print("Calculating referential distance maps...")
        num_threads = min(length(refnodes), Threads.nthreads())
        refnode_chunks = Iterators.partition([i for i in 1:length(refnodes)],
                                      length(refnodes) ÷ num_threads)


        refnode_tasks = map(refnode_chunks) do chunk
            Threads.@spawn process_refnodes_chunk(chunk, refnodes, maxlvl, lnks, nodlnkd)
        end
        refnode_results = fetch.(refnode_tasks)
        lvlref = vcat(refnode_results...)
        print("done\nStarting ring statistics...\n")
    
        chunks = Iterators.partition(nods, length(nods) ÷ Threads.nthreads())
        tasks = map(chunks) do chunk
            Threads.@spawn ring_statistics_single(numatoms,
                                                  nodlnkd,
                                                  lvlref,
                                                  lnks,
                                                  maxlvl,
                                                  mxpths,
                                                  progress,
                                                  chunk,
                                                  true,
                                                  verbosity)
        end

        results = fetch.(tasks)
        rs_s = [i[1] for i in results]
        ngfs = [i[2] for i in results]
        rings_s = [i[3] for i in results]

        rings = [Vector{Vector{Int64}}() for i in range(1,2*maxlvl)]
        for r in rings_s
            for (ct, ring_list_size_ct) in enumerate(r)
                append!(rings[ct], ring_list_size_ct)
            end
        end

        rs = sum(rs_s)
        ngf = sum(ngfs)
        sf = [i for i in range(1,length(rs_s[1]))]

        for element in rs ./ sf
            if element != floor(element)
                @warn "Non-integer found in ring stats. This means not
                all rings or extraneous rings have been found. Please check
                you have enough refnodes!" 
                break
            end
        end

        println("\nNumber of rings eliminated during search (even, odd): $(ngf[1]), $(ngf[2]) ",
            "$(round(
                (ngf[1] + ngf[2])/( sum(rs)+ngf[1]+ngf[2] )*100,
                sigdigits=4)
                )% of total")

        println("Ring statistics:")
        println("Ring Size |            Count ")
        println("----------|------------------")
        
        for i in 1:length(rs)
            printfmt("{:9d} | {:15.1f}\n", sf[i], (rs ./ sf)[i])
            if rs[i] == 0 && rs[i+1] == 0 && rs[max(1, i-1)] == 0 && i > 6
                break
            end
                
        end

        json_data = JSON.json(rings)
        file = open(outfile, "w")

        try
            write(file, json_data)
        finally
            close(file)
        end

        return rs ./ sf, ngf, rings
    end

end