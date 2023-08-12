
module rings

    using PyCall, JSON
    export ring_statistics

    function read_nodlnkd(filename)
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

    function dijkstra_nonwgt(nodsrc::Int64, lvlreq::Int64, lnks, nodlnkd)
        
        # initialise
        lvldist = lvlreq + 2 .+ zeros(Int64, length(lnks))
        lvldist[nodsrc] = 0
        quebgn = 0
        quend = 1
        queue = zeros(Int64, length(lnks))
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


    function strpth_record(nodcrt::Int64, lvlcrt::Int64, lvlprim::Int64, 
                           pths::Int64, strpth, strpthx,
                           lvldist, lnks, nodlnkd)

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


    function prime_ring(iodd::Int64, lvlprim::Int64, lnks, 
                        lvldist, lvlref, ringstat, 
                        nodlnkd, querng, 
                        srtpth1, srtpth2, 
                        prime1::Int64, prime2::Int64, 
                        nodsrc::Int64, ngf, rings, maxlvl::Int64)
        
        # cond = (nodsrc == 28509 && prime1 == 39689)
        cond = false
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
                        limit = zeros(Int64, length(lvlref), 2)
                        lr = length(lvlref)

                        limit[1:lr, 1] = [lvlref[i][nodmid] > maxlvl ? Int64(1e4) : lvlref[i][nodmid]+lvlprim-1 for i in 1:lr]
                        if cond
                            @show lvlprim
                        end
                        # limit[1:lr, 1] = [lvlref[i][nodmid]+lvlprim-1 for i in 1:lr]
                        limit[1:lr, 2] = [lvlref[i][nodmid] > maxlvl ? Int64(-1e4) : lvlref[i][nodmid]-lvlprim+1 for i in 1:lr]
                        # limit[1:lr, 2] = [lvlref[i][nodmid]-lvlprim+1 for i in 1:lr]

                        limit[lr, 1] = lvldist[nodmid] + lvlprim - 1
                        limit[lr, 2] = lvldist[nodmid] - lvlprim + 1
                        if cond
                            @show limit
                        end
                        goal_found = false
                        goal_found = pair_search(nodchk, nodmid, limit,
                                                 goal_found, lvlref, lvldist,
                                                 lnks, nodlnkd, maxlvl, cond)

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
                r = collect(Iterators.flatten(
                    [[nodsrc,],  
                    srtpth1[pth1, 1:Int64((rs-(2+iodd))/2)], 
                    primes, 
                    srtpth2[pth2, 1:Int64((rs-(2+iodd))/2)]]))
                push!(rings[rs], r)
                if cond
                    println("ring added")
                    @show r
                end
            end 
        @label next_ring
        end
    return ngf, rings, ringstat
    end


    function pair_search(nodcrt::Int64, nodgoal::Int64, limit,
                         goal_found::Bool, lvlref, lvldist,
                         lnks, nodlnkd, maxlvl, cond=false)
        # Input: goal_found, nodcrt, nodgoal, lnks(), lvldist(), nodlnkd(), lvlref(), limit()
        # Output: goal_found
        # Internal: lmtx()
        lr = length(lvlref)
        if nodcrt == nodgoal
            goal_found = true
            return goal_found
        end

        for lnkscrt in 1:lnks[nodcrt]
            nodprb = nodlnkd[nodcrt][lnkscrt]

            for iref in 1:lr-1
                if lvlref[iref][nodprb] >= maxlvl
                    # if cond
                    #     println("Refnode too far away ", lvlref[iref][nodprb], maxlvl)
                    #  end
                    continue
                end
                if ( lvlref[iref][nodprb] >= limit[iref, 1] || 
                    lvlref[iref][nodprb] <= limit[iref, 2] )
                    if cond
                        println("1st")
                        @show nodcrt, nodgoal, nodprb, lvldist[nodprb], lvldist[nodcrt], lvlref[iref][nodprb], limit[iref, :]
                    end
                    @goto next_search
                end
            end
            
            if ( lvldist[nodprb] >= limit[lr, 1] || 
                lvldist[nodprb] <= limit[lr, 2] )
                if cond
                    println("2nd")
                    @show nodcrt, nodgoal, nodprb, lvldist[nodprb], lvldist[nodcrt]
                end
                @goto next_search
            end
            
            lmtx = zeros(Int64, lr, 2)
            for iref in 1:length(lvlref)
                lmtx[iref, 1] = limit[iref, 1] - 1
                lmtx[iref, 2] = limit[iref, 2] + 1
            end

            if cond
                println("deeper")
            end
            goal_found = pair_search(nodprb, nodgoal, lmtx,
                                     goal_found, lvlref, lvldist,
                                     lnks, nodlnkd, maxlvl, cond)
            
            if goal_found
                if cond
                    println("Goal found!")
                end
                return goal_found
            end
            @label next_search
        end
    goal_found = false
    return goal_found
    end


    # Program itself - main function
    function ring_statistics_single(numatoms, nodlnkd, refnodes,
                                    maxlvl=12, mxpths=Int64(100),
                                    progress=true, nods=nothing,
                                    parallel=false)

        if (nods == nothing)
            nods = collect(range(1,numatoms))
        end

        # fixed array initialisations
        lnks = Int64[length(nodlnkd[m]) for m in 1:numatoms]
        ringstat = zeros(Float64, 2*maxlvl) # array for the ring stats
        L = maxlvl
        ngf = [0, 0]
        rings = [Vector{Vector{Int64}}() for i in range(1,2*maxlvl)]


        for (ct, nodsrc) in enumerate(nods)
            if progress && Threads.threadid()==1
                if nodsrc % floor(length(nods)/20) == 0
                    print(" ..$(Int64(floor(ct/length(nods) * 100)))%")
                    if nodsrc == numatoms-1
                        println()
                    end
                end
            end
            cond = false 
            # cond = !(nodsrc == 28509)
            if cond
                continue
            end

            lvlref = [dijkstra_nonwgt(i, maxlvl, lnks, nodlnkd) for i in refnodes]
            lvldist = dijkstra_nonwgt(nodsrc, maxlvl, lnks, nodlnkd)

            # find all the prime-mid-nodes
            primes = findall(x -> x<(L-maxlvl%2)/2 && x>=1, lvldist)
            # check each prime-mid-node
            for prime in primes
                # cond = (nodsrc == 28509 && prime in [39689])
                cond = false
                

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
                                if cond
                                    println("pushing ", [prime, strpth[pth1, :], strpth[pth2, :], nodsrc])
                                end
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
                                p1 = collect(Iterators.flatten([[prime1], strpth[pth1, :]]))
                                p2 = collect(Iterators.flatten([[prime2], strpth2[pth2, :]]))
                                inter = intersect(p1, p2)
                                if inter in [[0], []]
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

        if !parallel
            println("Number of rings elimated (even, odd): $(ngf[1]), $(ngf[2])")
            println("Ring statistics: \n")
            sf = [i for i in range(1,length(ringstat))]
            println(ringstat ./ sf)
            return ringstat ./ sf
        else
            return ringstat, ngf, rings
        end
    end


    function ring_statistics(numatoms, nodlnkd, refnodes,
                             maxlvl=12, mxpths=Int64(100), progress=true,
                             rings_out=true, outfile="rings_out.json")

        println("Running ring_statistics in parallel with $(Threads.nthreads()) threads")
        nods = collect(range(1,numatoms))
        chunks = Iterators.partition(nods, length(nods) ÷ Threads.nthreads())
        
        tasks = map(chunks) do chunk
            Threads.@spawn ring_statistics_single(numatoms, nodlnkd, refnodes,
                                                  maxlvl, mxpths, progress,
                                                  chunk, true)
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

        for element in rs
            if isa(element, Float64) && element != floor(element)
                @warn "Non-integer found in ring stats. This means not
                all rings or extraneous rings have been found. Please check
                you have enough refnodes!" 
                break
            end
        end

        println("Number of rings elimated during search (even, odd): $(ngf[1]), $(ngf[2])")
        println("$((ngf[1] + ngf[2])/(sum(rs) + ngf[1] + ngf[2]) * 100)% of total")

        println("Ring statistics: ")
        println(rs ./ sf)

        # Convert the data to JSON format
        json_data = JSON.json(rings)

        # Open the file for writing
        file = open(outfile, "w")

        try
            # Write JSON data to the file
            write(file, json_data)
        finally
            close(file)
        end

        return rs, ngf, rings
    end

end