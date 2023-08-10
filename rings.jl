using JuLIP
using ASE

# Functions

function dijkstra_nonwgt(nodsrc, lvlreq, lnks, nodlnkd, verbose=false)
    
    # initialise
    lvldist = lvlreq + 100 .+ zeros(Int64, length(lnks))
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


function strpth_record(nodcrt, lvlcrt, lvlprim, pths, strpth, strpthx, lvldist)
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
            pths, strpth, strpthx = strpth_record(nodprb, lvlprb, lvlprim, pths, strpth, strpthx, lvldist)
        end
    end
return pths, strpth, strpthx
end


function prime_ring(iodd, lvlprim, lnks, lvldist, lvlref, ringstat, nodlknkd, querng, srtpth1, srtpth2, prime1, prime2, nodsrc)
    global ngf, rings

    rngs = size(querng)[1]

    for irng in 1:rngs
        pth1 = querng[irng][1]
        pth2 = querng[irng][2]
        
        if pth1 > 0
            chkcount = 0
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

                    limit[1:5, 1] = [lvlref[i][nodmid]+lvlprim-1 for i in 1:5]
                    limit[1:5, 2] = [lvlref[i][nodmid]-lvlprim+1 for i in 1:5]
                    limit[5, 1] = lvldist[nodmid] + lvlprim - 1
                    limit[5, 2] = lvldist[nodmid] - lvlprim + 1

                    goal_found = false
                    goal_found = pair_search(nodchk, nodmid, limit, goal_found, lvlref, lvldist)

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
            r = collect(Iterators.flatten(
                [[nodsrc,], 
                srtpth1[pth1, 1:Int((rs-(2+iodd))/2)], 
                [prime1,], [prime2,], 
                srtpth2[pth2, 1:Int((rs-(2+iodd))/2)]]))
            push!(rings[rs], r)
        end 
    @label next_ring
    end
return 
end


function pair_search(nodcrt, nodgoal, limit, goal_found, lvlref, lvldist)
    # Input: goal_found, nodcrt, nodgoal, lnks(), lvldist(), nodlnkd(), lvlref(), limit()
    # Output: goal_found
    # Internal: lmtx()
    global lnks, nodlnkd

    if nodcrt == nodgoal
        goal_found = true
        return goal_found
    end

    for lnkscrt in 1:lnks[nodcrt]
        nodprb = nodlnkd[nodcrt][lnkscrt]

        for iref in 1:4
            if lvlref[iref][nodprb] >= limit[iref, 1] || lvlref[iref][nodprb] <= limit[iref, 2]
                @goto next_search
            end
        end
        
        if lvldist[nodprb] >= limit[5, 1] || lvldist[nodprb] <= limit[5, 2]
            @goto next_search
        end
        
        lmtx = zeros(Int, 5, 2)
        for iref in 1:5
            lmtx[iref, 1] = limit[iref, 1] - 1
            lmtx[iref, 2] = limit[iref, 2] + 1
        end

        goal_found = pair_search(nodprb, nodgoal, lmtx, goal_found, lvlref, lvldist)
        
        if goal_found
            return goal_found
        end
        @label next_search
    end
goal_found = false
return goal_found
end


# initialise - replace this for proper input later

at = Atoms(read_xyz("""aSi_500atom_11Ks_annealed.xyz"""))

pos = at.X
c = at.cell
inds = [1:size(at.X)[1];]

# calculate neighbourlist - should offload this to Ovito or ASE to reduce Julia requirements
nl = neighbourlist(at, 2.85)


# Actual Julia calculations start from here
nodlnkd = [nl.j[nl.i .== m] for m in 1:length(pos)]
lnks = [length(nodlnkd[m]) for m in 1:length(pos)]

# global parameters and arrays

# changeable
maxlvl = L = 12 # should be an adjustable argument
mxnodes = Int(1e4)
mxpths = Int(1e5)

# fixed array initialisations
limit = zeros(Int, 5, 2)
ngf = [0, 0]
ringstat = zeros(Float64, 2*maxlvl) # array for the ring stats
rings = [Vector{Vector{Int64}}() for i in range(1,20)] # array for the ring indices


# Program itself
for nodsrc in 1:length(pos)

    lvlref = [dijkstra_nonwgt(i, 30, lnks, nodlnkd) for i in [50, 150, 250, 400, 1]] # need to set these sensibly automatically
    lvldist = dijkstra_nonwgt(nodsrc, maxlvl, lnks, nodlnkd)

    # find all the prime-mid-nodes
    primes = findall(x -> x<(L-maxlvl%2)/2 && x>1, lvldist)

    # check each prime-mid-node
    for prime in primes

        # even ring case here
        if size(findall([lvldist[n] == lvldist[prime]-1 for n in nodlnkd[prime]]))[1] > 1

            # find shortest paths to source from this prime-mid-node
            pths = 0
            strpth = zeros(Int, mxpths, Int(maxlvl/2))
            strpthx = zeros(Int, Int(maxlvl/2))
            querng = []
            
            pths, strpth, strpthx = strpth_record(prime, lvldist[prime], lvldist[prime], pths, strpth, strpthx, lvldist)
            
            # check if there is a ring containing this prime-mid-node and source
            for pth1 in 1:pths
                for pth2 in pth1:pths
                    inter = intersect(strpth[pth1, :], strpth[pth2, :])
                    if intersect(strpth[pth1, :], strpth[pth2, :]) in [[0], []]
                        push!(querng, [pth1, pth2])
                    end
                end
            end

            # check if is prime ring
            prime_ring(0, lvldist[prime], lnks, lvldist, lvlref, ringstat, nodlnkd, querng, strpth, strpth, prime, prime, nodsrc)

        end

        # odd ring case here
        equal_lvls = findall([lvldist[n] == lvldist[prime] for n in nodlnkd[prime]])
        if size(equal_lvls)[1] > 0

            prime1 = prime
            # find shortest paths to source from this prime-mid-node
            pths = 0
            strpth = zeros(Int, mxpths, Int(maxlvl/2))
            strpthx = zeros(Int, Int(maxlvl/2))
            
            pths, strpth, strpthx = strpth_record(prime1, lvldist[prime1], lvldist[prime1], pths, strpth, strpthx, lvldist)
            for prime2 in nodlnkd[prime1][equal_lvls]
                if prime2 < prime1
                    continue
                end
                querng = []
                
                pths2 = 0
                strpth2 = zeros(Int, mxpths, Int(maxlvl/2))
                strpthx2 = zeros(Int, Int(maxlvl/2))
                pths2, strpth2, strpthx2 = strpth_record(prime2, lvldist[prime2], lvldist[prime2], pths2, strpth2, strpthx2, lvldist)

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
                prime_ring(1, lvldist[prime], lnks, lvldist, lvlref, ringstat, nodlnkd, querng, strpth, strpth2, prime1, prime2, nodsrc)
            end

        end
    end
end


println("Number of rings elimated (even, odd): $(ngf[1]), $(ngf[2])")

println("Ring statistics: \n")
sf = [i for i in range(1,length(ringstat))]
println(Int(ringstat) ./ sf)


