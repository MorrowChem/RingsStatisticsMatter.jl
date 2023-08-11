using DelimitedFiles
include("rings.jl")


file = open("nodlnkd.txt", "r")
nodlnkd = Vector{Vector{Int}}()

# Read each line from the file
for line in eachline(file)
    # Split the line into individual values
    values = split(line, ',')
    
    # Convert values to integers
    int_values = [parse(Int, v) for v in values]
    
    # Add the array of integers to the vector
    push!(nodlnkd, int_values)
end

# Close the file
close(file)


@show nodlnkd[1:3][:]

rs = rings.ring_statistics(4096, nodlnkd, [3920, 206, 1958, 2368, 1206])