module IOPiper
"""
Facilitates redirect of Julia stdout/stderr to IPython notebook output.
"""

const orig_stdin  = Ref{IO}()
const orig_stdout = Ref{IO}()
const orig_stderr = Ref{IO}()

function __init__()
    orig_stdin[]  = stdin
    orig_stdout[] = stdout
    orig_stderr[] = stderr
end

"""
If `d` ends with an incomplete UTF8-encoded character, return the number of trailing incomplete bytes.
Otherwise, return `0`.
Taken from IJulia.jl.
"""
function num_utf8_trailing(d::Vector{UInt8})
    i = length(d)
    # find last non-continuation byte in d:
    while i >= 1 && ((d[i] & 0xc0) == 0x80)
        i -= 1
    end
    i < 1 && return 0
    c = d[i]
    # compute number of expected UTF-8 bytes starting at i:
    n = c <= 0x7f ? 1 : c < 0xe0 ? 2 : c < 0xf0 ? 3 : 4
    nend = length(d) + 1 - i # num bytes from i to end
    return nend == n ? 0 : nend
end

function pipe_stream(sender::IO, receiver, buf::IO = IOBuffer())
    try
        while !eof(sender)
            nb = bytesavailable(sender)
            write(buf, read(sender, nb))

            # Taken from IJulia.send_stream:
            d = take!(buf)
            n = num_utf8_trailing(d)
            dextra = d[end-(n-1):end]
            resize!(d, length(d) - n)
            s = String(copy(d))

            write(buf, dextra)
            receiver(s)  # check isvalid(String, s)?
        end
    catch e
        if !isa(e, InterruptException)
            rethrow()
        end
        pipe_stream(sender, receiver, buf)
    end
end

const read_stdout = Ref{Base.PipeEndpoint}()
const read_stderr = Ref{Base.PipeEndpoint}()

function pipe_std_outputs(out_receiver, err_receiver)
    global readout_task
    global readerr_task
    read_stdout[], = redirect_stdout()
    readout_task = @async pipe_stream(read_stdout[], out_receiver)
    read_stderr[], = redirect_stderr()
    readerr_task = @async pipe_stream(read_stderr[], err_receiver)
end

function io_test()
    """Run this to see if redirect to notebook is working.
    Should see output Testing redircect if so"""
    println("Testing redirect")
end

end  # module