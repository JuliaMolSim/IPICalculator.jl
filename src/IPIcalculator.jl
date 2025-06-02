module IPIcalculator

# This is AtomsCalculators implementation for i-PI protocol.
# For reference see https://github.com/i-PI/i-PI/blob/master/drivers/py/driver.py

using AtomsBase
using AtomsCalculators
using Sockets
using StaticArrays
using Unitful
using UnitfulAtomic

export run_driver
export SocketServer

const hdrlen = 12  # number of characters in header messages

# Define context units for better conversions
const bohr    = Unitful.ContextUnits(u"bohr", u"Å")
const hartree = Unitful.ContextUnits(u"hartree", u"eV")

const energy_type   = typeof( zero(Float64) * hartree )
const pos_type      = typeof( zero( SVector{3, Float64} ) * bohr ) 
const force_el_type = typeof( zero( SVector{3, Float64} ) * (hartree/bohr) )
const virial_type   = typeof( zero( SMatrix{3, 3, Float64} ) * hartree )

function sendmsg(comm, message; nbytes=hdrlen)
    @debug "Sending message" message
    if length(message) == nbytes
        final_message = message
    elseif length(message) < nbytes
        l = nbytes - length(message)
        final_message = message * repeat(' ', l)
    else
        error("Message is too long")
    end
    write(comm, final_message)
end


function recvmsg(comm, nbytes=hdrlen)
    raw_message = read(comm, nbytes)
    if length(raw_message) == 0
        @info "Server was probably closed and did not send EXIT"
        return "EXIT"
    end
    @assert length(raw_message) == nbytes "recieved message did not have correct lenght"
    message = Char.(raw_message) |> String |> strip
    @debug "Recieved message" message
    return message
end

function recvinit(comm)
    @debug "Recieving INIT"
    bead_index = read(comm, Int32)
    nbytes = read(comm, Int32)
    raw_data = read(comm, nbytes)
    message = Char.(raw_data) |> String |> strip
    return (;
        :bead_index => bead_index,
        :message => message
    )
end

function send_init(comm)
    @debug "Sending INIT"
    write(comm, one(Int32))
    str = "ok"
    write(comm, sizeof(str))
    write(comm, str)
    return true
end


function recvposdata(comm)
    raw_cell = read(comm, 9*sizeof(Float64))
    raw_icell = read(comm, 9*sizeof(Float64)) # drop this (inverce cell)
    natoms = read(comm, Int32)
    raw_pos = read(comm, sizeof(Float64)*3*natoms)
    data_cell = reinterpret(pos_type, raw_cell)
    data_pos = reinterpret(pos_type, raw_pos)
    @debug "Position data recieved"
    return (;
        :cell => Tuple(Vector(data_cell)),
        :positions => Vector(data_pos)  # clean type a little bit
    )
end

function send_pos_data(comm, sys)
    @debug "Sending position data"
    box_tmp = reduce(vcat, cell_vectors(sys))
    box = (Float64 ∘ ustrip).(u"bohr", box_tmp)
    write(comm, box)
    write(comm, zeros(3,3) )  #inverse cell that is to be dropped
    l = length(sys)
    write(comm, Int32(l) )
    pos = map( 1:l ) do i 
        SVector{3, Float64}(ustrip.(u"bohr", position(sys, i)))
    end
    write(comm, pos)
    @debug "Position data sent"
    return true
end


function sendforce(comm, e::Number, forces::AbstractVector, virial::AbstractMatrix)
    etype = (eltype ∘ eltype)(forces)
    f_tmp = reinterpret(reshape, etype, forces)
    sendmsg(comm, "FORCEREADY")
    write(comm, (Float64 ∘ ustrip)(u"hartree", e) )
    write(comm, Int32( length(forces) ) )
    write(comm, (Float64 ∘ ustrip).(u"hartree/bohr", f_tmp) )
    write(comm, (Float64 ∘ ustrip).(u"hartree", virial) )

    # Send single byte at end to make sure we are alive
    write(comm, one(Int32) )
    write(comm, zero(UInt8) )
end


function recv_force(comm)
    sendmsg(comm, "GETFORCE")
    mess = recvmsg(comm)
    if mess == "FORCEREADY"
        @debug "Recieving forces"
        e = read(comm, Float64)
        n = read(comm, Int32)
        f_raw = read(comm, sizeof(Float64)*3*n)
        v_raw = read(comm, sizeof(Float64)*9)

        # Reading end message that is dropped
        i = read(comm, Int32)
        _ = read(comm, i)

        f = reinterpret(force_el_type, f_raw) |> Vector
        v = reinterpret(virial_type, v_raw)[1]
        return (
            energy = e * hartree,
            forces = f,
            virial = v
        )
    else
        error("Expected \"FORCEREADY\", but received \"$mess\"")
    end
end

"""
    run_driver(init_structure, calculator; address=ip"127.0.0.1", port=31415, unixsocket=nothing, basename="/tmp/ipi_" )

Connect i-PI driver to server at given `address` using `port` (default localhost and 31415).
If `unixsocket` is something, it will be used to connect to the server using Unix socket
at `basename*unixsocket`. If `unixsocket` is not given, the server will use network connection.

You need to give initial structure as i-PI protocol does not transfer atom symbols.
This means that, if you want to change the number of atoms or their symbols, you need
to lauch a new driver.

Calculator has major events logged at info level and all communication logged at debug level. 

You can use `ignore_virial=true` to ignore virial calculation. This is useful if your calculator does not
support virial calculation or when virial is not needed. In this case, the calculator will send zero virial to the server.

# Example
```julia
using IPIcalculator
using AtomsBase, AtomsCalculators

# Default for localhost and port 31415
run_driver(init_structure, calc)

# Custom port and address
run_driver(init_structure, calc; address="some.ip.address", port=12345)

# Custom Unix socket at /tmp/ipi_my_socket
run_driver(init_structure, calc; unixsocket="my_socket", basename="/tmp/ipi_")

# Calculation without virial (localhost and port 31415)
run_driver(init_structure, calc; ignore_virial=true)
```
"""
function run_driver(
        init_structure,
        calc;
        address=ip"127.0.0.1",
        port=31415,
        unixsocket=nothing,
        basename="/tmp/ipi_",
        ignore_virial=false
    )
    if isnothing( unixsocket )
        comm = connect(address, port)
    else 
        comm = connect(basename * unixsocket)
    end
    has_init = true  # we have init structure as an input
    has_data = false
    data = nothing

    pbc = periodicity(init_structure)
    masses = mass(init_structure, :)
    atom_species = species(init_structure, :)
    positions = position(init_structure, :)
    box = cell_vectors(init_structure)


    while true
        header = recvmsg(comm)

        if header == "STATUS"
            if !has_init
                sendmsg(comm, "NEEDINIT")
            elseif has_data
                sendmsg(comm, "HAVEDATA")
            else
                sendmsg(comm, "READY")
            end
        elseif header == "INIT"
            init = recvinit(comm)
            # we don't init anything for now
            has_init = true
            has_data = false
        elseif header == "POSDATA"
            pos = recvposdata(comm)
            positions = pos[:positions]
            box = pos[:cell]
            @assert length(atom_species) == length(positions) "received amount of position data does no match the atomic symbol data"
            system = FastSystem(box, pbc, positions, atom_species, masses)
            if ignore_virial
                tmp = AtomsCalculators.energy_forces(system, calc)
                data = (tmp..., virial=zero(virial_type))
            else
                data = AtomsCalculators.energy_forces_virial(system, calc)
            end
            has_data = true
        elseif header == "GETFORCE"
            sendforce(comm, data[:energy], data[:forces], data[:virial])
            has_data = false
        elseif header == "EXIT"
            @info "Server send EXIT signal - closing connection" 
            close(comm)
            break
        else
            close(comm)
            error("Message not recognised")
        end
        
    end
    return nothing
end



## Server specific part

"""
    SocketServer( port=31415, unixsocket=nothing, basename="/tmp/ipi_" )

Creates i-PI https://ipi-code.org/ server that works as an AtomsCalculators compatible calculator
once i-PI driver has been connected.

Calculator has major events logged at info level and all communication logged at debug level.

At the end of calculation you should call `close` on the calculator to send exit signal to the driver.

# Kwargs
- `basename="/tmp/ipi_"`    -  prefixed to unixsocket if given, otherwise ignored. This is used to create a Unix socket for communication.
- `port=31415`              -  network port the server is using
- `unixsocket=nothing`      -  if given will use `basename*unixsocket` as Unix socket for commucation.

# Example

```julia
using IPIcalculator

# Default port 31415
ipi_calculator = SocketServer()

# Custom port
ipi_calculator = SocketServer( port=12345 )

# Custom Unix socket at /tmp/ipi_mysocket
ipi_calculator = SocketServer( unixsocket="mysocket", basename="/tmp/ipi_" )

# At the end of calculation you should close the calculator
close(ipi_calculator)
```
"""
mutable struct SocketServer{TS, TC}
    server::TS
    sock::TC
    function SocketServer(; port=31415, unixsocket=nothing, basename="/tmp/ipi_" )
        server, sock = start_ipi_server( port=port, unixsocket=unixsocket, basename=basename)
        @info "Connection received"
        new{typeof(server), typeof(sock)}(server, sock)
    end
end

function start_ipi_server(; port=31415, unixsocket=nothing, basename="/tmp/ipi_", tries=5 )
    @info "Starting i-PI socket server"
    server = nothing
    if isnothing(unixsocket)
        server = listen(port)
    else
        server = listen( basename * unixsocket )
    end
    get_connection(server; tries=tries) # returns server, socket
end

function get_connection(server; tries=5)
    sock = accept(server)
    i = 1
    while isopen(sock) || i < tries
        sendmsg(sock, "STATUS")
        mess = recvmsg(sock)
        if mess == "NEEDINIT"
            sendmsg(sock, "INIT")
            send_init(sock)
            continue
        elseif mess == "READY"
            return server, sock
        else
            i += 1
            close(sock)
            sock = accept(server)
        end
    end
    error("Could not form a connection to a working i-PI driver")
end

function Base.close(ipi::SocketServer)
    if isopen(ipi.sock)
        sendmsg(ipi.sock, "EXIT")
    end
end


function AtomsCalculators.energy_forces_virial(sys, ipi::SocketServer; kwargs...)
    if ! isopen(ipi.sock)
        @debug "reconnecting to i-PI driver"
        _, sock = get_connection(ipi.server)
        ipi.sock = sock
    end
    sendmsg(ipi.sock, "POSDATA")
    send_pos_data(ipi.sock, sys)
    sendmsg(ipi.sock, "STATUS")
    mess = recvmsg(ipi.sock)
    if mess == "HAVEDATA"
        return recv_force(ipi.sock)
    else
        error("Expected \"HAVEDATA\", but received \"$mess\"")
    end
end


AtomsCalculators.@generate_interface function AtomsCalculators.potential_energy(sys, ipi::SocketServer; kwargs...)
    tmp = AtomsCalculators.energy_forces_virial(sys, ipi)
    return tmp.energy
end

AtomsCalculators.@generate_interface function AtomsCalculators.forces(sys, ipi::SocketServer; kwargs...)
    tmp = AtomsCalculators.energy_forces_virial(sys, ipi)
    return tmp.forces
end

AtomsCalculators.@generate_interface function AtomsCalculators.virial(sys, ipi::SocketServer; kwargs...)
    tmp = AtomsCalculators.energy_forces_virial(sys, ipi)
    return tmp.virial
end

function AtomsCalculators.energy_forces(sys, ipi::SocketServer; kwargs...)
    tmp = AtomsCalculators.energy_forces_virial(sys, ipi)
    return tmp
end



AtomsCalculators.energy_unit(::SocketServer) = hartree
AtomsCalculators.length_unit(::SocketServer) = bohr

AtomsCalculators.promote_force_type(::AtomsBase.AbstractSystem, ::SocketServer) = force_el_type

AtomsCalculators.zero_energy(sys, ::SocketServer) = zero(energy_type)
AtomsCalculators.zero_forces(sys, ::SocketServer) = zeros(force_el_type, length(sys))
AtomsCalculators.zero_virial(sys, ::SocketServer) = zeros(virial_type)

end # module IPIcalculator
