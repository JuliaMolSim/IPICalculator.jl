# IPICalculator.jl

Implements a server-client based calculator defined by [i-PI](https://github.com/i-pi/i-pi) .

Julia side server works as a [AtomsCalculators](https://github.com/JuliaMolSim/AtomsCalculators.jl) calculator and
allows using any program that implements [i-PI](https://github.com/i-pi/i-pi)
driver (client) or [Molssi Driver interface](https://molssi.org/software/mdi-2/) (MDI) engine
as an engine for calculations in Julia.

Julia side driver (client) allows using any AtomsCalculators combatible calculator in any program that implements either [i-PI](https://github.com/i-pi/i-pi) or [Molssi Driver interface](https://molssi.org/software/mdi-2/) server (note MDI calls server as driver).

IPIcalculator.jl implements only i-PI protocol not full MDI standard.

## Usage

Julia server can be started with

```julia
using IPICalculator

# network mode using port 33415
calc = SocketServer(port=33415)

# unixsocket mode with socket at /tmp/ipi_mysocket
calc = SocketServer( unixsocket="mysocket",  basename="/tmp/ipi_")

# perform calculations
# AtomsCalculators.energy_force_virial(sys, calc)
# AtomsCalculators.potential_energy(sys, calc)
# etc.   

# end by calling close
close(calc)
```

`SocketServer` implements [AtomsCalculators](https://github.com/JuliaMolSim/AtomsCalculators.jl) interface.


Julia driver can be started with

```julia
using IPICalculator

sys = # generate a system that sets up atom types for the calculator
calc = # generate a AtomsCalculators calculator

# network mode on localhost 
run_driver(sys, calc; port=33415)

# network mode on remote host
run_driver(sys, calc; address="some.ip.address", port=12345)

# unixsocket mode with socket at /tmp/ipi_mysocket
run_driver(sys, calc; unixsocket="mysocket",  basename="/tmp/ipi_")
```

Note, that you need to give driver [AtomsBase](https://github.com/JuliaMolSim/AtomsBase.jl) structure that defines atom types.
If you want to change either number or types of atoms, you need to create new driver.

i-PI inteface includes virial. When virial is not needed it can
be turned off, to allow calculators that don't support virial to work.

```julia
# Don't calculate virial and return zeros for virial instead
# using localhost and port 33415
run_driver(sys, calc; ignore_virial=true)
```

## Enable debug messages

You can enable communication logs by allowing debug level logs.
This can be done by calling (in Julia)

```julia
ENV["JULIA_DEBUG"] = IPICalculator 
```

## Using Julia based calculator from ASE

Fist start ASE [SocketIOCalculator](https://wiki.fysik.dtu.dk/ase/ase/calculators/socketio/socketio.html#ase.calculators.socketio.SocketIOCalculator) from Python

```python
from ase.calculators.socketio import SocketIOCalculator

# Using localhost and default port 33415
# Adjust as you see fit
calc = SocketIOCalculator()
```

Connect Julia driver to the socket

```julia
run_driver(sys, calc)
```

Do calculations in Python.

### Calling ASE calculators from Julia

You can use ASE calculators in Julia with using [SocketClient](https://wiki.fysik.dtu.dk/ase/ase/calculators/socketio/socketio.html#ase.calculators.socketio.SocketClient) from ASE and `SocketServer` from Julia.

Other options to call ASE from Julia include:
- [ASEconvert.jl](https://github.com/mfherbst/ASEconvert.jl) - AtomsCalculators combatible
- [Molly](https://github.com/JuliaMolSim/Molly.jl) has also its own way of calling [ASE calculators](https://juliamolsim.github.io/Molly.jl/stable/api/#Molly.ASECalculator)

The main difference is that using sockets Python and Julia are different processes, while ACEconvert and Molly start Python from Julia, so that Python is in the same prosess than Julia, which can cause problems sometimes.
