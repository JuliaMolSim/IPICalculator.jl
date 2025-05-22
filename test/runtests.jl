using Test

using AtomsBase
using AtomsCalculators
using AtomsCalculators.Testing
using AtomsCalculatorsUtilities.IPI
using AtomsCalculatorsUtilities.PairPotentials
using Unitful
using Base.Threads

@testset "IPI Calculator Tests" begin

    # Build a test system
    hydrogen = isolated_system([
        :H => [0.1, 0, 0.]u"Å",
        :H => [0, 0, 1.]u"Å",
        :H => [4., 0, 0.]u"Å",
        :H => [4., 1., 0.]u"Å"
    ])
    box = (
        [10.0, 0., 0.]u"Å",
        [0.0, 10., 0.]u"Å",
        [0.0, 0., 10.]u"Å",
    )
    pbc = (false, false, false)
    hydrogen = FlexibleSystem(hydrogen[:], cell_vectors=box, periodicity=pbc)

    V = SimplePairPotential(
        x-> (x-0.9)^2-1,
        1,
        1,
        2.0u"Å"
    )

    # We need to set up the calculator with asyncronous communication
    ipi_future = @spawn IPIserver(port=33415)
    sleep(1) # we need to yield to start the server

    ipi_driver = @spawn run_driver("127.0.0.1", V, hydrogen; port=33415)
    sleep(1) # we need to yield to connect to the server

    calc = fetch(ipi_future)

    # The tests itself
    test_energy_forces_virial(hydrogen, calc)

    # Check that the calculator returns the same energy, forces and virial as the potential
    @test AtomsCalculators.potential_energy(hydrogen, V) ≈ AtomsCalculators.potential_energy(hydrogen, calc)
    f_v   = AtomsCalculators.forces(hydrogen, V)
    f_ipi = AtomsCalculators.forces(hydrogen, calc)
    @test all( isapprox.(f_v, f_ipi) )
    @test AtomsCalculators.virial(hydrogen, V) ≈ AtomsCalculators.virial(hydrogen, calc)

    close(calc)
end
