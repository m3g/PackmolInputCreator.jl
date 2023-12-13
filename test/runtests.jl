using PackmolInputCreator
using TestItemRunner

@run_package_tests

@testitem "Consistency tests" begin
    using PDBTools
    using PackmolInputCreator
    using DelimitedFiles
    dir = PackmolInputCreator.project_dir*"/../test"

    # system with one solvent only, with constant density
    system = SolutionSystem(
        solute_pdbfile = "$dir/data/poly_h.pdb",
        solvent_pdbfile = "$dir/data/octane.pdb",
        cossolvent_pdbfile = "$dir/data/octane.pdb",
        density_table = fill(0.703, 10, 2)
    )

    @test find_molar_fraction(system; target_volume_fraction = 0.5) ≈ 0.5
    @test find_molar_fraction(system; target_volume_fraction = 0.0) ≈ 0.0 atol = 1e-5
    @test find_molar_fraction(system; target_volume_fraction = 1.0) ≈ 1.0 

    # True densities
    density_table = readdlm("$dir/data/density_table.dat", comments=true, comment_char='#')

end
