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
        density_table = hcat(0:0.1:1, 0.703*ones(11))
    )
    @test find_molar_fraction(system; target_volume_fraction = 0.5) ≈ 0.5 atol = 1e-5
    @test find_molar_fraction(system; target_volume_fraction = 0.0) ≈ 0.0 atol = 1e-5
    @test find_molar_fraction(system; target_volume_fraction = 1.0) ≈ 1.0 atol = 1e-5

    # system with ideal solution 
    system = SolutionSystem(
        solute_pdbfile = "$dir/data/poly_h.pdb",
        solvent_pdbfile = "$dir/data/water.pdb",
        cossolvent_pdbfile = "$dir/data/water.pdb",
        density_table = hcat([0.0 + 0.05*i for i in 0:20], [1.0 + 0.05*i for i in 0:20])
    )
    @test system.solvent_molar_mass ≈ 18.01 atol = 0.01
    @test system.cossolvent_molar_mass ≈ 18.01 atol = 0.01
    @test system.density_pure_solvent ≈ 1.0
    @test system.density_pure_cossolvent ≈ 2.0

    # Concentration conversions in this ideal system, from the molar fraction, x
    ρc = system.density_pure_cossolvent # g / mL
    ρw = system.density_pure_solvent # g / mL
    M = system.solvent_molar_mass # g / mol
    ρ(x) = ρc*x + ρw*(1-x) # g / mL
    vv(x) = (x / ρc) / (x / ρc + (1-x) / ρw) 
    v(x) = M / (1000*ρ(x)) # L/mol: volume of 1 mol (c + w) of solution
    mx(x) = x / v(x) # molarity of cossolute

    # Concentration conversions
    for x in (0.0, 0.2, 0.5, 0.7, 1.0)
        @test convert_concentration(system, x, "x" => "vv") ≈ vv(x)
        @test convert_concentration(system, x, "x" => "mol/L"; density = ρ(x)) ≈ mx(x)
        @test convert_concentration(system, vv(x), "vv" => "x") ≈ x 
        @test convert_concentration(system, vv(x), "vv" => "mol/L"; density = ρ(x)) ≈ mx(x) 
        @test convert_concentration(system, mx(x), "mol/L" => "x"; density = ρ(x)) ≈ x 
        @test convert_concentration(system, mx(x), "mol/L" => "vv"; density = ρ(x)) ≈ vv(x) 
    end

    # Molar fractions
    # voltar: precisa mesmo um processo iterativo para isso?
    @test find_molar_fraction(system; target_volume_fraction = 0.0) ≈ 0.0 atol = 1e-5
    @test find_molar_fraction(system; target_volume_fraction = 1.0) ≈ 1.0 atol = 1e-5
    @test find_molar_fraction(system; target_volume_fraction = 1.0) ≈ 1.0 atol = 1e-5
    @test find_molar_fraction(system; target_volume_fraction = 0.5) ≈ vv(0.5) atol = 1e-5

    # Ethanol-water mixture
    density_table = readdlm("$dir/data/water_ethanol.dat", comments=true, comment_char='#')
    system = SolutionSystem(
        solute_pdbfile = "$dir/data/poly_h.pdb",
        solvent_pdbfile = "$dir/data/water.pdb",
        cossolvent_pdbfile = "$dir/data/ethanol.pdb",
        density_table = density_table,
    )
    Mc = 1000 * system.density_pure_cossolvent / system.cossolvent_molar_mass # mol / L pure ethanol
    @test convert_concentration(system, 1.0, "x" => "mol/L"; density = system.density_pure_cossolvent) ≈ Mc

end
