module PackmolInputCreator

using PDBTools

export find_x
export convert_c
export write_input
export interpolate

# conversion factor from mL/mol to Å^3/molecule
const CMV = 1e24 / 6.02e23

# conversion factor from mol/L to molecules/Å^3
const CMC = 6.02e23 / 1e27

"""

Interpolate a value from a (x,y) data matrix

"""
function interpolate(x, ρ)
    i = findfirst(d -> d > x, ρ[:, 1])
    dρdx = (ρ[i, 2] - ρ[i-1, 2]) / (ρ[i, 1] - ρ[i-1, 1])
    d = ρ[i-1, 2] + dρdx * (x - ρ[i-1, 1])
    return d
end

"""

Find by bisection which is the molar fraction that is consistent
with a desired volume fraction (%), given a data table of densities
as a function of the molar fractions.

"""
function find_x(
    vv, cossolvent_mass, density_pure, densities::Matrix;
    tol=1e-5, maxit=1000
)

    xl = densities[1, 1]
    xr = densities[end, 1]
    if xl < xr
        increasing = true
    else
        increasing = false
    end
    x = (xr - xl) / 2
    xnew = convert_c(
        vv, "%vv" => "x", density=interpolate(x, densities),
        molar_mass=cossolvent_mass, density_pure=density_pure
    )
    it = 0
    while abs(xnew - x) > tol
        if xnew > x
            increasing ? xl = x : xr = x
        else
            increasing ? xr = x : xl = x
        end
        x = (xr + xl) / 2
        xnew = convert_c(
            vv, "%vv" => "x", density=interpolate(x, densities),
            molar_mass=cossolvent_mass, density_pure=density_pure
        )
        it += 1
        it > maxit && error("Maximum number of iterations achieved.")
    end
    return x
end

"""

Convert concentrations.

"""
function convert_c(
    cin, units;
    density=nothing, # solution
    density_pure=nothing, # cossolvent
    molar_mass=18.0, # cossolvent
    molar_mass_water=18.0,
    density_water=1.0
)

    # If the units didn't change, just return the input concentrations
    units[1] == units[2] && return cin

    ρ = density # density of the solution
    ρc = density_pure # density of the pure cossolvent
    Mw = molar_mass_water # molar mass of water
    Mc = molar_mass # molar mass of the cossolvent

    # nc and nw are the molar concentrations
    if units[1] == "%vv"
        if isnothing(ρc)
            error("Density of pure solvent is required to convert from %vv.")
        end
        vv = cin / 100
        nc = (ρc * vv / Mc)
        if units[2] == "x"
            if isnothing(ρ)
                error("Density of solution is required to convert to molar fraction.")
            end
            nw = (ρ - nc * Mc) / Mw
            return nc / (nc + nw)
        end
        if units[2] == "mol/L"
            return 1000 * nc
        end
    end

    if units[1] == "x"
        if isnothing(ρ)
            error("Density of solution is required to convert from molar fraction.")
        end
        x = cin
        nc = ρ / (Mc + Mw * (1 - x) / x)
        if units[2] == "mol/L"
            return 1000 * nc
        end
        if units[2] == "%vv"
            if isnothing(ρc)
                error("Density of pure solvent is required to convert to %vv.")
            end
            nw = nc * (1 - x) / x
            Vc = nc * Mc / ρc
            vv = 100 * Vc
            return vv
        end
    end

    if units[1] == "mol/L"
        if isnothing(ρ)
            error("Density of solution is required to convert from molarity.")
        end
        nc = cin / 1000
        nw = (ρ - nc * Mc) / Mw
        if units[2] == "x"
            return nc / (nc + nw)
        end
        if units[2] == "%vv"
            if isnothing(ρc)
                error("Density of pure solvent is required to convert to %vv.")
                return nothing
            end
            Vc = nc * Mc / ρc
            vv = 100 * Vc
            return vv
        end
    end

end

"""

Function that generates an input file for Packmol. By default, the concentrations is given in mol/L, but it can also be given in molar fraction "x" or volume percentage "%vv", using `cunit="x"` or `cunit="%vv"`. 

"""
function write_input(
    pdbfile::String, solvent_file::String, concentration::Real, box_side::Real;
    water_file="tip4p2005.pdb",
    density=1.0,
    density_pure_solvent=nothing,
    cunit="mol/L",
    packmol_input="box.inp",
    packmol_output="system.pdb")

    protein = readPDB(pdbfile)
    cossolvent = readPDB(solvent_file)

    # molar masses (g/mol)
    Mp = mass(protein)
    Mc = mass(cossolvent)
    Mw = 18.01528 # for water

    # aliases for clearer formulas
    ρ = density # of the solution
    ρc = density_pure_solvent # of the pure solvent

    # Convert concentration to mol/L
    cc_mol = convert_c(concentration, cunit => "mol/L", density=ρ, density_pure=ρc, molar_mass=Mc)
    c_vv = convert_c(concentration, cunit => "%vv", density=ρ, density_pure=ρc, molar_mass=Mc)
    c_x = convert_c(concentration, cunit => "x", density=ρ, density_pure=ρc, molar_mass=Mc)

    # Convert cossolvent concentration in molecules/Å³
    cc = CMC * cc_mol

    # Box volume (Å³)
    vbox = box_side^3

    # Solution volume (vbox - vprotein)
    vs = vbox - CMV * Mp / ρ

    # number of cossolvent molecules: cossolvent concentration × volume of the solution
    nc = round(Int, cc * vs)

    #
    # number of water molecules, obtained from the mass difference
    #
    # nc (molecules) / NA (mol) * Mc (g/mol) is the mass of the cossolvent molecules
    # ρ*vs is the total mass of the solution (density (g/L) × volume (L) 
    # ρ*vs - (nc/NA)*Mc is the mass of water in the solution (g)
    # (ρ*vs - (nc/NA)*Mc)/Mw is the number of mols of water in the solution (mol)
    # NA*(ρ*vs - (nc/NA)*Mc)/Mw is the number of water molcules
    # rearranging: nw = (NA*ρ*vs - nc*Mc)/Mw 
    # But since we have vs in Å³, we need the conversion vs = vs / 1e24, and we have
    # nw = (NA*ρ*vs/1e24 - nc*Mc)/Mw 
    # given that CMV = 1e24/NA, we have nw = (ρ*vs/CMV - nc*Mc)/Mw
    nw = round(Int, (ρ * vs / CMV - nc * Mc) / Mw)

    # Final density of the solution
    ρ = CMV * (Mc * nc + Mw * nw) / vs

    # Final cossolvent concentration (mol/L)
    cc_f = 1000 * (nc / vs) * CMV

    # Final water concentration (mol/L)
    cw_f = 1000 * (nw / vs) * CMV

    # Final recovered concentration in %vv
    vv = 100 * CMV * (nc * Mc / ρc) / vs

    println(
        """
        Summary:
        ========
        Target concentration = $cc_mol mol/L
                             = $c_vv %vv
                             = $c_x x
                             = $cc molecules/Å³

        Box volume = $vbox Å³
        Solution volume = $vs Å³   

        Density = $ρ g/mL
        Protein molar mass = $Mp g/mol
        Cossolvent molar mass = $Mc g/mol

        Number of cossolvent molecules = $nc molecules
        Number of water molecules = $nw molecules

        Final cossolvent concentration = $cc_f mol/L
        Final water concentration = $cw_f mol/L
                                  = $(CMC*cw_f) molecules/Å³
        Final solvent density = $ρ g/mL
        Final %vv concentration = $vv %
        Final molar fraction = $(nc/(nc+nw))
        """)

    l = box_side / 2
    open(packmol_input, "w") do io
        println(io,
            """
            tolerance 2.0
            output $packmol_output
            add_box_sides 1.0
            filetype pdb
            seed -1

            structure $pdbfile
              number 1
              center
              fixed 0. 0. 0. 0. 0. 0.
            end structure

            structure $water_file
              number $nw
              inside box -$l -$l -$l $l $l $l
            end structure
            """)
        if nc > 0
            println(io,
                """
                structure $solvent_file
                  number $nc
                  inside box -$l -$l -$l $l $l $l
                end structure
                """)
        end
    end
    println("Wrote file: $packmol_input")

end # function write_input

end # module

