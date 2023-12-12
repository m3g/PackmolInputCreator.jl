module PackmolInputCreator

using PDBTools

export find_molar_fraction
export convert_concentration
export write_input
export interpolate

# conversion factor from mL/mol to Å^3/molecule
const CMV = 1e24 / 6.02e23

# conversion factor from mol/L to molecules/Å^3
const CMC = 6.02e23 / 1e27

"""
    interpolate(x, ρ)

Obtain by interpolation the value of of ρ[:,2] that corresponds
to a the best estimate given a value of x corresponding to 
to the domain ρ[:,1].

"""
function interpolate(x, ρ)
    i = findfirst(d -> d > x, @view(ρ[:, 1]))
    dρdx = (ρ[i, 2] - ρ[i-1, 2]) / (ρ[i, 1] - ρ[i-1, 1])
    d = ρ[i-1, 2] + dρdx * (x - ρ[i-1, 1])
    return d
end

"""
    find_molar_fraction(;
        target_volume_fraction, # [0.0, 1.0] range
        densities::Matrix, # (x, density) data
        cossolvent_molar_mass, # g/mol  
        density_pure_cossolvent, # g/mL
        solvent_molar_mass = 18.0, # g/mol default: water
        density_pure_solvent = 1.0, # g/mL default: water
        tol=1e-5, maxit=1000, debug=false
    )

Find by bisection which is the molar fraction that is consistent
with a desired volume fraction (%), given a data table of densities
as a function of the molar fractions.

The input parameter are:

- `target_volume_fraction`: desired volume fraction 
- `cossolvent_molar_mass`: molar mass of the cossolvent
- `density_pure_cossolvent`: density of the pure solvent
- `densities`: a matrix with two columns, the first one with the
   molar fractions and the second one with the densities of the
   mixture.
- `solvent_molar_mass`: molar mass of the solvent (default: 18.0 - water)
- `density_pure_solvent`: density of the pure solvent (default: 1.0 - water)

Optional parameters:

- `tol`: tolerance for the bisection method
- `maxit`: maximum number of iterations
- `debug`: true/false to print debug information

"""
function find_molar_fraction(;
    target_volume_fraction,
    densities,
    cossolvent_molar_mass,
    density_pure_cossolvent,
    solvent_molar_mass = 18.0,
    density_pure_solvent = 1.0,
    tol=1e-5, maxit=1000, debug = false
)
    vv = target_volume_fraction 
    if !(0 < vv < 1)
        throw(ArgumentError("Target volume fraction must be in the [0,1] range."))
    end
    increasing = densities[1, 1] < densities[end, 1]
    if increasing
        xl, xr = densities[1, 1], densities[end, 1]
    else
        xl, xr = densities[end, 1], densities[1, 1]
    end
    x = (xr - xl) / 2
    xnew = convert_concentration(
        vv, "vv" => "x";
        density=interpolate(x, densities),
        cossolvent_molar_mass=cossolvent_molar_mass, 
        density_pure_cossolvent=density_pure_cossolvent,
        solvent_molar_mass=solvent_molar_mass,
        density_pure_solvent=density_pure_solvent,
    )
    it = 0
    debug && println("it = $it, x = $x, xnew = $xnew")
    while abs(xnew - x) > tol
        if xnew > x
            increasing ? xl = x : xr = x
        else
            increasing ? xr = x : xl = x
        end
        x = (xr + xl) / 2
        xnew = convert_concentration(
            vv, "vv" => "x", 
            density=interpolate(x, densities),
            cossolvent_molar_mass=cossolvent_molar_mass, 
            density_pure_cossolvent=density_pure_cossolvent,
            solvent_molar_mass=solvent_molar_mass,
            density_pure_solvent=density_pure_solvent,
        )
        it += 1
        debug && println("it = $it, x = $x, xnew = $xnew")
        it > maxit && error("Maximum number of iterations achieved.")
    end
    return x
end

"""
    convert_concentration(
        cin, units;
        density=nothing, # solution
        density_pure_cossolvent=nothing,
        cossolvent_molar_mass=nothing,
        solvent_molar_mass=18.0, # defaults to water
    )

Convert concentration from one unit to another. The input
concentration is given in `cin`, and the unit conversion 
is given by the pairs:

- `"mol/L"` => `"x"`: from molarity to molar fraction
- `"mol/L"` => `"vv"`: from molarity to volume fraction
- `"x"` => `"mol/L"`: from molar fraction to molarity
- `"x"` => `"vv"`: from molar fraction to volume fraction
- `"vv"` => `"mol/L"`: from volume fraction to molarity
- `"vv"` => `"x"`: from volume fraction to molar fraction

Depending on the conversion required, different combinations of 
densities as inputs are required. Appropriate error messages
are thrown if the required densities are not provided.

By default, the solvent is water, with a molar mass of 18.0 g/mol.

"""
function convert_concentration(
    cin, units;
    density=nothing, # of the solution
    density_pure_cossolvent=nothing,
    cossolvent_molar_mass=nothing,
    solvent_molar_mass=18.0, # default to water
)

    # If the units didn't change, just return the input concentrations
    units[1] == units[2] && return cin

    ρ = density # density of the solution
    ρc = density_pure_cossolvent # density of the pure cossolvent
    Mw = solvent_molar_mass # molar mass of solvent 
    Mc = cossolvent_molar_mass # molar mass of the cossolvent

    # nc and nw are the molar concentrations
    if units[1] == "vv"
        if isnothing(ρc)
            throw(ArgumentError("Density of pure cossolvent is required to convert from volume fraction."))
        end
        vv = cin
        nc = (ρc * vv / Mc)
        if units[2] == "x"
            if isnothing(ρ)
                throw(ArgumentError("Density of solution is required to convert to molar fraction."))
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
            throw(ArgumentError("Density of solution is required to convert from molar fraction."))
        end
        x = cin
        nc = ρ / (Mc + Mw * (1 - x) / x)
        if units[2] == "mol/L"
            return 1000 * nc
        end
        if units[2] == "vv"
            if isnothing(ρc)
                throw(ArgumentError("Density of pure cossolvent is required to convert to volume fraction."))
            end
            nw = nc * (1 - x) / x
            vv = nc * Mc / ρc
            return vv
        end
    end

    if units[1] == "mol/L"
        if isnothing(ρ)
            throw(ArgumentError("Density of solution is required to convert from molarity."))
        end
        nc = cin / 1000
        nw = (ρ - nc * Mc) / Mw
        if units[2] == "x"
            return nc / (nc + nw)
        end
        if units[2] == "vv"
            if isnothing(ρc)
                throw(ArgumentError("Density of pure cossolvent is required to convert to volume fraction."))
            end
            vv = nc * Mc / ρc
            return vv
        end
    end

end

"""
    write_input(;
        solute_pdbfile::String, 
        solvent_pdbfile::String,
        cossolvent_pdbfile::String,
        concentration::Real, 
        box_side::Real,
        density_solution=1.0,
        density_pure_solvent=nothing,
        density_pure_cossolvent=nothing,
        cunit="mol/L",
        packmol_input="box.inp",
        packmol_output="system.pdb"
    )

Function that generates an input file for Packmol. By default, 
the concentrations is given in `mol/L`, but it can also be 
given in molar fraction "x" or volume fraction "vv", using `cunit="x"` or `cunit="vv"`. 

Depending on the desired output concentration, the densities of the pure
components of the solution are required. 

"""
function write_input(;
    solute_pdbfile::String, 
    solvent_pdbfile::String,
    cossolvent_pdbfile::String,
    concentration::Real, 
    box_side::Real,
    density_solution=1.0,
    density_pure_solvent=nothing,
    density_pure_cossolvent=nothing,
    cunit="mol/L",
    packmol_input="box.inp",
    packmol_output="system.pdb"
)

    solute = readPDB(solute_pdbfile)
    solvent = readPDB(solvent_pdbfile)
    cossolvent = readPDB(cossolvent_pdbfile)

    # molar masses (g/mol)
    Mp = mass(solute) # Mp stands for "protein": solute
    Mc = mass(cossolvent) # Mc stands for "cossolvent"
    Mw = mass(solvent) # Mw is for "water": solvent

    # aliases for clearer formulas
    ρ = density_solution # of the solution
    ρc = density_pure_cossolvent # of the pure cossolvent

    # Convert concentration to mol/L
    cc_mol = convert_concentration(
        concentration, 
        cunit => "mol/L", 
        density=ρ, 
        density_pure_cossolvent=ρc, 
        cossolvent_molar_mass=Mc
    )
    c_vv = convert_concentration(
        concentration, 
        cunit => "vv", 
        density=ρ, 
        density_pure_cossolvent=ρc, 
        cossolvent_molar_mass=Mc
    )
    c_x = convert_concentration(
        concentration, 
        cunit => "x", 
        density=ρ, 
        density_pure_cossolvent=ρc, 
        cossolvent_molar_mass=Mc
    )

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

    # Final recovered concentration in vv
    vv = CMV * (nc * Mc / ρc) / vs

    println(
        """

        Summary:
        ========
        Target concentration = $cc_mol mol/L
                             = $c_vv volume fraction
                             = $c_x molar fraction
                             = $cc molecules/Å³

        Box volume = $vbox Å³
        Solution volume = $vs Å³   

        Density = $ρ g/mL
        Solute molar mass = $Mp g/mol
        Cossolvent molar mass = $Mc g/mol
        Solvent molar mass = $Mw g/mol

        Number of cossolvent ($cossolvent_pdbfile) molecules = $nc molecules
        Number of solvent ($solvent_pdbfile) molecules = $nw molecules

        Final cossolvent concentration = $cc_f mol/L
        Final solvent concentration = $cw_f mol/L
                                    = $(CMC*cw_f) molecules/Å³
                                    
        Final solvent density = $ρ g/mL
        Final volume fraction = $vv 
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

            structure $solute_pdbfile
              number 1
              center
              fixed 0. 0. 0. 0. 0. 0.
            end structure

            structure $solvent_pdbfile
              number $nw
              inside box -$l -$l -$l $l $l $l
            end structure
            """)
        if nc > 0
            println(io,
                """
                structure $cossolvent_pdbfile
                  number $nc
                  inside box -$l -$l -$l $l $l $l
                end structure
                """)
        end
    end
    println("Wrote file: $packmol_input")

end # function write_input

end # module

