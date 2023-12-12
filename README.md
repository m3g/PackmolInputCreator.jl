# PackmolInputCreator

This packages helps creating the input files for Packmol, by computing the number of 
molecules necessary to build a box with a desired concentration. 

## Installation

This is a Julia package. Install Julia, version 1.6 or higher, from the [official download site](https://julialang.org/downloads/). 

```julia
julia> using Pkg

julia> Pkg.add(url="https://github.com/m3g/PackmolInputCreator.jl")
```

## Example

```julia
#
# run with julia script.jl
#
using PDBTools
using PackmolInputCreator

# Density as a function of the molar fraction
# 1.0 means pure solvent, in the terminology used
# in this script. 0.0 means pure cossolvent. 
#     x (solvent)       ρ / g/mL
ρ = [ 0.0000 0.834 # 100% cyclooctane (pure cossolvent) 
      0.0155 0.82950177137763 
      0.0301 0.82719220548717
      0.0394 0.82573745713931
      0.0523 0.82373025397815
      0.0596 0.82260208345798
      0.0983 0.81671172734072
      0.1431 0.81006305253341
      0.1863 0.80383951557194
      0.2213 0.79887269483570
      0.2522 0.79454950559211
      0.2801 0.79070471811816
      0.3126 0.78625680868169
      0.3365 0.78304392274268
      0.3681 0.77881926482676
      0.3964 0.77505872004072
      0.4112 0.77323017810059
      0.5196 0.75922947979487
      0.5613 0.75400203132561
      0.5906 0.75038054304765
      0.6226 0.74646150275651
      0.6541 0.74264222538436
      0.6905 0.73827093000787
      0.7295 0.73363732026850
      0.7751 0.72828702684542
      0.8161 0.72351198260744
      0.8832 0.71588112250305
      0.9225 0.71151341183697
      0.9403 0.70953358268734
      0.9573 0.70768096611118
      0.9778 0.70542722961990 
      1.0000 0.703 # 100% octane (pure solvent)
]

# Directory where this script is hosted
script_dir = @__DIR__

# PDB files
data_dir="$script_dir/packmol"
println("Input directory: $data_dir")
solute_pdbfile = "$data_dir/poly_h.pdb"
solvent_pdbfile = "$data_dir/octane.pdb"
cossolvent_pdbfile = "$data_dir/cyclooctane.pdb"
box_size = 120.0

# Find to what molar fraction this volume fraction corresponds
x = find_molar_fraction(;
        target_volume_fraction = 0.5, # what we want 50%vv
        cossolvent_molar_mass = mass(readPDB(cossolvent_pdbfile)),
        density_pure_cossolvent = ρ[1,2], # g/mL
        solvent_molar_mass = mass(readPDB(solvent_pdbfile)),
        density_pure_solvent = ρ[end,2], # g/mL 
        densities = ρ,
)
println("Molar fraction = $x")

# Iterpolate to get density given molar fraction
density = interpolate(x,ρ)
println("Density = $density")

# Create input file
write_input(
    solute_pdbfile = "./packmol/poly_h.pdb",
    solvent_pdbfile = "./packmol/octane.pdb",
    cossolvent_pdbfile="./packmol/cyclooctane.pdb",
    concentration=0.5,
    box_side = 120.0,
    density_solution=density,
    density_pure_cossolvent=ρ[1,2],
    cunit="x",
)
```

