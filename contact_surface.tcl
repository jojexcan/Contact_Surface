# ==============================================================================================
# contact_surface.tcl
# ==============================================================================================
# Description:
#   Computes contact surface areas and estimated pseudo-binding affinity between  
#   two biomolecules (A and B) using a simplified chemical interaction model 
#   with three contact types:
#     - polar–polar (P:P)
#     - nonpolar–nonpolar (NP:NP)
#     - polar–nonpolar (P:NP and NP:P)
#
# Output:
#   per-frame contact surface areas decomposition and estimated surface interaction energy
#
# ==============================================================================================
#
# Author: Jorge Cantero
# Date: 29-06-2025
# version: 1.0
# Contact: Jorge Cantero
#
# ==============================================================================================

# ==============================================================================================
## User settings and configuration
# ==============================================================================================

# Atom selections for molecule A and B
set A "segname A and noh"
set B "segname B and noh"

# Reduce selections around 15 Å between A and B 
set selA "same residue as $A and within 15 of $B"
set selB "same residue as $B and within 15 of $A"

# Probe radius for SASA calculation (Å)
set probe_radius 1.4

# Frame range
set start 0
set stop -1 
set step 1

# Contact type weights (empirical, kcal/mol/Å²)
array set alpha {
    P:P    -0.015
    NP:NP  -0.025
    P:NP   -0.005
}

# Atom type definitions (polar and nonpolar)
set def_P  {name "N.*" "O.*" "P.*" "S.*"}
set def_NP "not ($def_P)"

# Output file
set outfile "contact_surface.dat"

# ==============================================================================================
## Contact Area function
# ==============================================================================================

# Compute SASA of a selection
proc sasa_area {sel probe} {
    return [measure sasa $probe $sel -points no]
}

# Conctact area function
proc contact_area {sel1 sel2 probe} {
    set sasa1 [sasa_area $sel1 $probe]
    set sasa2 [sasa_area $sel2 $probe]
    set sel12 [atomselect [molinfo top] "index [$sel1 get index] [$sel2 get index]"]
    set sasa12 [sasa_area $sel12 $probe]
    $sel12 delete
    return [expr {0.5 * ($sasa1 + $sasa2 - $sasa12)}]
}

# ==============================================================================================
## Progress bar function
# ==============================================================================================

proc progress_bar {current total} {
    clear
    set width 50
    set percent [expr {double($current + 1) / $total}]
    set filled [expr {int($width * $percent)}]
    set empty [expr {$width - $filled}]
    set bar "[string repeat "=" $filled][string repeat " " $empty]"
    puts -nonewline "Progress: \[$bar\] [format %.1f [expr {$percent * 100}]]%"
    flush stdout
    if {$current + 1 >= $total} {
        puts ""
    }
}

# ==============================================================================================
## Execution 
# ==============================================================================================
set mol [molinfo top]
set nframes [molinfo $mol get numframes]
if {$stop < 0} { set stop [expr {$nframes - 1}] }
set out [open $outfile "w"]
puts $out "#Frame	Polar-Polar	NoPolar-NoPolar	Polar-NoPolar	ContactSurface	ContactAffinity"

for {set i $start} {$i <= $stop} {incr i $step} {
    molinfo $mol set frame $i
    set selPA  [atomselect $mol "($selA) and ($def_P)"]
    set selNPA [atomselect $mol "($selA) and ($def_NP)"]
    set selPB  [atomselect $mol "($selB) and ($def_P)"]
    set selNPB [atomselect $mol "($selB) and ($def_NP)"]
    set area_PP   [contact_area $selPA $selPB $probe_radius]
    set area_NPNP [contact_area $selNPA $selNPB $probe_radius]
    set area_PNP1 [contact_area $selPA $selNPB $probe_radius]
    set area_PNP2 [contact_area $selNPA $selPB $probe_radius]
    set area_PNP  [expr {$area_PNP1 + $area_PNP2}]
    
    # Compute total surface
    set totalS [expr {
    $area_PP +
    $area_NPNP +
    $area_PNP
    }]

    # Estimate contact energy
    set totalG [expr {
        $alpha(P:P)*$area_PP +
        $alpha(NP:NP)*$area_NPNP +
        $alpha(P:NP)*$area_PNP
    }]
    

    # Write data to file
    puts $out "$i	$area_PP	$area_NPNP	$area_PNP	$totalS	$totalG"

    # Cleanup selections
    foreach sel [list $selPA $selNPA $selPB $selNPB] {
        $sel delete
    }
    progress_bar $i $nframes
}


close $out

puts "Contact surface Analysis complete"

# ==============================================================================================
