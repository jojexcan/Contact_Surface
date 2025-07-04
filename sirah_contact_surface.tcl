# ==============================================================================================
# sirah_contact_surface.tcl
# ==============================================================================================
# Description:
#   Computes contact surface areas and estimated pseudo-binding affinity between  
#   two biomolecules (A and B) using a simplified chemical interaction model 
#   with three contact types:
#     - polar–polar (P:P)
#     - nonpolar–nonpolar (NP:NP)
#     - polar–nonpolar (P:NP and NP:P)
# ==============================================================================================
# Author: Jorge Cantero
# Date: 29-06-2025
# version: 1.0
# Contact: jorgec@fq.edu.uy
# ==============================================================================================

# ==============================================================================================
## User settings and configuration
# ==============================================================================================

# Atom selections for molecule A and B
set selA "sirah_protein"
set selB "sirah_nucleic"

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
set def_P  {{charge >= 0.2 or charge <= -0.2} or name GN l2 K3 K4}
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
puts $out "Frame    Polar-Polar    NoPolar-NoPolar   Polar-NoPolar    ContactSurface    Affinity"

for {set i $start} {$i <= $stop} {incr i $step} {
    molinfo $mol set frame $i
    
    # Compute total contact surface
    set selAB [atomselect $mol "$selA and within 6.1 of $selB"]
    set selBA [atomselect $mol "$selB and within 6.1 of $selA"]
    set area_Total [contact_area $selAB $selBA $probe_radius]
    
    #Polar selection
    set selPA  [atomselect $mol "($selA) and ($def_P)"]
    set selPB  [atomselect $mol "($selB) and ($def_P)"]

    #NonPolar selection
    set selNPA [atomselect $mol "($selA) and ($def_NP)"]
    set selNPB [atomselect $mol "($selB) and ($def_NP)"]

    #Polar-polar interface
    #Polar A - Polar B
    set PAPB [atomselect $mol "index [$selPA get index] and within 6.1 of index [$selPB get index]"]
    set PBPA [atomselect $mol "index [$selPB get index] and within 6.1 of index [$selPA get index]"]

    #Surface Polar - Polar interface
    set area_PP   [contact_area $PAPB $PBPA $probe_radius]

    #NonPolar-NonPolar interface
    #NonPolar A - NonPolar B
    set NPA [atomselect $mol "index [$selNPA get index] and within 6.1 of index [$selNPB get index]"]
    set NPB [atomselect $mol "index [$selNPB get index] and within 6.1 of index [$selNPA get index]"]

    #Surface NonPolar - NonPolar interface
    set area_NPNP [contact_area $NPA $NPB $probe_radius]

    #Polar-Nonpolar interface
    #Polar A - NonPolar B
    set PANPB [atomselect $mol "index [$selPA get index] and within 6.1 of index [$selNPB get index]"]
    set NPBPA [atomselect $mol "index [$selNPB get index] and within 6.1 of index [$selPA get index]"]

    #Surface Polar A - NonPolar B
    set area_PNP1 [contact_area $PANPB $NPBPA $probe_radius]

    #NonPolar A - Polar B
    set NPAPB [atomselect $mol "index [$selNPA get index] and within 6.1 of index [$selPB get index]"]
    set PBNPA [atomselect $mol "index [$selPB get index] and within 6.1 of index [$selNPA get index]"]

    #Surface NonPolar A - Polar B
    set area_PNP2 [contact_area $NPAPB $PBNPA $probe_radius]

    #Surface Polar - NonPolar interface
    set area_PNP  [expr {$area_PNP1 + $area_PNP2}]
    
    # Estimate contact energy
    set totalG [expr {
        $alpha(P:P)*$area_PP +
        $alpha(NP:NP)*$area_NPNP +
        $alpha(P:NP)*$area_PNP
    }]
    

    # Write data to file
    puts $out "$i    $area_PP    $area_NPNP    $area_PNP    $area_Total    $totalG"

    # Cleanup selections
    foreach sel [list $selPA $selNPA $selPB $selNPB $PAPB $PBPA $NPA $NPB $PANPB $NPBPA $NPAPB $PBNPA] {
        $sel delete
    }
    progress_bar $i $nframes
}


close $out

puts "Contact surface Analysis complete"

# ==============================================================================================
