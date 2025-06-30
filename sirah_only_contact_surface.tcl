# ==============================================================================================
# sirah_only_contact_surface.tcl
# ==============================================================================================
# Description:
#   Computes contact surface areas between two biomolecules (A and B) using SASA method:
#                 Contact Surface (A and B) = 0.5 * [SASA_A + SASA_B - SASA_AuB]
# ==============================================================================================
# Author: Jorge Cantero
# Date: 30-06-2025
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
puts $out "Frame	ContactSurface"

for {set i $start} {$i <= $stop} {incr i $step} {
    molinfo $mol set frame $i
    
    # Compute total contact surface
    set selAB [atomselect $mol "$selA and within 4.45 of $selB"]
    set selBA [atomselect $mol "$selB and within 4.45 of $selA"]
    set area_Total [contact_area $selAB $selBA $probe_radius]
    
    # Write data to file
    puts $out "$i	$area_Total"

    # Cleanup selections
    foreach sel [list $selAB $selBA] {
        $sel delete
    }
    progress_bar $i $nframes
}

close $out

puts "Contact surface Analysis complete"

# ==============================================================================================
