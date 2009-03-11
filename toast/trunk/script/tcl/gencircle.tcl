# =============================================================================
# User interface for input of circle mesh parameters
# =============================================================================

# Default mesh parameters
set rad 25
set nring 10
set nsect 6
set nbnd 1

proc GenCircleDlg {} {
    global res rad nring nsect nbnd

    set dlg [toplevel .gencircledlg]
    wm title $dlg "Create circle mesh"

    set name [tixOptionName $dlg]
    option add *$name*TixControl*entry.width 7
    option add *$name*TixControl*label.width 14
    option add *$name*TixControl*label.anchor e

    frame $dlg.size -relief groove -bd 2
    label $dlg.size.l -text "Physical parameters" -anchor w
    pack  $dlg.size.l -side top -fill x
    tixControl $dlg.size.rad -label "Radius" -variable rad
    pack  $dlg.size.rad -side left
    pack  $dlg.size -side top -fill x -padx 5 -pady 5

    frame $dlg.grid -relief groove -bd 2
    label $dlg.grid.l -text "Mesh grid parameters" -anchor w
    pack  $dlg.grid.l -side top -fill x
    tixControl $dlg.grid.nring -label "Element rings" -variable nring \
	-integer yes -min 2
    tixControl $dlg.grid.nsect -label "Sectors" -variable nsect \
	-integer yes -min 4
    tixControl $dlg.grid.nbnd -label "Boundary layers" -variable nbnd \
	-integer yes -min 1
    pack  $dlg.grid.nring $dlg.grid.nsect $dlg.grid.nbnd -side top
    pack  $dlg.grid -side top -fill x -padx 5 -pady 5

    # Button bar ============================================
    tixButtonBox $dlg.bbar -orientation horizontal -pady 2 -relief flat
    $dlg.bbar add apply -text OK -underline 0 -width 6 \
	    -command [list GenCircleApply $dlg]
    $dlg.bbar add done -text Cancel -underline 0 -width 6 \
	    -command [list GenCircleCancel $dlg]
    pack $dlg.bbar -side bottom -fill x

    tkwait window $dlg
    return $res
}

proc GenCircleApply {w} {
    global res rad nring nsect nbnd
    $w.size.rad update
    $w.grid.nring update
    $w.grid.nsect update
    $w.grid.nbnd update
    GenCircleMesh $rad $nsect $nring $nbnd
    set res 1
    destroy $w
}

proc GenCircleCancel {w} {
    global res
    set res 0
    destroy $w
}