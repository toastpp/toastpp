# =============================================================================
# User interface for input of circle mesh parameters
# =============================================================================

# Default mesh parameters
set eltp 0
set cubetp 0
set xmin -10
set ymin -10
set zmin -10
set xmax 10
set ymax 10
set zmax 10
set xgrid 10
set ygrid 10
set zgrid 10

proc GenBrickDlg {} {
    global res eltp cubetp xmin ymin zmin xmax ymax zmax xgrid ygrid zgrid

    set dlg [toplevel .genbrickdlg]
    wm title $dlg "Create brick mesh"

    set name [tixOptionName $dlg]
    option add *$name*TixControl*entry.width 6
    option add *$name*TixControl*label.width 6
    option add *$name*TixControl*label.anchor e

    frame $dlg.size -relief groove -bd 2
    label $dlg.size.l -text "Physical dimensions" -anchor w
    pack $dlg.size.l -side top -fill x
    foreach grp {min max} {
	frame $dlg.size.$grp
	foreach dim {x y z} {
	    tixControl $dlg.size.$grp.$dim -label "$dim$grp" \
		    -variable "$dim$grp"
	    pack $dlg.size.$grp.$dim -side left
	}
	pack $dlg.size.$grp -side top -fill x
    }
    pack $dlg.size -side top -fill x -padx 5 -pady 5

    frame $dlg.grid -relief groove -bd 2
    label $dlg.grid.l -text "Cell grid dimensions" -anchor w
    set grp grid
    frame $dlg.grid.$grp
    foreach dim {x y z} {
	tixControl $dlg.grid.$grp.$dim -label $dim -variable $dim$grp
	pack $dlg.grid.$grp.$dim -side left
    }
    pack $dlg.grid.l $dlg.grid.$grp -side top -fill x -expand true
    pack $dlg.grid -side top -fill x -padx 5 -pady 5

    frame $dlg.eltp -relief groove -bd 2
    label $dlg.eltp.l -text "Element type" -anchor w
    radiobutton $dlg.eltp.1 -variable eltp -text "4-noded tetrahedra" -value 0
    radiobutton $dlg.eltp.2 -variable eltp -text "10-noded tetrahedra" -value 1
    radiobutton $dlg.eltp.3 -variable eltp -text "8-noded voxels" -value 2
    pack $dlg.eltp.l -fill x
    pack $dlg.eltp.1 $dlg.eltp.2 $dlg.eltp.3 -side top -anchor w
    pack $dlg.eltp -side top -fill x -padx 5 -pady 5

    frame $dlg.cubetp -relief groove -bd 2
    label $dlg.cubetp.l -text "Tetrahedra per element cube" -anchor w
    radiobutton $dlg.cubetp.1 -variable cubetp -text "5" -value 0
    radiobutton $dlg.cubetp.2 -variable cubetp -text "6" -value 1
    radiobutton $dlg.cubetp.3 -variable cubetp -text "24" -value 2
    pack $dlg.cubetp.l -fill x
    pack $dlg.cubetp.1 $dlg.cubetp.2 $dlg.cubetp.3 -side left -fill x
    pack $dlg.cubetp -side top -fill x -padx 5 -pady 5

    # Button bar ============================================
    tixButtonBox $dlg.bbar -orientation horizontal -pady 2 -relief flat
    $dlg.bbar add apply -text OK -underline 0 -width 6 \
	    -command [list GenBrickApply $dlg]
    $dlg.bbar add done -text Cancel -underline 0 -width 6 \
	    -command [list GenBrickCancel $dlg]
    pack $dlg.bbar -side bottom -fill x

    tkwait window $dlg
    return $res
}

proc GenBrickApply {w} {
    global res xmin ymin zmin xmax ymax zmax xgrid ygrid zgrid eltp cubetp
    foreach grp {min max} {
	foreach dim {x y z} {
	    $w.size.$grp.$dim update
	}
    }
    set grp grid
    foreach dim {x y z} {
	$w.grid.$grp.$dim update
    }
    GenBrickMesh $xmin $ymin $zmin $xmax $ymax $zmax $xgrid $ygrid $zgrid \
	    $eltp $cubetp
    set res 1
    destroy $w
}

proc GenBrickCancel {w} {
    global res
    set res 0
    destroy $w
}
