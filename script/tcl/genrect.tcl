# =============================================================================
# User interface for input of rectangle mesh parameters
# =============================================================================

# Default mesh parameters
set celltp 1
set xmin -10
set xmax 10
set ymin -10
set ymax 10
set xdim 10
set ydim 10

proc GenRectDlg {} {
    global res xmin ymin xmax ymax xdim ydim celltp

    set dlg [toplevel .genrectdlg]
    wm title $dlg "Create rectangular mesh"

    set name [tixOptionName $dlg]
    option add *$name*TixControl*entry.width 7
    option add *$name*TixControl*label.width 9
    option add *$name*TixControl*label.anchor e

    frame $dlg.size -relief groove -bd 2
    label $dlg.size.l -text "Physical parameters" -anchor w
    pack  $dlg.size.l -side top -fill x
    frame $dlg.size.x
    tixControl $dlg.size.x.min -label "xmin" -variable xmin
    tixControl $dlg.size.x.max -label "xmax" -variable xmax
    pack  $dlg.size.x.min $dlg.size.x.max -side left
    frame $dlg.size.y
    tixControl $dlg.size.y.min -label "ymin" -variable ymin
    tixControl $dlg.size.y.max -label "ymax" -variable ymax
    pack  $dlg.size.y.min $dlg.size.y.max -side left
    pack  $dlg.size.x $dlg.size.y -side top
    pack  $dlg.size -side top -fill x -padx 5 -pady 5

    frame $dlg.grid -relief groove -bd 2
    label $dlg.grid.l -text "Mesh grid parameters" -anchor w
    pack  $dlg.grid.l -side top -fill x
    tixControl $dlg.grid.x -label "x grid dim" -variable xdim \
	-integer yes -min 2
    tixControl $dlg.grid.y -label "y grid dim" -variable ydim \
	-integer yes -min 2
    pack $dlg.grid.x $dlg.grid.y -side top -fill x
    pack $dlg.grid -side top -fill x -padx 5 -pady 5

    frame $dlg.celltp -relief groove -bd 2
    label $dlg.celltp.l -text "Cell type" -anchor w
    pack  $dlg.celltp.l -side top -fill x
    radiobutton $dlg.celltp.1 -variable celltp -text "2 triangles" -value 0
    radiobutton $dlg.celltp.2 -variable celltp -text "4 triangles" -value 1
    radiobutton $dlg.celltp.3 -variable celltp -text "pixel" -value 2
    pack $dlg.celltp.1 $dlg.celltp.2 $dlg.celltp.3 -side left -fill x
    pack $dlg.celltp -side top -fill x -padx 5 -pady 5

    # Button bar ============================================
    tixButtonBox $dlg.bbar -orientation horizontal -pady 2 -relief flat
    $dlg.bbar add apply -text OK -underline 0 -width 6 \
	    -command [list GenRectApply $dlg]
    $dlg.bbar add done -text Cancel -underline 0 -width 6 \
	    -command [list GenRectCancel $dlg]
    pack $dlg.bbar -side bottom -fill x

    tkwait window $dlg
    return $res
}

proc GenRectApply {w} {
    global res xmin ymin xmax ymax xdim ydim celltp
    foreach dim {x y} {
	foreach grp {min max} {
	    $w.size.$dim.$grp update
	}
	$w.grid.$dim update
    }
    GenRectMesh $xmin $ymin $xmax $ymax $xdim $ydim $celltp
    set res 1
    destroy $w
}

proc GenRectCancel {w} {
    global res
    set res 0
    destroy $w
}
