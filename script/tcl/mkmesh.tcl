# ============================================================================
# Tcl interface for mesh generation routines

# Dialog for brick mesh parameters

proc GenBrickDlg {prmVar} {
    upvar $prmVar prm
    set dlg [toplevel .genbrickdlg]
    wm title $dlg "MeshGenerator: Brick parameters"

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
}
