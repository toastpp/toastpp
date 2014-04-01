wm title . "mkmesh_brick"

# ===========================================================================
# Globals

set name "mesh.msh"
set xmin "-10"
set ymin "-10"
set zmin "-10"
set xmax "10"
set ymax "10"
set zmax "10"
set xgrid 10
set ygrid 10
set zgrid 10
set eltp 0
set cubetp 0
set hmua 0.025
set hmus 2
set hn 1.4
set hreg 0
set pmua 0.025
set pmus 2
set pref 1.4
set preg 0
set pcnt "0 0 0"
set prad 1
set nds 0
set els 0
set bnds 0
set inds 0

# ===========================================================================
# scrolled listbox

proc Scroll_Set {scrollbar geoCmd offset size} {
    if {$offset != 0.0 || $size != 1.0} {
	eval $geoCmd     ;# Make sure it is visible
	$scrollbar set $offset $size
    } else {
	set manager [lindex $geoCmd 0]
	$manager forget $scrollbar ;# hide it
    }
}

proc Scrolled_Listbox { f args } {
    frame $f
    listbox $f.list \
	-xscrollcommand [list Scroll_Set $f.xscroll \
			[list grid $f.xscroll -row 1 -column 0 -sticky we]] \
	-yscrollcommand [list Scroll_Set $f.yscroll \
			[list grid $f.yscroll -row 0 -column 1 -sticky ns]]
    eval {$f.list configure} $args
    scrollbar $f.xscroll -orient horizontal -command [list $f.list xview]
    scrollbar $f.yscroll -orient vertical -command [list $f.list yview]
    grid $f.list $f.yscroll -sticky news
    grid $f.xscroll -sticky news
    grid rowconfigure $f 0 -weight 1
    grid columnconfigure $f 0 -weight 1
    return $f.list
}

# ===========================================================================
# parameter perturbation dialog

proc PertDlg { lb edit } {
    toplevel .f
    global pmua pmus pref preg pcnt prad
    wm title .f "Perturbation"
    frame .f.prm -relief groove -bd 2
    label .f.prm.l -text "Perturbation parameters" -anchor w
    frame .f.prm.val
    label .f.prm.val.tmua -text "MUA" -width 5 -anchor e
    entry .f.prm.val.emua -textvariable pmua -relief sunken -width 10
    label .f.prm.val.tmus -text "MUS" -width 5 -anchor e
    entry .f.prm.val.emus -textvariable pmus -relief sunken -width 10
    label .f.prm.val.tn   -text "N" -width 5 -anchor e
    entry .f.prm.val.en   -textvariable pref -relief sunken -width 10
    pack  .f.prm.val.tmua .f.prm.val.emua .f.prm.val.tmus .f.prm.val.emus \
	  .f.prm.val.tn .f.prm.val.en -side left
    frame .f.prm.reg
    label .f.prm.reg.l -text "Region"
    entry .f.prm.reg.reg -textvariable preg -relief sunken -width 10
    pack  .f.prm.reg.l .f.prm.reg.reg -side left -fill x
    pack  .f.prm.l .f.prm.val -fill x
    pack  .f.prm.reg -side right -fill x
    pack  .f.prm -side top -fill x -padx 4 -pady 2

    frame .f.geom -relief groove -bd 2
    label .f.geom.l -text "Perturbation dimensions" -anchor w
    frame .f.geom.val
    frame .f.geom.val.f1
    label .f.geom.val.f1.l1 -text "Centre \[x y z\]" -anchor e
    label .f.geom.val.f1.l2 -text "Radius" -anchor e
    pack  .f.geom.val.f1.l1 .f.geom.val.f1.l2 -side top -fill x
    frame .f.geom.val.f2
    entry .f.geom.val.f2.e1 -textvariable pcnt -relief sunken -width 10
    entry .f.geom.val.f2.e2 -textvariable prad -relief sunken -width 10
    pack  .f.geom.val.f2.e1 .f.geom.val.f2.e2 -side top -fill x
    pack  .f.geom.val.f1 .f.geom.val.f2 -side left
    pack  .f.geom.l .f.geom.val -fill x
    pack  .f.geom -side top -fill x -padx 4 -pady 2

    frame .f.bbar
    button .f.bbar.ok -text "OK" -command "AddPert $lb $edit; destroy .f"
    button .f.bbar.cancel -text "Cancel" -command { destroy .f }
    pack  .f.bbar.cancel .f.bbar.ok -side right -fill x
    pack  .f.bbar -side top -fill x -padx 4 -pady 2
}

# ===========================================================================

proc AddPert { lb edit } {
    global pmua pmus pref preg pcnt prad
    if {$edit} {
	$lb delete active
    }
    $lb insert active "mua=$pmua\; mus=$pmus\; n=$pref\; reg=$preg\; cnt=\[$pcnt\]\; rad=$prad"
}

# ===========================================================================

proc PertDelete { lb } {
    $lb delete active
}

# ===========================================================================
# Geometry dialog

proc GeomDlg {} {
    global f
    set f [frame .opt]
    frame $f.name -relief groove -bd 2
    label $f.name.l -text "Name" -anchor w
    entry $f.name.e -textvariable name -relief sunken
    pack $f.name.l $f.name.e -fill x
    pack $f.name -side top -fill x -padx 4 -pady 2

    frame $f.size -relief groove -bd 2
    label $f.size.l -text "Physical dimensions" -anchor w
    frame $f.size.val
    frame $f.size.val.b1
    label $f.size.val.b1.txmin -text "xmin" -width 5 -anchor e
    label $f.size.val.b1.txmax -text "xmax" -width 5 -anchor e
    pack  $f.size.val.b1.txmin $f.size.val.b1.txmax -side top
    frame $f.size.val.b2
    entry $f.size.val.b2.exmin -textvariable xmin -relief sunken -width 10
    entry $f.size.val.b2.exmax -textvariable xmax -relief sunken -width 10
    pack  $f.size.val.b2.exmin $f.size.val.b2.exmax -side top
    frame $f.size.val.b3
    label $f.size.val.b3.tymin -text "ymin" -width 5 -anchor e
    label $f.size.val.b3.tymax -text "ymax" -width 5 -anchor e
    pack  $f.size.val.b3.tymin $f.size.val.b3.tymax -side top
    frame $f.size.val.b4
    entry $f.size.val.b4.eymin -textvariable ymin -relief sunken -width 10
    entry $f.size.val.b4.eymax -textvariable ymax -relief sunken -width 10
    pack  $f.size.val.b4.eymin $f.size.val.b4.eymax -side top
    frame $f.size.val.b5
    label $f.size.val.b5.tzmin -text "zmin" -width 5 -anchor e
    label $f.size.val.b5.tzmax -text "zmax" -width 5 -anchor e
    pack  $f.size.val.b5.tzmin $f.size.val.b5.tzmax -side top
    frame $f.size.val.b6
    entry $f.size.val.b6.ezmin -textvariable zmin -relief sunken -width 10
    entry $f.size.val.b6.ezmax -textvariable zmax -relief sunken -width 10
    pack  $f.size.val.b6.ezmin $f.size.val.b6.ezmax -side top
    pack  $f.size.val.b1 $f.size.val.b2 $f.size.val.b3 $f.size.val.b4 \
	  $f.size.val.b5 $f.size.val.b6 -side left
    pack $f.size.l $f.size.val -fill x
    pack $f.size -side top -fill x -padx 4 -pady 2

    frame $f.grid -relief groove -bd 2
    label $f.grid.l -text "Element cell grid dimensions" -anchor w
    frame $f.grid.val
    label $f.grid.val.tx -text "x" -width 5 -anchor e
    entry $f.grid.val.ex -textvariable xgrid -relief sunken -width 10
    label $f.grid.val.ty -text "y" -width 5 -anchor e
    entry $f.grid.val.ey -textvariable ygrid -relief sunken -width 10
    label $f.grid.val.tz -text "z" -width 5 -anchor e
    entry $f.grid.val.ez -textvariable zgrid -relief sunken -width 10
    pack  $f.grid.val.tx $f.grid.val.ex $f.grid.val.ty $f.grid.val.ey \
	  $f.grid.val.tz $f.grid.val.ez -side left
    pack $f.grid.l $f.grid.val -fill x
    pack $f.grid -side top -fill x -padx 4 -pady 2

    frame $f.eltp -relief groove -bd 2
    label $f.eltp.l -text "Element type" -anchor w
    radiobutton $f.eltp.1 -variable eltp -text "4-noded tetrahedron" -value 0
    radiobutton $f.eltp.2 -variable eltp -text "10-noded tetrahedron" -value 1
    pack $f.eltp.l -fill x
    pack $f.eltp.1 $f.eltp.2 -side left -fill x
    pack $f.eltp -side top -fill x -padx 4 -pady 2

    frame $f.cubetp -relief groove -bd 2
    label $f.cubetp.l -text "Tetrahedra per element cube" -anchor w
    radiobutton $f.cubetp.1 -variable cubetp -text "5" -value 0
    radiobutton $f.cubetp.2 -variable cubetp -text "6" -value 1
    radiobutton $f.cubetp.3 -variable cubetp -text "24" -value 2
    pack $f.cubetp.l -fill x
    pack $f.cubetp.1 $f.cubetp.2 $f.cubetp.3 -side left -fill x
    pack $f.cubetp -side top -fill x -padx 4 -pady 2

    frame $f.bbar
    button $f.bbar.next -text "Next" -command { destroy .opt; ParamDlg }
    button $f.bbar.exit -text "Exit" -command exit
    pack $f.bbar.exit $f.bbar.next -side right -fill x
    pack $f.bbar -side top -fill x -padx 4 -pady 2

    pack $f -side top
}

# ===========================================================================
# Parameter dialog

proc ParamDlg {} {
    global f pvalid
    set f [frame .opt]
    frame $f.bg -relief groove -bd 2
    label $f.bg.l -text "Background parameters" -anchor w
    frame $f.bg.val
    label $f.bg.val.tmua -text "MUA" -width 5 -anchor e
    entry $f.bg.val.emua -textvariable hmua -relief sunken -width 10
    label $f.bg.val.tmus -text "MUS" -width 5 -anchor e
    entry $f.bg.val.emus -textvariable hmus -relief sunken -width 10
    label $f.bg.val.tn   -text "N" -width 5 -anchor e
    entry $f.bg.val.en   -textvariable hn -relief sunken -width 10
    pack  $f.bg.val.tmua $f.bg.val.emua $f.bg.val.tmus $f.bg.val.emus \
	  $f.bg.val.tn $f.bg.val.en -side left
    frame $f.bg.reg
    label $f.bg.reg.l -text "Region"
    entry $f.bg.reg.reg -textvariable hreg -relief sunken -width 10
    pack  $f.bg.reg.l $f.bg.reg.reg -side left -fill x
    pack  $f.bg.l $f.bg.val -fill x
    pack  $f.bg.reg -side right -fill x
    pack  $f.bg -side top -fill x -padx 4 -pady 2

    frame $f.blb -relief groove -bd 2
    label $f.blb.l -text "Perturbations" -anchor w
    set blobs [Scrolled_Listbox $f.blb.lb -height 5]
    pack  $f.blb.l $f.blb.lb -side top -fill x -expand true
    button $f.blb.add -text "Add" -command "PertDlg $blobs 0"
    button $f.blb.edit -text "Edit" -command "PertDlg $blobs 1"
    button $f.blb.delete -text "Delete" -command "PertDelete $blobs"
    pack  $f.blb.add $f.blb.edit $f.blb.delete -side left
    pack  $f.blb -side top -fill x -padx 4 -pady 2

    frame $f.bbar
    button $f.bbar.next -text "Next" -command { Execute; destroy .opt; \
						ResultDlg }
    button $f.bbar.back -text "Back" -command { destroy .opt; GeomDlg }
    button $f.bbar.exit -text "Exit" -command exit
    pack $f.bbar.exit $f.bbar.next $f.bbar.back -side right -fill x
    pack $f.bbar -side top -fill x -padx 4 -pady 2

    pack $f -side top
}

# ===========================================================================
# Result dialog

proc ResultDlg {} {
    global name els nds bnds inds
    global f
    set f [frame .opt]
    frame $f.data -relief groove -bd 2
    label $f.data.name -text "Wrote mesh to file $name" -anchor w
    label $f.data.els -text "$els elements" -anchor w
    label $f.data.nds -text "$nds nodes" -anchor w
    label $f.data.binds -text "($bnds boundary, $inds internal)" -anchor w
    pack $f.data.name $f.data.els $f.data.nds $f.data.binds -side top -fill x
    pack $f.data -side top -fill x -padx 4 -pady 2

    frame $f.bbar
    button $f.bbar.back -text "Restart" -command { destroy .opt; GeomDlg }
    button $f.bbar.exit -text "Exit" -command exit
    pack $f.bbar.exit $f.bbar.back -side right -fill x
    pack $f.bbar -side top -fill x -padx 4 -pady 2

    pack $f -side top
}

# ===========================================================================
# main

GeomDlg
