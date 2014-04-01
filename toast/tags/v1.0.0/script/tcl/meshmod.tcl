wm title . "meshmod"

set dimlist [list x y z]
#set meshmod_fname -
set meshmod_basisname []
set meshmod_dim 0
set meshmod_els 0
set meshmod_eltp -
set meshmod_nds 0
set meshmod_bnds 0
set meshmod_inds 0
set meshmod_bb -
set meshmod_sysbw -

set meshmod_meshsz 0
set meshmod_elsz 0

set meshmod_xscale 1
set meshmod_yscale 1
set meshmod_zscale 1

set meshmod_xshift 0
set meshmod_yshift 0
set meshmod_zshift 0

set meshmod_rotaxis x
set meshmod_rotangle 0

set meshmod_whichparam 0
set meshmod_regs [list Global]
set meshmod_p1 0.025
set meshmod_p2 2
set meshmod_p3 1.4

set meshmod_blobtype 0
set meshmod_blobcnt "0 0 0"
set meshmod_blobrad 1
set meshmod_blobreg 0

set meshmod_bndmode 0
set meshmod_optmode 0

set bitdim 400
set antialias 1

set rasterx 10
set rastery 10
set rasterz 10

# =============================================================================

proc ErrorBox {} {
    toplevel .ebox
    wm title .ebox "Error"
    label .ebox.l1 -text "TOAST library exception"
    message .ebox.msg -justify left -aspect 300 -anchor w -text \
	"An exception occurred in one of the TOAST libraries. See console for details. You can terminate the program or ignore the error. If you choose to continue, the program may become instable and abort or stop responding."
    pack .ebox.l1 .ebox.msg -side top -fill x -expand true

    # Button bar ============================================
    tixButtonBox .ebox.bbar -orientation horizontal -pady 2 -relief flat
    .ebox.bbar add term -text Terminate -underline 0 -width 8 \
	-command exit
    .ebox.bbar add cont -text Ignore -underline 0 -width 8 \
	-command "destroy .ebox"
    pack .ebox.bbar -side bottom -fill x

    tkwait window .ebox
    return $TCL_OK
}

# =============================================================================
# Status bar functions

#proc status:Echo {msg} {
#    .status.l config -text $msg
#    update
#    return TCL_OK
#}

#proc status:SetProgress {frac} {
#    .status.prog config -value $frac
#}


# =============================================================================
# ===========================  "Mesh" dialogs  ================================
# =============================================================================

# =============================================================================
# Create new circle mesh

proc NewCircleMesh {} {
    global meshmod_fname
    blt::busy hold .
    if {[GenCircleDlg] != 0} { ;# from gencircle.tcl
	set meshmod_fname "noname"
	UpdateInfo
	MeshAvailable 1
    }
    blt::busy release .
    return TCL_OK
}

# =============================================================================
# Create new rectangle mesh

proc NewRectangleMesh {} {
    global meshmod_fname
    blt::busy hold .
    if {[GenRectDlg] != 0} { ;# from genrect.tcl
	set meshmod_fname "noname"
	UpdateInfo
	MeshAvailable 1
    }
    blt::busy release .
    return TCL_OK
}

# =============================================================================
# Create new brick mesh

proc NewBrickMesh {} {
    global meshmod_fname
    blt::busy hold .
    if {[GenBrickDlg] != 0} { ;# from genbrick.tcl
	set meshmod_fname "noname"
	UpdateInfo
	MeshAvailable 1
    }
    blt::busy release .
    return TCL_OK
}

# =============================================================================
# "Open mesh" dialog

proc meshmod:OpenMesh {} {
    set dialog [tix filedialog tixExFileSelectDialog]
    $dialog config -command meshmod:DoOpenMesh -title "Open Mesh File"
    set dbox [$dialog subwidget fsbox]
    $dbox config -filetypes {
	{{*.msh *.opt *.mdf} {Mesh files}}
	{{*} {All files}}
    }
    $dbox subwidget types pick 0
    $dialog popup
}

proc MeshAvailable {what} {
    # enable/disable menu buttons
    set state [lindex [list disabled normal] $what]
    .mnu.mesh.menu entryconfigure "Save" -state $state
    .mnu.mesh.menu entryconfigure "Save As ..." -state $state
    .mnu.mesh.menu entryconfigure "Export surface nodes" -state $state
    .mnu.basis config -state $state
    .mnu.edit config -state $state
    .mnu.info config -state $state
}

proc meshmod:DoOpenMesh {file} {
    global meshmod_fname meshmod_basisfile
    set meshmod_fname $file
    set ccur [. cget -cursor]
    . config -cursor watch
    statusOut "Reading mesh ..."
    meshmod:Load
    UpdateInfo
    MeshAvailable 1

    # new mesh -> disable existing basis raster
    set meshmod_basisfile []
    foreach x {3 4 6} {
	.mnu.basis.menu entryconfigure $x -state disabled
    }

    statusReady
    . config -cursor $ccur
    meshmod:SysmatrixChanged
}

# =============================================================================
# "Save As" dialog

proc meshmod:SaveMeshAs {} {
    global meshmod_fname
    set dialog [tix filedialog tixExFileSelectDialog]
    $dialog config -command meshmod:DoSaveMeshAs -title "Save Mesh As"
    set dbox [$dialog subwidget fsbox]
    $dbox config -pattern $meshmod_fname
    $dialog popup
}

proc meshmod:DoSaveMeshAs {file} {
    global meshmod_fname
    set meshmod_fname $file
    meshmod:DoSaveMesh $file
    UpdateInfo
}

proc meshmod:DoSaveMesh {file} {
    puts "file=$file"
    if {$file == "noname"} {
	meshmod:SaveMeshAs
    } else {
	set ccur [. cget -cursor]
	. config -cursor watch
	statusOut "Saving mesh ..."
	meshmod:Save
	statusReady
	. config -cursor $ccur
    }
}

# =============================================================================
# ===========================  "Basis" dialogs  ===============================
# =============================================================================

# =============================================================================
# "New basis" dialog

proc meshmod:NewBasisDlg {} {
    global meshmod_dim dimlist
    global rasterx rastery rasterz

    toplevel .dlg
    wm title .dlg "Create basis raster"

    frame .dlg.0 -relief groove -bd 2
    label .dlg.0.l -text "Grid dimension" -anchor w
    pack  .dlg.0.l -side top -fill x
    for {set x 0} {$x < $meshmod_dim} {incr x} {
	set d [lindex $dimlist $x]
	tixControl .dlg.0.$x -label "$d-dim" -variable raster$d
	.dlg.0.$x subwidget entry config -width 5
	.dlg.0.$x subwidget label config -width 6
	pack .dlg.0.$x -side top -fill x -padx 4
    }
    checkbutton .dlg.0.cb -text "Square pixels" -anchor w
    if {$meshmod_dim == 3} {
	.dlg.0.cb config -text "Cubic voxels"
    }
    pack .dlg.0.cb -side top -fill x
    pack .dlg.0 -side top -fill x -padx 5 -pady 5

    # Button bar ============================================
    tixButtonBox .dlg.bbar -orientation horizontal -pady 2 -relief flat
    .dlg.bbar add apply -text OK -underline 0 -width 6 \
	    -command { meshmod:NewBasisApply .dlg }
    .dlg.bbar add done -text Cancel -underline 0 -width 6 \
	    -command {destroy .dlg}
    pack .dlg.bbar -side bottom -fill x
    return TCL_OK
}

proc meshmod:NewBasisApply {w} {
    global meshmod_dim rasterx rastery rasterz
    set ccur [$w cget -cursor]
    $w config -cursor watch
    statusOut "Generating basis raster ..."
    for {set x 0} {$x < $meshmod_dim} {incr x} {
	$w.0.$x update
    }
    meshmod:NewBasis $rasterx $rastery $rasterz
    statusReady
    $w config -cursor $ccur
    destroy $w
}

# ============================================================================
# "Open basis" dialog

proc meshmod:OpenBasis {} {
    set dialog [tix filedialog tixExFileSelectDialog]
    $dialog config -command meshmod:OpenBasisApply \
	    -title "Open Basis Raster File"
    set dbox [$dialog subwidget fsbox]
    $dbox config -filetypes {
	{{*.br} {Basis raster files}}
	{{*} {All files}}
    }
    $dbox subwidget types pick 0
    $dialog popup
    return TCL_OK
}

# =============================================================================
# "OpenBasisApply" function

proc meshmod:OpenBasisApply {file} {
    global meshmod_basisname
    set meshmod_basisname $file
    set ccur [. cget -cursor]
    . config -cursor watch
    statusOut "Reading basis ..."
    meshmod:LoadBasis $file

    # enable basis menu options
    foreach x {3 4 6} {
	.mnu.basis.menu entryconfigure $x -state normal
    }

    statusReady
    . config -cursor $ccur
    return TCL_OK
}

# =============================================================================
# "Save basis as" dialog

proc meshmod:SaveBasisAs {} {
    return TCL_OK
}

# =============================================================================
# "Save basis" function

proc meshmod:SaveBasis {file} {
    return TCL_OK
}

# =============================================================================
# "Filter basis" dialog

proc meshmod:FilterBasis {} {
    return TCL_OK
}

# =============================================================================
# "Transform" dialog

proc meshmod:Transform {} {
    global meshmod_dim meshmod_xscale meshmod_yscale meshmod_zscale \
	    meshmod_xshift meshmod_yshift meshmod_zshift \
	    meshmod_rotaxis meshmod_rotangle

    toplevel .dlg
    wm title .dlg "meshmod: Transform"

    set name [tixOptionName .dlg]
    option add *$name*TixControl*entry.width 10
    option add *$name*TixControl*label.width 10
    option add *$name*TixControl*label.anchor e
    option add *$name*TixControl*label.disabledForeground gray

    tixNoteBook .dlg.nb -ipadx 1 -ipady 1
    .dlg.nb subwidget nbframe config -tabpady 1
    .dlg.nb add scale -label "Scale" -underline 0
    .dlg.nb add translate -label "Translate" -underline 0
    .dlg.nb add rotate -label "Rotate" -underline 0
    pack .dlg.nb -expand true -fill both -padx 5 -pady 5 -side top

    # "Scale" tab ===========================================
    set f [.dlg.nb subwidget scale]
    frame $f.f
    pack  $f.f -side left -padx 2 -pady 2 -fill both -expand true
    tixControl $f.f.x -variable meshmod_xscale -min 0 -label "x-scale"
    tixControl $f.f.y -variable meshmod_yscale -min 0 -label "y-scale"
    tixControl $f.f.z -variable meshmod_zscale -min 0 -label "z-scale"
    if {$meshmod_dim < 3} {
	$f.f.z config -state disabled
    }
    pack  $f.f.x $f.f.y $f.f.z -side top -padx 20 -pady 2

    # "Translate" tab =======================================
    set f [.dlg.nb subwidget translate]
    frame $f.f
    pack  $f.f -side left -padx 2 -pady 2 -fill both -expand true
    tixControl $f.f.x -variable meshmod_xshift -label "x-ofs"
    tixControl $f.f.y -variable meshmod_yshift -label "y-ofs"
    tixControl $f.f.z -variable meshmod_zshift -label "z-ofs"
    if {$meshmod_dim < 3} {
	$f.f.z config -state disabled
    }
    pack  $f.f.x $f.f.y $f.f.z -side top -padx 20 -pady 2

    # "Rotate" tab ==========================================
    set f [.dlg.nb subwidget rotate]
    frame $f.f
    pack $f.f -side left -padx 2 -pady 2 -fill both -expand true
    tixSelect $f.f.axis -allowzero false -radio true -orient horizontal \
	    -label "Rotate around" -labelside top
    $f.f.axis add x -text "x-axis"
    $f.f.axis add y -text "y-axis"
    $f.f.axis add z -text "z-axis"
    $f.f.axis configure -variable meshmod_rotaxis
    tixControl $f.f.angle -variable meshmod_rotangle -min -360 -max 360 \
	    -label "angle \[deg\]"
    if {$meshmod_dim < 3} {
	$f.f.axis config -value z
	$f.f.axis subwidget x config -state disabled
	$f.f.axis subwidget y config -state disabled
    }
    pack  $f.f.axis $f.f.angle -side top -padx 20 -pady 2

    # Button bar ============================================
    tixButtonBox .dlg.bbar -orientation horizontal -pady 2 -relief flat
    .dlg.bbar add apply -text Apply -underline 0 -width 6 \
	    -command "meshmod:TransformApply .dlg"
    .dlg.bbar add done -text Done -underline 0 -width 6 \
	    -command "destroy .dlg"
    pack .dlg.bbar -side bottom -fill x
    return TCL_OK
}

proc meshmod:TransformApply {w} {
    set ccur [$w cget -cursor]
    $w config -cursor watch
    set curr_page [$w.nb raised]
    if {$curr_page == "scale"} {
	meshmod:ScaleApply $w
    } elseif {$curr_page == "translate"} {
	meshmod:TranslateApply $w
    } else {
	meshmod:RotateApply $w
    }
    $w config -cursor $ccur
    return TCL_OK
}

proc meshmod:ScaleApply {w} {
    statusOut "Scaling mesh ..."
    set f [$w.nb subwidget scale]
    foreach d {x y z} {
	$f.f.$d update
    }
    meshmod:Scale
    UpdateInfo
    statusReady
    return TCL_OK
}

proc meshmod:TranslateApply {w} {
    statusOut "Translating mesh ..."
    set f [$w.nb subwidget translate]
    foreach d {x y z} {
	$f.f.$d update
    }
    meshmod:Translate
    UpdateInfo
    statusReady
    return TCL_OK
}

proc meshmod:RotateApply {w} {
    statusOut "Rotating mesh ..."
    set f [$w.nb subwidget rotate]
    $f.f.angle update
    meshmod:Rotate
    UpdateInfo
    statusReady
    return TCL_OK
}

# =============================================================================
# "Parameter" dialog

proc meshmod:Parameter {} {
    global meshmod_whichparam meshmod_regs
    global meshmod_p1 meshmod_p2 meshmod_p3

    global do_mua do_mus do_n
    set do_mua 1
    set do_mus 1
    set do_n 1

    toplevel .dlg
    wm title .dlg "meshmod: Parameter"

    set name [tixOptionName .dlg]
    option add *$name*TixControl*entry.width 12
    option add *$name*TixControl*label.width 8
    option add *$name*TixControl*label.anchor e

    frame .dlg.f -relief groove -bd 2
    tixComboBox .dlg.f.reg -label Region -dropdown true -editable false \
	    -listwidth 30 -options { entry.width 12 }
    pack .dlg.f.reg -side top -fill x -expand true
    meshmod:ScanRegion
    foreach x $meshmod_regs {
	.dlg.f.reg insert end $x
    }
    tixSetSilent .dlg.f.reg [lindex $meshmod_regs 0]
    pack .dlg.f -side top -padx 5 -pady 5 -fill x \
	-expand true

    frame .dlg.p -relief groove -bd 2
    label .dlg.p.l -text "Apply parameter" -anchor w
    pack .dlg.p.l -side top -fill x

    frame .dlg.p.f1
    entry .dlg.p.f1.e -width 12 -textvariable meshmod_p1
    checkbutton .dlg.p.f1.c -text mua -width 5 -anchor w -variable do_mua \
	-command {toggle_entry $do_mua .dlg.p.f1.e}
    pack .dlg.p.f1.c .dlg.p.f1.e -side left
    frame .dlg.p.f2
    entry .dlg.p.f2.e -width 12 -textvariable meshmod_p2
    checkbutton .dlg.p.f2.c -text mus -width 5 -anchor w -variable do_mus \
	-command {toggle_entry $do_mus .dlg.p.f2.e}
    pack .dlg.p.f2.c .dlg.p.f2.e -side left
    frame .dlg.p.f3
    entry .dlg.p.f3.e -width 12 -textvariable meshmod_p3
    checkbutton .dlg.p.f3.c -text n -width 5 -anchor w -variable do_n \
	-command {toggle_entry $do_n .dlg.p.f3.e}
    pack .dlg.p.f3.c .dlg.p.f3.e -side left
    pack .dlg.p.f1 .dlg.p.f2 .dlg.p.f3 -side top -fill x

    pack .dlg.p -side top -padx 5 -pady 5 -fill x -expand true

    # Button bar ============================================
    tixButtonBox .dlg.bbar -orientation horizontal -pady 2 -relief flat
    .dlg.bbar add apply -text Apply -underline 0 -width 6 \
	-command "meshmod:ParameterApply .dlg"
    .dlg.bbar add done -text Done -underline 0 -width 6 \
	-command "destroy .dlg"
    pack .dlg.bbar -side bottom -fill x
    return TCL_OK
}

proc toggle_entry {s e} {
    if {$s == 1} {
	$e config -state normal -foreground black
    } else {
	$e config -state disabled -foreground gray
    }
}

proc toggle_shape {w v} {
    global meshmod_dim meshmod_blobrad
    if {$v == "rad"} {
	$w config -label "Radius"
	set meshmod_blobrad 1
    } elseif {$v == "ell"} {
	if {$meshmod_dim == 2} {
	    $w config -label "SM\[x y\] phi"
	    set meshmod_blobrad "1 1 0"
	} else {
	    $w config -label "SM\[x y z\] phi theta"
	    set meshmod_blobrad "1 1 1 0 0"
	}
    } else {
	if {$meshmod_dim == 2} {
	    $w config -label "Size \[x y\]"
	    set meshmod_blobrad "1 1"
	} else {
	    $w config -label "Size \[x y z\]"
	    set meshmod_blobrad "1 1 1"
	}
    }
}

proc meshmod:ParameterApply {w} {
    global meshmod_p1 meshmod_p2 meshmod_p3
    statusOut "Modifying parameters ..."
    foreach d {1 2 3} {
	if {[$w.p.f$d.e cget -state] == "normal"} {
	    set v meshmod_p$d
	    set p$d [subst $$v]
	} else {
	    set p$d -1
	}
    }
    set v [$w.f.reg cget -value]
    meshmod:RegionParam $v $p1 $p2 $p3
    UpdateInfo
    statusReady
    return TCL_OK
}

# =============================================================================
# "Insert blob" dialog

proc meshmod:AddBlob {} {
    global meshmod_dim meshmod_blobtype meshmod_blobcnt meshmod_blobrad

    global do_mua do_mus do_n do_reg
    set do_mua 1
    set do_mus 1
    set do_n 1
    set do_reg 1

    toplevel .dlg
    wm title .dlg "meshmod: Insert Blob"

    set name [tixOptionName .dlg]
    option add *$name*TixLabelEntry*entry.width 20
    option add *$name*TixLabelEntry*label.width 12
    option add *$name*TixLabelEntry*label.anchor e

    frame .dlg.cnt -relief groove -bd 2
    label .dlg.cnt.l -text "Blob geometry" -anchor w
    pack  .dlg.cnt.l -side top -fill x
    frame .dlg.cnt.shp
    label .dlg.cnt.shp.l -text Shape -width 12 -anchor e
    radiobutton .dlg.cnt.shp.r1 -variable meshmod_blobtype -value 0 -anchor w \
	    -command { toggle_shape .dlg.cnt.r rad }
    radiobutton .dlg.cnt.shp.r2 -variable meshmod_blobtype -value 1 -anchor w \
	    -command { toggle_shape .dlg.cnt.r rect }
    radiobutton .dlg.cnt.shp.r3 -variable meshmod_blobtype -value 2 -anchor w \
	    -command { toggle_shape .dlg.cnt.r ell }
    if {$meshmod_dim == 2} {
	.dlg.cnt.shp.r1 config -text "circle"
	.dlg.cnt.shp.r2 config -text "rectangle"
	.dlg.cnt.shp.r3 config -text "ellipse"
    } else {
	.dlg.cnt.shp.r1 config -text "sphere"
	.dlg.cnt.shp.r2 config -text "brick"
	.dlg.cnt.shp.r3 config -text "ellipsoid"
    }
    pack .dlg.cnt.shp.l .dlg.cnt.shp.r1 .dlg.cnt.shp.r2 .dlg.cnt.shp.r3 \
	-side left -fill x
    pack .dlg.cnt.shp -side top -fill x

    tixLabelEntry .dlg.cnt.c
    if {$meshmod_dim == 2} {
	.dlg.cnt.c config -label "Centre \[x y\]"
	set meshmod_blobcnt "0 0"
    } else {
	.dlg.cnt.c config -label "Centre \[x y z\]"
	set meshmod_blobcnt "0 0 0"
    }
    .dlg.cnt.c subwidget entry config -textvariable meshmod_blobcnt
    tixLabelEntry .dlg.cnt.r -label "Radius"
    .dlg.cnt.r subwidget entry config -textvariable meshmod_blobrad
    pack  .dlg.cnt.c .dlg.cnt.r -side top -fill x -expand true
    pack .dlg.cnt -side top -fill x -padx 5 -pady 5 -ipadx 2 -ipady 2

    frame .dlg.p -relief groove -bd 2
    label .dlg.p.l -text "Blob node attributes" -anchor w
    pack .dlg.p.l -side top -fill x

    frame .dlg.p.f1
    entry .dlg.p.f1.e -width 12 -textvariable meshmod_p1
    checkbutton .dlg.p.f1.c -text mua -width 6 -anchor w -variable do_mua \
	    -command {toggle_entry $do_mua .dlg.p.f1.e}
    pack .dlg.p.f1.c .dlg.p.f1.e -side left
    frame .dlg.p.f2
    entry .dlg.p.f2.e -width 12 -textvariable meshmod_p2
    checkbutton .dlg.p.f2.c -text mus -width 6 -anchor w -variable do_mus \
	    -command {toggle_entry $do_mus .dlg.p.f2.e}
    pack .dlg.p.f2.c .dlg.p.f2.e -side left
    frame .dlg.p.f3
    entry .dlg.p.f3.e -width 12 -textvariable meshmod_p3
    checkbutton .dlg.p.f3.c -text n -width 6 -anchor w -variable do_n \
	    -command {toggle_entry $do_n .dlg.p.f3.e}
    pack .dlg.p.f3.c .dlg.p.f3.e -side left
    pack .dlg.p.f1 .dlg.p.f2 .dlg.p.f3 -side top -fill x -padx 20

    frame .dlg.p.f4
    entry .dlg.p.f4.e -width 12 -textvariable meshmod_blobreg
    checkbutton .dlg.p.f4.c -text Region -width 6 -anchor w -variable do_reg \
	    -command {toggle_entry $do_reg .dlg.p.f4.e}
    pack .dlg.p.f4.c .dlg.p.f4.e -side left
    pack .dlg.p.f4 -side top -fill x -padx 20 -pady 6

    pack .dlg.p -side top -padx 5 -pady 5 -fill x -expand true

    # Button bar ============================================
    tixButtonBox .dlg.bbar -orientation horizontal -pady 2 -relief flat
    .dlg.bbar add apply -text Apply -underline 0 -width 6 \
	    -command "meshmod:AddBlobApply .dlg"
    .dlg.bbar add done -text Done -underline 0 -width 6 \
	    -command "destroy .dlg"
    pack .dlg.bbar -side bottom -fill x
    return TCL_OK
}

proc meshmod:AddBlobApply {w} {
    global meshmod_p1 meshmod_p2 meshmod_p3 meshmod_blobreg
    global meshmod_blobtype meshmod_blobcnt meshmod_blobrad
    statusOut "Inserting blob ..."

    foreach d {1 2 3} {
	if {[$w.p.f$d.e cget -state] == "normal"} {
	    set v meshmod_p$d
	    set p$d [subst $$v]
	} else {
	    set p$d -1
	}
    }
    if {[$w.p.f4.e cget -state] == "normal"} {
	set reg $meshmod_blobreg
    } else {
	set reg -1
    }
    meshmod:InsertBlob $meshmod_blobtype $meshmod_blobcnt $meshmod_blobrad \
	    $p1 $p2 $p3 $reg
    UpdateInfo
    statusReady
    return TCL_OK
}

# =============================================================================
# "Optimise" dialog

proc meshmod:OptimiseDlg {} {
    global meshmod_bndmode meshmod_optmode

    toplevel .dlg
    wm title .dlg "Optimise mesh"

    frame .dlg.f1 -relief groove -bd 2
    label .dlg.f1.l -text "Boundary nodes" -anchor w
    radiobutton .dlg.f1.r1 -text "Include in optimisation" \
	    -variable meshmod_bndmode -value 0 -anchor w
    radiobutton .dlg.f1.r2 -text "Sort to end of list" \
	    -variable meshmod_bndmode -value 1 -anchor w
    pack .dlg.f1.l .dlg.f1.r1 .dlg.f1.r2 -side top -fill x
    pack .dlg.f1 -side top -padx 5 -pady 5 -fill x -expand true

    frame .dlg.f2 -relief groove -bd 2
    label .dlg.f2.l -text "Optimisation method" -anchor w
    radiobutton .dlg.f2.r1 -text "Minimum bandwidth" \
	    -variable meshmod_optmode -value 0 -anchor w
    radiobutton .dlg.f2.r2 -text "Minimum degree" \
	    -variable meshmod_optmode -value 1 -anchor w
    radiobutton .dlg.f2.r3 -text "Tinney 2" -variable meshmod_optmode \
	    -value 2 -anchor w
    pack .dlg.f2.l .dlg.f2.r1 .dlg.f2.r2 .dlg.f2.r3 -side top -fill x
    pack .dlg.f2 -side top -padx 5 -pady 5 -expand true

    # Button bar ============================================
    tixButtonBox .dlg.bbar -orientation horizontal -pady 2 -relief flat
    .dlg.bbar add apply -text Apply -underline 0 -width 6 \
	    -command "meshmod:OptimiseApply .dlg"
    .dlg.bbar add done -text Done -underline 0 -width 6 \
	    -command "destroy .dlg"
    pack .dlg.bbar -side bottom -fill x
    return TCL_OK
}

proc meshmod:OptimiseApply {w} {
    global meshmod_bndmode meshmod_optmode
    set ccur [$w cget -cursor]
    $w config -cursor watch
    statusOut "Optimising mesh ..."
    meshmod:Optimise $meshmod_bndmode $meshmod_optmode 
    UpdateInfo
    statusReady
    $w config -cursor $ccur
    meshmod:SysmatrixChanged
    return TCL_OK
}

# =============================================================================
# "Mark boundary nodes" command

proc meshmod:MarkBoundary {} {
    statusOut "Scanning for boundary nodes ..."
    MarkMeshBoundary
    UpdateInfo
    statusReady
    return TCL_OK
}

# =============================================================================
# "System matrix info" dialog

proc meshmod:SysMatrix {} {
    global bitdim antialias

    set w [toplevel .smdlg]
    wm title $w "System Matrix Info"

    set labels1 [list "Pre" "Post"]
    set labels2 [list "Allocated entries" "Fill fraction"]

    for {set x 0} {$x < 2} {incr x} {
	frame $w.$x -relief groove -bd 2
	label $w.$x.l -text "[lindex $labels1 $x]-factorisation statistics" \
		-anchor w
	pack  $w.$x.l -side top -fill x
	for {set y 0} {$y < 2} {incr y} {
	    frame $w.$x.$y
	    label $w.$x.$y.l -text "[lindex $labels2 $y]:  " -width 16 \
		    -anchor e
	    label $w.$x.$y.v -text "-" -width 9 -anchor w -foreground blue
	    pack  $w.$x.$y.l $w.$x.$y.v -side left
	    pack  $w.$x.$y -side top -fill x
	}
	tixButtonBox $w.$x.bbar -orientation horizontal -pady 2 -relief flat
	$w.$x.bbar add disp -text "Display graph" \
		-command [list meshmod:DispSysmatrix $x]
	$w.$x.bbar add file -text "Write to file" \
		-command [list meshmod:SaveSysmatrix $x]
	pack $w.$x.bbar -side bottom -fill x
	pack $w.$x -side top -fill x -expand true -padx 5 -pady 5
    }

    frame $w.2 -relief groove -bd 2
    label $w.2.l -text "Graph options" -anchor w
    set c [tixControl $w.2.c -label "Max bitmap dimension" -min 4 \
	    -value $bitdim -variable bitdim]
    $c subwidget entry config -width 5
    checkbutton $w.2.cb -text "Anti-aliasing" -variable antialias -anchor w
    pack  $w.2.l $w.2.c $w.2.cb -side top -fill x
    pack  $w.2 -side top -fill x -expand true -padx 5 -pady 5

    # Button bar ============================================
    tixButtonBox $w.bbar -orientation horizontal -pady 2 -relief flat
    $w.bbar add done -text Done -underline 0 -width 6 \
	    -command "destroy $w"
    pack $w.bbar -side bottom -fill x
    meshmod:UpdateSysmatrixData $w
    return TCL_OK
}

proc meshmod:UpdateSysmatrixData {{w .}} {
    set ccur [$w cget -cursor]
    $w config -cursor watch
    statusOut "Generating system matrix ..."
    meshmod:UpdateSysmatrix
    if {$w != "."} {
	meshmod:GetSysmatrixData $w.0.0.v $w.0.1.v $w.1.0.v $w.1.1.v
    }
    statusReady
    $w config -cursor $ccur
    return TCL_OK
}

proc meshmod:DispSysmatrix {x} {
    global which_sysmatrix
    set which_sysmatrix $x
    set which [lindex [list sysmat factor] $x]
    meshmod:DoSaveSysmatrix "/tmp/$which.pgm" ;# temporary image file

    image create photo $which -file "/tmp/$which.pgm"
    if {[winfo exists .$which] == 0} {
	toplevel .$which
	label .$which.l -image $which
	pack  .$which.l
    } else {
	.$which.l config -image $which
    }
    #image delete $which
    return TCL_OK
}

proc meshmod:SaveSysmatrix {x} {
    global which_sysmatrix
    set which_sysmatrix $x
    set dialog [tix filedialog tixExFileSelectDialog]
    $dialog config -command {meshmod:DoSaveSysmatrix} \
	    -title "Save System Matrix"
    #set dbox [$dialog subwidget fsbox]
    #$dbox config -pattern $meshmod_fname
    $dialog popup
    return TCL_OK
}

proc meshmod:DoSaveSysmatrix {file} {
    global bitdim antialias which_sysmatrix
    if {[winfo exists .smdlg] == 1} {
	set w .smdlg
	$w.2.c update ;# bitmap dimension
    } else {
	set w .
    }
    set ccur [$w cget -cursor]
    $w config -cursor watch
    statusOut "Saving system matrix ..."
    if {$which_sysmatrix == 0} {
	meshmod:SaveSysmat $file $bitdim $antialias
    } else {
	meshmod:SaveFactor $file $bitdim $antialias
    }
    statusReady
    $w config -cursor $ccur
}

# =============================================================================
# "Element size" histogram

proc meshmod:SizeHisto {} {
    set w [toplevel .sizehisto]
    wm title $w "Element size histogram"
    global xgraph ygraph

    # Button bar ============================================
    tixButtonBox $w.bbar -orientation horizontal -pady 2 -relief flat
    $w.bbar add write -text "PS Output" -underline 0 \
	    -command "SizeHistoOutput $w.plot"
    $w.bbar add done -text Done -underline 0 -width 6 \
	    -command "destroy $w"
    pack $w.bbar -side bottom -fill x

    blt::barchart $w.plot
    $w.plot legend config -hide yes
    $w.plot axis config x -title "Element size"
    $w.plot axis config y -title "Element count"
    pack $w.plot -side top -fill both -expand true

    blt::vector create xgraph(100)
    blt::vector create ygraph(100)
    meshmod:GetSizeHisto 100
    set barwidth [expr ($xgraph(99)-$xgraph(0))*0.01]
    $w.plot element create elsize -xdata xgraph -ydata ygraph \
	-barwidth $barwidth -relief flat

    return TCL_OK
}

proc SizeHistoOutput {plot} {
    set dialog [tix filedialog tixExFileSelectDialog]
    $dialog config -command {DoSizeHistoOutput} -title "Output PS"
    set dbox [$dialog subwidget fsbox]
    $dbox config -filetypes {{{*.ps *.eps} {Postscript files}}}
    $dbox subwidget types pick 0
    $dialog popup
    return TCL_OK
}

proc DoSizeHistoOutput {file} {
    .sizehisto.plot postscript output $file -maxpect yes -decorations no
    return TCL_OK
}

# =============================================================================
# "Neighbour count" diagram

proc meshmod:NeighbourCount {} {
    set w .nbcount
    if [winfo exists $w] {
	raise $w
    } else {
	toplevel $w
	wm title $w "Neighbour count"

	# Button bar ============================================
	tixButtonBox $w.bbar -orientation horizontal -pady 2 -relief flat
	$w.bbar add write -text "PS Output" -underline 0 \
	    -command "NeighbourCountOutput $w.plot"
	$w.bbar add done -text Done -underline 0 -width 6 \
	    -command "destroy $w"
	pack $w.bbar -side bottom -fill x

	blt::barchart $w.plot
	$w.plot legend config -hide yes
	$w.plot axis config x -title "\# neighbours"
	$w.plot axis config y -title "\# nodes"
	$w.plot element create elsize -relief flat
	pack $w.plot -side top -fill both -expand true
    }
    blt::busy hold .
    statusOut "Computing neighbours ..."
    set step 1
    meshmod:GetNbCount xnbvec ynbvec step
    $w.plot element config elsize -xdata xnbvec -ydata ynbvec \
	-barwidth $step
    $w.plot axis config x -stepsize $step -subdivisions 1
    blt::busy release .
    statusReady

    return TCL_OK
}

proc NeighbourCountOutput {plot} {
    set dialog [tix filedialog tixExFileSelectDialog]
    $dialog config -command {DoNeighbourCountOutput} -title "Output PS"
    set dbox [$dialog subwidget fsbox]
    $dbox config -filetypes {{{*.ps *.eps} {Postscript files}}}
    $dbox subwidget types pick 0
    $dialog popup
    return TCL_OK
}

proc DoNeighbourCountOutput {file} {
    .nbcount.plot postscript output $file -maxpect yes -decorations no
    return TCL_OK
}

# =============================================================================
# Updates displays whenever system matrix structure has changed

proc meshmod:SysmatrixChanged {} {
    set updated 0
    if {[winfo exists .smdlg] == 1} {
	meshmod:UpdateSysmatrixData .smdlg
	set updated 1
    }
    if {[winfo exists .sysmat] == 1} {
	if {$updated == 0} {
	    meshmod:UpdateSysmatrixData
	    set updated 1
	}
	meshmod:DispSysmatrix 0
    }
    if {[winfo exists .factor] == 1} {
	if {$updated == 0} {
	    meshmod:UpdateSysmatrixData
	}
	meshmod:DispSysmatrix 1
    }
}

# =============================================================================
# Update info display

proc UpdateInfo {} {
    global meshmod_fname meshmod_dim meshmod_els meshmod_eltp meshmod_nds \
	   meshmod_bnds meshmod_inds meshmod_bb meshmod_meshsz meshmod_elsz \
	   meshmod_sysbw 
    
    UpdateParams
    if {$meshmod_fname == "noname"} {
	.info.0.v configure -text "\[None\]"
    } else {
	.info.0.v configure -text "$meshmod_fname"
    }
    .info.1.v configure -text "$meshmod_dim[]D"
    .info.2.v configure -text "$meshmod_els (type $meshmod_eltp)"
    .info.3.v configure -text \
	    "$meshmod_nds ($meshmod_bnds bnd, $meshmod_inds intl)"
    .info.4.v configure -text "$meshmod_bb"
    .info.5.v configure -text "$meshmod_meshsz"
    .info.6.v configure -text "$meshmod_elsz"
    .info.7.v configure -text "$meshmod_sysbw"

    if {[winfo exists .sizehisto] == 1} {
	blt::vector create xgraph(100)
	blt::vector create ygraph(100)
	meshmod:GetSizeHisto 100
	set barwidth [expr ($xgraph(99)-$xgraph(0))*0.01]
	.sizehisto.plot element config elsize -barwidth $barwidth
    }
    if {[winfo exists .nbcount] == 1} {
	meshmod:NeighbourCount
    }
}

# =============================================================================
# toplevel menu

lappend auto_path $env(TOAST_SCRIPT_PATH)
toastscheme

# ==== menu bar ====
MenuSetup .mnu
menubutton .mnu.mesh -text Mesh -underline 0 -menu .mnu.mesh.menu
menu .mnu.mesh.menu
.mnu.mesh.menu add cascade -label "New" -underline 0 -menu .mnu.mesh.menu.new
menu .mnu.mesh.menu.new
.mnu.mesh.menu.new add command -label "Circle ..." -underline 0 \
        -command NewCircleMesh
.mnu.mesh.menu.new add command -label "Rectangle ..." -underline 0 \
       -command NewRectangleMesh
.mnu.mesh.menu.new add command -label "Brick ..." -underline 0 \
	-command NewBrickMesh
.mnu.mesh.menu add command -label "Open ..." -underline 0 \
	-command meshmod:OpenMesh
.mnu.mesh.menu add command -label "Save" -underline 0 \
	-command {meshmod:DoSaveMesh $meshmod_fname} -state disabled
.mnu.mesh.menu add command -label "Save As ..." -underline 5 \
	-command meshmod:SaveMeshAs -state disabled
.mnu.mesh.menu add separator
.mnu.mesh.menu add command -label "Export surface nodes" \
        -command meshmod:ExportSurfNodes -state disabled
.mnu.mesh.menu add separator
.mnu.mesh.menu add command -label "Exit" -underline 1 -command exit

menubutton .mnu.basis -text Basis -underline 0 -menu .mnu.basis.menu \
	-state disabled
menu .mnu.basis.menu
.mnu.basis.menu add command -label "New ..." -underline 0 \
	-command meshmod:NewBasisDlg
.mnu.basis.menu add command -label "Open Basis Raster File ..." -underline 0 \
	-command meshmod:OpenBasis
.mnu.basis.menu add command -label "Save Basis Raster File" -underline 0 \
	-command {meshmod:SaveBasis $meshmod_basisname} -state disabled
.mnu.basis.menu add command -label "Save As ..." -underline 5 \
	-command meshmod:SaveBasisAs -state disabled
.mnu.basis.menu add separator
.mnu.basis.menu add command -label "Filter ..." -underline 0 \
	-command meshmod:FilterBasis -state disabled

menubutton .mnu.edit -text Edit -underline 0 -menu .mnu.edit.menu \
	-state disabled
menu .mnu.edit.menu
.mnu.edit.menu add command -label "Transform ..." -command meshmod:Transform
.mnu.edit.menu add command -label "Parameter ..." -command meshmod:Parameter
.mnu.edit.menu add command -label "Insert blob ..." -command meshmod:AddBlob
.mnu.edit.menu add separator
.mnu.edit.menu add command -label "Optimise ..." -command meshmod:OptimiseDlg
.mnu.edit.menu add command -label "Mark boundary nodes" \
        -command meshmod:MarkBoundary

menubutton .mnu.info -text Info -underline 0 -menu .mnu.info.menu \
	-state disabled
menu .mnu.info.menu
.mnu.info.menu add command -label "Element size distribution ..." \
	-command meshmod:SizeHisto
.mnu.info.menu add command -label "System matrix ..." \
	-command meshmod:SysMatrix
.mnu.info.menu add command -label "Neighbour count ..." \
        -command meshmod:NeighbourCount
.mnu.info.menu add command -label "Test mesh integrity" \
        -command meshmod:TestMesh
pack .mnu.mesh .mnu.basis .mnu.edit .mnu.info -side left -fill x

# ==== status bar ====
statusSetup

# ==== client area ====

frame .info -relief flat -bd 2

set infolabel [list "Mesh file" "Dimension" "Elements" "Nodes" "BBox" \
		  "Mesh size" "Element size" "Sysmat. SBW"]

for {set x 0} {$x < 8} {incr x} {
    set f [frame .info.$x]
    label $f.l -text "[lindex $infolabel $x]  " -width 14 -anchor e
    label $f.v -text "-" -width 35 -anchor w -foreground blue
    pack  $f.l -side left
    pack  $f.v -side left -fill x -expand true
    pack  $f -side top -fill x
}
pack .info -fill both -expand true

# load mesh, if a name is present
if {$meshmod_fname != "-"} {
    meshmod:DoOpenMesh $meshmod_fname
}
