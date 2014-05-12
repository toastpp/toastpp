wm title . "itdata"

set curset -1
set newset 0
set ndata 16
set ng 0
set graphlist {}
set dtype "unknown"
set linkdset 0
set upddset 0
set normalise 0
set enorm 0
set nmom 1
set laps 0.001
set transform ""
set enormstr [list "Chi-squared" "Unnormalised" "Normalised" \
	"Calculate" "Const. received" "Read"]

# =============================================================================

proc listInsert {idx setno type} {
    global nsets curset
    set lb [[.pw subwidget pane1].dlist subwidget listbox]
    $lb insert $idx "Set $setno \[$type\]"
    $lb activate $idx
    $lb selection clear 0 end
    $lb selection set $idx
    incr nsets
    set curset $setno
    selectData 0
    return TCL_OK
}

proc listDelete {idx} {
    global nsets
    set lb [[.pw subwidget pane1].dlist subwidget listbox]
    $lb delete $idx
    incr nsets -1
    return TCL_OK
}

# =============================================================================
# "Load definition file" dlg

proc dlgLoadDef {} {
    set dialog [tix filedialog tixExFileSelectDialog]
    $dialog config -command "doLoadDef" -title "Load definition"
    set dbox [$dialog subwidget fsbox]
    $dbox config -filetypes {{{*.def} {Definition files}}}
    $dbox subwidget types pick 0
    $dialog popup
    return TCL_OK
}

proc doLoadDef {file} {
    global curset nsets
    if {$nsets > 0} {
	if {[AskBox "This will replace all currently defined data sets. Continue?"] == 0} {
	    return TCL_OK
	}
    }
    femdataLoadEnv $file
    set dl [[.pw subwidget pane1].dlist subwidget listbox]
    $dl delete 0 end
    for {set i 0} {$i < $nsets} {incr i} {
	$dl insert end "Set $i \[A\]"
    }
    if {$nsets > 0} {
	$dl activate 0
	$dl selection set 0
	set curset -1
	selectData $dl
    }
    return TCL_OK
}

# =============================================================================
# "Insert data set" dialog

proc dlgInsertSet {} {

# Returns user's choice for the set number to be inserted.
# If the set already exists, the user wants to replace it.
# Return value -1 indicates cancel

    global iset setno

    set w [toplevel .dlgInsertSet]
    wm title $w "Insert new data set"

    blt::busy hold .
    focus $w
    bind $w <Destroy> { blt::busy release . }

    set iset 0
    while {[femdataDsetIndex $iset] >= 0} {incr iset}
    
    set c [tixControl $w.setno -label "Insert set no." -allowempty no \
	       -integer yes -min 0 -variable iset -value $iset]
    $c subwidget entry config -width 10
    pack $c -side top -padx 5 -pady 5

    # Button bar ============================================
    tixButtonBox $w.bbar -orientation horizontal -pady 2 -relief flat
    $w.bbar add apply -text Insert -underline 0 -width 6 \
	    -command [list applyInsertSet $w]
    $w.bbar add cancel -text Cancel -underline 0 -width 6 \
	    -command { set setno -1 }
    pack $w.bbar -side bottom -fill x

    tkwait variable setno
    destroy $w
    return $setno
}

proc dlgConfirmReplace {s} {
    set w [toplevel .dlgConfirmReplace]
    wm title $w "Warning"
    blt::busy hold .dlgInsertSet
    bind $w <Destroy> { blt::busy release .dlgInsertSet }

    global res
    message $w.msg -justify left -aspect 400 -anchor w -text \
	    "Data set $s exists. Replace?"
    pack $w.msg -side top -fill x -expand true -padx 10 -pady 5
    
    # Button bar ============================================
    tixButtonBox $w.bbar -orientation horizontal -pady 2 -relief flat
    $w.bbar add yesb -text Yes -underline 0 -width 6 -command { set res 1 }
    $w.bbar add nob -text No -underline 0 -width 6 -command { set res 0 }
    pack $w.bbar -side bottom -fill x

    tkwait variable res
    destroy $w
    return $res
}

proc applyInsertSet {w} {
    global setno
    $w.setno update
    set s [$w.setno cget -value]
    set idx [femdataDsetIndex $s]
    if {$idx < 0 || [dlgConfirmReplace $s] == 1} {
	set setno $s
    }
    return TCL_OK
}

# =============================================================================
# "Mesh and QM" dialog

proc dlgMesh {} {
    global meshfile qmfile bc solver

    set w [toplevel .dlgMesh]
    wm title $w "Load mesh and QM"

    blt::busy hold .
    focus $w
    bind $w <Destroy> { blt::busy release . }

    tixFileEntry $w.mesh -dialogtype tixExFileSelectDialog -labelside top \
	    -label "Mesh file" -variable meshfile \
	    -activatecmd [list configFileDlg_Mesh $w.mesh]
    $w.mesh subwidget entry config -width 30
    tixFileEntry $w.qm -dialogtype tixExFileSelectDialog -labelside top \
	    -label "QM file" -variable qmfile \
	    -activatecmd [list configFileDlg_QM $w.qm]
    $w.qm subwidget entry config -width 30

    set bndc [tixOptionMenu $w.bndc -label "BC" -labelside left \
	    -dynamicgeometry false]
    $bndc subwidget label config -width 8
    $bndc add command robin -label Robin
    $bndc add command dirichlet -label Dirichlet
    $bndc config -value $bc -variable bc

    set slvr [tixOptionMenu $w.slvr -label "Solver" -labelside left \
            -dynamicgeometry false]
    $slvr subwidget label config -width 8
    $slvr add command cholesky -label "Cholesky"
    $slvr add command cg -label "Conjugate gradient"
    $slvr config -value $solver -variable solver

    pack $w.mesh $w.qm $w.bndc $w.slvr -side top -padx 10 -pady 4 -fill x

    # Button bar ============================================
    set bb [tixButtonBox $w.bbar -orientation horizontal -pady 2 -relief flat]
    $bb add apply -text Apply -underline 0 -width 6 \
	    -command [list applyMesh $w]
    $bb add cancel -text Cancel -underline 0 -width 6 \
	    -command [list destroy $w]
    pack $bb -side bottom -fill x

    return TCL_OK
}

proc configFileDlg_Mesh {w} {
    $w filedialog config -title "Load mesh file"
    set fsb [$w filedialog subwidget fsbox]
    $fsb config -filetypes {
	{{*.msh *.opt} {Mesh files}}
	{{*} {All files}}
    }
    $fsb subwidget types pick 0
    return TCL_OK
}

proc configFileDlg_QM {w} {
    $w filedialog config -title "Load QM file"
    set fsb [$w filedialog subwidget fsbox]
    $fsb config -filetypes {
	{{*.qm} {QM files}}
	{{*} {All files}}
    }
    $fsb subwidget types pick 0
    return TCL_OK
}

proc configFileDlg_Data {w} {
    $w filedialog config -title "Load data"
    set fsb [$w filedialog subwidget fsbox]
    $fsb config -filetypes {
	{{*.fem *.dat} {Data files}}
	{{*} {All files}}
    }
    $fsb subwidget types pick 0
    return TCL_OK
}

proc applyMesh {w} {
    global meshfile qmfile bc solver
    $w.mesh update
    $w.qm update
    femdataLoadQMMesh $meshfile $bc $qmfile $solver
    destroy $w
    return TCL_OK
}

# =============================================================================
# "Source spec" dialog

proc dlgSource {} {
    global srctp srcprf srcwdt

    set w [toplevel .dlgSource]
    wm title $w "Source specification"

    blt::busy hold .
    focus $w
    bind $w <Destroy> { blt::busy release . }

    set sctp [tixOptionMenu $w.sctp -label "Type" -labelside left \
	    -dynamicgeometry false]
    $sctp subwidget label config -width 8
    $sctp add command isotropic -label Isotropic
    $sctp add command neumann -label Neumann
    $sctp add command explicit -label Explicit
    $sctp config -value $srctp -variable srctp

    set scpr [tixOptionMenu $w.scpr -label "Profile" -labelside left \
	    -dynamicgeometry false]
    $scpr subwidget label config -width 8
    $scpr add command point -label Point
    $scpr add command gaussian -label Gaussian
    $scpr config -value $srcprf -variable srcprf

    set scwd [tixControl $w.scwd -allowempty false -disablecallback true \
		  -integer false -label Width -labelside left -min 0 \
		  -value $srcwdt -variable srcwdt]
    $scwd subwidget label config -width 8
    $scwd subwidget entry config -width 12

    pack $w.sctp $w.scpr $w.scwd -side top -padx 10 -pady 4 -fill x

    # Button bar ============================================
    set bb [tixButtonBox $w.bbar -orientation horizontal -pady 2 -relief flat]
    $bb add apply -text Apply -underline 0 -width 6 \
	    -command [list applySource $w]
    $bb add cancel -text Cancel -underline 0 -width 6 \
	    -command [list destroy $w]
    pack $bb -side bottom -fill x

    return TCL_OK
}

proc applySource {w} {
    global srctp srcprf srcwdt

    $w.scwd update
    femdataSetSourceSpec $srctp $srcprf $srcwdt
    destroy $w
    return TCL_OK
}

# =============================================================================
# "Parameter" dialog

proc dlgParams {} {
    global meshmod_whichparam meshmod_regs
    global meshmod_p1 meshmod_p2 meshmod_p3

    global do_mua do_mus do_n
    set do_mua 1
    set do_mus 1
    set do_n 1

    toplevel .dlg
    wm title .dlg "Parameters"

    set name [tixOptionName .dlg]
    option add *$name*TixControl*entry.width 12
    option add *$name*TixControl*label.width 8
    option add *$name*TixControl*label.anchor e

    frame .dlg.f -relief groove -bd 2
    tixComboBox .dlg.f.reg -label Region -dropdown true -editable false \
	    -listwidth 30 -options { entry.width 12 }
    pack .dlg.f.reg -side top -fill x -expand true
    ScanRegions meshmod_regs
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
	-command "applyParameters .dlg"
    .dlg.bbar add done -text Done -underline 0 -width 6 \
	-command "destroy .dlg"
    pack .dlg.bbar -side bottom -fill x
    return TCL_OK
}

proc toggle_entry {s e} {
    if {$s == 1} {
	$e config -state normal
    } else {
	$e config -state disabled
    }
}

proc applyParameters {w} {
    global meshmod_p1 meshmod_p2 meshmod_p3
    foreach d {1 2 3} {
	if {[$w.p.f$d.e cget -state] == "normal"} {
	    set v meshmod_p$d
	    set p$d [subst $$v]
	} else {
	    set p$d -1
	}
    }
    set v [$w.f.reg cget -value]
    femdataSetRegparam $v $p1 $p2 $p3
    puts "Set region $v to $p1 $p2 $p3"
    return TCL_OK
}

# =============================================================================
# Execute femdata

proc femdataGo {} {
    global nsets ndata
    blt::busy hold .
    statusOut "FEM solver running ..."
    update
    femdataExec
    femdataEvalTransforms
    updateAllGraphs
    statusReady
    blt::busy release .
    return TCL_OK
}

# =============================================================================

proc dlgNewSet {type} {
    global qmfile
    if {$type == "Static" && $qmfile == ""} {
	ErrorNotice "Cannot load static data sets without a mesh and QM file. Use \"Solver | Mesh and QM\"."
	return TCL_OK
    }
    set setno [dlgInsertSet]
    if {$setno < 0} {return TCL_OK}
    set idx [femdataDsetIndex $setno]
    if {$idx >= 0} {
	femdataDsetDelete $setno
	listDelete $idx
    }
    set idx [femdataDsetInsert $setno $type]
    if {$idx < 0} {
	return TCL_ERROR
    } else {
	if {[dlgEdit$type $setno] == 1} {
	    listInsert $idx $setno [string index $type 0]
	} else {
	    femdataDsetDelete $setno
	}
	return TCL_OK
    }
}

# =============================================================================

proc DelSet {} {
    global curset nsets
    set idx [femdataDsetDelete $curset]
    if {$idx >= 0} {
	listDelete $idx
	if {$idx >= $nsets} {
	    set idx [expr $nsets - 1]
	}
	if {$idx >= 0} {
	    set lb [[.pw subwidget pane1].dlist subwidget listbox]
	    set curset [femdataDsetSetno $idx]
	    $lb activate $idx
	    $lb selection clear 0 end
	    $lb selection set $idx
	} else {
	    set curset -1;
	}
	selectData 0
	return TCL_OK
    } else {
	return TCL_ERROR
    }
}

# =============================================================================

proc dlgEditSet {} {
    global curset

    if {$curset < 0} {
	ErrorNotice "No data set available!"
	return 0
    }

    set type [femdataDsetType $curset]
    dlgEdit$type $curset
}

# =============================================================================

proc dlgEditActive {s} {
    global res dtype normalise nmom laps dfile enorm enormstr
    set dfile "unknown"

    set w [toplevel .dlgEditActive]
    wm title $w "Edit Active Set $s"

    blt::busy hold .
    focus $w
    bind $w <Destroy> { blt::busy release . }

    # Button bar ============================================
    tixButtonBox $w.bbar -orientation horizontal -pady 2 -relief flat
    $w.bbar add apply -text Apply -underline 0 -width 6 \
	    -command [list applyEditActive $w $s res]
    $w.bbar add cancel -text Cancel -underline 0 -width 6 \
	    -command { set res 0 }
    pack $w.bbar -side bottom -fill x

    set d [frame $w.data -relief groove -bd 2]
    pack $d -side top -fill x -expand true -padx 5 -pady 5

    label $d.l -text "Data"
    pack $d.l -side top

    set tp [tixOptionMenu $d.type -label Type \
	    -labelside left -dynamicgeometry false -variable dtype]
    $tp subwidget label config -width 6
    $tp add command logE -label "Integrated intensity"
    $tp add command Moment -label "Moment"
    $tp add command CMoment -label "Central moment"
    $tp add command Laplace -label "Laplace transform"
    $tp add command MellinLap -label "Mellin-Laplace"

    frame $d.norm
    label $d.norm.l -width 6
    checkbutton $d.norm.c -text "Normalise with intensity" \
	    -variable normalise -anchor w
    pack $d.norm.l $d.norm.c -side left

    set prm [frame $d.prm]
    tixControl $prm.nmom -allowempty false -label Nmom -labelside left \
	    -integer true -min 1 -variable nmom
    $prm.nmom subwidget label config -width 6
    $prm.nmom subwidget entry config -width 5
    $prm.nmom update
    tixControl $prm.laps -allowempty false -label S -labelside left \
	    -variable laps
    $prm.laps subwidget label config -width 6 -anchor e
    $prm.laps subwidget entry config -width 6
    $prm.laps update
    pack $prm.nmom $prm.laps -side left

    tixFileEntry $d.file -dialogtype tixExFileSelectDialog -label File \
	    -variable dfile -activatecmd [list configFileDlg_Data $d.file] 
    $d.file subwidget label config -width 6

    pack $d.type $d.norm $d.prm -side top -fill x -padx 5
    pack $d.file -side top -fill x -padx 5 -pady 3

    set d [frame $w.noise -relief groove -bd 2]
    pack $d -side top -fill x -expand true -padx 5 -pady 5

    label $d.l -text "Noise and SD"
    pack $d.l -side top

    set enm [tixOptionMenu $d.enm -label Enorm -labelside left \
	    -dynamicgeometry false -variable enorm]
    $enm subwidget label config -width 6
    for {set i 0} {$i < 6} {incr i} {
	$enm add command $i -label [lindex $enormstr $i]
    }

    pack $d.enm -side top -fill x -padx 5

    set tmp 0
    femdataGetDset $s dtype tmp normalise dfile enorm
    set prms [split $tmp " "]
    maskPrmEntries $prm $dtype
    setPrmEntries $prm $dtype $prms
    $tp config -command callbackDtype

    tkwait variable res
    destroy $w
    return $res
}

proc callbackDtype {dt} {
    maskPrmEntries .dlgEditActive.data.prm $dt 
    return TCL_OK
}

proc maskPrmEntries {w dt} {
    switch $dt {
      CMoment   -
      Moment    {$w.nmom config -state normal; $w.laps config -state disabled}
      Laplace   {$w.nmom config -state disabled; $w.laps config -state normal}
      MellinLap {$w.nmom config -state normal; $w.laps config -state normal}
      default   {$w.nmom config -state disabled; $w.laps config -state disabled}
    }
    return TCL_OK
}

proc setPrmEntries {w dt prms} {
    global nmom laps
    switch $dt {
	CMoment   -
	Moment    { set nmom [lindex $prms 0] }
	Laplace   { set laps [lindex $prms 0] }
	MellinLap { set nmom [lindex $prms 0]; set laps [lindex $prms 1] }
    }
    return TCL_OK
}

proc applyEditActive {w s res} {
    global dtype nmom laps normalise dfile enorm
    upvar $res r
    .dlgEditActive.data.prm.nmom update
    .dlgEditActive.data.prm.laps update
    .dlgEditActive.data.file update
    femdataSetDset $s $dtype $nmom $laps $normalise $dfile $enorm
    selectData 0
    set r 1
    return TCL_OK
}

# =============================================================================

proc dlgEditTransform {s} {
    global res transform

    set w [toplevel .dlgEditTransform]
    wm title $w "Edit Transform Set $s"

    blt::busy hold .
    focus $w
    bind $w <Destroy> { blt::busy release . }

    # Button bar ============================================
    tixButtonBox $w.bbar -orientation horizontal -pady 2 -relief flat
    $w.bbar add apply -text Apply -underline 0 -width 6 \
	    -command [list applyEditTransform $w $s res]
    $w.bbar add cancel -text Cancel -underline 0 -width 6 \
	    -command { set res 0 }
    pack $w.bbar -side bottom -fill x

    tixFileEntry $w.file -dialogtype tixExFileSelectDialog -labelside top \
	-label "Output file" -variable datafile \
	-activatecmd [list configFileDlg_Data $w.file] -options {
	    entry.width 30
	}
    # $w.file subwidget entry config -width 30

    tixLabelEntry $w.trans -label Transform -labelside top -options {
	entry.width 30
    }
    $w.trans subwidget entry config -textvariable transform
    set transform [femdataGetTransform $s]
    message $w.doc -justify left -aspect 1000 -anchor w -text \
    "Use \"sn\" to specify set n\nExamples:\n(s0-s1)/s0\ns5*0.1-log(s2)"
    pack $w.file $w.trans $w.doc -side top -padx 10 -pady 5 -fill x \
	-expand true

    tkwait variable res
    destroy $w
    return $res
}

proc applyEditTransform {w s res} {
    global transform datafile
    upvar $res r
    femdataSetTransform $s $transform $datafile
    selectData 0
    set r 1
    femdataEvalTransforms
    updateAllGraphs
    return TCL_OK
}

# =============================================================================

proc dlgEditStatic {s} {
    global res

    set w [toplevel .dlgEditStatic]
    wm title $w "Edit Static Set $s"

    blt::busy hold .
    focus $w
    bind $w <Destroy> { blt::busy release . }

    # Button bar ============================================
    tixButtonBox $w.bbar -orientation horizontal -pady 2 -relief flat
    $w.bbar add apply -text Load -underline 0 -width 6 \
	    -command [list applyEditStatic $w $s res]
    $w.bbar add cancel -text Cancel -underline 0 -width 6 \
	    -command { set res 0 }
    pack $w.bbar -side bottom -fill x

    tixFileEntry $w.file -dialogtype tixExFileSelectDialog -labelside top \
	    -label "Data file" -variable datafile \
	    -activatecmd [list configFileDlg_Data $w.file]
    $w.file subwidget entry config -width 30
    pack $w.file -side top -fill x -padx 10 -pady 4

    Line $w.line1

    label $w.setno -text "Set \# in file (multiset files only)"
    tixControl $w.fileidx -label "Set \#" -min 1 -value 1
    pack $w.setno $w.fileidx -side top -fill x

    tkwait variable res
    destroy $w
    return $res
}

proc applyEditStatic {w s res} {
    global datafile
    upvar $res r
    $w.fileidx update
    set fileidx [$w.fileidx cget -value]
    set idx [femdataDsetLoad $s $datafile $fileidx]
    if {$idx < 0} {
	ErrorNotice "The length of this data set differs from the current QM specification. Load failed."
	set r 0
    } else {
	set r 1
	selectData 0
	femdataEvalTransforms
	updateAllGraphs
    }
    return TCL_OK
}

# =============================================================================
# Select data set

proc selectData {lb} {
    global curset linkdset upddset enormstr
    set p2 [.pw subwidget pane2]

    if {$lb != 0} {
	set idx [$lb index active]
	set sel [femdataDsetSetno $idx]
	if {$sel == $curset} {
	    return TCL_OK
	}
	set curset $sel
    }

    if {$curset < 0} {
	$p2.l config -text "Set -"
	set linkdset 0
	set upddset 0
	$p2.mode.link config -state disabled
	$p2.mode.upd config -state disabled
	for {set x 0} {$x < 3} {incr x} {
	    $p2.data.$x.l config -text ""
	    $p2.data.$x.v config -text ""
	}
	return TCL_OK
    }

    set dsettype [femdataDsetType $curset]
    set datalabel [list "" "" ""]
    $p2.l config -text "Set $curset \[$dsettype\]"
    set linkdset [femdataGetDsetLinkState $curset]
    set upddset [femdataGetDsetUpdState $curset]
    set dfile [file tail [femdataGetDfile $curset]]
    if {$dfile == ""} {set dfile "<none>"}
    $p2.data.0.v config -text $dfile
    switch $dsettype {
	Active {
	    set datalabel [list "File" "Type" "Mode"]
	    $p2.mode.link config -state normal
	    $p2.mode.upd config -state normal
	    femdataGetDset $curset tp prm mod dfile enm
	    set prms [split $prm " "]
	    switch $tp {
		logE       { set prmline "" }
		CMoment    -
		Moment     { set prmline "(m=[lindex $prms 0])" }
		LTransform { set prmline "(s=[lindex $prms 0])" }
		MellinLap  { set prmline "(m=[lindex $prms 0], s=[lindex $prms 1])" }
		default    { set prmline $prm }
	    }
	    set m [list "not normalised" "normalised"]
	    set modline [lindex $m $mod]
	    set enmline [lindex $enormstr $enm]
	    $p2.data.1.v config -text [concat $tp $prmline]
	    $p2.data.2.v config -text $modline
	    $p2.noise.0.v config -text $enmline
	}
	Static {
	    set datalabel [list "File" "" "" ""]
	    $p2.mode.link config -state disabled
	    $p2.mode.upd config -state disabled
	    $p2.data.1.v config -text ""
	    $p2.data.2.v config -text ""
	}
	Transform {
	    set datalabel [list "File" "Trans." "" ""]
	    $p2.mode.link config -state disabled
	    $p2.mode.upd config -state normal
	    $p2.data.1.v config -text [femdataGetTransform $curset]
	    $p2.data.2.v config -text ""
	}
    }
    for {set x 0} {$x < 3} {incr x} {
	$p2.data.$x.l config -text [lindex $datalabel $x]
    }

    return TCL_OK
}

# =============================================================================
# Graph window functions

# Add new graph window

proc dataGraph {ng} {
    global graphlist
    set w [toplevel .graph$ng]
    wm title $w "Graph $ng"
    set gf [toastGraph $w.g]
    toastGraph_setCallbackAdd $gf dataGraph_addSet
    toastGraph_setCallbackDel $gf dataGraph_delSet
    toastGraph_setCallbackSet $gf dataGraph_setSet
    bind $gf <Destroy> [list dataGraph_Del $w $ng]
    lappend graphlist $ng

    # Add a few special items to the graph menu
    menubutton $gf.mnu.src -text Source -underline 0 -menu $gf.mnu.src.menu
    set m [menu $gf.mnu.src.menu]
    $m add command -label "Plot mode ..." -underline 0 \
	    -command [list dlgPlotMode $gf]
    tixControl $gf.mnu.csrc -allowempty false -integer true -min -1 -value 0 \
	    -command [list callbackSource $gf] -options { entry.width 4 }
    pack $gf.mnu.src $gf.mnu.csrc -side left

    upvar #0 plotdata$ng pdata
    set pdata(against) 0
    set pdata(mode) 0
}

# Kill graph window

proc dataGraph_Del {w n} {
    global graphlist
    destroy $w
    set idx [lsearch $graphlist $n]
    if {$idx >= 0} {
	set graphlist [lreplace $graphlist $idx $idx]
    }
}

# Callback for adding a new data set

proc dataGraph_addSet {w s xvec yvec} {
    upvar #0 $xvec xv
    upvar #0 $yvec yv
    if {[femdataDsetIndex $s] < 0} {
	ErrorNotice "Set $s does not exist!"
	set xv ""
	set yv ""
    } else {
	regexp {\.graph([0-9]+).*} $w dummy n ;# get graph id
	set id [join [list $n "x" $s] {}]
	upvar #0 plotdata$n pdata

	blt::vector create absc$id
	blt::vector create dv$id

	# find out which source we need
	$w.mnu.csrc update
	set q [$w.mnu.csrc cget -value]

	# read the data vectors
	femdataGetData $s $q $pdata(against) $pdata(mode) absc$id dv$id
	set xv absc$id
	set yv dv$id
    }
}

# Callback for deleting a data set

proc dataGraph_delSet {w s} {
    regexp {\.graph([0-9]+).*} $w dummy n ;# get graph id
    set id [join [list $n "x" $s] {}]
    blt::vector destroy dv$id
    blt::vector destroy absc$id
}

# Callback for selecting a data set

proc dataGraph_setSet {w s} {
}

# Callback for selecting a source

proc callbackSource {w q} {
    regexp {\.graph([0-9]+).*} $w dummy n ;# get graph id
    updateGraph $n
}

# Update a graph's data

proc updateGraph {n} {
    upvar #0 plotdata$n pdata
    .graph$n.g.mnu.csrc update
    set q [.graph$n.g.mnu.csrc cget -value]
    foreach el [.graph$n.g.g element names] {
	regexp {line([0-9]+)} $el line s
	set id [join [list $n "x" $s] {}]
	femdataGetData $s $q $pdata(against) $pdata(mode) absc$id dv$id
    }
    return TCL_OK
}

# Update the data for all graphs

proc updateAllGraphs {} {
    global graphlist
    foreach n $graphlist { updateGraph $n }
    return TCL_OK
}

# Dialog for modifying the abscissa mode

proc dlgPlotMode {w} {
    regexp {\.graph([0-9]+).*} $w dummy n ;# get graph id
    upvar #0 plotdata$n pdata
    set dlgname .toastgraph_mode
    if {[winfo exists $dlgname] == 0} {
	set dlg [toplevel $dlgname]
	wm title $dlg "Plot mode"

	set f [frame $dlg.bbar]
	button $f.apply -text Apply
	button $f.done -text Done -command [list destroy $dlg]
	pack $f.apply $f.done -side left -padx 5 -pady 5
	pack $f -side bottom

	set f [frame $dlg.1 -relief groove -bd 2]
	pack $f -side top -fill x
	set mlist {"Detector number" "Optode spacing" "Detector x-position" \
		"Detector y-position" "Detector z-position" \
		"Detector azimuth angle"}
	tixComboBox $f.ag -dropdown yes -editable no \
		-label "Plot data against" -labelside top -options {
	    entry.width 20
	    slistbox.scrollbar auto
	}
	for {set x 0} {$x < 6} {incr x} {
	    $f.ag insert end [lindex $mlist $x]
	}
	set mlist {"Absolute" "Relative to source"}
	tixComboBox $f.md -dropdown yes -editable no \
		-label "Mode" -labelside top -options {
	    entry.width 20
	    slistbox.scrollbar auto
	}
	for {set x 0} {$x < 6} {incr x} {
	    $f.md insert end [lindex $mlist $x]
	}
	pack $f.ag $f.md -side top
    } else {
	set dlg $dlgname
	raise $dlg
    }
    $dlg.bbar.apply config -command [list doPlotMode $dlg $n]
    $dlg.1.ag pick $pdata(against)
    $dlg.1.md pick $pdata(mode)
}

proc doPlotMode {dlg n} {
    global plotdata$n
    set against [join [list plotdata $n (against)] {}]
    set mode [join [list plotdata $n (mode)] {}]
    set $against [$dlg.1.ag subwidget listbox index active]
    set $mode [$dlg.1.md subwidget listbox index active]
    updateGraph $n
}

# =============================================================================
# main window

lappend auto_path $env(TOAST_SCRIPT_PATH)
toastscheme
set im [image create photo -file [file join $toastBitmapPath itdata.gif]]
toplevel .icon
pack [button .icon.b -image $im -relief flat -bd 0 -highlightthickness 0 -padx 0 -pady 0] -padx 0 -pady 0 -ipadx 0 -ipady 0
wm iconwindow . .icon
wm iconmask . @[file join $toastBitmapPath itdata.mask]


# ==== menu bar ====
MenuSetup .mnu
menubutton .mnu.file -text File -underline 0 -menu .mnu.file.menu
menu .mnu.file.menu
.mnu.file.menu add command -label "Load definition ..." -underline 5 \
	-command dlgLoadDef
.mnu.file.menu add separator
.mnu.file.menu add command -label "Exit" -underline 1 -command exit

menubutton .mnu.data -text Data -underline 0 -menu .mnu.data.menu
menu .mnu.data.menu
.mnu.data.menu add command -label "New active ..." -underline 4 \
	-command {dlgNewSet Active}
.mnu.data.menu add command -label "New transform ..." -underline 4 \
	-command {dlgNewSet Transform}
.mnu.data.menu add command -label "Load static ..." -underline 0 \
	-command {dlgNewSet Static}
.mnu.data.menu add command -label "Edit ..." -underline 0 \
	-command dlgEditSet
.mnu.data.menu add command -label "Delete" -underline 0 \
	-command DelSet

menubutton .mnu.cnfg -text Solver -underline 0 -menu .mnu.cnfg.menu
menu .mnu.cnfg.menu
.mnu.cnfg.menu add command -label "Mesh and QM ..." -underline 0 \
	-command dlgMesh
.mnu.cnfg.menu add command -label "Source specs ..." -underline 0 \
        -command dlgSource
.mnu.cnfg.menu add command -label "Parameters ..." -underline 0 \
	-command dlgParams
.mnu.cnfg.menu add separator
.mnu.cnfg.menu add command -label "Run" -underline 0 -command femdataGo

menubutton .mnu.win -text Window -underline 0 -menu .mnu.win.menu
menu .mnu.win.menu
.mnu.win.menu add command -label "New data graph ..." -underline 9 \
	-command {dataGraph [incr ng]}
pack .mnu.file .mnu.data .mnu.cnfg .mnu.win -side left -fill x

# ==== status bar ====
statusSetup

# ==== client area ====
set m [tixPanedWindow .pw -dynamicgeometry true -orient horizontal \
	-panerelief flat -panebd 1]
pack $m -side left -fill both -expand true
$m add pane1 -size 75
$m add pane2
set p1 [$m subwidget pane1]
set p2 [$m subwidget pane2]

label $p1.l -text Data
pack $p1.l -side top
set dl [tixScrolledListBox $p1.dlist -scrollbar auto]
$dl config -browsecmd "selectData [$dl subwidget listbox]" \
	-command dlgEditActive
for {set i 0} {$i < $nsets} {incr i} {
    $dl subwidget listbox insert end "Set $i \[A\]"
}
pack $dl -side top -fill both -expand true -padx 5 -pady 5

label $p2.l -text "Set -"
pack $p2.l -side top

set f [frame $p2.mode]
checkbutton $p2.mode.link -text "Link to solver" -variable linkdset \
    -command {femdataSetDsetLinkState $curset $linkdset} -anchor w \
    -state disabled
checkbutton $p2.mode.upd -text "Update file(s)" -variable upddset \
    -command {femdataSetDsetUpdState $curset $upddset} -anchor w \
    -state disabled
pack $p2.mode.link $p2.mode.upd -side left
pack $p2.mode -side top

set d [frame $p2.data -relief groove -bd 2]
pack $d -side top -padx 5 -pady 5 -fill x

label $d.l -text Data
pack $d.l -side top

for {set x 0} {$x < 3} {incr x} {
    set f [frame $d.$x]
    label $f.l -width 7 -anchor e
    label $f.v -width 30 -anchor w -foreground blue
    pack $f.l -side left -padx 2
    pack $f.v -side left -fill x -padx 2 -expand true
    pack $f -side top -fill x
}

set d [frame $p2.noise -relief groove -bd 2]
pack $d -side top -padx 5 -pady 5 -fill x

label $d.l -text Noise
pack $d.l -side top

set datalabel [list "Enorm" "Scale"]
for {set x 0} {$x < 2} {incr x} {
    set f [frame $d.$x]
    label $f.l -text [lindex $datalabel $x] -width 7 -anchor e
    label $f.v -text "-" -width 30 -anchor w -foreground blue
    pack $f.l -side left -padx 2
    pack $f.v -side left -fill x -padx 2 -expand true
    pack $f -side top -fill x
}

if {$nsets > 0} {
    $dl subwidget listbox activate 0
    $dl subwidget listbox selection set 0
    selectData [$dl subwidget listbox]
}
focus [$dl subwidget listbox]
update

femdataGo

