wm title . "Toast 'im!"

set imgparam(fixaspect) 1
set imgparam(applymode) 0
set imgparam(3dmode) 0
set tmp_3dmode 0
set imgparam(gridx) 10
set imgparam(gridy) 10
set imgparam(gridz) 10
set imgparam(nimg) 0
set imgparam(cimg) 0
set imgparam(mbbox) 0
set imgparam(sbbox) 0
set imgparam(dim) 2
set imgparam(imgmin) 0
set imgparam(imgmax) 0
set imgparam(totmin) 0
set imgparam(totmax) 0
set imgparam(imgfile) 0
set imgparam(mesh) 0
set imgparam(logscale) 0
set imgparam(imgfmt) 0
set imgparam(orient3d) 0

# contour parameters
set cntrparam(draw) 0
set cntrparam(ncntr) 16
set cntrparam(col) white
set cntrparam(noimg) 0
set cntrparam(thickline) 0

set colourmode(rgb) 0
set colourmode(pal) ""

set loadparam(reraster) 1
set loadparam(rescale) 1

#image tags
set rimgtag 0
set rimgxytag 0
set rimgxztag 0
set rimgyztag 0

set transform ""
set ntransform 0

# convert window pixel coordinates to image pixel coordinates,
# by removing border offset

proc AdjustScreenCoord {x y} {
    upvar $x lx
    upvar $y ly
    incr lx -2
    incr ly -2
    if {$lx < 0} {set lx 0}
    if {$ly < 0} {set ly 0}
    if {$lx >= [image width rimg]} {
	set lx [image width rimg]
	incr lx -1
    }
    if {$ly >= [image height rimg]} {
	set ly [image height rimg]
	incr ly -1
    }
}

# =============================================================================
# "Load image" dlg

proc dlgLoadImg {} {
    global rc
    set dialog [tix filedialog tixExFileSelectDialog]
    $dialog config -command "doLoadImg" -title "Load image"
    set dbox [$dialog subwidget fsbox]
    $dbox config -directory $rc(pwd) -filetypes {{{*.nim} {Nodal image files}}}
    $dbox subwidget types pick 0
    $dialog popup
    return TCL_OK
}

proc doLoadImg {file} {
    global imgparam loadparam mbbox sbbox
    blt::busy hold .
    puts "toastim: Loading image $file"

    set range(totmin) $imgparam(totmin)
    set range(totmax) $imgparam(totmax)
    set range(imgmin) $imgparam(imgmin)
    set range(imgmax) $imgparam(imgmax)
    set range(reqmin) $imgparam(imgmin)
    set range(reqmax) $imgparam(imgmax)
    if {$loadparam(rescale) == 0} {
	toastimGetRange range
    }

    # Scan image file, store image parameters in local variable
    set err [toastimOpenImage $file loc_imgparam]
    if $err {
	set msg [list "ok" \
		     "Image file not found" \
		     "Invalid image format" \
		     "Error parsing image header" \
		     "File contains no images" \
		     "Image size not defined in header" \
		     "Mesh file not defined in image header" \
		     "Error loading mesh"]
	ErrorNotice "Error during opening of image file [file tail $file]: \
                     [lindex $msg $err]. Load terminated."
	return TCL_OK
    }

    # Copy params to global array
    foreach field [array names loc_imgparam] {
	set old_imgparam($field) $imgparam($field)
	set imgparam($field) $loc_imgparam($field)
    }

    if {$loadparam(reraster)} {

	# init bounding box variables
	toastimGetBB mesh mbbox
	toastimGetBB sel sbbox
	toastimGet3Dparam imgparam

    } else {

	foreach d {x y z} {
	    set imgparam(grid$d) $old_imgparam(grid$d)
#	    set imgparam(plane$d) $old_imgparam(plane$d)
	}
	set imgparam(mbbox) $old_imgparam(mbbox)
	set imgparam(sbbox) $old_imgparam(sbbox)
	set imgparam(orient3d) $old_imgqparam(orient3d)
	toastimSet3Dparam imgparam
	UpdateCanvas
    }

    # Update image list
    set dl [[.pw subwidget pane1].dlist subwidget listbox]
    $dl delete 0 end
    for {set i 1} {$i <= $imgparam(nimg)} {incr i} {
	$dl insert end "$i"
    }
    $dl selection set $imgparam(cimg)
    $dl see $imgparam(cimg)

    # Sample image
    statusOut "Pre-sampling image ..."
    set err [toastimResampleImage imgparam sbbox]
    if $err {
	set msg [list "ok" "Invalid image format"]
	ErrorNotice "Error during sampling: [lindex $msg $err]. \
                     Load terminated."
	return TCL_OK
    }

    if {$loadparam(rescale) == 0} {
	set rangeprm(apply) 1
	set rangeprm(mode) 2
	set rangeprm(vmin) $range(reqmin)
	set rangeprm(vmax) $range(reqmax)
	toastimSetRange rangeprm
    }

    # Rasterize image
    statusOut "Rasterizing image ..."
    set cimg $imgparam(cimg)
    toastimRasterImage $cimg
    UpdateCanvas
    incr cimg

    # unlock menu and init image selector
    .mnu.image config -state normal
    .mnu.list config -state normal
    .mnu.anls config -state normal
    .mnu.file.menu entryconfig 2 -state normal
    #.mnu.cimg config -disablecallback yes
    #.mnu.cimg config -min 1 -max $imgparam(nimg) -value $cimg
    #.mnu.cimg config -disablecallback no

    # set 3D controls
    if {$imgparam(dim) < 3} {
	.mnu.3d config -state disabled
	hide3dcontrols
    } else {
	.mnu.3d config -state normal
	show3dcontrols $imgparam(3dmode)
    }

    # bring up the image window if necessary
    showImage
    ImageResized
    Refresh

    statusReady
    blt::busy release .
    return TCL_OK
}

# =============================================================================
# Image list operations
# =============================================================================

# =============================================================================
# Delete list entry

proc delListEntry {idx} {
    global imgparam
    toastimDeleteEntry imgparam $idx
    set dl [[.pw subwidget pane1].dlist subwidget listbox]
    $dl delete $idx
    $dl selection set $imgparam(cimg)
    set i $imgparam(cimg)
    incr i
    set imgparam(cimg) -1
    callbackChangeImage $i
}

proc delListBefore {idx} {
    global imgparam
    set dl [[.pw subwidget pane1].dlist subwidget listbox]
    while {$imgparam(cimg) > 0} {
	toastimDeleteEntry imgparam 0
	$dl delete 0
    }
    $dl selection set 0
    set imgparam(cimg) -1
    callbackChangeImage 1
}

proc delListAfter {idx} {
    global imgparam
    set dl [[.pw subwidget pane1].dlist subwidget listbox]
    while {$imgparam(cimg) < $imgparam(nimg)-1} {
	toastimDeleteEntry imgparam [expr $imgparam(nimg)-1]
	$dl delete end
    }
    $dl selection set end
    set imgparam(cimg) -1
    callbackChangeImage $imgparam(nimg)
}

proc UpdateCanvas {} {
    global imgparam
    if {$imgparam(3dmode) == 5} {
	set imgparam(3dmode) 0
	set imgparam(orient3d) z
	toastimSet3Dparam imgparam
	toastimMapImage rimgxy
	set imgparam(orient3d) y
	toastimSet3Dparam imgparam
	toastimMapImage rimgxz
	set imgparam(orient3d) x
	toastimSet3Dparam imgparam
	toastimMapImage rimgyz
	set imgparam(3dmode) 5
	toastimSet3Dparam imgparam
	set w [image width rimgyz]
	set h [image height rimgxy]
	toastimDrawCSCube rimg $w $h
	multi3dRedrawCSIndicators
    } else {
	toastimMapImage rimg
    }
}

proc multi3dRedrawCSIndicators {} {
    global imgparam rimgxytag rimgxztag rimgyztag
    set c .imgwin.c

    set x [expr ($imgparam(planex)*1.0) / $imgparam(gridx)]
    set y [expr ($imgparam(planey)*1.0) / $imgparam(gridy)]
    set z [expr ($imgparam(planez)*1.0) / $imgparam(gridz)]

    set w [image width rimgxy]
    set h [image height rimgxy]
    set crd [$c coords $rimgxytag]
    set ix [expr $x * $w]
    set iy [expr (1.0 - $y) * $h - 1]

    $c coords xycshorz [lindex $crd 0] [expr [lindex $crd 1] + $iy] \
	[expr [lindex $crd 0] + $w] [expr [lindex $crd 1] + $iy]
    $c coords xycsvert [expr [lindex $crd 0] + $ix] [lindex $crd 1] \
	[expr [lindex $crd 0] + $ix] [expr [lindex $crd 1] + $h]

    set w [image width rimgxz]
    set h [image height rimgxz]
    set crd [$c coords $rimgxztag]
    set ix [expr $x * $w]
    set iz [expr (1.0 - $z) * $h - 1]
    $c coords xzcshorz [lindex $crd 0] [expr [lindex $crd 1] + $iz] \
	[expr [lindex $crd 0] + $w] [expr [lindex $crd 1] + $iz]
    $c coords xzcsvert [expr [lindex $crd 0] + $ix] [lindex $crd 1] \
	[expr [lindex $crd 0] + $ix] [expr [lindex $crd 1] + $h]

    set w [image width rimgyz]
    set h [image height rimgyz]
    set crd [$c coords $rimgyztag]
    set iy [expr $y * $w]
    set iz [expr (1.0 - $z) * $h - 1]
    $c coords yzcshorz [lindex $crd 0] [expr [lindex $crd 1] + $iz] \
	[expr [lindex $crd 0] + $w] [expr [lindex $crd 1] + $iz]
    $c coords yzcsvert [expr [lindex $crd 0] + $iy] [lindex $crd 1] \
	[expr [lindex $crd 0] + $iy] [expr [lindex $crd 1] + $h]
}

# =============================================================================
# Show image window

proc showImage {} {
    global imgparam rimgtag
    if [winfo exists .imgwin] {
	raise .imgwin
    } else {
	set geom [wm geometry .]
	regexp {([0-9]+)x([0-9]+)[+-]([0-9]+)[+-]([0-9]+)} $geom \
	    match width height x y
	incr x $width
	set w [toplevel .imgwin]
	wm geometry $w "+$x+$y"
	if {$imgparam(3dmode) == 5} {
	    set cw [image width rimgxy]
	    incr cw [image width rimgyz]
	    set ch [image height rimgxy]
	    incr ch [image height rimgxz]
	} else {
	    set cw [image width rimg]
	    set ch [image height rimg]
	}
	pack [canvas $w.c -cursor crosshair -width $cw \
		  -height $ch -closeenough 5 \
		  -relief flat -borderwidth 0] \
	    -padx 0 -pady 0 -ipadx 0 -ipady 0
	set rimgtag [$w.c create image 0 0 -anchor nw -image rimg]
	wm resizable $w false false
    }
    wm title .imgwin [file tail $imgparam(imgfile)]

    # set mouse bindings
    if {$imgparam(dim) == 2} {
	bind .imgwin <Motion> {callback_motion_2d %x %y}
    } elseif {$imgparam(3dmode) == 0} {
	bind .imgwin <Motion> {callback_motion_3d %x %y}
    } else {
	bind .imgwin <Motion> {}
    }

    ImageChanged
}

proc ImageResized {} {
    global imgparam rimgtag rimgxytag rimgxztag rimgyztag
    set c .imgwin.c

    if [winfo exists .imgwin] {
	if {$imgparam(3dmode) == 5} {
	    set cw [image width rimgxy]
	    incr cw [image width rimgyz]
	    set ch [image height rimgxy]
	    incr ch [image height rimgxz]

	    set yzw [image width rimgyz]
	    set xyh [image height rimgxy]
	    $c coords $rimgtag 0 0

	    set crds [$c coords $rimgxytag]
	    $c move $rimgxytag [expr $yzw - [lindex $crds 0]] 0
	    set crds0 [$c coords xyaxes]
	    $c move xyaxes [expr $yzw - [lindex $crds0 0] + 5] \
		[expr [image height rimgxy] - [lindex $crds0 1] - 25]

	    set crds [$c coords $rimgyztag]
	    $c move $rimgyztag 0 [expr $xyh - [lindex $crds 1]]
	    set crds0 [$c coords yzaxes]
	    $c move yzaxes [expr 0 - [lindex $crds0 0] + 5] \
		[expr $xyh + [image height rimgyz] - [lindex $crds0 1] - 25]

	    set crds [$c coords $rimgxztag]
	    $c move $rimgxztag [expr $yzw - [lindex $crds 0]] \
		[expr $xyh - [lindex $crds 1]]
	    set crds0 [$c coords xzaxes]
	    $c move xzaxes [expr $yzw - [lindex $crds0 0] + 5] \
		[expr $xyh + [image height rimgxz] - [lindex $crds0 1] - 25]
	    
	    multi3dRedrawCSIndicators

	} else {
	    set cw [image width rimg]
	    set ch [image height rimg]
	}

	.imgwin.c config -width $cw -height $ch
	.imgwin.c xview moveto 0
	.imgwin.c yview moveto 0
    }
    ImageChanged
}

proc ImageChanged {} {
    winProfile_update
}

# =============================================================================
# "Add transform" dialog

proc dlgAddTransform {} {
    global res transform

    set dlg .transformdlg
    toplevel $dlg
    wm title $dlg "Add transform image"
	
    # Button bar ============================================
    tixButtonBox $dlg.bbar -orientation horizontal -pady 2 -relief flat
    $dlg.bbar add apply -text Apply -underline 0 -width 6 \
	-command [list applyAddTransform $dlg res]
    $dlg.bbar add cancel -text Cancel -underline 0 -width 6 \
	-command { set res 0 }
    pack $dlg.bbar -side bottom -fill x
    
    tixLabelEntry $dlg.trans -label Transform -labelside top -options {
	entry.width 30
    }
    $dlg.trans subwidget entry config -textvariable transform
    #set transform [femdataGetTransform $s]
    message $dlg.doc -justify left -aspect 1000 -anchor w -text \
	"Use \"sn\" to specify set n\nExamples:\n(s0-s1)/s0\ns5*0.1-log(s2)"
    pack $dlg.trans $dlg.doc -side top -padx 10 -pady 5 -fill x \
	-expand true
    
    tkwait variable res
    destroy $dlg
    return $res
}

proc applyAddTransform {dlg res} {
    global transform ntransform imgparam
    upvar $res r
    toastimNewTransform $transform imgparam
    incr ntransform

    set dl [[.pw subwidget pane1].dlist subwidget listbox]
    $dl insert end "T$ntransform"

    set r 1
    return TCL_OK
}

# =============================================================================
# "Resample image" dialog

proc dlgResampleImg {} {
    global imgparam sbbox
    set dlg .resampledlg
    if [winfo exists $dlg] {
	raise $dlg
    } else {
	toplevel $dlg
	wm title $dlg "Resample image"
	wm resizable $dlg false false

	# button bar
	set f [frame $dlg.bbar]
	button $f.apply -text Apply -command [list doResampleImg]
	button $f.reset -text "Reset BB" -command [list callbackResetBB $dlg]
	button $f.done -text Done -command [list destroy $dlg]
	pack $f.apply $f.reset $f.done -side left -padx 5
	pack $f -side bottom -pady 10

	set name [tixOptionName $dlg]
	option add *$name*TixControl*entry.width 10
	option add *$name*TixControl*label.width 6
	option add *$name*TixControl*label.anchor e

	# image bbox selector
	set f [frame $dlg.sbbox]
	pack $f -side top -fill x -expand yes -padx 10 -pady 10
	label $f.l -text "Image bounding box" -anchor w
	pack $f.l -side top -fill x
	foreach d {x y z} {
	    set f1 [frame $f.$d]
	    foreach m {min max} {
		tixControl $f1.$m -label "$d $m"
	    }
	    pack $f1.min $f1.max -side left -padx 2 -pady 1
	    pack $f1
	}

	Line $dlg.line1

	# image size selector
	set f [frame $dlg.size]
	pack $f -side top -fill x -expand yes -padx 10 -pady 10
	label $f.l -text "Image size" -anchor w
	pack $f.l -side top -fill x
	pack [frame $f.dim] -side top -fill x
	foreach d {x y z} {
	    tixControl $f.dim.$d -label $d -integer yes -options {
		entry.width 6
		label.width 3
	    }
	    pack $f.dim.$d -side left -padx 2 -pady 1
	}
	checkbutton $f.aspect -text "Maintain aspect ratio" \
	    -variable imgparam(fixaspect)
	pack $f.aspect -side top -fill x

	Line $dlg.line2

	initResampleImg
	foreach d {x y z} {
	    foreach m {min max} {
		$dlg.sbbox.$d.$m config -command \
		    [list callbackResampleAspect $dlg bb $d]
	    }
	    $dlg.size.dim.$d config -command \
		[list callbackResampleAspect $dlg sz $d]
	}
	$dlg.size.aspect config -command \
	    [list callbackResampleAspect $dlg sz x 0]
	
	callbackResampleAspect $dlg bb x 0
    }
    
    return TCL_OK
}

proc initResampleImg {} {
    global imgparam mbbox sbbox
    set dlg .resampledlg
    if {[winfo exists $dlg] == 0} {return TCL_OK}
    foreach d {x y z} {
	foreach m {min max} {
#	    $dlg.sbbox.$d.$m config -min $mbbox(min$d) -max $mbbox(max$d) \
		-value $sbbox($m$d)
	    $dlg.sbbox.$d.$m config -value $sbbox($m$d)
	}
	$dlg.size.dim.$d config -value $imgparam(grid$d)
    }
    if {$imgparam(dim) < 3} {
	$dlg.sbbox.z.min config -state disabled
	$dlg.sbbox.z.max config -state disabled
	$dlg.size.dim.z config -state disabled
    } else {
	$dlg.sbbox.z.min config -state normal
	$dlg.sbbox.z.max config -state normal
	$dlg.size.dim.z config -state normal
    }
}

proc callbackResetBB {dlg} {
    global sbbox mbbox
    foreach d {x y z} {
	foreach m {min max} {
	    set sbbox($m$d) $mbbox($m$d)
	    $dlg.sbbox.$d.$m config -value $sbbox($m$d)
	}
    }
}

proc callbackResampleAspect {dlg who dim v} {
    global imgparam sbbox
    if {$imgparam(fixaspect) == 0} {
	return TCL_OK
    }
    foreach d {x y z} {
	set rold($d) [expr $sbbox(max$d) - $sbbox(min$d)]
	foreach m {min max} {
	    $dlg.sbbox.$dim.$m update
	    set sbbox($m$d) [$dlg.sbbox.$d.$m cget -value]
	}
	$dlg.size.dim.$d update
	set grid($d) [$dlg.size.dim.$d cget -value]
	set rnew($d) [expr $sbbox(max$d) - $sbbox(min$d)]
    }

    if {$who == "bb"} {
	# we are modifying the bounding box
	if {$rold($dim) > 0} {
	    set scale [expr $rnew($dim) / $rold($dim)]
	    set gnew  [expr round ($grid($dim) * $scale)]
	    if {$grid($dim) != $gnew} {
		$dlg.size.dim.$dim config -value $gnew
	    }
	}
    } elseif {$who == "sz"} {
	# we are modifying the image size
	if {$rnew($dim) > 0} {
	    foreach d {x y z} {
		if {$d != $dim} {
		    set scale [expr $rnew($d) / $rnew($dim) * $grid($dim)]
		    set gnew  [expr round ($scale)]
		    if {$grid($d) != $gnew} {
			$dlg.size.dim.$d config -value $gnew
		    }
		}
	    }
	}
    }
}

proc doResampleImg {} {
    global sbbox imgparam
    set dlg .resampledlg
    if [winfo exists $dlg] {
	foreach d {x y z} {
	    foreach m {min max} {$dlg.sbbox.$d.$m update}
	    $dlg.size.dim.$d update
	}
	foreach d {x y z} {
	    foreach m {min max} {
		set sbbox($m$d) [$dlg.sbbox.$d.$m cget -value]
	    }
	    set imgparam(grid$d) [$dlg.size.dim.$d cget -value]
	}
    }
    toastimResampleImage imgparam sbbox
    toastimRasterImage $imgparam(cimg)
    UpdateCanvas
    ImageResized
    Refresh
}

# =============================================================================
# "Contours" dialog

proc dlgContours {} {
    global cntrparam
    set dlg .contourdlg
    if [winfo exists $dlg] {
	raise $dlg
    } else {
	toplevel $dlg
	wm title $dlg "Edit contours"
	wm resizable $dlg false false

	# button bar
	set f [frame $dlg.bbar]
	button $f.apply -text Apply -command [list doContours $dlg]
	button $f.done -text Done -command [list destroy $dlg]
	pack $f.apply $f.done -side left -padx 5
	pack $f -side bottom -pady 10

	set f [frame $dlg.f1]
	pack $f -side top -fill x -expand yes -padx 10 -pady 10
	checkbutton $f.con -text "Draw contours" -anchor w \
	    -variable cntrparam(draw)
	checkbutton $f.noimg -text "Blackout image" -anchor w \
	    -variable cntrparam(noimg)
	checkbutton $f.thickline -text "Thick lines" -anchor w \
	    -variable cntrparam(thickline)
	pack $f.con $f.noimg $f.thickline -side top -fill x

	Line $dlg.line1

	set f [frame $dlg.f2]
	pack $f -side top -fill x -expand yes -padx 10 -pady 10
	grid [label $f.l2 -text "\# Levels"] -row 0 -column 0 -sticky e

	tixControl $f.ncntr -integer yes \
	    -min 2 -variable cntrparam(ncntr) -options {
		entry.width 6
	    }
	grid $f.ncntr -row 0 -column 1 -sticky w

	grid [label $f.l3 -text Colour] -row 1 -column 0 -sticky e
	grid [colorChooser $f.col -options {entry.width 8}] -row 1 -column 1 \
	    -sticky w
	colorChooser_setColor $f.col $cntrparam(col)

	Line $dlg.line2
    }
}

proc doContours {dlg} {
    global imgparam cntrparam
    $dlg.f2.ncntr update
    set cntrparam(col) [colorChooser_getColor $dlg.f2.col]
    toastimSetContours cntrparam
    toastimRasterImage $imgparam(cimg)
    UpdateCanvas
}

# =============================================================================
# "Contrast" dialog

proc dlgContrast {} {
    global imgparam
    set dlg .contrastdlg
    if [winfo exists $dlg] {
	raise $dlg
    } else {
	toplevel $dlg
	wm title $dlg "Image contrast"
	wm resizable $dlg false false

	# button bar
	set f [frame $dlg.bbar]
	button $f.resglob -text "Reset global" \
	    -command [list doContrast $dlg 1]
	button $f.resimg -text "Reset individual" \
	    -command [list doContrast $dlg 0]
	button $f.apply -text Apply -command [list doContrast $dlg 2]
	button $f.done -text Done -command [list destroy $dlg]
	pack $f.resglob $f.resimg $f.apply $f.done -side left -padx 5
	pack $f -side bottom -pady 10

	set f [frame $dlg.f1]
	pack $f -side top -fill x -expand yes -padx 10 -pady 10
	for {set r 0} {$r < 2} {incr r} {
	    grid [label $f.l$r -text "[lindex [list Total Image] $r] range"] \
		-row $r -column 0 -sticky e -padx 5
	    set what [lindex [list tot img] $r]
	    grid [label $f.min$what] -row $r -column 1 -sticky w
	    grid [label $f.dsh$what -text "to"] -row $r -column 2 -padx 5
	    grid [label $f.max$what] -row $r -column 3 -sticky w
	}
	grid [label $f.l2 -text "Set range"] -row 2 -column 0 -sticky e -padx 5
	grid [entry $f.minreq -width 9] -row 2 -column 1
	grid [label $f.dshset -text "to"] -row 2 -column 2 -padx 5
	grid [entry $f.maxreq -width 9] -row 2 -column 3

	# the image range barchart
	set f [frame $dlg.f2]
	pack $f -side top -fill x -expand yes -padx 10 -pady 0
	set bc [blt::barchart $f.bc -invertxy yes -height 1.3i -width 4.5i \
		    -barwidth 1 -plotpadx 0 -plotpady 0]
	$bc legend config -hide yes
	$bc axis config x -majorticks {0 1 2} -command callbackRangeLabels
	$bc element create bg_lo -xdata {0 1 2} -ydata {0 0 0} \
	    -barwidth 1 -relief solid -foreground white
	$bc element create req -xdata {0} -ydata {0} -relief solid \
	    -barwidth 1 -foreground green
	$bc element create img -xdata {1} -ydata {0} -relief solid \
	    -barwidth 1 -foreground grey50
	$bc element create tot -xdata {2} -ydata {0} -relief solid \
	    -barwidth 1 -foreground grey70
	$bc element create bg_hi -xdata {0 1 2} -ydata {1 1 1} -relief solid \
	    -barwidth 1 -foreground white
	pack $bc -side top

	set f [frame $dlg.f3]
	pack $f -side top -fill x -expand yes -padx 10 -pady 10
	radiobutton $f.single -text "Apply to current image" -anchor w \
	    -variable imgparam(applymode) -value 0
	radiobutton $f.all -text "Apply to all images" -anchor w \
	    -variable imgparam(applymode) -value 1
	pack $f.single $f.all -side top -fill x
    }
    dlgContrast_update
}

proc callbackRangeLabels {w v} {
    return [lindex [list User Image Total] $v]
}

proc dlgContrast_update {} {
    global imgparam
    set dlg .contrastdlg
    if {[winfo exists $dlg] == 0} {return}

    set data(totmin) 0
    set data(totmax) 0
    set data(imgmin) 0
    set data(imgmax) 0
    set data(reqmin) 0
    set data(reqmax) 0
    toastimGetRange data

    set f $dlg.f1
    for {set r 0} {$r < 2} {incr r} {
	set what [lindex [list tot img] $r]
	foreach c {1 3} {
	    set m [lindex [list x min x max] $c]
	    $f.$m$what config -text $data($what$m)
	}
    }
    $f.minreq delete 0 end
    $f.maxreq delete 0 end
    $f.minreq insert 0 $data(reqmin)
    $f.maxreq insert 0 $data(reqmax)

    # set up the barchart
    set bc $dlg.f2.bc
    set bmin $data(totmin)
    if {$data(reqmin) < $data(totmin)} { set bmin $data(reqmin) }
    set bmax $data(totmax)
    if {$data(reqmax) > $data(totmax)} { set bmax $data(reqmax) }
    $bc axis config y -min $bmin -max $bmax
    $bc element config bg_lo -ydata {$data(reqmin) $data(imgmin) $data(totmin)}
    $bc element config bg_hi -ydata {$bmax $bmax $bmax}
    $bc element config tot -ydata {$data(totmax)}
    $bc element config img -ydata {$data(imgmax)}
    $bc element config req -ydata {$data(reqmax)}
}

proc doContrast {dlg mode} {
    global imgparam
    set rangeprm(apply) $imgparam(applymode)
    set rangeprm(mode) $mode
    set rangeprm(vmin) [$dlg.f1.minreq get]
    set rangeprm(vmax) [$dlg.f1.maxreq get]
    toastimSetRange rangeprm
    toastimRasterImage $imgparam(cimg)
    UpdateCanvas
    Refresh
}

# =============================================================================
# "Colour" dialog

proc dlgColour {} {
    global colourmode
    set dlg .colourdlg
    if [winfo exists $dlg] {
	raise $dlg
    } else {
	toplevel $dlg
	wm title $dlg "Colour"
	wm resizable $dlg false false

	# button bar
	set f [frame $dlg.bbar]
	button $f.apply -text Apply -command [list doColour $dlg]
	button $f.done -text Done -command [list destroy $dlg]
	pack $f.apply $f.done -side left -padx 5
	pack $f -side bottom -pady 10

	set f [frame $dlg.f1]
	pack $f -side top -fill x -expand yes -padx 10 -pady 10
	radiobutton $f.r1 -text "Linear greyscale" -anchor w \
	    -variable colourmode(rgb) -value 0
	radiobutton $f.r2 -text "RGB colours from file" -anchor w \
	    -variable colourmode(rgb) -value 1
	pack $f.r1 $f.r2 -side top -fill x

	tixFileEntry $f.fe -dialogtype tixExFileSelectDialog -labelside top \
	    -label "Palette file" -variable colourmode(pal) \
	    -activatecmd [list callbackPaletteFile $f.fe]
	$f.fe subwidget entry config -width 30
	pack $f.fe -side top -padx 10 -pady 10 -fill x
    }
}

proc callbackPaletteFile {w} {
    global rc
    $w filedialog config -title "Load palette"
    set fsb [$w filedialog subwidget fsbox]
    if {[string length $rc(scaledir)]} {
	$fsb config -directory $rc(scaledir)
    }
    $fsb config -filetypes {{{*.pal} {Palette files}}}
    $fsb subwidget types pick 0
}

proc doColour {w} {
    global imgparam colourmode
    toastimSetColour colourmode
    toastimRasterImage $imgparam(cimg)
    UpdateCanvas
}

# =============================================================================
# "Mapping" dialog

proc dlgMapping {} {
    global imgparam

    set dlg .mappingdlg
    if [winfo exists $dlg] {
	raise $dlg
    } else {
	toplevel $dlg
	wm title $dlg "Mapping"
	wm resizable $dlg false false

	# button bar
	set f [frame $dlg.bbar]
	button $f.apply -text Apply -command [list doMapping $dlg]
	button $f.done -text Done -command [list destroy $dlg]
	pack $f.apply $f.done -side left -padx 5
	pack $f -side bottom -pady 10

	# linear/log
	set f [frame $dlg.f1]
	pack $f -side top -fill x -expand yes -padx 10 -pady 10
	radiobutton $f.r1 -text Linear -anchor w -value 0 \
	    -variable imgparam(logscale)
	radiobutton $f.r2 -text Logarithmic -anchor w -value 1\
	    -variable imgparam(logscale)
	pack $f.r1 $f.r2 -side top -fill x
    }
}

proc doMapping {w} {
    global imgparam
    toastimSetMapping imgparam
    toastimRasterImage $imgparam(cimg)
    UpdateCanvas
}

# =============================================================================
# "Profile" window

proc win_profile {} {
    global profparam imgw imgh
    set w .profilewin
    if [winfo exists $w] {
	raise $w
    } else {
	toplevel $w
	wm title $w "Profile"
	bind $w <Destroy> { DrawProfile 0 }

	set ornt [tixOptionMenu $w.ornt -label Orientation -labelside left \
		      -dynamicgeometry false]
	$ornt subwidget label config -width 10
	$ornt add command horz -label Horizontal
	$ornt add command vert -label Vertical
	$ornt add command oblq -label Oblique
	$ornt config -command callback_ProfileOrient
	set profparam(orient) horz
	$ornt config -variable profparam(orient)
	pack $ornt -side bottom

	set gf [toastGraph $w.prf]
	$gf.g legend config -hide yes
	toastGraph_setCallbackAdd $gf profileGraph_addSet
	toastGraphSetAdd $gf
	toastGraph_setCallbackAdd $gf ""
	$gf.g element config line0 -linewidth 1 -symbol ""
	bind $gf <Destroy> [list destroy $w]

	set imgw [image width rimg]
	set imgh [image height rimg]

	winProfile_update
    }
}

proc profileGraph_addSet {w s xvec yvec} {
    upvar #0 $xvec xv
    upvar #0 $yvec yv
    set xv [blt::vector create profilex]
    set yv [blt::vector create profiley]
}

proc DrawProfile {{draw 1}} {
    global profparam
    set c .imgwin.c
    if {[winfo exists $c] == 0} { return TCL_OK }
    if {$draw == 0} {
	$c delete profile
    } elseif {[$c find withtag profile] == ""} {
	set p [$c create line $profparam(x0) $profparam(y0) $profparam(x1) \
		   $profparam(y1) -arrow last -fill green]
	$c addtag profile withtag $p
	$c bind $p <Button-1> {callback_ClickProfile %x %y}
	$c bind $p <B1-Motion> {callback_DragProfile %x %y}
    } else {
	$c coords profile $profparam(x0) $profparam(y0) $profparam(x1) \
		   $profparam(y1)
    }
    toastimCalcProfile profparam profilex profiley
}

proc callback_ClickProfile {x y} {
    global profparam
    AdjustScreenCoord x y
    if {$profparam(orient) == "oblq"} {
	set dx [expr $profparam(x0) - $x]
	set dy [expr $profparam(y0) - $y]
	set d0 [expr sqrt ($dx * $dx + $dy * $dy)]
	set dx [expr $profparam(x1) - $x]
	set dy [expr $profparam(y1) - $y]
	set d1 [expr sqrt ($dx * $dx + $dy * $dy)]

	if {$d0 < $d1} {
	    set profparam(x0) $x
	    set profparam(y0) $y
	    set profparam(dragpoint) 0
	} else {
	    set profparam(x1) $x
	    set profparam(y1) $y
	    set profparam(dragpoint) 1
	}
    }
}

proc callback_DragProfile {x y} {
    global profparam
    AdjustScreenCoord x y
    switch $profparam(orient) {
	horz {set profparam(y0) $y; set profparam(y1) $y}
	vert {set profparam(x0) $x; set profparam(x1) $x}
	oblq {
	    if {$profparam(dragpoint) == 0} {
		set profparam(x0) $x; set profparam(y0) $y
	    } else {
		set profparam(x1) $x; set profparam(y1) $y
	    }
	}
    }
    DrawProfile
}

proc callback_ProfileOrient {which} {
    global profparam
    switch $which {
	horz {set profparam(x0) 0; set profparam(x1) [image width rimg]; \
	      set profparam(y0) [expr [image height rimg]/2]; \
	      set profparam(y1) $profparam(y0)}
	vert {set profparam(x0) [expr [image width rimg]/2]; \
	      set profparam(x1) $profparam(x0); \
	      set profparam(y0) [image height rimg]; set profparam(y1) 0}
	oblq {set profparam(x0) 0; set profparam(x1) [image width rimg]; \
	      set profparam(y0) 0; set profparam(y1) [image height rimg]}
    }
    DrawProfile
}

proc winProfile_update {} {
    global imgparam profparam range imgw imgh
    if {[winfo exists .profilewin] == 0} { return TCL_OK }
    if {$imgparam(dim) == 2 || $imgparam(3dmode) == 0} {
	if {[image width rimg] != $imgw || [image height rimg] != $imgh} {
	    set profparam(x0) [expr ($profparam(x0)*[image width rimg])/$imgw]
	    set profparam(x1) [expr ($profparam(x1)*[image width rimg])/$imgw]
	    set profparam(y0) [expr ($profparam(y0)*[image height rimg])/$imgh]
	    set profparam(y1) [expr ($profparam(y1)*[image height rimg])/$imgh]
	    set imgw [image width rimg]
	    set imgh [image height rimg]
	}
	DrawProfile
    } else {
	DrawProfile 0
    }
    toastimGetRange range
    set prf .profilewin.prf.g
    $prf axis config y -min $range(imgmin) -max $range(imgmax)
}

# =============================================================================
# "Single cross section" dialog

proc dlg_single3d {} {
    set dlg .single3d
    if [winfo exists $dlg] {
	raise $dlg
    } else {
	toplevel $dlg
	wm title $dlg "3-D: Single cross section"
	wm resizable $dlg false false

	# button bar
	set f [frame $dlg.bbar]
	button $f.done -text Done -command [list destroy $dlg]
	pack $f.done -side left -padx 5
	pack $f -side bottom -pady 10

	set f [frame $dlg.f1]
	pack $f -side top -fill x -expand yes -padx 10 -pady 10
	set ont [tixSelect $f.orient -label "Cross-section orientation" \
		     -labelside top -radio yes -allowzero no \
		     -orient horizontal -command callback_single3d_orient]
	$ont add z -text xy-plane
	$ont add y -text xz-plane
	$ont add x -text yz-plane
	pack $ont -side top

	set f [frame $dlg.f2]
	pack $f -side top -fill x -expand yes -padx 10 -pady 10
	set info [label $f.info -anchor w]
	set cs [scale $f.cs -orient horizontal -from 1 \
		    -showvalue no -command callback_single3d_plane]
	pack $info $cs -side top -fill x
    }
    single3d_update
}

proc callback_single3d_orient {d s} {
    # change cross section orientation
    global imgparam
    if {$s == 0} { return TCL_OK }
    if {$d != $imgparam(orient3d)} {
	set imgparam(orient3d) $d
	set maxpln $imgparam(grid$d)
	set pln $imgparam(plane$d)
	.single3d.f2.cs config -to $maxpln
	.single3d.f2.cs set [incr pln]
	single3d_update_info $d
	toastimSet3Dparam imgparam
	UpdateCanvas
	ImageResized
    }
}

proc callback_single3d_plane {v} {
    # change cross section plane
    global imgparam
    incr v -1
    set d $imgparam(orient3d)
    if {$v == $imgparam(plane$d)} { return TCL_OK }
    set imgparam(plane$d) $v
    single3d_update_info $d
    toastimSet3Dparam imgparam
    UpdateCanvas
    ImageChanged
}

proc single3d_update {} {
    # update dialog controls from data
    global imgparam
    set d $imgparam(orient3d)
    set pln $imgparam(plane$d)
    set maxpln $imgparam(grid$d)
    .single3d.f1.orient config -value $imgparam(orient3d)
    .single3d.f2.cs config -to $maxpln
    .single3d.f2.cs set [incr pln]
    single3d_update_info $d
}

proc single3d_update_info {d} {
    global imgparam sbbox
    set pln $imgparam(plane$d)
    incr pln
    set maxpln $imgparam(grid$d)
    set level [expr ($sbbox(max$d) - $sbbox(min$d))/($maxpln - 1)*($pln - 1) \
		   + $sbbox(min$d)]
    .single3d.f2.info config -text \
	"Slice $pln of $maxpln ($d = [format "%+0.3g" $level])"
}

# =============================================================================
# "Multiple cross sections" dialog

proc dlg_multi3d {} {
    set dlg .multi3d
    if [winfo exists $dlg] {
	raise $dlg
    } else {
	toplevel $dlg
	wm title $dlg "3-D Multiple cross sections"
	wm resizable $dlg false false

	# the plane sliders
	set f [frame $dlg.f1]
	pack $f -side top -fill x -expand yes -padx 10 -pady 10
	foreach d {z y x} {
	    pack [label $f.info$d -anchor w] -side top -fill x
	    pack [scale $f.cs$d -orient horizontal -length 200 -from 1 \
		-showvalue no -command [list callback_multiplane $d 0]] \
		-side top -fill x
	}
    }
    multi3d_update
}

proc callback_multiplane {which setslider v} {
    global imgparam
    incr v -1
    if {$v < 0} {
	set v 0
    } elseif {$v >= $imgparam(grid$which)} {
	set v [expr $imgparam(grid$which)-1]
    }
    if {$v == $imgparam(plane$which)} { return TCL_OK }
    set imgparam(plane$which) $v
    multi3d_update_info setslider
    toastimSet3Dparam imgparam
    UpdateCanvas
}

proc multi3d_update {} {
    global imgparam
    foreach which {x y z} {
	set pln $imgparam(plane$which)
	incr pln
	.multi3d.f1.cs$which config -to $imgparam(grid$which)
	.multi3d.f1.cs$which set $pln
    }
    multi3d_update_info 0
}

proc multi3d_update_info {setslider} {
    global imgparam sbbox
    for {set x 0} {$x < 3} {incr x} {
	set plane [lindex {xy xz yz} $x]
	set d [lindex {z y x} $x]
	set pln $imgparam(plane$d)
	incr pln
	set maxpln $imgparam(grid$d)
	set level [expr ($sbbox(max$d) - $sbbox(min$d))/($maxpln - 1) * \
		       ($pln - 1) + $sbbox(min$d)]
	.multi3d.f1.info$d config -text \
	    "$plane slice $pln of $maxpln ($d = [format "%+0.3g" $level])"
	if {$setslider != 0} {.multi3d.f1.cs$d set $pln}
    }
}

# =============================================================================
# "Orthogonal cross sections" dialog

proc dlg_ortho3d {} {
    global imgparam
    set dlg .ortho3d
    if [winfo exists $dlg] {
	raise $dlg
    } else {
	toplevel $dlg
	wm title $dlg "3-D: Orthogonal cross sections"
	wm resizable $dlg false false

	# button bar
	set f [frame $dlg.bbar]
	button $f.done -text Done -command [list destroy $dlg]
	pack $f.done -side left -padx 5
	pack $f -side bottom -pady 10

	# the plane sliders
	set f [frame $dlg.f1]
	pack $f -side top -fill x -expand yes -padx 10 -pady 10
	foreach d {z y x} {
	    pack [label $f.info$d -anchor w] -side top -fill x
	    pack [scale $f.cs$d -orient horizontal -length 200 -from 1 \
		-showvalue no -command [list callback_orthoplane $d]] \
		-side top -fill x
	}

	Line $dlg.line1

	# pack rotation angle sliders
	pack [frame_3d_viewangle $dlg.rot] -side top -fill x -expand yes \
	    -padx 10 -pady 10

	Line $dlg.line2

	# pack transparency slider
	pack [frame_transparency $dlg.trans] -side top -fill x -expand yes \
	    -padx 10 -pady 10
    }
    ortho3d_update
}

proc callback_orthoplane {which v} {
    global imgparam
    incr v -1
    if {$v == $imgparam(plane$which)} { return TCL_OK }
    set imgparam(plane$which) $v
    ortho3d_update_info
    toastimSet3Dparam imgparam
    UpdateCanvas
}

proc ortho3d_update {} {
    global imgparam
    foreach which {x y z} {
	set pln $imgparam(plane$which)
	incr pln
	.ortho3d.f1.cs$which config -to $imgparam(grid$which)
	.ortho3d.f1.cs$which set $pln
    }
    ortho3d_update_info
}

proc ortho3d_update_info {} {
    global imgparam sbbox
    for {set x 0} {$x < 3} {incr x} {
	set plane [lindex {xy xz yz} $x]
	set d [lindex {z y x} $x]
	set pln $imgparam(plane$d)
	incr pln
	set maxpln $imgparam(grid$d)
	set level [expr ($sbbox(max$d) - $sbbox(min$d))/($maxpln - 1) * \
		       ($pln - 1) + $sbbox(min$d)]
	.ortho3d.f1.info$d config -text \
	    "$plane slice $pln of $maxpln ($d = [format "%+0.3g" $level])"
    }
}

# =============================================================================
# "Multislice stack" dialog

proc dlg_mslice3d {} {
    set dlg .mslice3d
    if [winfo exists $dlg] {
	raise $dlg
    } else {
	toplevel $dlg
	wm title $dlg "3-D: Multislice stack"
	wm resizable $dlg false false

	# button bar
	set f [frame $dlg.bbar]
	button $f.done -text Done -command [list destroy $dlg]
	pack $f.done -side left -padx 5
	pack $f -side bottom -pady 10

	# pack rotation angle sliders
	pack [frame_3d_viewangle $dlg.rot] -side top -fill x -expand yes \
	    -padx 10 -pady 10

	Line $dlg.line1

	# pack transparency slider
	pack [frame_transparency $dlg.trans] -side top -fill x -expand yes \
	    -padx 10 -pady 10
    }
}

# =============================================================================
# 3-D view rotation angle package

proc frame_3d_viewangle {f} {
    global imgparam
    frame $f
    pack [label $f.title -text "View angles"] -side top
    pack [label $f.lphi -anchor w -text "Phi = $imgparam(rphi)"] \
	-side top -fill x
    pack [scale $f.phi -orient horizontal -length 200 -from 0 -to 360 \
	-showvalue no -command [list callback_3d_viewangle $f rphi]] \
	-side top -fill x
    $f.phi set $imgparam(rphi)
    pack [label $f.ltheta -anchor w -text "Theta = $imgparam(rtheta)"] \
	-side top -fill x
    pack [scale $f.theta -orient horizontal -length 200 -from -90 -to 90 \
	-showvalue no -command [list callback_3d_viewangle $f rtheta]] \
	-side top -fill x
    $f.theta set $imgparam(rtheta)
    return $f
}

proc callback_3d_viewangle {f which v} {
    global imgparam
    if {$v == $imgparam($which)} { return TCL_OK }
    set imgparam($which) $v
    $f.lphi config -text "Phi = $imgparam(rphi)"
    $f.ltheta config -text "Theta = $imgparam(rtheta)"
    toastimSet3Dparam imgparam
    UpdateCanvas
}

# =============================================================================
# Transparency package

proc frame_transparency {f} {
    global imgparam
    frame $f
    pack [label $f.ltrans -anchor w -text \
	      "Transparency = $imgparam(transparency)"] -side top -fill x
    pack [scale $f.trans -orient horizontal -length 200 -from 0 -to 100 \
	-showvalue no -command [list callback_transparency $f]] \
	-side top -fill x
    $f.trans set $imgparam(transparency)
    return $f
}

proc callback_transparency {f v} {
    global imgparam
    if {$v == $imgparam(transparency)} { return TCL_OK }
    set imgparam(transparency) $v
    $f.ltrans config -text "Transparency = $imgparam(transparency)"
    toastimSet3Dparam imgparam
    UpdateCanvas
}

# =============================================================================
# 3D mode change notification

proc callbackChange3Dmode {mode} {
    global imgparam rimgtag rimgxytag rimgxztag rimgyztag
    if {$mode == $imgparam(3dmode)} { return TCL_OK }
    set c .imgwin.c

    if {$imgparam(3dmode) == 5} {
	$c delete multi3dtag
	image delete rimgxy
	image delete rimgxz
	image delete rimgyz
	set rimgtag [$c create image 0 0 -anchor nw -image rimg]
    }

    set imgparam(3dmode) $mode
    show3dcontrols $mode

    # set mouse bindings
    if {$imgparam(3dmode) == 0} {
	bind .imgwin <Motion> {callback_motion_3d %x %y}
    } else {
	bind .imgwin <Motion> {}
    }

    toastimSet3Dparam imgparam

    if {$imgparam(3dmode) == 5} {
	$c delete rimg
	image create photo rimgxy
	image create photo rimgxz
	image create photo rimgyz
	set rimgxytag [$c create image 0 0 -anchor nw -image rimgxy]
	set rimgxztag [$c create image 0 0 -anchor nw -image rimgxz]
	set rimgyztag [$c create image 0 0 -anchor nw -image rimgyz]

	# coordinate indicators
	$c addtag xyaxes withtag [$c create line \
				      0 -20 0 0 20 0 -arrow both -fill green]
	$c addtag xyaxes withtag [$c create text \
				      20 0 -anchor w -fill green -text x]
	$c addtag xyaxes withtag [$c create text \
				      0 -20 -anchor s -fill green -text y]
	$c addtag xzaxes withtag [$c create line \
				      0 -20 0 0 20 0 -arrow both -fill green]
	$c addtag xzaxes withtag [$c create text \
				      20 0 -anchor w -fill green -text x]
	$c addtag xzaxes withtag [$c create text \
				      0 -20 -anchor s -fill green -text z]
	$c addtag yzaxes withtag [$c create line \
				      0 -20 0 0 20 0 -arrow both -fill green]
	$c addtag yzaxes withtag [$c create text \
				      20 0 -anchor w -fill green -text y]
	$c addtag yzaxes withtag [$c create text \
				      0 -20 -anchor s -fill green -text z]

	# cross section indicator lines
	$c addtag xycshorz withtag [$c create line 0 0 100 0 -fill orange]
	$c addtag xycsvert withtag [$c create line 0 0 0 100 -fill orange]
	$c addtag xzcshorz withtag [$c create line 0 0 100 0 -fill orange]
	$c addtag xzcsvert withtag [$c create line 0 0 0 100 -fill orange]
	$c addtag yzcshorz withtag [$c create line 0 0 100 0 -fill orange]
	$c addtag yzcsvert withtag [$c create line 0 0 0 100 -fill orange]

	$c addtag xyelements withtag $rimgxytag
	$c addtag xyelements withtag xyaxes
	$c addtag xyelements withtag xycshorz
	$c addtag xyelements withtag xycsvert
	$c addtag xzelements withtag $rimgxztag
	$c addtag xzelements withtag xzaxes
	$c addtag xzelements withtag xzcshorz
	$c addtag xzelements withtag xzcsvert
	$c addtag yzelements withtag $rimgyztag
	$c addtag yzelements withtag yzaxes
	$c addtag yzelements withtag yzcshorz
	$c addtag yzelements withtag yzcsvert

	$c bind xyelements <Button-1> {callback_DragCS %x %y xy}
	$c bind xzelements <Button-1> {callback_DragCS %x %y xz}
	$c bind yzelements <Button-1> {callback_DragCS %x %y yz}
	$c bind xyelements <B1-Motion> {callback_DragCS %x %y xy}
	$c bind xzelements <B1-Motion> {callback_DragCS %x %y xz}
	$c bind yzelements <B1-Motion> {callback_DragCS %x %y yz}

	$c addtag multi3dtag all
    }

    UpdateCanvas
    Refresh
    ImageResized
}

proc callback_DragCS {x y plane} {
    global imgparam rimgxytag rimgxztag rimgyztag

    AdjustScreenCoord x y
    set c .imgwin.c
    set x0 [image width rimgyz]
    set y0 [image height rimgxy]

    switch $plane {
	xy {
	    set crd [$c coords $rimgxytag]
	    set w [image width rimgxy]
	    set h [image height rimgxy]
	}
	xz {
	    set crd [$c coords $rimgxztag]
	    set w [image width rimgxz]
	    set h [image height rimgxz]
	}
	yz {
	    set crd [$c coords $rimgyztag]
	    set w [image width rimgyz]
	    set h [image height rimgyz]
	}
    }

    set x [expr ($x - [lindex $crd 0]) / $w]
    set y [expr ($y - [lindex $crd 1]) / $h]
    if {$x < 0} {
	set x 0
    } elseif {$x > 1} {
	set x 1
    }
    if {$y < 0} {
	set y 0
    } elseif {$y > 1} {
	set y 1
    }

    switch $plane {
	xy {
	    callback_multiplane x 1 [expr int ($x * $imgparam(gridx))+1]
	    callback_multiplane y 1 [expr int ((1 - $y) * $imgparam(gridy))+1]
	}
	xz {
	    callback_multiplane x 1 [expr int ($x * $imgparam(gridx))+1]
	    callback_multiplane z 1 [expr int ((1 - $y) * $imgparam(gridz))+1]
	}
	yz {
	    callback_multiplane y 1 [expr int ($x * $imgparam(gridy))+1]
	    callback_multiplane z 1 [expr int ((1 - $y) * $imgparam(gridz))+1]
	}
    }
}

proc hide3dcontrols {} {
    set dlgs [list single3d ortho3d]
    foreach d $dlgs {
	if [winfo exists .$d] { destroy .$d }
    }
}

proc show3dcontrols {mode} {
    set dlgs [list single3d ortho3d mslice3d ortho3d ortho3d multi3d]
    set cdlg [lindex $dlgs $mode]
    foreach d $dlgs {
	if {$d == $cdlg} {
	    dlg_$d
	} elseif [winfo exists .$d] { destroy .$d }
    }
}

# =============================================================================
# image change notification

proc callbackPickImage {lb} {
    if {$lb != 0} {
	set idx [$lb index active]
	incr idx
	callbackChangeImage $idx
    }
    return TCL_OK
}

proc callbackChangeImage {imgno} {
    global imgparam
    incr imgno -1
    if { $imgno != $imgparam(cimg) && $imgno < $imgparam(nimg) } {
	set imgparam(cimg) $imgno
	toastimRasterImage $imgno
	UpdateCanvas
	Refresh
    }
    return TCL_OK
}   

# =============================================================================
# Refresh all windows after something has changed

proc Refresh {} {
    global imgparam

    # grid raster string
    set grid "$imgparam(gridx) x $imgparam(gridy)"
    if {$imgparam(dim) > 2} {
	set grid [concat $grid "x $imgparam(gridz)"]
    }

    # refresh info window
    set p [.pw subwidget pane2]
    $p.f1.0.v config -text [file tail $imgparam(imgfile)]
    $p.f1.1.v config -text $imgparam(imgfmt)
    $p.f1.2.v config -text $imgparam(dim)
    $p.f1.3.v config -text $imgparam(nimg)
    $p.f1.4.v config -text [file tail $imgparam(mesh)]
    $p.f1.5.v config -text $imgparam(mbbox)
    $p.f1.6.v config -text $imgparam(sbbox)
    $p.f1.7.v config -text $grid
    $p.f1.8.v config -text "$imgparam(totmin)  -  $imgparam(totmax)"

    toastimGetRange range
    $p.f2.0.v config -text "$range(imgmin)  -  $range(imgmax)"
    $p.f2.1.v config -text "$range(reqmin)  -  $range(reqmax)"

    # resample dialog
    if [winfo exists .resampledlg] { initResampleImg }

    # 3d dialog
    if [winfo exists .single3d] { single3d_update }
    if [winfo exists .multi3d] { multi3d_update }

    # contrast dialog
    dlgContrast_update

    # profile window
    winProfile_update
}

# =============================================================================

proc callback_motion_2d {x y} {
    global imgparam sbbox rc
    AdjustScreenCoord x y
    set x [expr $x / $rc(init_zoom_2d)]
    set y [expr $imgparam(gridy) - $y / $rc(init_zoom_2d)]
    set xcoord [expr ($sbbox(maxx) - $sbbox(minx)) * $x / $imgparam(gridx) \
		    + $sbbox(minx)]
    set ycoord [expr ($sbbox(maxy) - $sbbox(miny)) * $y / $imgparam(gridy) \
		    + $sbbox(miny)]
    set p [.pw subwidget pane2]
    $p.f2.2.v config -text \
	"x=[format "%0.2f" $xcoord], y=[format "%0.2f" $ycoord]"
    $p.f2.3.v config -text [toastimGetValueFromPixel $x $y]
}

proc callback_motion_3d {x y} {
    AdjustScreenCoord x y
    set p [.pw subwidget pane2]
    $p.f2.2.v config -text "x=$x, y=$y"
}

# =============================================================================
# main window

# ==== menu bar ====
MenuSetup .mnu
menubutton .mnu.file -text File -underline 0 -menu .mnu.file.menu
menu .mnu.file.menu
.mnu.file.menu add command -label "Load image ..." -underline 0 \
    -command dlgLoadImg
.mnu.file.menu add command -label "Export image ..." -underline 0 \
    -state disabled \
    -command {toastExportPixmapDialog $colourmode(rgb) toastimWritePixmap}
.mnu.file.menu add command -label "Export image data ..." -underline 1 \
    -command {toastimWritePixmapData}
.mnu.file.menu add separator
.mnu.file.menu add command -label "Exit" -underline 1 -command exit

menubutton .mnu.list -text List -underline 0 -state disabled \
    -menu .mnu.list.menu
set m [menu .mnu.list.menu]
$m add cascade -label "Delete" -menu $m.del

set sm [menu $m.del]
$sm add command -label "Current" -underline 0 \
    -command {delListEntry $imgparam(cimg)}
$sm add command -label "All before current" -underline 0 \
    -command {delListBefore $imgparam(cimg)}
$sm add command -label "All after current" -underline 0 \
    -command {delListAfter $imgparam(cimg)}

$m add cascade -label "Add" -menu $m.add

set sm [menu $m.add]
$sm add command -label "Transform ..." -underline 0 \
    -command dlgAddTransform

menubutton .mnu.image -text Image -underline 0 -state disabled \
    -menu .mnu.image.menu
menu .mnu.image.menu
.mnu.image.menu add command -label "Window" -underline 0 \
    -command showImage
.mnu.image.menu add separator
.mnu.image.menu add command -label "Resample ..." -underline 0 \
    -command dlgResampleImg
.mnu.image.menu add command -label "Mapping ..." -underline 0 \
    -command dlgMapping
.mnu.image.menu add command -label "Contrast ..." -underline 0 \
    -command dlgContrast
.mnu.image.menu add command -label "Colour ..." -underline 1 \
    -command dlgColour
.mnu.image.menu add command -label "Contours ..." -underline 2 \
    -command dlgContours

menubutton .mnu.3d -text "3-D" -underline 0 -state disabled -menu .mnu.3d.menu
set m [menu .mnu.3d.menu]
$m add cascade -label "Mode" -menu $m.mode
$m add separator
$m add command -label "Controls ..." -underline 0 \
    -command {show3dcontrols $imgparam(3dmode)}

set sm [menu $m.mode]
$sm add radiobutton -label "Single cross section" -value 0 \
    -variable tmp_3dmode -command [list callbackChange3Dmode 0]
$sm add radiobutton -label "Multiple cross sections" -value 5 \
    -variable tmp_3dmode -command [list callbackChange3Dmode 5]
$sm add radiobutton -label "Orthogonal cross sections" -value 1 \
    -variable tmp_3dmode -command [list callbackChange3Dmode 1]
$sm add radiobutton -label "Multislice stack" -value 2 \
    -variable tmp_3dmode -command [list callbackChange3Dmode 2]
$sm add radiobutton -label "Surface" -value 3 \
    -variable tmp_3dmode -command [list callbackChange3Dmode 3]
$sm add radiobutton -label "Mapped surface" -value 4 \
    -variable tmp_3dmode -command [list callbackChange3Dmode 4]

menubutton .mnu.anls -text Analysis -underline 0 -state disabled \
    -menu .mnu.anls.menu
set m [menu .mnu.anls.menu]
$m add command -label "Profile ..." -underline 0 \
    -command win_profile

menubutton .mnu.help -text Help -underline 0 -menu .mnu.help.menu
set m [menu .mnu.help.menu]
$m add command -label "Contents ..." -underline 0 \
    -command [list htmldoc toastim index.html]

#tixControl .mnu.cimg -label "Image" -integer yes -min 0 -max 0 \
#    -command [list callbackChangeImage] -options { entry.width 5 }

pack .mnu.file .mnu.list .mnu.image .mnu.3d .mnu.anls .mnu.help -side left
#pack .mnu.cimg -side right -padx 5

# ==== status bar ====
statusSetup

# ==== client area ====

set m [tixPanedWindow .pw -dynamicgeometry true -orient horizontal \
	   -panerelief flat -panebd 1]
pack $m -fill both -expand true
$m add pane2 -size 350 -expand 1
$m add pane1 -size 75 -expand 0
set p1 [$m subwidget pane1]
set p2 [$m subwidget pane2]

set f [frame $p2.f1]
pack $f -side top -fill x -expand yes -padx 10 -pady 10
label $f.l -text "Image file information" -anchor w
pack  $f.l -side top -fill x

set infolabel [list "Image file" "Image format" "Dimension" "Image count" \
	"Mesh file" "Mesh BBox" "Image BBox" "Raster grid" "Total range"]
for {set x 0} {$x < 9} {incr x} {
    set f1 [frame $f.$x]
    label $f1.l -text "[lindex $infolabel $x]  " -width 15 -anchor e
    label $f1.v -text "-" -width 35 -anchor w -foreground blue
    pack  $f1.l -side left
    pack  $f1.v -side left -fill x -expand true
    pack $f1 -side top -fill x
}

Line $p2.line1

set f [frame $p2.f2]
pack $f -fill x -expand yes -padx 10 -pady 10
label $f.l -text "Current image information" -anchor w
pack  $f.l -side top -fill x

set infolabel [list "Data range" "Greyscale range" "Cursor pos" "Cursor value"]
for {set x 0} {$x < 4} {incr x} {
    set f2 [frame $f.$x]
    label $f2.l -text "[lindex $infolabel $x]  " -width 15 -anchor e
    label $f2.v -text "-" -width 35 -anchor w -foreground blue
    pack  $f2.l -side left
    pack  $f2.v -side left -fill x -expand true
    pack  $f2 -side top -fill x
}

# === create the image list box ===

set dl [tixScrolledListBox $p1.dlist -scrollbar auto]
pack $dl -side top -fill both -expand true -padx 5 -pady 5
$dl config -browsecmd "callbackPickImage [$dl subwidget listbox]"

# ==== create an image handle
image create photo rimg

# ==== default palette ====
set colourmode(pal) [file join $rc(scaledir) $rc(default_scale)]

# ==== load dialog customisation ===
set fdlg [tix filedialog tixExFileSelectDialog]
set f [frame $fdlg.custom]
checkbutton $f.reraster -text "Reset bounding box and grid" \
    -variable loadparam(reraster)
checkbutton $f.rescale -text "Rescale image" \
    -variable loadparam(rescale)
pack $f.reraster -side top -anchor w
pack $f.rescale -side top -anchor w
pack $f -side bottom -fill x -padx 10 -pady 5

# ==== read initial image file ====
if {[string length $imgfile] > 0} {
    doLoadImg $imgfile
}
