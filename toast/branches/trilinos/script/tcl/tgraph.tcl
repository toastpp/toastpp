# =============================================================================
# TOAST v.15                                          (c) Martin Schweiger 1999
# tgraph (Tcl script)
# Interface to the BLT graph widget
# Implements a generic customizable graph window
# =============================================================================

# The following functions are exported:
# toastGraph {w}
#    Generate a new graph in frame w
#
# toastGraph_setCallbackAdd {w func}
#    Install callback function func to be called when addition of a new
#    data set (element) is requested. w is the graph window passed to
#    toastGraph
#    func should have the following interface: callback {w s xvec yvec}
#    where w is the graph frame, s is the number for the new set, and xvec
#    and yvec are the names of the BLT vectors to receive the data for the
#    new set
#    If func = "" then the existing add method is disabled.
#    If this function is not invoked then the interactive addition of sets
#    is disabled
#
# toastGraph_setCallbackDel {w func}
#    Install callback function func to be called when the user request
#    to delete an element. w is the graph window passed to toastGraph
#    func should have the following interface: callback {w s}
#    where w is the graph frame and s is the number of the set to be deleted
#    If this function is not invoked then the interactive deletion of sets
#    is disabled
#
# toastGraph_setCallbackSet {w func}
#    Install callback function func to be called when the user selects a
#    different element. w is the graph window passed to toastGraph
#    func should have the following interface: callback {w s}
#    where w is the graph frame and s is the newly selected set.
#    If this function is not invoked then the interactive change of sets
#    is disabled

proc toastGraph {w} {
    # This is a unique identifier for the graph
    regsub -all {\.} $w x id

    # The list of sets (elements) in the graph
    global slist$id
    set slist$id {}

    # Create the plot's main frame
    set gf [frame $w]
    pack $gf -fill both -expand yes

    # Create the graph menu
    MenuSetup $gf.mnu
    menubutton $gf.mnu.graph -text Graph -underline 0 -menu $gf.mnu.graph.menu
    set m [menu $gf.mnu.graph.menu]
    $m add command -label "Configure ..." -underline 1 \
	    -command [list toastConfigGraph $gf.g]
    $m add command -label "Save Postscript ..." -underline 0 \
	    -command [list toastGraphPSOut $gf.g]
    $m add command -label "Export data ..." -underline 0 \
	    -command [list toastExportData $gf.g]
    $m add command -label "Close" -underline 0 -command [list destroy $w]

    menubutton $gf.mnu.axis -text Axis -underline 0 -menu $gf.mnu.axis.menu
    set m [menu $gf.mnu.axis.menu]
    $m add cascade -label "Configure" -menu $m.conf
    set sm [menu $m.conf]
    $sm add command -label "x-axis 1" \
	    -command [list toastConfigGraphAxis $gf.g x 1]
    $sm add command -label "y-axis 1" \
	    -command [list toastConfigGraphAxis $gf.g y 1]
    $sm add command -label "x-axis 2" \
	    -command [list toastConfigGraphAxis $gf.g x 2]
    $sm add command -label "y-axis 2" \
	    -command [list toastConfigGraphAxis $gf.g y 2]

    menubutton $gf.mnu.set -text Set -underline 0 -menu $gf.mnu.set.menu
    set m [menu $gf.mnu.set.menu]
    $m add command -label Add -command [list toastGraphSetAdd $gf] \
	-state disabled
    $m add command -label Delete -command [list toastGraphSetDel $gf $id] \
	-state disabled
    $m add command -label "Configure ..." \
	    -command [list toastConfigGraphSet $gf]
    tixControl $w.mnu.cset -allowempty false -integer true -min 0 -value 0 \
	-options { entry.width 4 } -command [list toastGraphSetSet $w $id] \
	-state disabled

    pack $gf.mnu.graph $gf.mnu.axis $gf.mnu.set $gf.mnu.cset -side left -fill x

    blt::graph $gf.g -plotbackground white -plotrelief flat
    pack $gf.g -side top -fill both -expand yes

    return $gf
}

# ============================================================================
# Dialog for general graph customisation

proc toastConfigGraph {g} {
    global relieftype ttl titlefont ta
    set dname .toastgraph_config
    if {[winfo exists $dname] == 0} {
	set w [toplevel $dname]
	wm title $w "Configure graph"

	set f [frame $w.bbar]
	button $f.apply -text Apply
	button $f.done -text Done -command [list destroy $w]
	pack $f.apply $f.done -side left -padx 5 -pady 5
	pack $f -side bottom

	set f [frame $w.0 -relief groove -bd 2]
	pack $f -padx 5 -pady 5 -fill x
	tixLabelEntry $f.lb -label "Title" -labelside top
	$f.lb subwidget entry config -textvariable ttl
	button $f.fb -text "Font ..." -command [list selectFont titlefont]
	colorChooser $f.tcc -label "Colour" -options { entry.width 8 }
	pack $f.lb -side top -fill x -padx 2 -pady 2
	pack $f.fb $f.tcc -side left -fill x -padx 2 -pady 2

	set f [frame $w.1 -relief groove -bd 2]
	pack $f -padx 5 -pady 5 -fill x

	tixComboBox $f.rs -dropdown yes -editable no -history no \
		-label "Plot relief" -labelside left -options {
	    label.width 13
	    entry.width 10
	    slistbox.scrollbar auto
	    slistbox.width 100
	    slistbox.height 100
	}
	for {set x 0} {$x < 6} {incr x} {
	    $f.rs insert end [lindex $relieftype $x]
	}
	pack $f.rs -side top -fill x -padx 2

	pack [colorChooser $f.pcc -label "Plot colour" -options {
	    label.width 13
	    entry.width 10
	}] -side top -fill x -padx 2
	pack [colorChooser $f.fcc -label "Frame colour" -options {
	    label.width 13
	    entry.width 10
	}] -side top -fill x -padx 2

	checkbutton $f.trans -text "transpose axes" -variable ta
	pack $f.trans -side top

    } else {
	set w $dname
	raise $w
    }
    $w.bbar.apply config -command  [list applyConfigGraph $w $g]
    $w.1.rs config -value [$g cget -plotrelief]
    set ttl [$g cget -title]
    set titlefont [$g cget -font]
    colorChooser_setColor $w.0.tcc [$g cget -foreground]
    set ta [$g cget -invertxy]
    colorChooser_setColor $w.1.pcc [$g cget -plotbackground]
    colorChooser_setColor $w.1.fcc [$g cget -background]
    return TCL_OK
}

proc applyConfigGraph {w g} {
    global ttl titlefont global ta
    $g config -title $ttl \
	    -font $titlefont \
	    -foreground [colorChooser_getColor $w.0.tcc] \
	    -plotbackground [colorChooser_getColor $w.1.pcc] \
	    -background [colorChooser_getColor $w.1.fcc] \
	    -plotrelief [$w.1.rs cget -value] \
	    -invertxy $ta
    return TCL_OK
}

# ============================================================================
# Configure axis a n of graph g, where 'a' is x or y, and 'n' is 1 or 2

proc toastConfigGraphAxis {g a n} {
    global axislabel axisfont rangemode rangemin rangemax axisiv axissd \
	    axislw showtick tickfont ticklength
    set dname .toastgraph_configaxis
    if {[winfo exists $dname] == 0} {
	set w [toplevel $dname]
	set f [frame $w.bbar]
	button $f.apply -text Apply
	button $f.done -text Done -command [list destroy $w]
	pack $f.apply $f.done -side left -padx 5 -pady 5
	pack $f -side bottom

	set f [frame $w.1 -relief groove -bd 2]
	pack $f -padx 5 -pady 5 -fill x
	tixLabelEntry $f.lb -label "Axis label" -labelside top
	$f.lb subwidget entry config -textvariable axislabel
	button $f.fb -text "Font ..." -command [list selectFont axisfont]
	colorChooser $f.tcc -label "Colour" -options { entry.width 8 }
	pack $f.lb -side top -fill x -padx 2 -pady 2
	pack $f.fb $f.tcc -side left -fill x -padx 2 -pady 2

	set f [frame $w.3 -relief groove -bd 2]
	pack $f -padx 5 -pady 5 -fill x
	pack [frame $f.0] -side top -fill x
	pack [label $f.0.l -text Ticks -width 12] -side left
	pack [checkbutton $f.0.c -text Show -variable showtick] \
		-side left -fill x

	frame $f.1
	pack [tixControl $f.1.iv -label "Interval" -min 0 -variable axisiv \
		-integer false -options {
	    entry.width 10
	    label.width 12
	}] -side left
	frame $f.2
	pack [tixControl $f.2.sd -label "Subdivisions" -min 1 \
		-variable axissd -integer true -options {
	    entry.width 10
	    label.width 12
	}] -side left
	frame $f.3
	pack [tixControl $f.3.lw -label "Line width" -min 0 \
		-variable axislw -integer true -options {
	    entry.width 10
	    label.width 12
	}] -side left
	frame $f.4
	pack [tixControl $f.4.tl -label "Tick length" -min 0 \
		-variable ticklength -integer true -options {
	    entry.width 10
	    label.width 12
	}] -side left
	pack $f.1 $f.2 $f.3 $f.4 -side top -fill x

	button $f.fb -text "Font ..." -command [list selectFont tickfont]
	colorChooser $f.tcc -label "Colour" -options { entry.width 8 }
	pack $f.fb $f.tcc -side left -fill x -padx 2 -pady 2

	set f [frame $w.2 -relief groove -bd 2]
	pack $f -padx 5 -pady 5 -fill x
	label $f.l -text "Range" -anchor w
	pack $f.l -side top -fill x
	radiobutton $f.r1 -text "Auto (tight)" -variable rangemode -value 0 \
	    -command [list changeRangeMode $w rangemode] -anchor w
	radiobutton $f.r2 -text "Auto (loose)" -variable rangemode -value 1 \
	    -command [list changeRangeMode $w rangemode] -anchor w
	pack $f.r1 $f.r2 -side top -fill x
	set fs [frame $f.1]
	radiobutton $fs.r3 -text "Manual" -variable rangemode -value 2 \
	    -command [list changeRangeMode $w rangemode] -anchor w
	tixLabelEntry $fs.min -label from -options {
	    entry.width 6
	}
	$fs.min subwidget entry config -textvariable rangemin
	tixLabelEntry $fs.max -label to -options {
	    entry.width 6
	}
	$fs.max subwidget entry config -textvariable rangemax
	pack $fs.r3 $fs.min $fs.max -side left
	pack $fs -side top -fill x
	tixOptionMenu $f.log -label "Scale" -labelside left -options {
	    label.width 12
	}
	$f.log add command 0 -label Linear
	$f.log add command 1 -label Logarithmic
	pack $f.log -side top -fill x

    } else {
	set w $dname
	raise $w
    }
    wm title $w "Configure $a axis $n"
    $w.bbar.apply config -command  [list applyConfigGraphAxis $w $g $a $n]
    if {$n == 1} {set n {}}
    set ax [join [list $a $n] {}]
    set axislabel [$g axis cget $ax -title]
    set axisfont [$g axis cget $ax -titlefont]
    set tickfont [$g axis cget $ax -tickfont]
    colorChooser_setColor $w.1.tcc [$g axis cget $ax -titlecolor]
    colorChooser_setColor $w.3.tcc [$g axis cget $ax -color]
    if {[$g axis cget $ax -min] == "" || [$g axis cget $ax -max] == ""} {
	if {[$g axis cget $ax -loose] == 0} {
	    set rangemode 0
	} else {
	    set rangemode 1
	}
    } else {
	set rangemode 2
    }
    set rangemin [$g axis cget $ax -min]
    set rangemax [$g axis cget $ax -max]
    changeRangeMode $w rangemode
    set showtick [$g axis cget $ax -showticks]
    set axisiv [$g axis cget $ax -stepsize]
    set axissd [$g axis cget $ax -subdivisions]
    set axislw [$g axis cget $ax -linewidth]
    set ticklength [$g axis cget $ax -ticklength]
    $w.2.log config -value [$g axis cget $ax -logscale]
    return TCL_OK
}

proc applyConfigGraphAxis {w g a n} {
    global axislabel axisfont rangemode rangemin rangemax axisiv axissd \
	    axislw showtick tickfont ticklength
    if {$n == 1} {set n {}}
    $w.3.1.iv update
    $w.3.2.sd update
    $w.3.3.lw update
    $w.3.4.tl update
    set ax [join [list $a $n] {}]
    set scl [$w.2.log cget -value]
    $g axis config $ax \
	-title $axislabel \
	-titlefont $axisfont \
	-titlecolor [colorChooser_getColor $w.1.tcc] \
	-color [colorChooser_getColor $w.3.tcc] \
	-logscale $scl \
	-stepsize $axisiv \
	-subdivisions $axissd \
	-linewidth $axislw \
	-tickfont $tickfont \
	-showticks $showtick \
	-ticklength $ticklength
    switch $rangemode {
	0 { $g axis config $ax -min "" -max "" -loose no }
	1 { $g axis config $ax -min "" -max "" -loose yes }
	2 { $g axis config $ax -min $rangemin -max $rangemax }
    }
    return TCL_OK
}

proc changeRangeMode {w rangevar} {
    upvar $rangevar md
    if {$md == 2} {
	$w.2.1.min config -state normal
	$w.2.1.max config -state normal
    } else {
	$w.2.1.min config -state disabled
	$w.2.1.max config -state disabled
    }
}

# ============================================================================
# File dialog for Postscript output of the graph

proc toastGraphPSOut {g} {
    set dialog [tix filedialog tixExFileSelectDialog]
    $dialog config -command [list doGraphPSOut $g] -title "Output PS"
    set dbox [$dialog subwidget fsbox]
    $dbox config -filetypes {{{*.ps *.eps} {Postscript files}}}
    $dbox subwidget types pick 0
    $dialog popup
    return TCL_OK
}

proc doGraphPSOut {g file} {
    $g postscript output $file -maxpect yes -decorations no
}

# ============================================================================
# File dialog for exporting data

proc toastExportData {g} {
    set dialog [tix filedialog tixExFileSelectDialog]
    $dialog config -command [list doExportData $g] -title "Export data"
    set dbox [$dialog subwidget fsbox]
    $dbox config -filetypes {{{*.dat} {Data files}}}
    $dbox subwidget types pick 0
    $dialog popup
    return TCL_OK
}

proc doExportData {g file} {
    puts "graph = $g"
    set cid [open $file w]
    foreach x [$g element names] {
	set xd [$g element cget $x -xdata]
	$xd variable local_xd
	set yd [$g element cget $x -ydata]
	$yd variable local_yd
	set len [$xd length]
	for {set i 0} {$i < $len} {incr i} {
	    puts $cid "$local_xd($i) $local_yd($i)"
	}
    }
    close $cid
}

# ============================================================================
# Add new graph element

proc toastGraph_setCallbackAdd {w func} {
    regsub -all {\.} $w x id
    global toastgraph_prm$id
    if {$func == ""} {
	$w.mnu.set.menu entryconfig 1 -state disabled
    } else {
	set callback [join [list "toastgraph_prm" $id "(addfunc)"] {}]
	set $callback $func
	$w.mnu.set.menu entryconfig 1 -state normal
    }
}

proc toastGraphSetAdd {w} {
    regsub -all {\.} $w x id
    global toastGraph_xvec toastGraph_yvec
    upvar #0 slist$id sl
    set s [$w.mnu.cset cget -value]
    if {[lsearch $sl $s] >= 0} {
	ErrorNotice "Set is already present"
	return TCL_OK
    }

    # Call the user-supplied callback function for adding a set
    set callback toastgraph_prm$id
    upvar #0 $callback cb
    set xvec ""
    set yvec ""
    $cb(addfunc) $w $s toastGraph_xvec toastGraph_yvec
    if {$toastGraph_xvec != "" && $toastGraph_yvec != ""} {
	lappend sl $s
	global default_color_set
	set cols $default_color_set
	set ncols [llength $cols]
	set idx [expr $s % $ncols]
	set col [lindex $cols $idx] 
	$w.g element create line$s -xdata $toastGraph_xvec \
	    -ydata $toastGraph_yvec -pixels 4 -linewidth 0 -label "Set $s" \
	    -color $col -trace increasing
	$w.g element bind line$s <Double-1> [list toastConfigGraphSet $w $s]
    }
}

# ============================================================================
# Delete existing graph element

proc toastGraph_setCallbackDel {w func} {
    regsub -all {\.} $w x id
    global toastgraph_prm$id
    set callback [join [list "toastgraph_prm" $id "(delfunc)"] {}]
    set $callback $func
    $w.mnu.set.menu entryconfig 2 -state normal
}

proc toastGraphSetDel {w id} {
    upvar #0 slist$id sl
    set s [$w.mnu.cset cget -value]
    set idx [lsearch $sl $s]
    if {$idx < 0} {
	ErrorNotice "Graph does not contain set $s!"
	return TCL_OK
    }
    $w.g element delete line$s
    set sl [lreplace $sl $idx $idx]

    # Call the user-supplied callback function for deleting a set
    set callback toastgraph_prm$id
    upvar #0 $callback cb
    $cb(delfunc) $w $s
}

# ============================================================================
# Change current graph element

proc toastGraph_setCallbackSet {w func} {
    regsub -all {\.} $w x id
    global toastgraph_prm$id
    set callback [join [list "toastgraph_prm" $id "(chngfunc)"] {}]
    set $callback $func
    $w.mnu.cset config -state normal
}

proc toastGraphSetSet {w id s} {
    # If the config dialog is open, reset for new set
    if [winfo exists .toastgraph_configel] {
	toastConfigGraphSet $w $s
    }

    # Call the user-supplied callback function for changing a set
    set callback toastgraph_prm$id
    upvar #0 $callback cb
    $cb(chngfunc) $w $s
}
   
# ============================================================================
# Configure existing graph element

proc toastConfigGraphSet {gf {s -1}} {
    global linestyle linecode symbstyle smooth lw ss ow lb
    set g $gf.g
    set dname .toastgraph_configel
    if {[winfo exists $dname] == 0} {
	set w [toplevel $dname]

	set f [frame $w.bbar]
	button $f.apply -text Apply
	button $f.done -text Done -command [list destroy $w]
	pack $f.apply $f.done -side left -padx 5 -pady 5
	pack $f -side bottom

	set f [frame $w.0 -relief groove -bd 2]
	pack $f -padx 5 -pady 5 -fill x
	tixLabelEntry  $f.lb -label Label -labelside top
	$f.lb subwidget entry config -textvariable lb
	pack $f.lb -side top -fill x -padx 2 -pady 2

	set f [frame $w.1 -relief groove -bd 2]
	pack $f -padx 5 -pady 5 -fill x

	tixComboBox $f.ls -dropdown yes -editable no -history no \
		-label "Line style" -labelside left -options {
	    entry.width 10
	    label.width 12
	    slistbox.scrollbar auto
	    slistbox.width 100
	    slistbox.height 150
	}
	for {set x 0} {$x < 7} {incr x} {
	    $f.ls insert end [lindex $linestyle $x]
	}

	tixComboBox $f.sm -dropdown yes -editable no -history no \
		-label "Smoothing" -labelside left -options {
	    entry.width 10
	    label.width 12
	    slistbox.scrollbar auto
	    slistbox.width 100
	    slistbox.height 100
	}
	for {set x 0} {$x < 4} {incr x} {
	    $f.sm insert end [lindex $smooth $x]
	}

	frame $f.1
	tixControl $f.1.lw -label "Line width" -min 0 -variable lw \
		-integer true -options {
	    entry.width 5
	    label.width 12
	}
	pack $f.1.lw -side left
	pack $f.ls $f.sm $f.1 -side top -fill x

	colorChooser $f.lcc -label "Line colour" -options {
	    entry.width 10
	    label.width 12
	}
	pack $f.lcc -side top

	set f [frame $w.2 -relief groove -bd 2]
	pack $f -padx 5 -pady 5 -fill x
	tixComboBox $f.sy -dropdown yes -editable no -history no \
		-label "Symbol style" -labelside left -options {
	    entry.width 10
	    label.width 12
	    slistbox.scrollbar auto
	}
	for {set x 0} {$x < 9} {incr x} {
	    $f.sy insert end [lindex $symbstyle $x]
	}
	frame $f.1
	tixControl $f.1.ss -label "Symbol size" -min 0 -variable ss \
		-integer true -options {
	    entry.width 5
	    label.width 12
	}
	pack $f.1.ss -side left
	frame $f.2
	tixControl $f.2.ow -label "Outline width" -min 0 -variable ow \
		-integer true -options {
	    entry.width 5
	    label.width 12
	}
	pack $f.2.ow -side left
	pack $f.sy $f.1 $f.2 -side top -fill x

	colorChooser $f.fcc -label "Fill colour" -options {
	    entry.width 10
	    label.width 12
	}
	pack $f.fcc -side top

	colorChooser $f.occ -label "Outline colour" -options {
	    entry.width 10
	    label.width 12
	}
	pack $f.occ -side top
    } else {
	set w $dname
	raise $w
    }
    if {$s < 0} {
	set s [$gf.mnu.cset cget -value]
    }
    wm title $w "Configure Set $s"
    $w.bbar.apply config -command  [list applyConfigGraphSet $w $g $s]

    # set dialog entries to current element style
    set lb [$g element cget line$s -label]
    set idx [lsearch $linecode [$g element cget line$s -dashes]]
    if {$idx < 0} {set idx 0}
    $w.1.ls config -value [lindex $linestyle $idx]
    set sm [$g element cget line$s -smooth]
    $w.1.sm config -value $sm
    set lw [$g element cget line$s -linewidth]
    set lc [$g element cget line$s -color]
    colorChooser_setColor $w.1.lcc $lc
    set sy [$g element cget line$s -symbol]
    $w.2.sy config -value $sy
    set ss [$g element cget line$s -pixels]
    set ow [$g element cget line$s -outlinewidth]
    set fc [$g element cget line$s -fill]
    if {$fc == "defcolor"} {set fc $lc}
    colorChooser_setColor $w.2.fcc $fc
    set oc [$g element cget line$s -outline]
    if {$oc == "defcolor"} {set oc $lc}
    colorChooser_setColor $w.2.occ $oc
    return TCL_OK
}

proc applyConfigGraphSet {w g s} {
    global linestyle linecode lw ss ow lb
    $w.1.1.lw update
    $w.2.1.ss update
    $w.2.2.ow update
    set dsh [$w.1.ls cget -value]
    set sy [$w.2.sy cget -value]
    if {$sy == "<none>"} {set sy ""}
    $g element config line$s \
	    -linewidth $lw -pixels $ss -outlinewidth $ow \
	    -color [colorChooser_getColor $w.1.lcc] \
	    -fill [colorChooser_getColor $w.2.fcc] \
	    -outline [colorChooser_getColor $w.2.occ] \
	    -symbol $sy \
	    -dashes [lindex $linecode [lsearch $linestyle $dsh]] \
	    -smooth [$w.1.sm cget -value] \
	    -label $lb

    $g axis config y -logscale [$g axis cget y -logscale]
    # we need this to force the graph to update

    return TCL_OK
}

# some line styles and their codes
set linestyle {"solid" "long dash" "medium dash" "short dash" "dot" \
	"dash dot" "dash dot dot"}
set linecode {"" {16 8} {10 5} {6 3} {2 2} {10 3 2 3} {10 3 2 3 2 3}}

# symbol styles
set symbstyle {"<none>" "circle" "square" "diamond" "triangle" "plus" \
	"cross" "splus" "scross"}

# smoothing styles
set smooth {"linear" "step" "quadratic" "natural"}

# relief types
set relieftype {"flat" "solid" "raised" "sunken" "groove" "ridge"}

# some explicit colors
set default_color_set {"black" "red" "blue" "green" "yellow" "cyan" "magenta" \
			   "orange" "darkred" "darkblue" "darkgreen" "white"}
