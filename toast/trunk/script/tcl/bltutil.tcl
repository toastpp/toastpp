# ============================================================================
# Some utility routines in support of the BLT Tcl extension package
# ============================================================================

# ============================================================================
# set the font attribute for widget 'wg' from user dialog input
# return value is the new font, or "" if user cancels

proc selectFont {fntVar} {
    upvar $fntVar fnt
    global fntsel
    set fntsel(ok) 0
    set fntsel(fmly) ""
    set fntsel(size) 0
    set stdsize {6 8 10 12 14 18 24}

    if {[winfo exists .fntsel] == 0} {
	set w [toplevel .fntsel]

	set f [frame $w.bbar]
	button $f.apply -text Select -command {set fntsel(ok) 1}
	button $f.done -text Cancel -command {set fntsel(ok) 0}
	pack $f.apply $f.done -side left -padx 5 -pady 5
	pack $f -side bottom

	set f [frame $w.sample]
	message $f.l -aspect 1000 -borderwidth 5 -text \
"ABCDEFGHIJKLMNOPQRSTUVWXYZ
abcdefghijklmnopqrstuvwxyz
0123456789!$%^&*()_+"
	pack $f.l -side top -padx 10 -pady 10 -expand no
	pack $f -side top -fill x

	set f [frame $w.1]
	pack $f -side top
	tixComboBox $f.fmly -dropdown yes -editable no -history no \
		-label "Family" -labelside left \
		-variable fntsel(fmly) \
		-command [list setFont $w.sample.l] \
		-options {
	    entry.width 14
	    slistbox.scrollbar auto
	    slistbox.width 100
	    slistbox.height 100
	}

	tixComboBox $f.size -dropdown yes -editable yes -history no \
		-label "Size" -labelside left \
		-variable fntsel(size) \
		-command [list setFont $w.sample.l] \
		-options {
	    entry.width 6
	    slistbox.scrollbar auto
	    slistbox.width 20
	    slistbox.height 100
	}
	for {set x 0} {$x < 7} {incr x} {
	    $f.size insert end [lindex $stdsize $x]
	}
	pack $f.fmly $f.size -side left

    } else {
	set w .fntsel
    }
    wm title $w "Select Font"
    set allfml [font families]
    set nfml [llength $allfml]
    array set fontspec [font actual $fnt]
    $w.1.fmly subwidget listbox delete 0 end
    for {set x 0} {$x < $nfml} {incr x} {
	$w.1.fmly insert end [lindex $allfml $x]
    }
    $w.1.fmly config -value $fontspec(-family)
    $w.1.size config -value $fontspec(-size)
    $w.sample.l config -font [array get fontspec]

    tkwait variable fntsel(ok)
    if {$fntsel(ok) == 1} {
	set fontspec(-family) [$w.1.fmly cget -value]
	set fontspec(-size) [$w.1.size cget -value]
    }
    set fnt [array get fontspec]
    destroy $w
    return TCL_OK
}

proc setFont {w f} {
    global fntsel
    $w config -font [list -family $fntsel(fmly) -size $fntsel(size)]
}

# =============================================================================
# Color chooser widget

proc colorChooser {w args} {
    global toastBitmapPath default_color_set

    eval {tixComboBox $w} $args -fancy yes -editable yes \
	    -tickbitmap @[file join $toastBitmapPath dots.xbm]
    $w config -command [list colorChooser_cmd $w]
    set colors $default_color_set
    set ncol [llength $colors]
    for {set x 0} {$x < $ncol} {incr x} {
	$w insert end [lindex $colors $x]
    }
    destroy [$w subwidget cross]
    set tickb [$w subwidget tick]
    $tickb config -command [list colorChooser_dlg $w]
    return $w
}

proc colorChooser_getColor {w} {
    return [$w cget -value]
}

proc colorChooser_setColor {w col} {
    $w config -value $col
    colorChooser_cmd $w $col
}

proc colorChooser_cmd {w col} {
    set col [string trim $col]

    # Check to make sure that the color is valid
    if [catch {set color [winfo rgb . $col]} ] {
	$w config -value [$w subwidget tick cget -background]
	return
    }

    set rgb [winfo rgb $w $col]
    if {[expr [lindex $rgb 0] + [lindex $rgb 1] + [lindex $rgb 2]] > 98304} {
	set icol black
    } else {
	set icol white
    }
    $w subwidget tick config -background $col -activebackground $col \
	    -foreground $icol -activeforeground $icol
}

proc colorChooser_dlg {w} {
    set col [tk_chooseColor -initialcolor [colorChooser_getColor $w]]
    if {$col != ""} {
	colorChooser_setColor $w $col
    }
}

set default_color_set {"black" "red" "blue" "green" "yellow" "cyan" "magenta" \
			   "orange" "darkred" "darkblue" "darkgreen" "white"}
