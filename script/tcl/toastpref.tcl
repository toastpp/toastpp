# =============================================================================
# Default colour and font schemes for TOAST applications
# Initialise environment toastenv(browser,docs)

proc toastscheme {} {
    global env toastenv tixOption
    global toastBitmapPath

    # tix resetoptions TixGray 12Point
    # option add *disabledForeground gray
    set toastBitmapPath [file join $env(TOAST_SCRIPT_PATH) bitmaps]

    # reverse the Tix italic scale font (sorry, I don't like it)
    # option add *Scale.font $tixOption(bold_font)

    # Set the help browser application
    if {[llength [array get env TOAST_BROWSER]] > 0} {
	set toastenv(browser) $env(TOAST_BROWSER)
    } elseif {[llength [array get env BROWSER]] > 0} {
	set toastenv(browser) $env(BROWSER)
    } else {
	set toastenv(browser) netscape
    }

    #set the path to html documentation
    if {[llength [array get env TOAST_DOCS_PATH]] > 0} {
	set toastenv(docs) $env(TOAST_DOCS_PATH)
    } else {
	set toastenv(docs) [file join $env(TOASTDIR) docs]
    }

    return TCL_OK
}

# =============================================================================
# Read toast resource file <name> and store the results in the resource array
# rcvar. Resource files are searched in these locations:
#
#     $TOASTDIR/toastrc.d/<name>rc
#     $HOME/toastrc.d/<name>rc
#     $PWD/<name>rc
#
# Entries in $TOASTDIR/toastrc.d/<file> are overridden by
# $HOME/toastrc.d/<file>, and those by $PWD/<file>
#
# Resource files are ASCII with one entry per line. Entries are of the form
#
#     <variable> = <value>
#
# The rcvar resource will have fields corresponding to the <variables>'s
# defined in the file, and corresponding <value>'s. The should be set to
# default values before the call to be used if not specified in the file.

proc toastReadResource {name rcvar} {
    global env
    upvar $rcvar rc

    # some default entries
    set rc(pwd) $env(PWD)
    set rc(home) $env(HOME)

    # read from resource scripts in order
    set rcname [join [list $name "rc"] ""]
    foreach dir [list [file join $env(TOASTDIR) "toastrc.d"] \
		      [file join $env(HOME) "toastrc.d"] "."] {
	set rcfile [file join $dir $rcname]
	if {[catch {open $rcfile r} fileId] == 0} {
	    while {[gets $fileId line] >= 0} {
		# remove comments (starting with '#')
		set line [string trim [lindex [split $line \#] 0]]
		# subsitute environment variables
		set line [substEnv $line]
		# split into field and value
		if [string length $line] {
		    set cat [string trim [lindex [split $line "="] 0]]
		    set val [string trim [lindex [split $line "="] 1]]
		    set rc($cat) $val
		}
	    }
	    close $fileId
	}
    }
}

# =============================================================================
# Status bar management

proc statusSetup {} {
    frame .status -relief raised -bd 2
    label .status.l -width 10 -anchor w -text Ready
    pack  .status.l -side left -fill x -expand true -padx 2 -pady 2
    pack  .status -side bottom -fill x
}

proc statusOut {msg} {
    .status.l config -text $msg
    update
    return TCL_OK
}

proc statusReady {} {
    statusOut "Ready"
    return TCL_OK
}

# =============================================================================
# Menu bar management

proc MenuSetup {m} {
    frame $m -relief raised -bd 2
    pack  $m -side top -fill x
    return TCL_OK
}

# =============================================================================
# A horizontal line (separator in dialog boxes)

proc Line {w} {
    frame $w -relief sunken -bd 1 -height 2
    pack  $w -side top -fill x -expand true -padx 10
}

# =============================================================================
# Error boxes

proc ErrorNotice {msg} {
    set w [toplevel .ebox]
    blt::busy hold .
    wm title $w "Error"
    message $w.msg -justify left -aspect 400 -anchor w -text $msg
    pack $w.msg -side top -fill x -expand true -padx 10 -pady 5
    button $w.b -text "Continue" -command [list destroy $w]
    pack $w.b -side top -padx 10 -pady 5
    tkwait window $w
    blt::busy release .
    return TCL_OK
}

proc AskBox {msg} {
    global res
    set w [toplevel .ebox]
    blt::busy hold .
    wm title $w "Warning"
    message $w.msg -justify left -aspect 400 -anchor w -text $msg
    pack $w.msg -side top -fill x -expand true -padx 10 -pady 5
    set f [frame $w.b]
    button $f.y -text "Yes" -command { set res 1 }
    button $f.n -text "No" -command { set res 0 }
    pack $f.y $f.n -side left -padx 10 -pady 5
    pack $f -side top
    tkwait variable res
    destroy $w
    blt::busy release .
    return $res
}

# =============================================================================
# An auxiliary routine to substitute environment variables in a string

proc substEnv {str} {
    global env
    while {[regexp {[^\$]*\$\(([^\)]+)\).*} $str match name] == 1} {
	if {[llength [array get env $name]] > 0} {
	    set val $env($name)
	} else {
	    set val ""
	}
	regsub {\$\([^\)]+\)} $str $val str
    }
    return $str
}
# =============================================================================
# HTML documentation (application app, document doc, e.g. index.html#toc1)
# This assumes documenation to be stored in the default location
# $toastenv(docs)/$app/$doc

proc htmldoc {app doc} {
    global toastenv
    return [exec $toastenv(browser) [file join $toastenv(docs) $app $doc] &]
}
