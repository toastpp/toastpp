# Examples from Brent B. Welch
# Practical Programming in Tcl and Tk (second edition)

#
# Example 33-1
# Procedures to help build dialogs.
#

proc Dialog_Create {top title args} {
	global dialog
	if [winfo exists $top] {
		switch -- [wm state $top] {
			normal {
				# Raise a buried window
				raise $top
			}
			withdrawn -
			iconified {
				# Open and restore geometry
				wm deiconify $top
				catch {wm geometry $top $dialog(geo,$top)}
			}
		}
		return 0
	} else {
		eval {toplevel $top} $args
		wm title $top $title
		return 1
	}
}
proc Dialog_Wait {top varName {focus {}}} {
	upvar $varName var

	# Poke the variable if the user nukes the window
	bind $top <Destroy> [list set $varName $var]

	# Grab focus for the dialog
	if {[string length $focus] == 0} {
		set focus $top
	}
	set old [focus -displayof $top]
	focus $focus
	catch {tkwait visibility $top}
	catch {grab $top}

	# Wait for the dialog to complete
	tkwait variable $varName
	catch {grab release $top}
	focus $old
}
proc Dialog_Dismiss {top} {
	global dialog
	# Save current size and position
	catch {
		# window may have been deleted
		set dialog(geo,$top) [wm geometry $top]
		wm withdraw $top
	}
}
#
# Example 36-5
# Font selection dialog.
#

proc Font_Select {{top .fontsel}} {
	global font
        set font(ok) 0

	# Create the menus for the menu bar

	toplevel $top -class Fontsel -bd 10
	set menubar [menu $top.menubar]
	$top config -menu $menubar
	foreach x {File Font Size Format} {
		set menu [menu $menubar.[string tolower $x]]
		$menubar add cascade -menu $menu -label $x
	}
	$menubar.file add command -label Reset -command FontReset
	$menubar.file add command -label OK \
		-command {set font(ok) 1}
	$menubar.file add command -label Cancel \
		-command {set font(ok) 0}

	# Build a cascaded family menu if there are lots of fonts

	set allfonts [font families]
	set numfonts [llength $allfonts]
	set limit 20
	if {$numfonts < $limit} {
		# Just a single level menu
		foreach family $allfonts {
			$menubar.font add radio -label $family \
				-variable font(-family) \
				-value $family \
				-command FontUpdate
		}
	} else {
		set c 0 ; set l 0
		foreach family $allfonts {
			if {$l == 0} {
				$menubar.font add cascade -label $family... \
					-menu $menubar.font.$c
				set m [menu $menubar.font.$c]
				incr c
			}
			$m add radio -label $family \
				-variable font(-family) \
				-value $family \
				-command FontUpdate
			set l [expr ($l +1) % $limit]
		}
	}

	# Complete the other menus

	foreach size {7 8 10 12 14 18 24 36 72} {
		$menubar.size add radio -label $size \
			-variable font(-size) \
			-value $size \
			-command FontUpdate
	}
	$menubar.size add command -label Other... \
			-command [list FontSetSize $top]
	$menubar.format add check -label Bold \
			-variable font(-weight) \
			-onvalue bold -offvalue normal \
			-command FontUpdate
	$menubar.format add check -label Italic \
			-variable font(-slant) \
			-onvalue italic -offvalue roman \
			-command FontUpdate
	$menubar.format add check -label underline \
		-variable font(-underline) \
		-command FontUpdate
	$menubar.format add check -label overstrike \
		-variable font(-overstrike) \
		-command FontUpdate

	FontReset

	# This label displays the current font
	label $top.font -textvar font(name) -bd 5

	# A message displays a string in the font.
	message $top.msg -aspect 1000 \
					-borderwidth 10 -font fontsel \
					-text  "
ABCDEFGHIJKLMNOPQRSTUVWXYZ
abcdefghijklmnopqrstuvwxyz
0123456789
!@#$%^&*()_+-=[]{};:\"'` ~,.<>/?\\|
"

	# Lay out the dialog

	pack $top.font $top.msg  -side top
	set f [frame $top.buttons]
	button $f.ok -text Ok -command {set font(ok) 1}
	button $f.cancel -text Cancel -command {set font(ok) 0}
	pack $f.ok $f.cancel -padx 10 -side left
	pack $f -side top

	# Dialog_Wait is defined in Example 33-1 on page 437
	Dialog_Wait $top font(ok)
	destroy $top
	if {$font(ok)} {
		return [array get font -*]
	} else {
		return {}
	}
}
proc FontUpdate { } {
	global font

	# The elements of font that have a leading - are
	# used directly in the font configuration command.

	eval {font configure fontsel} [array get font -*]
	FontSet
}
proc FontReset {} {
	catch {font delete fontsel}
	font create fontsel
	FontSet
}
proc FontSet {} {
	global font

	# The name is the font configuration information
	# with a line break so it looks nicer

	set font(name) [font actual fontsel]
	regsub -- "-slant" $font(name) "\n-slant" font(name)

	# Save the actual parameters after any font substitutions

	array set font [font actual fontsel]
}
proc FontSetSize {top} {

	# Add an entry to enter a specific size.

	set f [frame $top.size -borderwidth 10]
	pack $f -side top -fill x
	label $f.msg -text "Size:"
	entry $f.entry -textvariable font(-size)
	bind $f.entry <Return> FontUpdate
	pack $f.msg -side left
	pack $f.entry -side top -fill x
}


