# =============================================================================
# Definition of some useful "common" dialog boxes
# =============================================================================

# =============================================================================
# Dialog for writing a pixmap out to a file. Optionally filters for GIF and
# PS output can be supplied.
# The following fields in the global "rc" array are examined:
#   ppm_to_gif_filt: Filter program for ppm->gif conversion (stdin->stdout)
#   ppm_to_ps_filt:  Filter program for ppm->ps conversion (stdin->stdout)
# If any of these is not defined then the corresponding output format does not
# appear in the dialog.
# Arguments:
#   colour: 0 indicates greyscale image, otherwise colour
#   ppmcmd: Name of the command that writes the ppm image to file. The command
#           must accept the file name as its only argument
# =============================================================================

proc toastExportPixmapDialog {colour ppmcmd} {

    global rc
    set dlg .toastexportpixmapdialog
    if [winfo exists $dlg] { destroy $dlg }
    toplevel $dlg
    wm title $dlg "Export pixmap image"

    set f [tixExFileSelectBox $dlg.fsb -command [list DoExport $dlg $ppmcmd]]
    $f subwidget cancel config -command [list destroy $dlg]
    pack $f -side top -fill both -expand yes

    # Build up the file types list
    set tlist [list [list {*} {All files}]]
    if {$colour} {
	lappend tlist [list {*.ppm} {Portable pixmap files}]
    } else {
	lappend tlist [list {*.pgm} {Portable graymap files}]
    }
    if {[catch {set tmp $rc(ppm_to_gif_filt)}] == 0} {
	lappend tlist [list {*.gif} {GIF files}]
    }
    if {[catch {set tmp $rc(ppm_to_ps_filt)}] == 0} {
	lappend tlist [list {*.ps *.eps} {Postscript files}]
    }
    $dlg.fsb config -filetypes $tlist
    $dlg.fsb subwidget types pick 1
}

proc DoExport {dlg ppmcmd file} {
    global rc
    switch [file extension $file] {
	".gif" {
	    set tmpfile [join [list $file "$"] ""]
	    $ppmcmd $tmpfile
	    if [catch {exec $rc(ppm_to_gif_filt) < $tmpfile > $file \
			   2> /dev/null}] {
		ErrorNotice "GIF filter failed ($rc(ppm_to_gif_filt))"
	    }
	    file delete $tmpfile
	}
	".ps" -
	".eps" {
	    set tmpfile [join [list $file "$"] ""]
	    $ppmcmd $tmpfile
	    if [catch {exec $rc(ppm_to_ps_filt) < $tmpfile > $file \
			   2> /dev/null}] {
		ErrorNotice "Postscript filter failed ($rc(ppm_to_ps_filt))"
	    }
	    file delete $tmpfile
	}
	default {
	    $ppmcmd $file
	}
    }
    destroy $dlg
}
