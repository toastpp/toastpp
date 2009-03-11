wm title . "lem - linear elasticity module"

set runprm(meshfile) "\[undefined\]"
set runprm(matfile) "\[none\]"
set runprm(dispfile) "\[none\]"
set runprm(bforcefile) "\[none\]"

set whichdsp x

# =============================================================================

proc loadMeshDlg {w} {
    $w filedialog config -title "Load mesh file"
    set fsb [$w filedialog subwidget fsbox]
    $fsb config -filetypes {
	{{*.msh *.opt} {Mesh files}}
	{{*} {All files}}
    }
    $fsb subwidget types pick 0
    return TCL_OK
}

proc loadMatDlg {w} {
    $w filedialog config -title "Load material file"
    set fsb [$w filedialog subwidget fsbox]
    $fsb config -filetypes {
	{{*.dat} {Data files}}
	{{*} {All files}}
    }
    $fsb subwidget types pick 0
    return TCL_OK
}

proc loadDispDlg {w} {
    $w filedialog config -title "Load displacement file"
    set fsb [$w filedialog subwidget fsbox]
    $fsb config -filetypes {
	{{*.dat} {Data files}}
	{{*} {All files}}
    }
    $fsb subwidget types pick 0
    return TCL_OK
}

proc loadBforceDlg {w} {
    $w filedialog config -title "Load boundary force file"
    set fsb [$w filedialog subwidget fsbox]
    $fsb config -filetypes {
	{{*.dat} {Data files}}
	{{*} {All files}}
    }
    $fsb subwidget types pick 0
    return TCL_OK
}

# This is called when a slider is changed.
proc setAngle {axis value} {
    global xAngle yAngle zAngle

    switch -exact $axis {
	x {.renderwin.f1.rw setXrot $value}
	y {.renderwin.f1.rw setYrot $value}
	z {.renderwin.f1.rw setZrot $value}
    }
    .renderwin.f2.s$axis configure -label "$axis axis $value"
}

proc setDist {value} {
    global Dist

    .renderwin.f1.rw setDist $value
    .renderwin.f2.dst configure -label "dist $value"
}

proc showDisp {axis} {
    .renderwin.f1.rw showDisp $axis
}

# =============================================================================
# runSolver

proc runSolver {} {
    global runprm
    .renderwin.f1.rw Solve runprm
    return TCL_OK
}

# =============================================================================
# main window

lappend auto_path $env(TOAST_SCRIPT_PATH)
toastscheme

# ==== menu bar ====
MenuSetup .mnu
menubutton .mnu.file -text File -underline 0 -menu .mnu.file.menu
menu .mnu.file.menu
.mnu.file.menu add command -label "Exit" -underline 1 -command exit

pack .mnu.file -side left -fill x

frame .info -relief flat -bd 2

tixFileEntry .info.mesh -dialogtype tixExFileSelectDialog -labelside top \
    -label "Mesh file" -variable runprm(meshfile) \
    -activatecmd [list loadMeshDlg .info.mesh]
.info.mesh subwidget entry config -width 30

tixFileEntry .info.mat -dialogtype tixExFileSelectDialog -labelside top \
	-label "Material file" -variable runprm(matfile) \
	-activatecmd [list loadMatDlg .info.mat]
.info.mat subwidget entry config -width 30

tixFileEntry .info.disp -dialogtype tixExFileSelectDialog -labelside top \
    -label "Displacement file" -variable runprm(dispfile) \
    -activatecmd [list loadDispDlg .info.mesh]
.info.disp subwidget entry config -width 30

tixFileEntry .info.bforce -dialogtype tixExFileSelectDialog -labelside top \
    -label "Boundary force file" -variable runprm(bforcefile) \
    -activatecmd [list loadBforceDlg .info.mesh]
.info.bforce subwidget entry config -width 30

pack .info.mesh .info.mat .info.disp .info.bforce -side top -fill x -expand true
pack .info -fill both -expand true

set f [frame .bbar]
button $f.run -text Run -command runSolver
button $f.quit -text Quit -command exit
pack $f.run $f.quit -side left -padx 5
pack $f -side bottom -pady 10

set w [toplevel .renderwin]

set f [frame $w.f1]
togl $f.rw -width 600 -height 600 -rgba true -double true -depth true \
    -ident render
pack $f.rw -side left -padx 3 -pady 3 -fill both -expand t
pack $f -fill both -expand t -side top

set f [frame $w.f2]
scale $f.sx -label {x rotation} -from 0 -to 360 -command {setAngle x} \
    -orient horizontal -showvalue no -length 200
scale $f.sy -label {y rotation} -from 0 -to 360 -command {setAngle y} \
    -orient horizontal -showvalue no -length 200
scale $f.sz -label {z rotation} -from 0 -to 360 -command {setAngle z} \
    -orient horizontal -showvalue no -length 200
scale $f.dst -label {camera distance} -from 10 -to 500 -command {setDist} \
    -orient horizontal -showvalue no -length 200
pack $f.sx $f.sy $f.sz $f.dst -side top
pack $f -side left -fill x -expand f

set f [frame $w.f3]
radiobutton $f.dspx -text {x displacement} -command {showDisp x} \
    -variable whichdsp -value x
radiobutton $f.dspy -text {y displacement} -command {showDisp y} \
    -variable whichdsp -value y
radiobutton $f.dspz -text {z displacement} -command {showDisp z} \
    -variable whichdsp -value z
pack $f.dspx $f.dspy $f.dspz -side top
pack $f -side right

