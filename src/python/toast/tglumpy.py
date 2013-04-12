import os
import numpy as np
import OpenGL.GL as gl
import glumpy as gp
from glumpy import figure, Trackball
from numpy import matrix
from scipy import sparse
from scipy.sparse import linalg
from numpy.random import rand
from types import *
import mesh

import pdb

class Mesh3D(object):
    def __init__(self,hmesh,nim,cm):
        nlist,elist,perm = mesh.SurfData (hmesh)
        bb = mesh.MeshBB(hmesh)
        bbmin = bb[:,0]
        bbmax = bb[:,1]
        bbcnt = (bbmin+bbmax)/2
        scale = 2/np.max(bbmax-bbmin)
        nlen = nlist.shape[0]
        self.elen = elist.shape[0]
        self.indices = np.zeros((self.elen,3), dtype=np.int32)
        self.vertices = np.zeros((self.elen*3,3), dtype=np.float32)
        self.normals = np.zeros((self.elen*3,3), dtype=np.float32)
        for xi in range(self.elen):
            d1 = nlist[elist[xi,1],:]-nlist[elist[xi,0],:]
            d2 = nlist[elist[xi,2],:]-nlist[elist[xi,0],:]
            nml = np.cross(d1,d2)
            nml = nml/np.linalg.norm(nml)
            for yi in range(3):
                ns = elist[xi,yi]
                nt = xi*3+yi
                self.indices[xi,yi] = nt
                for zi in range(3):
                    self.vertices[nt,zi] = (nlist[ns,zi]-bbcnt[zi])*scale
                    self.normals[nt,zi] = nml[zi]
        if nim == None:
            self.values = None
        else:
            self.values = np.zeros((self.elen*3,3), dtype=np.float32)
            nmin = np.min(nim)
            nmax = np.max(nim)
            if nmin == nmax:
                nmin = nmin-0.5
                nmax = nmax+0.5
            for xi in range(self.elen):
                for yi in range(3):
                    ns = elist[xi,yi]
                    nt = xi*3+yi
                    v = (nim[perm[ns]]-nmin)/(nmax-nmin)
                    col = cm.get_color(v)
                    self.values[nt,:] = col._get_rgb()
    def draw(self):
        gl.glEnableClientState(gl.GL_VERTEX_ARRAY)
        gl.glEnableClientState(gl.GL_NORMAL_ARRAY)
        gl.glVertexPointerf(self.vertices)
        gl.glNormalPointerf(self.normals)
        if self.values != None:
            gl.glEnableClientState(gl.GL_COLOR_ARRAY)
            gl.glColorPointer(3, gl.GL_FLOAT, 0, self.values)
        gl.glDrawElements(gl.GL_TRIANGLES, self.elen*3, gl.GL_UNSIGNED_INT, self.indices)
        gl.glDisableClientState(gl.GL_VERTEX_ARRAY)
        gl.glDisableClientState(gl.GL_NORMAL_ARRAY)
        gl.glDisableClientState(gl.GL_COLOR_ARRAY)
    def drawwire(self):
        gl.glEnableClientState(gl.GL_VERTEX_ARRAY)
        gl.glVertexPointerf(self.vertices)
        gl.glDrawElements(gl.GL_TRIANGLES, self.elen*3, gl.GL_UNSIGNED_INT, self.indices)
        gl.glDisableClientState(gl.GL_VERTEX_ARRAY)

    
class Mesh2D(object):
    def __init__(self,hmesh,nim,cm):
        nlist,elist,eltp = mesh.Data (hmesh)
        bb = mesh.BB(hmesh)
        bbmin = bb[:,0]
        bbmax = bb[:,1]
        bbcnt = (bbmin+bbmax)/2
        scale = 1.5/np.max(bbmax-bbmin)
        nlen = nlist.shape[0]
        elen = elist.shape[0]
        vertices = np.zeros((elen*3),dtype=[('position','f4',3)])
        #for i in range(elen):
        #    for j in range(3):
        #        vertices[i*3+j] = nlist[elist[i,j]
        # TODO!!!

        self.elen = elist.shape[0]
        self.indices = np.zeros((self.elen,3), dtype=np.int32)
        self.vertices = np.zeros((nlen,3), dtype=np.float32)
        self.values = np.zeros((nlen,3), dtype=np.float32)
        if nim == None:
            nim = np.zeros((mesh.NodeCount(hmesh),1), dtype=np.float32)
        nmin = np.min(nim)
        nmax = np.max(nim)
        if nmin == nmax:
            nmin = nmin-0.5
            nmax = nmax+0.5
        for xi in range(nlen):
            for zi in range(2):
                self.vertices[xi,zi] = (nlist[xi,zi]-bbcnt[zi])*scale
            self.vertices[xi,2] = 1
            v = (nim[xi]-nmin)/(nmax-nmin)
            col = cm.get_color(v)
            self.values[xi,:] = col._get_rgb()
        for xi in range(self.elen):
            for yi in range(3):
                ns = elist[xi,yi]
                self.indices[xi,yi] = ns
    def draw(self):
        gl.glEnableClientState(gl.GL_VERTEX_ARRAY)
        gl.glEnableClientState(gl.GL_COLOR_ARRAY)
        gl.glVertexPointerf(self.vertices)
        gl.glColorPointer(3, gl.GL_FLOAT, 0, self.values)
        gl.glDrawElements(gl.GL_TRIANGLES, self.elen*3, gl.GL_UNSIGNED_INT, self.indices)
        gl.glDisableClientState(gl.GL_VERTEX_ARRAY)
        gl.glDisableClientState(gl.GL_COLOR_ARRAY)
    def drawwire(self):
        gl.glEnableClientState(gl.GL_VERTEX_ARRAY)
        gl.glVertexPointerf(self.vertices)
        gl.glDrawElements(gl.GL_TRIANGLES, self.elen*3, gl.GL_UNSIGNED_INT, self.indices)
        gl.glDisableClientState(gl.GL_VERTEX_ARRAY)

    
#if __name__ == '__main__':
def ShowMesh3D(hmesh,nim,col,cmap,lighting,mode):

    cm = gp.colormap.Grey
    if cmap=='Hot':
        cm = gp.colormap.Hot
    elif cmap=='Fire':
        cm = gp.colormap.Fire
    elif cmap=='Ice':
        cm = gp.colormap.Ice
    elif cmap=='IceAndFire':
        cm = gp.colormap.IceAndFire

    wire = True
    fill = True
    if mode=='Wire':
        fill = False
    elif mode=='Fill':
        wire = False

    mesh = []
    if type(hmesh) is list:
        for i in range(len(hmesh)):
            mesh.append(Mesh3D(hmesh[i],nim,cm))
    else:
        mesh.append(Mesh3D(hmesh,nim,cm))
    
    window = gp.Window(800,800)
    trackball = gp.Trackball(0,0,2)

    @window.event
    def on_draw():
        gl.glClearColor(0,0,0,1)
        window.clear()
        trackball.push()
        if lighting==True:
            gl.glEnable (gl.GL_LIGHTING)
        gl.glEnable(gl.GL_DEPTH_TEST)
        gl.glEnable(gl.GL_BLEND)
        gl.glBlendFunc (gl.GL_SRC_ALPHA, gl.GL_ONE_MINUS_SRC_ALPHA)
        #I.shader.bind(I.texture,I._lut)
        for i in range(len(mesh)):
            if type(col) is list:
                c = col[i]
            else:
                c = col
            gl.glColor4f(c[0],c[1],c[2],c[3])
            if fill==True:
                mesh[i].draw()
            if wire==True:
                gl.glPolygonMode (gl.GL_FRONT, gl.GL_LINE)
                gl.glEnable(gl.GL_POLYGON_OFFSET_LINE)
                gl.glPolygonOffset (-1,0)
                mesh[i].drawwire()
                gl.glPolygonMode (gl.GL_FRONT, gl.GL_FILL)
                gl.glDisable(gl.GL_POLYGON_OFFSET_LINE)

        #I.shader.unbind()
        trackball.pop()

    @window.event
    def on_init():
        gl.glEnable (gl.GL_LIGHT0)
        gl.glLightfv (gl.GL_LIGHT0, gl.GL_DIFFUSE, (1.0, 0.7, 0.5, 1))
        gl.glLightfv (gl.GL_LIGHT0, gl.GL_AMBIENT, (0.2, 0.2, 0.2, 1))
        gl.glLightfv (gl.GL_LIGHT0, gl.GL_SPECULAR,(1.0, 1.0, 1.0, 1))
        #gl.glLightfv (gl.GL_LIGHT0, gl.GL_POSITION,(-1.0, 2.0, -1.0, 0.0))
        gl.glLightfv (gl.GL_LIGHT0, gl.GL_POSITION,(-0.5, -0.2, -1, 0))
        gl.glEnable (gl.GL_BLEND)
        #gl.glColorMaterial(gl.GL_FRONT_AND_BACK, gl.GL_AMBIENT_AND_DIFFUSE)
        gl.glColorMaterial(gl.GL_FRONT_AND_BACK, gl.GL_DIFFUSE)
        gl.glMaterialfv(gl.GL_FRONT, gl.GL_SHININESS, 50.0);
        gl.glEnable (gl.GL_COLOR_MATERIAL)
        gl.glBlendFunc (gl.GL_SRC_ALPHA, gl.GL_ONE_MINUS_SRC_ALPHA)
        gl.glPolygonMode (gl.GL_FRONT, gl.GL_FILL)
        gl.glFrontFace (gl.GL_CW)
        gl.glEnable (gl.GL_CULL_FACE)
        gl.glShadeModel (gl.GL_SMOOTH)
        
    @window.event
    def on_mouse_drag(x, y, dx, dy, button):
        trackball.drag_to(x,y,dx,dy)
        window.draw()

    @window.event
    def on_mouse_scroll(x, y, dx, dy):
        trackball.zoom_to(x,y,dx,dy)
        window.draw()

    window.mainloop()


def ShowMesh2D(hmesh,nim,cmap,mode):

    cm = gp.colormap.Grey
    if cmap=='Hot':
        cm = gp.colormap.Hot
    elif cmap=='Fire':
        cm = gp.colormap.Fire
    elif cmap=='Ice':
        cm = gp.colormap.Ice
    elif cmap=='IceAndFire':
        cm = gp.colormap.IceAndFire

    wire = True
    fill = True
    if mode=='Wire':
        fill = False
    elif mode=='Fill':
        wire = False
    
    mesh = Mesh2D(hmesh,nim,cm)

    fig = figure(size=(800,800))
    trackball = Trackball(0,0,2)
    fig.push (mesh)
    fig.show()

#    @window.event
#    def on_draw():
#        gl.glClearColor(0,0,0,1)
#        window.clear()
#        trackball.push()
#        gl.glDisable(gl.GL_DEPTH_TEST)
#        gl.glEnable(gl.GL_BLEND)
#        gl.glBlendFunc (gl.GL_SRC_ALPHA, gl.GL_ONE_MINUS_SRC_ALPHA)
#        #I.shader.bind(I.texture,I._lut)
#        mesh.draw()
#        if wire==True:
#            gl.glPolygonMode (gl.GL_FRONT, gl.GL_LINE)
#            gl.glColor4f(0,1,0,1)
#            mesh.drawwire()
#            gl.glPolygonMode (gl.GL_FRONT, gl.GL_FILL)
#
#       #I.shader.unbind()
#        trackball.pop()
#
#    @window.event
#    def on_init():
#        gl.glEnable (gl.GL_BLEND)
#        gl.glColorMaterial(gl.GL_FRONT_AND_BACK, gl.GL_AMBIENT_AND_DIFFUSE)
#        gl.glEnable (gl.GL_COLOR_MATERIAL)
#        gl.glBlendFunc (gl.GL_SRC_ALPHA, gl.GL_ONE_MINUS_SRC_ALPHA)
#        gl.glPolygonMode (gl.GL_FRONT, gl.GL_FILL)
#        gl.glFrontFace (gl.GL_CCW)
#        gl.glEnable (gl.GL_CULL_FACE)
#        gl.glShadeModel (gl.GL_SMOOTH)
#        
#    #@window.event
#    #def on_mouse_drag(x, y, dx, dy, button):
#    #    trackball.drag_to(x,y,dx,dy)
#    #    window.draw()
#
#    @window.event
#    def on_mouse_scroll(x, y, dx, dy):
#        trackball.zoom_to(x,y,dx,dy)
#        window.draw()
#
#    window.mainloop()


