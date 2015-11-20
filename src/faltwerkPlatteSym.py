# -*- coding: utf-8 -*-
'''
Created on 22.02.2013

@author: heller
'''
import sys
import math
import numpy as np

#from feap import *

NODE_COUNTER = 1
ELEM_COUNTER = 1
MATE_COUNTER = 0

class BoundaryInput:
   def __init__(self, cmd, comment='', params=()):
      self.cmd = cmd
      self.comment = comment
      self.params = params
      
   def __str__(self):
      s = "%s   ** %s\n" % (self.cmd, self.comment)
      for line in self.params:
         for val in line:
            s+= "%s," % str(val)
         s = s[:-1] + "\n"
         
      return s

class Header:
   def __init__(self, name, params, solv=2):
      self.name = name
      self.params = params
      self.solver = solv
      
   def __str__(self):
      s = "feap   ***%s***\n" % self.name
      for p in self.params:
         s += "%s," % str(p)
      s = s[:-1] + "\n" # cut last comma, add newline
      s += "\nsolv\n%i\n\n" % self.solver
      return s
      

class Footer:
   def __init__(self, tie=True, macr=(), macrBlanks=0, link=()):
      self.tie = tie
      self.macr = macr
      self.link = link
      self.macrBlanks = macrBlanks
      
   def __str__(self):
      s = "\n\nend\n"
      if self.tie: s+= "tie,nopr\ntie\ntie,prin\n"
      
      if len(self.link)>0:
         s+= "link\n"
         for l in self.link:
            s+= "%s\n" % str(l)
      
      if len(self.macr)>0:
         s+= "\nmacr\n"
         for m in self.macr:
            s+= "%s\n" % str(m)
         s+= "end\n"
         s+= self.macrBlanks*"\n"
      
      s+= "inte\nstop\n\n"
      return s

class Material:
   def __init__(self, elmt, params):
      global MATE_COUNTER
      MATE_COUNTER += 1
      self.id = MATE_COUNTER
      self.elmt = elmt
      self.params = params
      
   def __str__(self):
      s = "%i,%i\n" % (self.id, self.elmt)
      for line in self.params:
         for val in line:
            s += "%s," % str(val)
         s = s[:-1] + "\n" # cut last comma, add newline
      
      return s

class Element:
   def __init__(self, node_a, node_b, node_c, node_d, matn=1):
      global ELEM_COUNTER
      ELEM_COUNTER += 1
      self.id = ELEM_COUNTER
      self.node_a = node_a.id
      self.node_b = node_b.id
      self.node_c = node_c.id
      self.node_d = node_d.id
      self.matn = matn
      
   def __str__(self):
      return "%i,%i,%i,%i,%i,%i,0\n" % (self.id, self.matn, self.node_a, self.node_b, self.node_c, self.node_d)
   
   def SetMaterial(self, matn):
      self.matn = matn
      
class Block2d:
   def __init__(self, numNodes, rInc, sInc, node1, elmt1, matn, coords, comment=''):
      global ELEM_COUNTER
      global NODE_COUNTER
      
      if len(coords) != numNodes:
         raise Exception("Block error: specified number of nodes doesn't match number of nodal coordinates")
      
      ELEM_COUNTER += rInc*sInc
      NODE_COUNTER += (rInc+1)*(sInc+1)
      
      self.numNodes = numNodes
      self.rInc = rInc
      self.sInc = sInc
      self.node1 = node1
      self.elmt1 = elmt1
      self.matn = matn
      self.coords = coords
      self.comment = comment
      
      # determine orientation of a block by
      # finding the plane on which the first,
      # second and last nodes A, B, D lie.
      
      # determine vectors
      A, B, D = coords[0], coords[1], coords[3]
      AB = B[0]-A[0], B[1]-A[1], B[2]-A[2]
      AD = D[0]-A[0], D[1]-A[1], D[2]-A[2]
      
      # determine plane normal as cross product AB x AD
      self.normal = np.cross(AB, AD) # this is the local z-Axis of the block
      
      # normalize normal vector to unit normal vector
      n2 = math.sqrt(self.normal[0]**2 + self.normal[1]**2 + self.normal[2]**2)
      self.normal = (self.normal[0]/n2, self.normal[1]/n2, self.normal[2]/n2)
      
      print "Block Unit Normal Vector:", self.normal
      
   def __str__(self):
      s = 'bloc%s\n' % (" **" + self.comment if self.comment != '' else '')
      s+='%i,%i,%i,%i,%i,%i\n' % (self.numNodes, self.rInc, self.sInc, self.node1, self.elmt1, self.matn)
      i = 1
      for node in self.coords:
         s+= '%i,%f,%f' % (i, node[0], node[1])
         if len(node)>2: s+= ',%f' % node[2]
         s+='\n' 
         i+= 1
      s+='\n'
      return s

class Node:
   def __init__(self, x=0., y=0., z=0.):
      global NODE_COUNTER
      NODE_COUNTER += 1
      self.id = NODE_COUNTER
      self.x = x
      self.y = y
      self.z = z
      
   def __str__(self):
      return "%.1f, %.1f, %.1f" % (self.x,self.y,self.z)
   
   
def main():
   global NODE_COUNTER
   global ELEM_COUNTER
   
   
   #########################
   #### OPTIONS
   #########################
   
   filename = "C:/FG/feap/exe/rve/microShell/Faltwerke/iPlatteSym"
   
   #rand = "Einspannung"
   #rand = "Navier"
   rand = "NavierFullZ"
   
   # parameters
   li = 0      # linear geometry?  0=lin,1=nili (moderate),2=nili (finite) 
   t = 0.1     # shell thickness / DECKSCHICHT
   tSteg=0.3
   e = 7000    # Young's modulus
   v = 0.34    # Poisson ratio
   f = -0.005  # loading force
   
   # material - choose one, comment the other(s)
   
#    matn = 1 # linear-elastisch
#    mateCard = (e,v,0)
#    nmnbnz = '000'
    
   matn = 4 # elasto-plastisch
   y0 = 10.
   xk = 0.
   nmnbnz = '002'
   mateCard = (e,v,y0,0,xk,0,0)
   
   
   # linear elements on the outer edges
   linmatn = 1 # linear-elastisch
   linmateCard = (e,v,0)
   linnmnbnz = '000'
   
   lincells = 0 # per side!
   
   
   
   #########################
   #### END OPTIONS
   #########################
   
   if lincells>0:
      filename += "2Mat"
      assert (matn != 1) # don't use two different linear materials!
   
   # generate simplified "honeycomb" mesh
   
   #segments = []
   materials = []
   boundaries = []
   
   #additionalElements = []
   
   # create mesh
   
   # if closeEdges is True, there are two more outer shells per direction
   nx = 16 # nr. of vertically aligned shells in y direction (x=const.)
   ny = 16 # nr. of vertically aligned shells in x direction (y=const.)
   
   if not (nx%2==0 and ny%2==0):
      print "nx, ny müssen gerade Zahlen sein!"
      return 
   
   # hx/hy müssen Vielfache von 2 sein, damit man Stege in der Mitte eines "RVEs" mit der Deckschicht noch erwischt, außerdem
   # falls nx > 0 muss hy Vielf. von 6 sein, falls ny > 0 muss hy Vielf. von 6 sein !
   hx = 16 # elements between two vertical shells (x direction) // total elements in x-direction if nx <= 1
   hy = 16 # elements between two vertical shells (y direction) // total elements in y-direction if ny <= 1
   
   hz = 8 # elements in thickness direction (mind. 2 wegen Mittenlagerung!!)
   
   closeEdges = False # set to True to add vertical shells on the outer edges of the system 
   
   # global dimensions               [kN,cm]
   lx = 2.*nx #2.*nx
   ly = 2.*ny #2.*ny
   h  = 0.8
   
   # element lengths
   #dx = lx / (hx*(nx+1))
   #dy = ly / (hy*(ny+1))
   #dz = h / hz
   
   eps = (lx+ly+h)*1.e-4 # small length parameter

   # ob. Deckschicht und Waben/unt. Deckschicht trennen, lässt sich dann besser plotten (zB. in ParaView)
   materials.append(Material(35, (
                  (1,0,0,0),
                  (1,2,t,-0.5*t,matn,li,5,nmnbnz,1,3),
                  (0,),
                  mateCard,
                  (0,t)     )))
   
   materials.append(Material(35, (
                  (1,0,0,0),
                  (1,2,t,-0.5*t,matn,li,5,nmnbnz,1,3),
                  (0,),
                  mateCard,
                  (0,t)     )))
   
   # Stege ggf. mit anderer Dicke
   materials.append(Material(35, (
                  (1,0,0,0),
                  (1,2,tSteg,-0.5*tSteg,matn,li,5,nmnbnz,1,3),
                  (0,),
                  mateCard,
                  (0,tSteg)     )))
   
   
   if lincells>0:
      # dasselbe für lineares Material
      materials.append(Material(35, (
                     (1,0,0,0),
                     (1,2,t,-0.5*t,linmatn,li,5,linnmnbnz,1,3),
                     (0,),
                     linmateCard,
                     (0,t)     )))
      
      materials.append(Material(35, (
                     (1,0,0,0),
                     (1,2,t,-0.5*t,linmatn,li,5,linnmnbnz,1,3),
                     (0,),
                     linmateCard,
                     (0,t)     )))
      
      # Stege ggf. mit anderer Dicke
      materials.append(Material(35, (
                     (1,0,0,0),
                     (1,2,tSteg,-0.5*tSteg,linmatn,li,5,linnmnbnz,1,3),
                     (0,),
                     linmateCard,
                     (0,tSteg)     )))
   
   header = Header("Schalen Sandwich Mikroelement", ('','','',3,6,4), solv=4)

   
   # top and bottom plates
   top = True
   btm = True
   
   if lincells > 0 and nx/2>lincells and ny/2>lincells:
      # regulaer / plastischer Block
      sizeX = (hx*(nx/2-lincells))
      linSizeX = lincells*hx
      
      lx_inner = lx/2. * (nx/2-lincells)/(nx/2)
      lx_outer = lx/2. - lx_inner
      
      
      sizeY = (hy*(ny/2-lincells))
      linSizeY = lincells*hy
      
      ly_inner = ly/2. * (ny/2-lincells)/(ny/2)
      ly_outer = ly/2. - ly_inner
      
      
      if top:
         topBloc = Block2d(4, sizeX, sizeY, NODE_COUNTER, ELEM_COUNTER, 2, 
                          ((0,0,h/2.), (lx_inner,0,h/2.), (lx_inner,ly_inner,h/2.), (0,ly_inner,h/2.) ), "Deckschicht oben (regulaer)")
         
         # lineare Blöcke am Rand
         linTopBlocX = Block2d(4, sizeX, linSizeY, NODE_COUNTER, ELEM_COUNTER, 5, 
                          ((0,ly_inner,h/2.), (lx_inner,ly_inner,h/2.), (lx_inner,ly/2.,h/2.), (0,ly/2.,h/2.) ), "Deckschicht oben (linear) X")
         linTopBlocY = Block2d(4, linSizeX, sizeY, NODE_COUNTER, ELEM_COUNTER, 5, 
                          ((lx_inner,0,h/2.), (lx/2.,0,h/2.), (lx/2.,ly_inner,h/2.), (lx_inner,ly_inner,h/2.) ), "Deckschicht oben (linear) Y")
         linTopBlocXY = Block2d(4, linSizeX, linSizeY, NODE_COUNTER, ELEM_COUNTER, 5,
                          ((lx_inner,ly_inner,h/2.), (lx/2., ly_inner, h/2.), (lx/2., ly/2., h/2.), (lx_inner, ly/2., h/2.) ), "Deckschicht oben (linear) XY")
         
      if btm:
         btmBloc = Block2d(4, sizeX, sizeY, NODE_COUNTER, ELEM_COUNTER, 1, 
                          ((0,0,-h/2.), (lx_inner,0,-h/2.), (lx_inner,ly_inner,-h/2.), (0,ly_inner,-h/2.) ), "Deckschicht unten (regulaer)")
         
         # lineare Blöcke am Rand
         linBtmBlocX = Block2d(4, sizeX, linSizeY, NODE_COUNTER, ELEM_COUNTER, 4, 
                          ((0,ly_inner,-h/2.), (lx_inner,ly_inner,-h/2.), (lx_inner,ly/2.,-h/2.), (0,ly/2.,-h/2.) ), "Deckschicht unten (linear) X")
         linBtmBlocY = Block2d(4, linSizeX, sizeY, NODE_COUNTER, ELEM_COUNTER, 4, 
                          ((lx_inner,0,-h/2.), (lx/2.,0,-h/2.), (lx/2.,ly_inner,-h/2.), (lx_inner,ly_inner,-h/2.) ), "Deckschicht unten (linear) Y")
         linBtmBlocXY = Block2d(4, linSizeX, linSizeY, NODE_COUNTER, ELEM_COUNTER, 4,
                          ((lx_inner,ly_inner,-h/2.), (lx/2., ly_inner, -h/2.), (lx/2., ly/2., -h/2.), (lx_inner, ly/2., -h/2.) ), "Deckschicht unten (linear) XY")
         
         
   else:
      if top: topBloc = Block2d(4, hx*nx/2 if nx > 0 else hx, hy*ny/2 if ny > 0 else hy, NODE_COUNTER, ELEM_COUNTER, 2, 
                       ((0,0,h/2.), (lx/2.,0,h/2.), (lx/2.,ly/2.,h/2.), (0,ly/2.,h/2.) ), "Deckschicht oben")
      if btm: btmBloc = Block2d(4, hx*nx/2 if nx > 0 else hx, hy*ny/2 if ny > 0 else hy, NODE_COUNTER, ELEM_COUNTER, 1, 
                    ((0,0,-h/2.), (lx/2.,0,-h/2.), (lx/2.,ly/2.,-h/2.), (0,ly/2.,-h/2.) ), "Deckschicht unten")
   
   # node to plot in TPLO: Mid node (n+1)/2 of top layer, but node counter is already n+1 as it denotes the number of the *next* node to be inserted
   # with symmetry: tploNode = 1 (Origin)
   tploNode = 1
   
   # Draufsicht
   #footer = Footer(tie=True, macr=('plot,rot1,0', 'plot,rot3,0', 'plot,mesh', 'plot,axis', 'plot,boun,,1', 'plot,boun,,2', 'plot,boun,,3', 'plot,load'), macrBlanks=0 )
   
   # schräge Ansicht
   footer = Footer(tie=True, macr=('prop', 'plot,rot1,-60', 'plot,rot3,30', 'plot,mesh', 'plot,axis', 'plot,boun,,1', 'plot,boun,,2', 'plot,boun,,3',
                                   'plot,boun,,4', 'plot,boun,,5', 'plot,node,,-%i' % tploNode, 'tplo,init,%i,-3,1' % tploNode, 'tplo,mark,2,4,1', 'stre,node',
                                   'parv,init,10'), macrBlanks=1 )
   
   vxBlocs = []
   vyBlocs = []
   xLines = [] # speichere x/y Koordinaten des Gitters, dies sind (zumindest oben und unten) Verschneidungskanten wo 6.DOF freigegeben werden muss
   yLines = []
   
   # vertikal ausgerichtete Blöcke, konstantes X
   for ix in range(nx//2-lincells):
      deltaX = 1.*  lx / nx
      curX = deltaX/2. + deltaX*ix
      
      if lincells == 0:
         vxBlocs.append(  Block2d(4, hy*ny/2 if ny>0 else hy, hz, NODE_COUNTER, ELEM_COUNTER, 3, 
                       ((curX,0,-h/2), (curX,ly/2,-h/2), (curX,ly/2,h/2), (curX,0,h/2) ), "x=%.2f Schale" % curX))
      else:
         vxBlocs.append(  Block2d(4, sizeY, hz, NODE_COUNTER, ELEM_COUNTER, 3, 
                       ((curX,0,-h/2), (curX,ly_inner,-h/2), (curX,ly_inner,h/2), (curX,0,h/2) ), "x=%.2f Schale (reg.)" % curX))
         vxBlocs.append(  Block2d(4, linSizeY, hz, NODE_COUNTER, ELEM_COUNTER, 6, 
                       ((curX,ly_inner,-h/2), (curX,ly/2.,-h/2), (curX,ly/2.,h/2), (curX,ly_inner,h/2) ), "x=%.2f Schale (lin.)" % curX))
   
      xLines.append(curX)

   for ix in range(nx//2-lincells, nx//2):
      deltaX = 1.*  lx / nx
      curX = deltaX/2. + deltaX*ix
      
      vxBlocs.append(   Block2d(4, hy*ny/2 if ny>0 else hy, hz, NODE_COUNTER, ELEM_COUNTER, 6, 
                       ((curX,0,-h/2), (curX,ly/2.,-h/2), (curX,ly/2.,h/2), (curX,0,h/2) ), "x=%.2f Schale (lin.)" % curX))
      
      xLines.append(curX)
      
#    if nx>0 and closeEdges:
#       curX = lx/2.
#       vxBlocs.append(  Block2d(4, hy*ny/2 if ny>0 else hy, hz, NODE_COUNTER, ELEM_COUNTER, 3, 
#                     ((curX,0,-h/2), (curX,ly/2,-h/2), (curX,ly/2,h/2), (curX,0,h/2) ), "x=%.2f Schale (Abschluss)" % curX))
#       xLines.append(curX)
   
   # vertikal ausgerichtete Blöcke, konstantes Y
   for iy in range(ny//2-lincells):
      deltaY = 1.*  ly / ny
      curY = deltaY/2. + deltaY*iy
      
      if lincells == 0:
         vyBlocs.append( Block2d(4, hx*nx/2 if nx>0 else hx, hz, NODE_COUNTER, ELEM_COUNTER, 3, 
                       ((0,curY,-h/2), (lx/2,curY,-h/2), (lx/2,curY,h/2), (0,curY,h/2) ), "y=%.2f Schale" % curY))
      else:
         vyBlocs.append( Block2d(4, sizeX, hz, NODE_COUNTER, ELEM_COUNTER, 3, 
                       ((0,curY,-h/2), (lx_inner,curY,-h/2), (lx_inner,curY,h/2), (0,curY,h/2) ), "y=%.2f Schale (reg.)" % curY))
         vyBlocs.append( Block2d(4, linSizeX, hz, NODE_COUNTER, ELEM_COUNTER, 6, 
                       ((lx_inner,curY,-h/2), (lx/2,curY,-h/2), (lx/2,curY,h/2), (lx_inner,curY,h/2) ), "y=%.2f Schale (lin.)" % curY))
         
      yLines.append(curY)
      
   for iy in range(ny//2-lincells, ny//2):
      deltaY = 1.*  ly / ny
      curY = deltaY/2. + deltaY*iy
      
      vyBlocs.append(   Block2d(4, hx*nx/2 if nx>0 else hx, hz, NODE_COUNTER, ELEM_COUNTER, 6, 
                       ((0,curY,-h/2), (lx/2,curY,-h/2), (lx/2,curY,h/2), (0,curY,h/2) ), "y=%.2f Schale (lin.)" % curY))
      
      yLines.append(curY)
   
#    if ny>0 and closeEdges:
#       curY = ly/2.
#       vyBlocs.append( Block2d(4, hx*nx/2 if nx>0 else hx, hz, NODE_COUNTER, ELEM_COUNTER, 3, 
#                     ((0,curY,-h/2), (lx/2,curY,-h/2), (lx/2,curY,h/2), (0,curY,h/2) ), "y=%.2f Schale" % curY))
#       yLines.append(curY)
    
   
   boundaries.append(BoundaryInput('vbou', '6. FHG komplett sperren', ((-eps,lx/2.+eps,-eps,ly/2.+eps,-h/2.-eps,h/2.+eps,0,0,0,0,0,1),)))
   
   print xLines
   print yLines
   
   
   # Freigabe 6.DOF entlang "Gitternetz" (Verschneidung mit vertikalen Schalen) an den Deckschichten oben und unten.
   edgeLines = []
   for xVal in xLines:
      print "EdgeLine: x=",xVal
      edgeLines.append( (xVal,0.,-h/2.,xVal,ly/2.,-h/2.,0,0,0,1000,1,0) )
      edgeLines.append( (0,0,0,0,0,0) )
      edgeLines.append( (xVal,0.,h/2.,xVal,ly/2.,h/2.,0,0,0,1000,1,0) )
      edgeLines.append( (0,0,0,0,0,0) )
      
   for yVal in yLines:
      print "EdgeLine: y=",yVal
      edgeLines.append( (0.,yVal,-h/2.,lx/2.,yVal,-h/2.,0,0,0,1000,1,0) )
      edgeLines.append( (0,0,0,0,0,0) )
      edgeLines.append( (0.,yVal,h/2.,lx/2.,yVal,h/2.,0,0,0,1000,1,0) )
      edgeLines.append( (0,0,0,0,0,0) )
      
   # Freigabe 6.DOF entlang Linien in z-Richtung, an denen sich die x-/y-Schalen treffen (Knotenpunkte des Gitters in der Draufsicht)
   edgeLinesZ = []
   for xVal in xLines:
      for yVal in yLines:
         print "EdgeLine Z: x=",xVal,"y=",yVal
         edgeLinesZ.append( (xVal,yVal,-h/2.,xVal,yVal,h/2.,0,0,0,1000,1,0) )
         edgeLinesZ.append( (0,0,0,0,0,0) )
    
   
   if len(edgeLines)>0:
      boundaries.append(BoundaryInput('edge', 'An Verschneidungen wieder freigeben', edgeLines))
   
   if len(edgeLinesZ)>0:
      boundaries.append(BoundaryInput('edge', 'An Verschneidungen wieder freigeben', edgeLinesZ))
   
   
   if rand=="Navier":
      boundaries.append(BoundaryInput('edge', 'z-Lagerung Rand',
                                   ((0., +ly/2., 0., +lx/2., +ly/2., 0., 0, 0, 0, 1000, 0, 0), (0,0,1,0,0,0),
                                   (lx/2.,0.,0.,lx/2.,ly/2.,0.,0,0,0,1000,0,0) , (0,0,1,0,0,0) )))
      
   elif rand=="NavierFullZ":
      boundaries.append(BoundaryInput('ebou', 'z-Lagerung Rand komplette Hoehe',
                                   ( (1, +lx/2., 0,0,1,0,0,0),
                                     (2, +ly/2., 0,0,1,0,0,0) )))
      filename = filename + "NavierZ"
      
   elif rand=="Einspannung":
      boundaries.append(BoundaryInput('ebou', 'Einspannung Rand',
                                   ( (1, +lx/2., 1,1,1,1,1,1),
                                     (2, +ly/2., 1,1,1,1,1,1) )))
      filename = filename + "Einsp"
   
   #boundaries.append(BoundaryInput('edge', 'Symmetriebed.',
   #                                ((0., 0., +lx/2., 0., 0., 0, 0, 0, 1000, 0, 0), (0,1,0,1,0,0),
   #                                (0., 0., 0., +ly/2., 0., 0, 0, 0, 1000, 0, 0), (1,0,0,0,1,0) )))
   
   
   #============================================================================
   # Symmetrie-Randbedingungen
   #============================================================================
   #
   # am Symmetrierand müssen nur an den Kanten der oberen und unteren
   # Deckschichten zusätzliche Rotationen gesperrt werden.
   # Die Zellwände sind so orientiert, dass die lokale z-Achse immer der
   # Symmetrieachse entspricht, es muss also der (ohnehin schon gesperrte)
   # 6. DOF gesperrt werden.
   #
   # Die Verschiebung normal zum Symmetrierand ist an der kompletten
   # Symmetrieebene gesperrt.
   # 
   #============================================================================
   
   boundaries.append(BoundaryInput('ebou', 'Symmetriebed. Verschiebungen',
                                   ((1, 0., 1,0,0,0,0,0),
                                    (2, 0., 0,1,0,0,0,0) )))
   
   boundaries.append(BoundaryInput('edge', 'Symmetriebed. Rot. X=0 oben',
                                   ((0., 0., +h/2., 0., +ly/2., +h/2., 0, 0, 0, 1000, 0, 0 ),
                                    (0, 0, 0, 0, 1, 0) )))
   boundaries.append(BoundaryInput('edge', 'Symmetriebed. Rot. X=0 unten',
                                   ((0., 0., -h/2., 0., +ly/2., -h/2., 0, 0, 0, 1000, 0, 0 ),
                                    (0, 0, 0, 0, 1, 0) )))
   boundaries.append(BoundaryInput('edge', 'Symmetriebed. Rot. Y=0 oben',
                                   ((0., 0., +h/2., +lx/2., 0., +h/2., 0, 0, 0, 1000, 0, 0 ),
                                    (0, 0, 0, 1, 0, 0) )))
   boundaries.append(BoundaryInput('edge', 'Symmetriebed. Rot. Y=0 unten',
                                   ((0., 0., -h/2., +lx/2., 0., -h/2., 0, 0, 0, 1000, 0, 0 ),
                                    (0, 0, 0, 1, 0, 0) )))
                      
   #boundaries.append(BoundaryInput('ebou', 'Symmetriebed.',
   #                                ((1, 0., 1,0,0,0,1,1),
   #                                (2, 0., 0,1,0,1,0,1) )))
   
   
   
#    if nx == 0:
#       boundaries.append(BoundaryInput('poin', 'x-/y-Lagerung Ecke',
#                                    ((+lx/2., +ly/2.-ly/ny/2., 0., 10),
#                                    (0,0,0), (1,1,0,0,0,0))))
#    else:
#       boundaries.append(BoundaryInput('poin', 'x-/y-Lagerung Ecke',
#                                    ((+lx/2.-lx/nx/2., +ly/2., 0., 10),
#                                    (0,0,0), (1,1,0,0,0,0))))
   
   with open(filename,'w') as file:
      file.write(str(header))
      
      if top: file.write(str(topBloc))
      if btm: file.write(str(btmBloc))
      
      if top and lincells>0:
         file.write(str(linTopBlocX))
         file.write(str(linTopBlocY))
         file.write(str(linTopBlocXY))
      if btm and lincells>0:
         file.write(str(linBtmBlocX))
         file.write(str(linBtmBlocY))
         file.write(str(linBtmBlocXY))
      
      for bloc in vxBlocs:
         file.write(str(bloc))
      for bloc in vyBlocs:
         file.write(str(bloc))
         
      file.write('\n\n')
      
      #file.write('epsq\n')
      #file.write('3,8,2,0,0\n0.001,0,0,0,0,0,0,0')
         
      file.write('\n\n')
      for b in boundaries:
         file.write(str(b)+'\n')
         
#       file.write('eloa\n')
#       file.write('%f,%f,%f,%f,%f,%f,0,0,0,1000\n' % (-lx/6., -ly/2., h/2., -lx/6., ly/2., h/2.))
#       file.write('1,1,0,0,0,0,%f\n' % f)
#       file.write('%f,%f,%f,%f,%f,%f,0,0,0,1000\n' % (lx/6., -ly/2., h/2., lx/6., ly/2., h/2.))
#       file.write('1,1,0,0,0,0,%f\n' % f)
#       file.write('\n')
      
      file.write('aloa\n')
      file.write('%f,%f,%f,%f,%f,%f\n' % ( -lx/2. - eps, lx/2. + eps, -ly/2. - eps, ly/2. + eps, h/2.-eps, h/2.+eps))
      file.write('%i,%i\n' % (3, 4))
      file.write('%i,%f\n' % (1, f))
      file.write('\n')
      
      if len(materials)>0:
         for mate in materials:
            file.write('\n')
            file.write('mate\n')
            file.write(str(mate))
            
      file.write(str(footer))
      
   print "File written:", filename
   sys.exit(0)
   
if __name__=="__main__":
   main()