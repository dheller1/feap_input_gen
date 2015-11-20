# -*- coding: utf-8 -*-
'''
Created on 22.02.2013

@author: heller
'''
import sys

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
      if self.tie: s+= "tie\n"
      
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
   
   
   # generate simplified "honeycomb" mesh
   
   segments = []
   materials = []
   boundaries = []
   
   additionalElements = []
   
   # create mesh
   
   # if closeEdges is True, there are two more outer shells per direction
   nx = 0 # nr. of vertically aligned shells in y direction (x=const.)
   ny = 12 # nr. of vertically aligned shells in x direction (y=const.)
   
   # hx/hy müssen Vielfache von 2 sein, damit man Stege in der Mitte eines "RVEs" mit der Deckschicht noch erwischt, außerdem
   # falls nx > 0 muss hy Vielf. von 6 sein, falls ny > 0 muss hx Vielf. von 6 sein !
   hx = 24 # elements between two vertical shells (x direction) // total elements in x-direction if nx <= 1
   hy = 2 # elements between two vertical shells (y direction) // total elements in y-direction if ny <= 1
   
   hz = 2 # elements in thickness direction
   
   closeEdges = False # set to True to add vertical shells on the outer edges of the system 
   
   # global dimensions               [kN,cm]
   lx = 100.
   ly = 6.
   h  = 3.
   
   # element lengths
   #dx = lx / (hx*(nx+1))
   #dy = ly / (hy*(ny+1))
   #dz = h / hz
   
   # parameters
   li = 0      # linear geometry?  0=lin,1=nili (moderate),2=nili (finite) 
   t = 0.1     # shell thickness
   e = 7000.   # Young's modulus
   v = 0.#34    # Poisson ratio
   f = -0.1    # loading force
   
   y0 = 24.    # initial yield stress (only for nonlinear material)
   xk = 100.   # linear hardening modulus (only for nonlinear material)
   
   eps = (lx+ly+h)*1.e-4 # small length parameter
   
   mat1 = Material(35, (
                  (1,0,0,0),
                  (1,2,t,-0.5*t,1,li,5,'000',1,3),
                  (0,),
                  (e,v,0),
                  (0,t)     ))
   
#    mat1 = Material(35, (
#                   (1,0,0,0),
#                   (1,2,t,-0.5*t,4,li,5,'002',1,3),
#                   (0,),
#                   (e,v,y0,0,xk,0,0),
#                   (0,t)     ))
   
   materials.append(mat1)
   
   header = Header("Schalen Sandwich Mikroelement", ('','','',3,6,4), solv=4)
   
   tploNode = 963
   footer = Footer(tie=True, macr=('prop','plot,rot1,-80', 'plot,rot3,15', 'plot,mesh', 'plot,axis', 'plot,boun,,1', 'plot,boun,,2', 'plot,boun,,3', 'plot,load',
                                   'plot,node,,-%i' % tploNode, 'tplo,init,%i,-3,1' % tploNode, 'tplo,mark,2,4,1'), macrBlanks=1 )
   
   filename = "C:\\FG\\feap\\exe\\rve\\microShell\\Faltwerke\\iBalken"
   
   # top and bottom plates
   top = True
   btm = True
   
   if top: topBloc = Block2d(4, hx*nx if nx > 0 else hx, hy*ny if ny > 0 else ny, NODE_COUNTER, ELEM_COUNTER, 1, 
                    ((-lx/2.,-ly/2.,h/2.), (lx/2.,-ly/2.,h/2.), (lx/2.,ly/2.,h/2.), (-lx/2.,ly/2.,h/2.) ), "Deckschicht oben")
   if btm: btmBloc = Block2d(4, hx*nx if nx > 0 else hx, hy*ny if ny > 0 else ny, NODE_COUNTER, ELEM_COUNTER, 1, 
                    ((-lx/2.,-ly/2.,-h/2.), (lx/2.,-ly/2.,-h/2.), (lx/2.,ly/2.,-h/2.), (-lx/2.,ly/2.,-h/2.) ), "Deckschicht unten")
   
   vxBlocs = []
   vyBlocs = []
   xLines = [] # speichere x/y Koordinaten des Gitters, dies sind (zumindest oben und unten) Verschneidungskanten wo 6.DOF freigegeben werden muss
   yLines = []
   
   # vertikal ausgerichtete Blöcke, konstantes X
   for ix in range(nx):
      deltaX = 1.*  lx / nx
      curX = -lx/2. + deltaX/2. + deltaX*ix
      
      vxBlocs.append(  Block2d(4, hy*ny if ny>0 else hy, hz, NODE_COUNTER, ELEM_COUNTER, 1, 
                    ((curX,-ly/2,-h/2), (curX,ly/2,-h/2), (curX,ly/2,h/2), (curX,-ly/2,h/2) ), "x=%.2f Schale" % curX))
      xLines.append(curX)
      
   if nx>0 and closeEdges:
      curX = -lx/2.
      vxBlocs.append(  Block2d(4, hy*ny if ny>0 else hy, hz, NODE_COUNTER, ELEM_COUNTER, 1, 
                    ((curX,-ly/2,-h/2), (curX,ly/2,-h/2), (curX,ly/2,h/2), (curX,-ly/2,h/2) ), "x=%.2f Schale (Abschluss)" % curX))
      xLines.append(curX)
      
      curX = lx/2.
      vxBlocs.append(  Block2d(4, hy*ny if ny>0 else hy, hz, NODE_COUNTER, ELEM_COUNTER, 1, 
                    ((curX,-ly/2,-h/2), (curX,ly/2,-h/2), (curX,ly/2,h/2), (curX,-ly/2,h/2) ), "x=%.2f Schale (Abschluss)" % curX))
      xLines.append(curX)
   
   # vertikal ausgerichtete Blöcke, konstantes Y
   for iy in range(ny):
      deltaY = 1.*  ly / ny
      curY = -ly/2. + deltaY/2. + deltaY*iy
       
      vyBlocs.append( Block2d(4, hx*nx if nx>0 else hx, hz, NODE_COUNTER, ELEM_COUNTER, 1, 
                    ((-lx/2,curY,-h/2), (lx/2,curY,-h/2), (lx/2,curY,h/2), (-lx/2,curY,h/2) ), "y=%.2f Schale" % curY))
      yLines.append(curY)
   
   if ny>0 and closeEdges:
      curY = -ly/2.
      vyBlocs.append( Block2d(4, hx*nx if nx>0 else hx, hz, NODE_COUNTER, ELEM_COUNTER, 1, 
                    ((-lx/2,curY,-h/2), (lx/2,curY,-h/2), (lx/2,curY,h/2), (-lx/2,curY,h/2) ), "y=%.2f Schale" % curY))
      yLines.append(curY)
      
      curY = ly/2.
      vyBlocs.append( Block2d(4, hx*nx if nx>0 else hx, hz, NODE_COUNTER, ELEM_COUNTER, 1, 
                    ((-lx/2,curY,-h/2), (lx/2,curY,-h/2), (lx/2,curY,h/2), (-lx/2,curY,h/2) ), "y=%.2f Schale" % curY))
      yLines.append(curY)
    
   
   boundaries.append(BoundaryInput('vbou', '6. FHG komplett sperren', ((-lx/2.-eps,lx/2.+eps,-ly/2.-eps,ly/2.+eps,-h/2.-eps,h/2.+eps,0,0,0,0,0,1),)))
   
   print xLines
   print yLines
   
   
   # Freigabe 6.DOF entlang "Gitternetz" (Verschneidung mit vertikalen Schalen) an den Deckschichten oben und unten.
   edgeLines = []
   for xVal in xLines:
      print "EdgeLine: x=",xVal
      edgeLines.append( (xVal,-ly/2.,-h/2.,xVal,ly/2.,-h/2.,0,0,0,1000,1,0) )
      edgeLines.append( (0,0,0,0,0,0) )
      edgeLines.append( (xVal,-ly/2.,h/2.,xVal,ly/2.,h/2.,0,0,0,1000,1,0) )
      edgeLines.append( (0,0,0,0,0,0) )
      
   for yVal in yLines:
      print "EdgeLine: y=",yVal
      edgeLines.append( (-lx/2.,yVal,-h/2.,lx/2.,yVal,-h/2.,0,0,0,1000,1,0) )
      edgeLines.append( (0,0,0,0,0,0) )
      edgeLines.append( (-lx/2.,yVal,h/2.,lx/2.,yVal,h/2.,0,0,0,1000,1,0) )
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
   
   boundaries.append(BoundaryInput('edge', 'z-Lagerung Rand',
                                   ((-lx/2., -ly/2., 0., -lx/2., ly/2., 0., 0, 0, 0, 1000, 0, 0), (0,0,1,0,0,0),
                                    (lx/2.,-ly/2.,0.,lx/2.,ly/2.,0.,0,0,0,1000,0,0) , (0,0,1,0,0,0) )))
   
   
   boundaries.append(BoundaryInput('poin', 'x-/y-Lagerung Ecke',
                                   ((-lx/2., -ly/2.+ly/ny/2.,0., 10),
                                   (0,0,0), (1,1,0,0,0,0),
                                   (lx/2.,-ly/2.+ly/ny/2.,0., 10),
                                   (0,0,0), (0,1,0,0,0,0))))
   
   with open(filename,'w') as file:
      file.write(str(header))
      
      if top: file.write(str(topBloc))
      if btm: file.write(str(btmBloc))
      
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
         
      file.write('eloa\n')
      file.write('%f,%f,%f,%f,%f,%f,0,0,0,1000\n' % (-lx/6., -ly/2., h/2., -lx/6., ly/2., h/2.))
      file.write('1,1,0,0,0,0,%f\n' % f)
      file.write('%f,%f,%f,%f,%f,%f,0,0,0,1000\n' % (lx/6., -ly/2., h/2., lx/6., ly/2., h/2.))
      file.write('1,1,0,0,0,0,%f\n' % f)
      file.write('\n')
      
      if len(materials)>0:
         for mate in materials:
            file.write('\n')
            file.write('mate\n')
            file.write(str(mate))
            
      file.write(str(footer))
      
      print "Written %s." % filename
      
   sys.exit(0)
   
if __name__=="__main__":
   main()