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
   def __init__(self, name, params, solv=4, nopr=False):
      self.name = name
      self.params = params
      self.solver = solv
      self.nopr = nopr
      
   def __str__(self):
      s = "feap   ***%s***\n" % self.name
      for p in self.params:
         s += "%s," % str(p)
      s = s[:-1] + "\n" # cut last comma, add newline
      if self.nopr: s+= "\nnopr\n"
      s += "\nsolv\n%i\n\n" % self.solver
      return s
      

class Footer:
   def __init__(self, tie=True, macr=(), macrBlanks=0, link=(), noInteStop=False):
      self.tie = tie
      self.macr = macr
      self.link = link
      self.macrBlanks = macrBlanks
      self.noInteStop = noInteStop
      
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
      
      if not self.noInteStop: s+= "inte\nstop\n\n"
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
   global MATE_COUNTER
   
   #mode = "einzeln"
   mode = "FE2"
   
   # generate simplified "honeycomb" mesh
   
   segments = []
   materials = []
   boundaries = []
   
   additionalElements = []
   
   # create mesh
   
   # Vorsicht: Um Mittelknoten lagern zu können, müssen die Summen (nx+hx) und (ny+hy) gerade Zahlen (z.B. 2) sein!
   
   nx = 0 # nr. of vertically aligned shells in y direction (x=const.)
   ny = 1 # nr. of vertically aligned shells in x direction (y=const.)

   hx = 2 # elements between two vertical shells (x direction)
   hy = 2 # elements between two vertical shells (y direction)
   
   # gerade Zahl wg. Mittenknotenlagerung
   hz = 2 # elements in thickness direction
   
   closeEdges = False # set to True to add vertical shells on the outer edges of the system 
   
   # global dimensions               [kN,cm]
   # !!!unbedingt mit Dezimalpunkt eingeben!!!
   #lx = 0.5
   ly = 0.5*ny  # ly_Mikro = ly/ny|_Makro
   lx = ly 
   h  = 3.
   
   # element lengths
#    dx = lx / (hx*(nx+1))
#    dy = ly / (hy*(ny+1))
#    dz = h / hz
   
   # parameters
   li = 0       # linear geometry?  0=lin,1=nili (moderate),2=nili (finite) 
   t = 0.1      # shell thickness
   e = 7000.    # Young's modulus
   v = 0.       # Poisson ratio
   #f = -0.05   # loading force
   
   eps = (lx+ly+h)*1.e-4 # small length parameter
   
   # CHOOSE MATERIAL
   
   # LINEAR MATERIAL
   
   matn = 1
   mateCard = (e,v,0)
   nmnbnz = '000'
   if mode == "FE2" and li == 0:
      filename = "C:/FG/feap/exe/rve/microShell/balkenLin/iFE2BalkenRVELin/irves8_01"
   elif mode == "einzeln" and li == 0:
      filename = "C:/FG/feap/exe/rve/microShell/iFE2BalkenRVELin_einzeln"
   elif mode == "FE2" and li>0:
      filename = "C:/FG/feap/exe/rve/microShell/balkenGeomNL/iFE2BalkenRVEGeomNL/irves8_01"
   elif mode == "einzeln" and li>0:
      filename = "C:/FG/feap/exe/rve/microShell/iFE2BalkenRVEGeomNL_einzeln"
   
   # NONLINEAR MATERIAL
   
#    y0 = 24.
#    xk = 100.
#    matn = 4
#    mateCard = (e,v,y0,0,xk,0,0)
#    nmnbnz = '002'
#    if mode == "FE2" and li == 0:
#       filename = "C:/FG/feap/exe/rve/microShell/balkenPhysNL/iFE2BalkenRVEPhysNL/irves8_01"
#    elif mode == "einzeln" and li == 0:
#       filename = "C:/FG/feap/exe/rve/microShell/balkenPhysNL/iFE2BalkenRVEPhysNL_einzeln"
#    if mode == "FE2" and li>0:
#       filename = "C:/FG/feap/exe/rve/microShell/balkenFullNL/iFE2BalkenRVEFullNL/irves8_01"
#    elif mode == "einzeln" and li>0:
#       filename = "C:/FG/feap/exe/rve/microShell/balkenFullNL/iFE2BalkenRVEFullNL_einzeln"
   
   # END MATERIAL
   
   materials.append(Material(35, (
                  (1,0,0,0),
                  (1,4,t,-0.5*t,matn,li,5,nmnbnz,1,3),
                  (0,),
                  mateCard,
                  (0,t)     )))
   
   if li==0 and matn == 1: caption = "Linear"
   elif li>0 and matn == 1: caption = "Geom. Nichtlin. li=%i" % li
   elif li ==0 and matn == 4: caption = "Phys. Nichtlin. mate=%i" % matn
   else: caption ="Geom. Nichtlin. li=%i, Phys. Nichtlin. mate=%i" % (li, matn)
   header = Header("Schalen Sandwich Mikroelement %s" % caption, ('','','',3,6,4), nopr=(mode=="FE2"))
   #footer = Footer(tie=True, macr=('plot,rot1,-80', 'plot,rot3,15', 'plot,mesh', 'plot,axis', 'plot,boun', 'prop', 'dt,,0.01'), macrBlanks=1,
   #                link=('6,1,%f,%f,1,1,0,1,1,0' % (-lx/2.,lx/2.), '6,2,%f,%f,1,1,0,1,1,0' % (-ly/2.,ly/2.) ))
   
   footer = Footer(tie=True, macr=(), noInteStop=True, link=('6,1,%f,%f,1,1,0,1,1,1' % (-lx/2.,lx/2.), '6,2,%f,%f,1,1,0,1,1,1' % (-ly/2.,ly/2.) ))
   
   
   if mode == "FE2":
      footer = str(footer) + """
   
   
batc,tang
nopr
prop,,1
rest,,0,
epsq 
loop,,2
tang,,1 
next    
sigq
end,,0,
1,1,0,10000,1

stop,tang


batc,updh
nopr
prop,,1
rest,,0,
updh,,2
end,,0,
1,1,0,10000,1

stop,updh


batc,micr
nopr
prop,,1
rest,,0
plot,pers
plot,mesh 
plot,axis
end,,0



inte
stop,micr

"""
   
   elif mode == "einzeln":
      footer = str(footer) + """

macr
epsq,,1
plot,rot1,-65
plot,rot3,55
plot,mesh
plot,axis
plot,boun
end
inte
stop

"""
   
   
   # top and bottom plates
   top = True
   btm = True
   
#    if top: topBloc = Block2d(4, hx*(nx+1), hy*(ny+1), NODE_COUNTER, ELEM_COUNTER, 1, 
#                     ((-lx/2.,-ly/2.,h/2.), (lx/2.,-ly/2.,h/2.), (lx/2.,ly/2.,h/2.), (-lx/2.,ly/2.,h/2.) ), "Deckschicht oben")
#    if btm: btmBloc = Block2d(4, hx*(nx+1), hy*(ny+1), NODE_COUNTER, ELEM_COUNTER, 1, 
#                     ((-lx/2.,-ly/2.,-h/2.), (lx/2.,-ly/2.,-h/2.), (lx/2.,ly/2.,-h/2.), (-lx/2.,ly/2.,-h/2.) ), "Deckschicht unten")
   if top: topBloc = Block2d(4, hx*nx if nx > 0 else hx, hy*ny if ny > 0 else hy, NODE_COUNTER, ELEM_COUNTER, 1, 
                    ((-lx/2.,-ly/2.,h/2.), (lx/2.,-ly/2.,h/2.), (lx/2.,ly/2.,h/2.), (-lx/2.,ly/2.,h/2.) ), "Deckschicht oben")
   if btm: btmBloc = Block2d(4, hx*nx if nx > 0 else hx, hy*ny if ny > 0 else hy, NODE_COUNTER, ELEM_COUNTER, 1, 
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
   
   # vertikal ausgerichtete Blöcke, konstantes Y
   for iy in range(ny):
      deltaY = 1.*  ly / ny
      curY = -ly/2. + deltaY/2. + deltaY*iy
       
      vyBlocs.append( Block2d(4, hx*nx if nx>0 else hx, hz, NODE_COUNTER, ELEM_COUNTER, 1, 
                    ((-lx/2,curY,-h/2), (lx/2,curY,-h/2), (lx/2,curY,h/2), (-lx/2,curY,h/2) ), "y=%.2f Schale" % curY))
      yLines.append(curY)
   
#    # vertikal ausgerichtete Bl�cke, konstantes X
#    for ix in range(1,nx+1):
#       deltaX = 1.*lx / (nx+1)
#       curX = -lx/2. + deltaX*ix
#       
#       vxBlocs.append(  Block2d(4, hy*(ny+1), hz, NODE_COUNTER, ELEM_COUNTER, 1, 
#                     ((curX,-ly/2,-h/2), (curX,ly/2,-h/2), (curX,ly/2,h/2), (curX,-ly/2,h/2) ), "x=%.2f Schale" % curX))
#       xLines.append(curX)
#    
#    # vertikal ausgerichtete Bl�cke, konstantes Y
#    for iy in range(1,ny+1):
#       deltaY = 1.*ly / (ny+1)
#       curY = -ly/2. + deltaY*iy
#        
#       vyBlocs.append( Block2d(4, hx*(nx+1), hz, NODE_COUNTER, ELEM_COUNTER, 1, 
#                     ((-lx/2,curY,-h/2), (lx/2,curY,-h/2), (lx/2,curY,h/2), (-lx/2,curY,h/2) ), "y=%.2f Schale" % curY))
#       yLines.append(curY)
    
   
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
   
   
#    boundaries.append(BoundaryInput('edge', 'Randlagerung',
#                                    ((-lx/2.,-ly/2.,-h/2.,lx/2.,-ly/2.,-h/2.,0,0,0,1000,1,0), (1,1,1,1,1,1),
#                                    (-lx/2.,-ly/2.,-h/2.,-lx/2.,ly/2.,-h/2.,0,0,0,1000,1,0), (1,1,1,1,1,1),
#                                    (-lx/2.,ly/2.,-h/2.,lx/2.,ly/2.,-h/2.,0,0,0,1000,1,0), (1,1,1,1,1,1),
#                                    (lx/2.,-ly/2.,-h/2.,lx/2.,ly/2.,-h/2.,0,0,0,1000,1,0), (1,1,1,1,1,1))))

   boundaries.append(BoundaryInput('ebou', 'x-y Lagerung an Aussenkanten',
                                   ((1,-lx/2.,1,1,0), (1,lx/2.,1,1,0), (2,-ly/2.,1,1,0), (2,ly/2.,1,1,0)) ))
   boundaries.append(BoundaryInput('poin', 'z-Lagerung Mittelknoten', ((0.,0.,0.,20.),(0,0,0,0,0,0),(0,0,1,0,0,0)) ))
   
   with open(filename,'w+') as file:
      file.write(str(header))
      
      if top: file.write(str(topBloc))
      if btm: file.write(str(btmBloc))
      
      for bloc in vxBlocs:
         file.write(str(bloc))
      for bloc in vyBlocs:
         file.write(str(bloc))
         
      file.write('\n\n')
      
      if mode == "einzeln":
         file.write('epsq\n')
         file.write('3,8,2,1,1\n0,0,0,0,0,0,0.02,0')
         
      file.write('\n\n')
      for b in boundaries:
         file.write(str(b)+'\n')
         
#       file.write('aloa\n')
#       file.write('%f,%f,%f,%f,%f,%f\n' % (-lx/2.,lx/2.,-ly/2.,ly/2.,h/2., h/2.))
#       file.write('3,4\n')
#       file.write('1,%s\n' % f ) # load
      
      
#      file.write('poin\n')
#      file.write('5,0,0,0,-%f,0,0,0\n' % f )
      
      if len(materials)>0:
         for mate in materials:
            file.write('\n')
            file.write('mate\n')
            file.write(str(mate))
            
      file.write(str(footer))
      
      
   print "%s written." % filename.replace('/','\\')

   sys.exit(0)
   
if __name__=="__main__":
   main()