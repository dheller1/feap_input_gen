# -*- coding: utf-8 -*-
'''
Created on 22.02.2013

@author: heller
'''
import sys, math, os
from math import sqrt

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
      if self.tie: s+= "tie,nopr\ntie\ntie,prin\n"
      #if self.tie: s+= "tie\n"
      
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
      self.id = ELEM_COUNTER
      ELEM_COUNTER += 1
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
      self.id = NODE_COUNTER
      NODE_COUNTER += 1
      self.x = x
      self.y = y
      self.z = z
      
   def __str__(self):
      return "%i,0,%f, %f, %f\n" % (self.id, self.x,self.y,self.z)
   
   
class HoneycombCell:
   def __init__(self, x=0., y=0., z=0., lx=2., ly=2., h=1., hx=1, hz=1, sectPt=(1./3, -1./3)):
      """ (x,y) == bottom left corner """
      global NODE_COUNTER
      global ELEM_COUNTER
      
      self.blocks = []
      self.edgeLines = []
      self.sectCorners = []
      
      mid_x,mid_y = (x+(1.+sectPt[0])*lx/2, y+(1.+sectPt[1])*ly/2) # mid node coordinates
      tr_x,tr_y = (x+lx, y+ly) # top right corner
      br_x,br_y = (x+lx, y)    # bottom right corner
      bl_x,bl_y = (x, y)       # bottom left corner
      z1, z2 = z-h/2, z+h/2
      
      self.blocks.append( Block2d(4, hx, hz, NODE_COUNTER, ELEM_COUNTER, 1, 
            ((mid_x,mid_y,z1),(tr_x,tr_y,z1),(tr_x,tr_y,z2),(mid_x,mid_y,z2)), "Steg I"))
      
      self.blocks.append( Block2d(4, hx, hz, NODE_COUNTER, ELEM_COUNTER, 1, 
            ((mid_x,mid_y,z1),(br_x,br_y,z1),(br_x,br_y,z2),(mid_x,mid_y,z2)), "Steg II"))
      
      self.blocks.append( Block2d(4, hx, hz, NODE_COUNTER, ELEM_COUNTER, 1, 
            ((mid_x,mid_y,z1),(bl_x,bl_y,z1),(bl_x,bl_y,z2),(mid_x,mid_y,z2)), "Steg III"))
      
      # free 6th DOF at mid-node intersection
      self.edgeLines.append(  (mid_x,mid_y,z1,mid_x,mid_y,z2,0,0,0,1000,1,0) )
      self.edgeLines.append(  (0,0,0,0,0,0)  )
      
      # save "touched" corner nodes (x,y) where intersections occur when combined with other honeycomb segments
      self.sectCorners.extend( ((tr_x, tr_y), (br_x,br_y), (bl_x,bl_y)) )
      
      
#       
# class OverLaySimple:
#    def __init__(self, x=0., y=0., z=0., lx=10., ly=10., h=1., olx=1., oly=1.):
#       """ (x,y) == center """
#       global NODE_COUNTER
#       global ELEM_COUNTER
#         
#       self.blocks = []
#       self.edgeLines = []
#         
#       tr_x,tr_y = (x+lx/2., y+ly/2.) # top right corner
#       tl_x,tl_y = (x-lx/2., y+ly/2.) # top left corner
#       br_x,br_y = (x+lx/2., y-ly/2.) # bottom right corner
#       bl_x,bl_y = (x-lx/2., y-ly/2.) # bottom left corner
#       z1, z2 = z-h/2., z+h/2.
#         
#       self.blocks.append( Block2d(4, olx, oly, NODE_COUNTER, ELEM_COUNTER, 1, 
#             ((bl_x,bl_y,z2),(br_x,br_y,z2),(tr_x,tr_y,z2),(tl_x,tl_y,z2)), "Deckschicht oben"))
#         
#       self.blocks.append( Block2d(4, olx, oly, NODE_COUNTER, ELEM_COUNTER, 1, 
#             ((bl_x,bl_y,z1),(br_x,br_y,z1),(tr_x,tr_y,z1),(tl_x,tl_y,z1)), "Deckschicht unten"))
#        
# #       # free 6th DOF at mid-node intersection
# #       self.edgeLines.append(  (mid_x,mid_y,z1,mid_x,mid_y,z2,0,0,0,1000,1,0) )
# #       self.edgeLines.append(  (0,0,0,0,0,0)  )
# #        
#       # save "touched" corner nodes (x,y) where intersections might occur when combined with other honeycomb segments
#       self.sectCorners.extend( ((tr_x, tr_y), (br_x,br_y), (bl_x,bl_y)) )


class OverLay:
   def __init__(self, x=0., y=0., z=0., lx=2., ly=2., h=1., sectPt=(1./3, -1./3)):
      """ (x,y) == bottom left corner """
      global NODE_COUNTER
      global ELEM_COUNTER
        
      #self.blocks = []
      self.edgeLines = []
        
      tr_x,tr_y = (x+lx, y+ly) # top right corner
      tl_x,tl_y = (x, y+ly) # top left corner
      br_x,br_y = (x+lx, y) # bottom right corner
      bl_x,bl_y = (x, y) # bottom left corner
      mid_x,mid_y = (x+lx/2.+sectPt[0]/(lx/2.), y+ly/2.+sectPt[1]/(ly/2.))
      z1, z2 = z-h/2., z+h/2.
      
      nodes = self.nodes = [ Node(bl_x, bl_y, z2), #0
                Node(bl_x, bl_y, z1), #1
                Node(br_x, br_y, z1), #2
                Node(br_x, br_y, z2), #3
                Node(tl_x, tl_y, z2), #4
                Node(tl_x, tl_y, z1), #5
                Node(tr_x, tr_y, z1), #6
                Node(tr_x, tr_y, z2), #7
                Node(mid_x, mid_y, z2), #8
                Node(mid_x, mid_y, z1), #9
                Node(mid_x, bl_y, z2), #10
                Node(br_x, mid_y, z2), #11
                Node(mid_x, tl_y, z2), #12
                Node(bl_x, mid_y, z2), #13
                Node(mid_x, bl_y, z1), #14
                Node(br_x, mid_y, z1), #15
                Node(mid_x, tl_y, z1), #16
                Node(bl_x, mid_y, z1), #17
               ]
      
      self.elems = [ Element( nodes[0], nodes[10], nodes[8], nodes[13], ), # A
                     Element( nodes[10], nodes[3], nodes[11], nodes[8], ), # B
                     Element( nodes[8], nodes[11], nodes[7], nodes[12], ), # C
                     Element( nodes[13], nodes[8], nodes[12], nodes[4], ), # D
                     Element( nodes[1], nodes[14], nodes[9], nodes[17], ), # E
                     Element( nodes[14], nodes[2], nodes[15], nodes[9], ), # F
                     Element( nodes[9], nodes[15], nodes[6], nodes[16], ), # G
                     Element( nodes[17], nodes[9], nodes[16], nodes[5], ), # H
               ]
      
class OverLayInter:
   def __init__(self, x=0., y=0., z=0., lx=2., ly=2., h=1., sectPt=(1./3, -1./3)):
      """ (x,y) == bottom left corner """
      global NODE_COUNTER
      global ELEM_COUNTER
        
      #self.blocks = []
      self.edgeLines = []
        
      tr_x,tr_y = (x+lx, y+ly) # top right corner
      tl_x,tl_y = (x, y+ly) # top left corner
      br_x,br_y = (x+lx, y) # bottom right corner
      bl_x,bl_y = (x, y) # bottom left corner
      mid_x,mid_y = (x+lx/2.+sectPt[0]/(lx/2.), y+ly/2.+sectPt[1]/(ly/2.)) # mid node
      
      
      int1_x,int1_y = (bl_x+mid_x)/2, (bl_y+mid_y)/2 # steg 1 internal midpoint
      int2_x,int2_y = (br_x+mid_x)/2, (br_y+mid_y)/2 # steg 2 internal midpoint
      int3_x,int3_y = (tr_x+mid_x)/2, (tr_y+mid_y)/2 # steg 3 internal midpoint
      
      z1, z2 = z-h/2., z+h/2.
      
      topNodes = []
      btmNodes = []
      topElems = []
      btmElems = []
      
      self.sectNodes = []
      
      for i in range(2): 
         if i == 0: z = z2
         elif i == 1: z = z1
            
         nds = [ # ROW 1 (y = tl_y)
             Node(tl_x,   tl_y, z), #0 
             Node(int1_x, tl_y, z), #1 
             Node(mid_x,  tl_y, z), #2 
             Node(int2_x, tl_y, z), #3 
             Node(tr_x,   tl_y, z), #4 
             
             # ROW 2 (y = int3_y)
             Node(tl_x,   int3_y, z), #5 
             Node(int1_x, int3_y, z), #6 
             Node(mid_x,  int3_y, z), #7 
             Node(int2_x, int3_y, z), #8 
             Node(tr_x,   int3_y, z), #9 
             
             # ROW 3 (y = mid_y)
             Node(tl_x,   mid_y, z), #10
             Node(int1_x, mid_y, z), #11
             Node(mid_x,  mid_y, z), #12
             Node(int2_x, mid_y, z), #13
             Node(tr_x,   mid_y, z), #14
             
             # ROW 4 (y = int1_y)
             Node(tl_x,   int1_y, z), #15
             Node(int1_x, int1_y, z), #16
             Node(mid_x,  int1_y, z), #17
             Node(int2_x, int1_y, z), #18
             Node(tr_x,   int1_y, z), #19
             
             # ROW 5 (y = bl_y)
             Node(tl_x,   bl_y, z), #20
             Node(int1_x, bl_y, z), #21
             Node(mid_x,  bl_y, z), #22
             Node(int2_x, bl_y, z), #23
             Node(tr_x,   bl_y, z)] #24
         
         self.sectNodes.extend((nds[8],nds[16],nds[18],)) # innere Stegknoten müssen später stets im 6.DOF freigegeben werden
      
         els = [Element( nds[0], nds[5], nds[6], nds[1], ), # A
             Element( nds[1], nds[6], nds[7], nds[2], ), # B
             Element( nds[2], nds[7], nds[8], nds[3], ), # C
             Element( nds[3], nds[8], nds[9], nds[4], ), # D
             
             Element( nds[5], nds[10], nds[11], nds[6], ), # E
             Element( nds[6], nds[11], nds[12], nds[7], ), # F
             Element( nds[7], nds[12], nds[13], nds[8], ), # G
             Element( nds[8], nds[13], nds[14], nds[9], ), # H
             
             Element( nds[10], nds[15], nds[16], nds[11], ), # I
             Element( nds[11], nds[16], nds[17], nds[12], ), # J
             Element( nds[12], nds[17], nds[18], nds[13], ), # K
             Element( nds[13], nds[18], nds[19], nds[14], ), # L
             
             Element( nds[15], nds[20], nds[21], nds[16], ), # M
             Element( nds[16], nds[21], nds[22], nds[17], ), # N
             Element( nds[17], nds[22], nds[23], nds[18], ), # O
             Element( nds[18], nds[23], nds[24], nds[19], )] # P
         
         if i == 0:
            topNodes = nds
            topElems = els
         elif i == 1:
            btmNodes = nds
            btmElems = els
      
      self.nodes = [] 
      self.elems = []
      self.nodes.extend(topNodes)
      self.nodes.extend(btmNodes)
      
      self.elems.extend(topElems)
      self.elems.extend(btmElems)

class OverLayHx:
   def __init__(self, x=0., y=0., z=0., lx=2., ly=2., h=1., hx=2, sectPt=(1./3, -1./3)):
      """ (x,y) == bottom left corner """
      global NODE_COUNTER
      global ELEM_COUNTER
        
      #self.blocks = []
      self.edgeLines = []
        
      tr_x,tr_y = (x+lx, y+ly) # top right corner
      tl_x,tl_y = (x, y+ly) # top left corner
      br_x,br_y = (x+lx, y) # bottom right corner
      bl_x,bl_y = (x, y) # bottom left corner
      mid_x,mid_y = (x+lx/2.+sectPt[0]/(lx/2.), y+ly/2.+sectPt[1]/(ly/2.)) # mid node
      
      z1, z2 = z-h/2., z+h/2.
      
      w = 2 * hx + 1 # number of nodes in both x- and y-directions.
      
      # global lengths
      width_a = mid_x - bl_x # length from mid point to left side of the cell
      width_b = br_x - mid_x # length from mid point to right side of the cell
      height_a = mid_y - bl_y # length from mid point to bottom side of the cell
      height_b = tl_y - mid_y # length from mid point to top side of the cell  
      
      # generate coordinate grid
      cx = [] # coordinates cx[0], ... , cx[w-1] for nodes in column j
      cy = [] # coordinates cy[0], ... , cy[w-1] for nodes in row i
      
      for i in range(w):
         if(i<=hx): # bottom / left of mid node
            cx.append ( bl_x + width_a * (1.*i/hx) )
            cy.append ( bl_y + height_a * (1.*i/hx) )
         else: # top / right of mid node
            cx.append ( mid_x + width_b * (1.*i/hx - 1) )
            cy.append ( mid_y + height_b * (1.*i/hx - 1) )
            
      # generate nodes and elements
      topNodes = []
      btmNodes = []
      topElems = []
      btmElems = []
      
      self.sectNodes = []
      
      # loop top and bottom face layers
      for k in range(2): 
         if k == 0:
            matn = 2
            z = z2
         elif k == 1:
            matn = 1
            z = z1
            
         # generate nodes
         nds = []
         for i in range(w): # row
            for j in range(w): # col
               nds.append( Node( cx[j], cy[w-1-i], z ))  # y-coordinates are numbered starting from the bottom, but node numbering starts at the top, so cy[w-1-i]
         
         # free inner intersection nodes' 6th DOF
         for i in range(2, w//2+1):
            self.sectNodes.append( nds[i*w - i] )
         for i in range(w//2 + 2, w):
            self.sectNodes.append( nds[i*w - i] )
            self.sectNodes.append( nds[i*w - w + i - 1])
      
         # generate elements
         els = []
         for i in range(w-1): # row
            for j in range(w-1): # col
               n = i*w + j # number of top left node of the element
               els.append( Element( nds[n], nds[n+w], nds[n+w+1], nds[n+1], matn ) )
         
         if k == 0:
            topNodes = nds
            topElems = els
         elif k == 1:
            btmNodes = nds
            btmElems = els
      
      self.nodes = [] 
      self.elems = []
      self.nodes.extend(topNodes)
      self.nodes.extend(btmNodes)
      
      self.elems.extend(topElems)
      self.elems.extend(btmElems)

class HoneycombMesh:
   def __init__(self, lx=5., ly=5., h=2., cell_lx=2., cell_ly=2., hx=2, hz=2, origin=(0.,0.,0.), sectPt=(1./3, -1./3)):
      self.cells = []
      self.blocks = []
      self.edgeLines = []
      self.bounLines = []
      self.elems = []
      self.nodes = []
      
      self.sectNodes = {}
      
      # determine global measurements
      nx = int(math.ceil( lx/cell_lx ))
      ny = int(math.ceil( ly/cell_ly ))
      lx = cell_lx*nx
      ly = cell_ly*ny
      
      z1, z2 = origin[2]-h/2., origin[2]+h/2.
      
      epsz = h*1.e-4
      
      #######
      # create mesh
      # complete cells (blocks) first, then the overlays (coor/elem) - to not confuse FEAP with interfering coor/bloc commands
      #######
      for i in range(nx):
         for j in range(ny):
            x = origin[0] - lx/2. + i*cell_lx
            y = origin[1] - ly/2. + j*cell_ly
            
            # create cells
            hcc = HoneycombCell( x, y, origin[2], cell_lx, cell_ly, h, hx, hz, sectPt )
            self.cells.append(hcc)
            self.blocks.extend(hcc.blocks)
            self.edgeLines.extend(hcc.edgeLines)
            
            for pt in hcc.sectCorners:
               if pt in self.sectNodes: self.sectNodes[pt]+=1
               else: self.sectNodes[pt] = 1
               
      for i in range(nx):
         for j in range(ny):
            x = origin[0] - lx/2. + i*cell_lx
            y = origin[1] - ly/2. + j*cell_ly
            
            # create overlay
#             if hx == 1:
#                ol = OverLay( x, y, origin[2], cell_lx, cell_ly, h, sectPt )
#             elif hx == 2:
#                ol = OverLayInter( x, y, origin[2], cell_lx, cell_ly, h, sectPt )
#             elif hx > 1:
            ol = OverLayHx( x, y, origin[2], cell_lx, cell_ly, h, hx, sectPt )
            for nd in ol.sectNodes:
               #self.edgeLines.append( (nd.x,nd.y,nd.z-epsz,nd.x,nd.y,nd.z+epsz,0,0,0,1000,1,0) )
               #self.edgeLines.append( ((0,0,0,0,0,0) ))
               
               # USE BOUN INSTEAD OF EDGE FOR FASTER INPUT FILE PARSING
               self.bounLines.append( (nd.id, 0, 0, 0, 0, 0, 0, 0) )
            self.elems.extend(ol.elems)
            self.nodes.extend(ol.nodes)
               

      
      #self.blocks.extend(ol.blocks)
      
      # release 6th DOF on intersection nodes vertically
      
      
      for pt in self.sectNodes.keys():
         (x,y) = pt
         #### Internal stiffener/stiffener intersection for all z
         if self.sectNodes[pt] >= 2:
            self.edgeLines.append( (x,y,z1,x,y,z2,0,0,0,1000,1,0) )
            self.edgeLines.append( ((0,0,0,0,0,0) ))
         #### Boundary stiffener/face layer intersection for z = +h/2 and -h/2 only
         else:
            ############
            # despite only needing to free two nodes (x,y,z1), (x,y,z2) here,
            # edge is used instead of poin, as poin doesn't support SETting bcs instead of ADDing them.
            ############
            self.edgeLines.append( (x,y,z1-epsz,x,y,z1+epsz,0,0,0,1000,1,0) )
            self.edgeLines.append( ((0,0,0,0,0,0) ))
            self.edgeLines.append( (x,y,z2-epsz,x,y,z2+epsz,0,0,0,1000,1,0) )
            self.edgeLines.append( ((0,0,0,0,0,0) ))
            
      print "Generated honeycomb mesh with %ix%i cells with lengths lx=%f, ly=%f, h=%f." % (nx, ny, cell_lx, cell_ly, h)
      print "Global mesh has dimensions of lx=%f, ly=%f." % (lx, ly)
      print "Mesh contains %i blocks and %i edge commands." % (len(self.blocks), len(self.edgeLines)//2)

 
 


#################################################### START
 
 
 
 
   
def main():
   global NODE_COUNTER
   global ELEM_COUNTER
   
   
   #########################
   #### OPTIONS
   #########################
   
   filename = os.path.join("C:\\","FG","feap","exe","rve","microShell","Faltwerke","iPlatteWabe")
   
   ### WÄHLE RANDBEDINGUNG UND BELASTUNGSART
   
   #rand = ""
   rand = 'Einspannung'          # allseitig 
   #rand = 'NavierLinksRechts'    # links und rechts nur Mittelfläche
   #rand = 'Navier'               # allseitig nur Mittelfläche        
   #rand = 'EinspannungRechts'    # gesamtes Profil
   #rand = 'Punktlagerung'        # in den Ecken
   #rand = 'RVE'
   
   #last = 'Torsion'
   #last = '4PktBiegung'
   last = 'Flächenlast'
   #last = 'RVEeinzeln'
   #last = 'RVEkopp'
   
   if last=='RVEkopp': filename = os.path.join("C:\\","FG","feap","exe","rve","microShell","platteWabeLin","FE2PlatteWabeLinRVE","irves8_01")
   
   # material parameters
   li = 0      # linear geometry?  0=lin,1=nili (moderate),2=nili (finite)
   e = 7000    # Young's modulus
   v = 0.34    # Poisson ratio
   f = -0.005  # loading force
   
   # material - choose one, comment the other(s)
   
   matn = 1 # linear-elastisch
   mateCard = (e,v,0)
   nmnbnz = '440'
   
   #matn = 4 # elasto-plastisch
   #y0 = 10.
   #xk = 40.
   #nmnbnz = '002'
   #mateCard = (e,v,y0,0,xk,0,0)
   
   
   
   
   
   #########################
   #### END OPTIONS
   #########################
   
   # generate simplified "honeycomb" mesh
   
   #segments = []
   materials = []
   boundaries = []
   elems = []
   #blocks = []
   
   #additionalElements = []
   
   # create mesh
   
   # Inkrementierung
   hz = 1        # Anzahl Elemente Hoehe
   hx = 4        # Anzahl Elemente Breite
   # Abmessungen global (gerade zahlen!)
   a = 4.         # LÃ¤nge Deckschicht
   b = 4.         # Breite Deckschicht
   # Abmessungen Zelle
   lx = 2.
   ly = 2.
   h = 0.5     
   t = 0.1
   # Koordinaten global
   x1 = -a/2.    # links
   x2 = a/2.     # rechts
   y1 = -b/2.    # unten
   y2 = b/2.     # oben
   z1 = -h/2.
   z2 = h/2.
   # Koordinaten Schnittpunkt Zelle Dreibein: werden vorgegeben!
   #xs = lx/6.
   #ys = -ly/6.
   
   
   
   shmMode = 0 # 0 = files, 2 = shared memory history files
   
   nx = int(a/lx)  # Anzahl Zellen x-Richtung
   ny = int(b/ly)  # Anzahl Zellen y-Richtung
   
   # Knoten mit Last
   #f = (x+1)*(y/2.+1)                            # Lastfall 1
   #f = (x+1)*(y+1)+(x+1)*(z/2.+1)                # Lastfall 2
   #f = (x+1)*(y+1)+(x+1)*(z+1)+(x+1)*(y/2.+1)    # Lastfall 3
   
   
   eps = (t/4.)*1.e-4 # small length parameter
   

   ######## HEADER/FOOTER/MATERIALS
   
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
   
   tploNode = 1
   
   if rand=='RVE' and last=='RVEeinzeln':
      footer = Footer(tie=True, macr=('epsq,,1','plot,rot1,-65','plot,rot2,0', 'plot,rot3,20', 'plot,mesh', 'plot,axis', 'plot,boun,,1', 'plot,boun,,2', 'plot,boun,,3',
                                      'plot,boun,,4', 'plot,boun,,5'), link=('6,1,%f,%f,1,1,0,1,1,1' % (-a/2.,a/2.), '6,2,%f,%f,1,1,0,1,1,1' % (-b/2.,b/2.) ), macrBlanks=1)
   elif rand=='RVE' and last=='RVEkopp':
      footer = Footer(tie=True, noInteStop=True, macr=(), link=('6,1,%f,%f,1,1,0,1,1,1' % (-a/2.,a/2.), '6,2,%f,%f,1,1,0,1,1,1' % (-b/2.,b/2.)) )
      footer = str(footer) + """
   
   
batc,tang
nopr
prop,,1
rest,,0,%i
epsq 
loop,,2
tang,,1 
next
sigq
end,,0,%i
1,1,0,10000,1

stop,tang


batc,updh
nopr
prop,,1
rest,,0,%i
updh,,2
end,,0,%i
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

""" % (shmMode, shmMode, shmMode, shmMode)
   else:
      footer = Footer(tie=True, macr=('prop', 'plot,rot1,-65','plot,rot2,0', 'plot,rot3,20', 'plot,mesh', 'plot,axis', 'plot,boun,,1', 'plot,boun,,2', 'plot,boun,,3',
                                   'plot,boun,,4', 'plot,boun,,5', 'plot,node,,-%i' % tploNode, 'tplo,init,%i,-3,1' % tploNode, 'tplo,mark,2,4,1'), macrBlanks=1 )
   

 
   ######## MESH: NODES/ELEMENTS/BLOCKS
   
   edgeLines = []
   bounLines = []
   
   mesh = HoneycombMesh(a,b,h,lx,ly,hx,hz)
   blocks = mesh.blocks
   elems = mesh.elems
   nodes = mesh.nodes
   edgeLines.extend(mesh.edgeLines)
   bounLines.extend(mesh.bounLines)
   
   ######## BOUNDARY INPUTS

   # Rotationsfreiheitsgrad um die 3-Achse aller Knoten sperren   
   boundaries.append(BoundaryInput('vbou', '6. FHG komplett sperren', ((-a-eps,a+eps,-b-eps,b+eps,-h-eps,h+eps,0,0,0,0,0,1),)))  # sperrt dofs im defin. volumen
   # Rotationen um die 3-Achse auf Verschneidungskanten wieder freigeben
   if len(bounLines) > 0:
      boundaries.append(BoundaryInput('boun', 'An Verschneidungen 6. FHG freigeben (Deckschicht/Kern innen fuer hx>1)', bounLines))
   
   boundaries.append(BoundaryInput('edge', 'An  allen Verschneidungen 6. FHG wieder freigeben', edgeLines))
   
   

   if rand=="EinspannungRechts":
      boundaries.append(BoundaryInput('ebou', 'Einspannung rechts', ((1, a/2., 1,1,1,1,1,1),)))
      
   elif rand=="NavierLinksRechts":
      boundaries.append(BoundaryInput('edge', 'Auflager rechts', ( (x2,y1,0,x2,y2,0,0,0,0,1000,0,0), (0,0,1,0,0,0))))
      boundaries.append(BoundaryInput('edge', 'Auflager links', ( (x1,y1,0,x1,y2,0,0,0,0,1000,0,0), (0,0,1,0,0,0)))) 
      boundaries.append(BoundaryInput('poin', 'Starrk�rperbewegungen unterdr�cken - linker Randknoten', ((x1, y1, 0), (0,0,0,0,0), (1,1,0,0,0,0))))
      boundaries.append(BoundaryInput('poin', 'Starrk�rperbewegungen unterdr�cken - rechter Randknoten', ((x2, y1, 0), (0,0,0,0,0), (0,1,0,0,0,0))))
      
   elif rand=="Punktlagerung":
      boundaries.append(BoundaryInput('poin', 'Randknoten links unten', ((x1, y1, 0), (0,0,0,0,0), (1,1,1,0,0,0))))
      boundaries.append(BoundaryInput('poin', 'Randknoten links oben', ((x1, y2, 0), (0,0,0,0,0), (0,0,1,0,0,0))))
      boundaries.append(BoundaryInput('poin', 'Randknoten rechts unten', ((x2, y1, 0), (0,0,0,0,0), (0,0,1,0,0,0))))
      boundaries.append(BoundaryInput('poin', 'Randknoten rechts oben', ((x2, y2, 0), (0,0,0,0,0), (0,0,1,0,0,0))))
     
   elif rand=="Einspannung":
      boundaries.append(BoundaryInput('ebou', 'Einspannung rechts', ((1, a/2., 1,1,1,1,1,1),)))
      boundaries.append(BoundaryInput('ebou', 'Einspannung links', ((1, -a/2., 1,1,1,1,1,1),)))
      boundaries.append(BoundaryInput('ebou', 'Einspannung oben', ((2, b/2., 1,1,1,1,1,1),)))
      boundaries.append(BoundaryInput('ebou', 'Einspannung unten', ((2, -b/2., 1,1,1,1,1,1),)))
   
   elif rand=="Navier":
      boundaries.append(BoundaryInput('edge', 'Auflager rechts', ( (x2,y1,0,x2,y2,0,0,0,0,1000,0,0), (0,0,1,0,0,0))))
      boundaries.append(BoundaryInput('edge', 'Auflager links', ( (x1,y1,0,x1,y2,0,0,0,0,1000,0,0), (0,0,1,0,0,0))))
      boundaries.append(BoundaryInput('edge', 'Auflager oben', ( (x1,y2,0,x2,y2,0,0,0,0,1000,0,0), (0,0,1,0,0,0))))
      boundaries.append(BoundaryInput('edge', 'Auflager unten', ( (x1,y1,0,x2,y1,0,0,0,0,1000,0,0), (0,0,1,0,0,0))))
      boundaries.append(BoundaryInput('poin', 'Starrk�rperbewegungen unterdr�cken - linker Randknoten', ((x1, y1, 0), (0,0,0,0,0), (1,1,0,0,0,0))))
      boundaries.append(BoundaryInput('poin', 'Starrk�rperbewegungen unterdr�cken - rechter Randknoten', ((x2, y1, 0), (0,0,0,0,0), (0,1,0,0,0,0))))
  
   elif rand=="RVE":
      boundaries.append(BoundaryInput('ebou', 'RVE: Inplane-Lagerung Aussenkanten',((1,-a/2.,1,1,0), (1,a/2.,1,1,0), (2,-b/2.,1,1,0), (2,b/2.,1,1,0)) ))
      if nx%2==0 and ny%2==0: # z-Lagerung Mittelknoten
         boundaries.append(BoundaryInput('poin', 'RVE: z-Lagerung Mittelknoten',((0., 0., 0., 1.), (0,0,0,0,0), (0,0,1,0,0,0))))
      elif nx%2==1 and ny%2==1:# z-Lagerung verschobener Mittelknoten
         boundaries.append(BoundaryInput('poin', 'RVE: z-Lagerung Mittelknoten',((1./3.,-1./3., 0., 1.), (0,0,0,0,0), (0,0,1,0,0,0))))
   
   
   header = Header("Honigwabe mit Schalenelementen", (NODE_COUNTER-1,ELEM_COUNTER-1,len(materials),3,6,4), solv=4 if rand!='RVE' else 2, nopr=(last=="RVEkopp"))
 
 
################ BARRIER ####################################
   
   with open(filename,'w') as file:
      file.write(str(header))
      
      for bloc in blocks:
         file.write(str(bloc))
         
      if len(nodes)>0:
         file.write('coor\n')
         for nd in nodes:
            file.write(str(nd))
         file.write('\n')
      
      if len(elems)>0:
         file.write('elem\n')
         for elem in elems:
            file.write(str(elem))
         file.write('\n')
      
      file.write('\n\n')
      file.write('\n\n')
      
      for b in boundaries:
         file.write(str(b)+'\n')
         
      
      # Torsion Kragtr�ger
      if last=="Torsion":
         file.write('eloa\n')
         file.write('%f,%f,%f,%f,%f,%f,%f,%f,%f,%i,%i,%i\n' % ( x1, y1, z2, x1, y2, z2, 0. ,0., 0., 1000, 1, 0))
         file.write('%i,%i,%f,%f,%f,%f,%f\n' % (1,1,0,0,-f,0,0))
         file.write('%f,%f,%f,%f,%f,%f,%f,%f,%f,%i,%i,%i\n' % ( x1, y1, z1, x1, y2, z1, 0. ,0., 0., 1000, 1, 0))
         file.write('%i,%i,%f,%f,%f,%f,%f\n' % (1,1,0,0,f,0,0))
         file.write('\n')
      
      # 4-Punkt-Biegung
      elif last=="4PktBiegung":
         file.write('eloa\n')
         file.write('%f,%f,%f,%f,%f,%f,%f,%f,%f,%i,%i,%i\n' % ( x2/3., y1, z2, x2/3., y2, z2, 0. ,0., 0., 1000, 1, 0))
         file.write('%i,%i,%f,%f,%f,%f,%f\n' % (1,1,0,0,0,0,f))
         file.write('%f,%f,%f,%f,%f,%f,%f,%f,%f,%i,%i,%i\n' % ( x1/3., y1, z2, x1/3., y2, z2, 0. ,0., 0., 1000, 1, 0))
         file.write('%i,%i,%f,%f,%f,%f,%f\n' % (1,1,0,0,0,0,f))
         file.write('\n')
      
      # Plattenbelastung
      elif last=="Flächenlast":
         file.write('aloa\n')
         file.write('%f,%f,%f,%f,%f,%f\n' % ( x1 - eps, x2 + eps, y1 - eps, y2 + eps, z2 -eps, z2 +eps))
         file.write('%i,%i\n' % (3, 4))
         file.write('%i,%f\n' % (1, f))
         file.write('\n')
         
      # RVE-"Last" => aufgebrachte Schalenverzerrung
      elif last=='RVEeinzeln':
         file.write('epsq\n')
         file.write('3,8,2,1,1\n0,0,0,0.01,0,0,0,0\n')
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