
# Enter your Python code here

## custom_pcell ##

import pya
import math
import sys
import numpy as np
import math
import numpy as np

def _fresnel(R0, s, num_pts, n_iter=8):
    """ Fresnel integral using a series expansion """
    t = np.linspace(0,s/(np.sqrt(2)*R0), num_pts)
    x = np.zeros(num_pts)
    y = np.zeros(num_pts)

    for n in range(0,n_iter):
      x += (-1)**n * t**(4*n+1)/(np.math.factorial(2*n) * (4*n+1))
      y += (-1)**n * t**(4*n+3)/(np.math.factorial(2*n+1) * (4*n+3))

    return np.array([np.sqrt(2)*R0*x, np.sqrt(2)*R0*y])

def _rotate_points(points, angle = 45, center = (0, 0)):
    """ Rotates points around a centerpoint defined by ``center``.  ``points``
    may be input as either single points [1,2] or array-like[N][2], and will
    return in kind.
    Parameters
    ----------
    points : array-like[N][2]
        Coordinates of the element to be rotated.
    angle : int or float
        Angle to rotate the points.
    center : array-like[2]
        Centerpoint of rotation.
    Returns
    -------
    A new set of points that are rotated around ``center``.
    """
    if angle == 0:
         return points
    angle = angle * math.pi/180
    ca = np.cos(angle)
    sa = np.sin(angle)
    sa = np.array((-sa, sa))
    c0 = np.array(center)
    if np.asarray(points).ndim == 2:
        return (points - c0)*ca + (points - c0)[:,::-1]*sa + c0
    if np.asarray(points).ndim == 1:
        return (points - c0)*ca + (points - c0)[::-1]*sa + c0



def euler(radius = 3, angle = 90, p = 1.0, use_eff = False, num_pts = 720):
    """ Create an Euler bend (also known as "racetrack" or "clothoid" curves)
    that adiabatically transitions from straight to curved.  By default,
    `radius` corresponds to the minimum radius of curvature of the bend.
    However, if `use_eff` is set to True, `radius` corresponds to the effective
    radius of curvature (making the curve a drop-in replacement for an arc). If
    p < 1.0, will create a "partial euler" curve as described in Vogelbacher et.
    al. https://dx.doi.org/10.1364/oe.27.031394
    Parameters
    ----------
    radius : int or float
        Minimum radius of curvature
    angle : int or float
        Total angle of curve
    p : float
        Proportion of curve that is an Euler curve
    use_eff : bool
        If False: `radius` corresponds to minimum radius of curvature of the bend
        If True: The curve will be scaled such that the endpoints match an arc
        with parameters `radius` and `angle`
    num_pts : int
        Number of points used per 360 degrees
    Returns
    -------
    Path
        A Path object with the specified Euler curve
    """
    if (p < 0) or (p > 1):
        raise ValueError('[PHIDL] euler() requires argument `p` be between 0 and 1')
    if p == 0:
        P = arc(radius = radius, angle = angle, num_pts = num_pts)
        P.info['Reff'] = radius
        P.info['Rmin'] = radius
        return P

    if angle < 0:
        mirror = True
        angle = np.abs(angle)
    else:
        mirror = False

    # What is R0
    R0 = 1
    alpha = np.radians(angle)

    # Curvature at (Xp,Yp) 
    Rp = R0 / (np.sqrt(p*alpha))
    print(Rp)
    

    sp = R0 * np.sqrt(p*alpha)
    s0 = 2*sp + Rp*alpha*(1-p)
    num_pts = abs(int(num_pts*angle/360))
    num_pts_euler = int(np.round(sp/(s0/2)*num_pts))
    num_pts_arc = num_pts - num_pts_euler

    xbend1, ybend1 = _fresnel(R0, sp, num_pts_euler)
    xp, yp = xbend1[-1], ybend1[-1]

    dx = xp - Rp*np.sin(p*alpha/2)
    dy = yp - Rp*(1-np.cos(p*alpha/2))

    s = np.linspace(sp, s0/2, num_pts_arc)
    xbend2 = Rp*np.sin((s-sp)/Rp + p*alpha/2) + dx
    ybend2 = Rp*(1 - np.cos((s-sp)/Rp + p*alpha/2)) + dy

    x = np.concatenate([xbend1, xbend2[1:]])
    y = np.concatenate([ybend1, ybend2[1:]])
    points1 = np.array([x,y]).T
    points2 = np.flipud(np.array([x,-y]).T)

    points2 = _rotate_points(points2, angle-180)
    points2 += -points2[0,:] + points1[-1,:]

    points = np.concatenate([points1[:-1],points2])

    # Find y-axis intersection point to compute Reff
    start_angle = 180*(angle<0)
    end_angle = start_angle + angle
    dy = np.tan(np.radians(end_angle-90)) * points[-1][0]
    Reff = points[-1][1] - dy
    Rmin = Rp

    # Fix degenerate condition at angle == 180
    if np.abs(180-angle) < 1e-3:
        Reff = points[-1][1]/2

    # Scale curve to either match Reff or Rmin
    if use_eff == True:
        scale = radius/Reff
    else:
        scale = radius/Rmin
    points *= scale
    
    return points
euler_bend = euler(radius = 3, angle = 90, p = 1.0, use_eff = False, num_pts = 720)

for i in range(len(euler_bend)):
    print(tuple(euler_bend[i]))



"""
This sample PCell implements a library called "MyLib" with a single PCell that
draws a sprial. It demonstrates the basic implementation techniques for a PCell 
and how to use the "guiding shape" feature to implement a handle for the inner
and out radii of the spiral.

NOTE: after changing the code, the macro needs to be rerun to install the new
implementation. The macro is also set to "auto run" to install the PCell 
when KLayout is run.
"""

def _fresnel(R0, s, num_pts, n_iter=8):
    """ Fresnel integral using a series expansion """
    t = np.linspace(0,s/(np.sqrt(2)*R0), num_pts)
    x = np.zeros(num_pts)
    y = np.zeros(num_pts)

    for n in range(0,n_iter):
      x += (-1)**n * t**(4*n+1)/(np.math.factorial(2*n) * (4*n+1))
      y += (-1)**n * t**(4*n+3)/(np.math.factorial(2*n+1) * (4*n+3))

    return np.array([np.sqrt(2)*R0*x, np.sqrt(2)*R0*y])

def _rotate_points(points, angle = 45, center = (0, 0)):
    """ Rotates points around a centerpoint defined by ``center``.  ``points``
    may be input as either single points [1,2] or array-like[N][2], and will
    return in kind.
    Parameters
    ----------
    points : array-like[N][2]
        Coordinates of the element to be rotated.
    angle : int or float
        Angle to rotate the points.
    center : array-like[2]
        Centerpoint of rotation.
    Returns
    -------
    A new set of points that are rotated around ``center``.
    """
    if angle == 0:
         return points
    angle = angle * math.pi/180
    ca = np.cos(angle)
    sa = np.sin(angle)
    sa = np.array((-sa, sa))
    c0 = np.array(center)
    if np.asarray(points).ndim == 2:
        return (points - c0)*ca + (points - c0)[:,::-1]*sa + c0
    if np.asarray(points).ndim == 1:
        return (points - c0)*ca + (points - c0)[::-1]*sa + c0



def euler(radius = 3, angle = 90, p = 1.0, use_eff = False, num_pts = 720):
    """ Create an Euler bend (also known as "racetrack" or "clothoid" curves)
    that adiabatically transitions from straight to curved.  By default,
    `radius` corresponds to the minimum radius of curvature of the bend.
    However, if `use_eff` is set to True, `radius` corresponds to the effective
    radius of curvature (making the curve a drop-in replacement for an arc). If
    p < 1.0, will create a "partial euler" curve as described in Vogelbacher et.
    al. https://dx.doi.org/10.1364/oe.27.031394
    Parameters
    ----------
    radius : int or float
        Minimum radius of curvature
    angle : int or float
        Total angle of curve
    p : float
        Proportion of curve that is an Euler curve
    use_eff : bool
        If False: `radius` corresponds to minimum radius of curvature of the bend
        If True: The curve will be scaled such that the endpoints match an arc
        with parameters `radius` and `angle`
    num_pts : int
        Number of points used per 360 degrees
    Returns
    -------
    Path
        A Path object with the specified Euler curve
    """
    if (p < 0) or (p > 1):
        raise ValueError('[PHIDL] euler() requires argument `p` be between 0 and 1')
    if p == 0:
        P = arc(radius = radius, angle = angle, num_pts = num_pts)
        P.info['Reff'] = radius
        P.info['Rmin'] = radius
        return P

    if angle < 0:
        mirror = True
        angle = np.abs(angle)
    else:
        mirror = False

    # What is R0
    R0 = 1
    alpha = np.radians(angle)

    # Curvature at (Xp,Yp) 
    Rp = R0 / (np.sqrt(p*alpha))
    print(Rp)
    

    sp = R0 * np.sqrt(p*alpha)
    s0 = 2*sp + Rp*alpha*(1-p)
    num_pts = abs(int(num_pts*angle/360))
    num_pts_euler = int(np.round(sp/(s0/2)*num_pts))
    num_pts_arc = num_pts - num_pts_euler

    xbend1, ybend1 = _fresnel(R0, sp, num_pts_euler)
    xp, yp = xbend1[-1], ybend1[-1]

    dx = xp - Rp*np.sin(p*alpha/2)
    dy = yp - Rp*(1-np.cos(p*alpha/2))

    s = np.linspace(sp, s0/2, num_pts_arc)
    xbend2 = Rp*np.sin((s-sp)/Rp + p*alpha/2) + dx
    ybend2 = Rp*(1 - np.cos((s-sp)/Rp + p*alpha/2)) + dy

    x = np.concatenate([xbend1, xbend2[1:]])
    y = np.concatenate([ybend1, ybend2[1:]])
    points1 = np.array([x,y]).T
    points2 = np.flipud(np.array([x,-y]).T)

    points2 = _rotate_points(points2, angle-180)
    points2 += -points2[0,:] + points1[-1,:]

    points = np.concatenate([points1[:-1],points2])

    # Find y-axis intersection point to compute Reff
    start_angle = 180*(angle<0)
    end_angle = start_angle + angle
    dy = np.tan(np.radians(end_angle-90)) * points[-1][0]
    Reff = points[-1][1] - dy
    Rmin = Rp

    # Fix degenerate condition at angle == 180
    if np.abs(180-angle) < 1e-3:
        Reff = points[-1][1]/2

    # Scale curve to either match Reff or Rmin
    if use_eff == True:
        scale = radius/Reff
    else:
        scale = radius/Rmin
    points *= scale
    
    return points
euler_bend = euler(radius = 3, angle = 90, p = 1.0, use_eff = False, num_pts = 720)

for i in range(len(euler_bend)):
    print(tuple(euler_bend[i]))



class Euler_90_bend(pya.PCellDeclarationHelper):
  """
  The PCell declaration for the sprial
  """

  def __init__(self):

    # Important: initialize the super class
    super(Euler_90_bend, self).__init__()

    # declare the parameters
    self.param("l", self.TypeLayer, "Layer", default = pya.LayerInfo(1, 0))
    self.param("n", self.TypeInt, "Number of points", default = 64)
    self.param("radius", self.TypeDouble, "Radius", default = 10)
    self.param("angle", self.TypeDouble, "Angle", default = 90)
    self.param("width", self.TypeDouble, "Width", default = 1)
    self.param("p_parameter", self.TypeDouble, "P parameter", default = 1)

     
    # this hidden parameter is used to determine whether the radii have changed
    # or the "inner/outer" handles have been moved
    # self.param("inner_mem", self.TypeDouble, "Inner_mem", default = 0.0, hidden = True)
    # self.param("outer_mem", self.TypeDouble, "outer_mem", default = 0.0, hidden = True)

  def display_text_impl(self):
    # Provide a descriptive text for the cell
    return "Euler_90_bend(L=" + str(self.l) + ",Outer radius=" + ('%.3f' % self.radius) + ")"
    
    # n must be larger or equal than 4
    if self.n <= 4:
      self.n = 4
      
    # n must be larger or equal than 4
   # if self.p_parameter <= 4:
     # self.np_parameter= 4
  
  def can_create_from_shape_impl(self):
    # Implement the "Create PCell from shape" protocol: we can use any shape which 
    # has a finite bounding box
    return self.shape.is_box() or self.shape.is_polygon() or self.shape.is_path()
  
  #def parameters_from_shape_impl(self):
    # Implement the "Create PCell from shape" protocol: we set inner/outr_r and l 
    #from the shape's bounding box width and layer
    #self.inner_r = self.shape.bbox().width() * self.layout.dbu / 4
    #self.outer_r = self.shape.bbox().width() * self.layout.dbu / 2
    #self.l = self.layout.get_info(self.layer)
  
  def transformation_from_shape_impl(self):
    # Implement the "Create PCell from shape" protocol: we use the center of the shape's
    # bounding box to determine the transformation
    return pya.Trans(self.shape.bbox().center())
  
  def produce_impl(self):


    #euler_points = euler(radius = self.radius, angle = self.angle, p = self.p_parameter, num_pts = self.n)
    euler_bend = euler(self.radius, self.angle, self.p_parameter , self.n)
    print(euler_bend)
    pts = []
    for i in range(len(euler_bend)):

        pts.append(pya.Point.from_dpoint(pya.DPoint(euler_bend[i][0]/ self.layout.dbu,euler_bend[i][1]/ self.layout.dbu)))
#        pts.append(pya.Point.from_dpoint(pya.DPoint(euler_points[i][0]/self.layout.dbu , euler_points[i][1]/self.layout.dbu )))


    self.cell.shapes(self.l_layer).insert(pya.Path(pts,self.width/self.layout.dbu))




class Spiral(pya.PCellDeclarationHelper):
  """
  The PCell declaration for the sprial
  """

  def __init__(self):

    # Important: initialize the super class
    super(Spiral, self).__init__()

    # declare the parameters
    self.param("l", self.TypeLayer, "Layer", default = pya.LayerInfo(1, 0))
    self.param("n", self.TypeInt, "Number of points", default = 64)
    self.param("inner_r", self.TypeDouble, "Inner Radius", default = 1)
    self.param("outer_r", self.TypeDouble, "Outer Radius", default = 10)
    self.param("width", self.TypeDouble, "Width", default = 1)
    self.param("spacing", self.TypeDouble, "Spacing", default = 1)
    self.param("inner_handle", self.TypeShape, "", default = pya.DPoint(0, 0))
    self.param("outer_handle", self.TypeShape, "", default = pya.DPoint(0, 0))
     
    # this hidden parameter is used to determine whether the radii have changed
    # or the "inner/outer" handles have been moved
    self.param("inner_mem", self.TypeDouble, "Inner_mem", default = 0.0, hidden = True)
    self.param("outer_mem", self.TypeDouble, "outer_mem", default = 0.0, hidden = True)

  def display_text_impl(self):
    # Provide a descriptive text for the cell
    return "Spiral(L=" + str(self.l) + ",Outer radius=" + ('%.3f' % self.outer_r) + ")"
  
  def coerce_parameters_impl(self):
  
    # We employ coerce_parameters_impl to decide whether the handle or the 
    # numeric parameter has changed (by comparing against the effective 
    # radii *_handle_radius) and set *_mem to the effective radius. We also update the 
    # numerical value or the shape, depending on which on has not changed.
    inner_handle_radius = None
    outer_handle_radius = None
    if isinstance(self.inner_handle, pya.DPoint): 
      # compute distance in micron
      inner_handle_radius = self.inner_handle.distance(pya.DPoint(0, 0))
    if isinstance(self.outer_handle, pya.DPoint): 
      # compute distance in micron
      outer_handle_radius = self.outer_handle.distance(pya.DPoint(0, 0))
    if abs(self.inner_r-self.inner_mem) < 1e-6 and abs(self.outer_r-self.outer_mem) < 1e-6:
      self.inner_mem = inner_handle_radius
      self.inner_r = inner_handle_radius
      self.outer_mem = outer_handle_radius
      self.outer_r = outer_handle_radius
    else:
      self.inner_mem = self.inner_r
      self.inner_handle = pya.DPoint(-self.inner_r, 0)
      self.outer_mem = self.outer_r
      self.outer_handle = pya.DPoint(-self.outer_r, 0)
    
    # n must be larger or equal than 4
    if self.n <= 4:
      self.n = 4
  
  def can_create_from_shape_impl(self):
    # Implement the "Create PCell from shape" protocol: we can use any shape which 
    # has a finite bounding box
    return self.shape.is_box() or self.shape.is_polygon() or self.shape.is_path()
  
  def parameters_from_shape_impl(self):
    # Implement the "Create PCell from shape" protocol: we set inner/outr_r and l 
    #from the shape's bounding box width and layer
    self.inner_r = self.shape.bbox().width() * self.layout.dbu / 4
    self.outer_r = self.shape.bbox().width() * self.layout.dbu / 2
    self.l = self.layout.get_info(self.layer)
  
  def transformation_from_shape_impl(self):
    # Implement the "Create PCell from shape" protocol: we use the center of the shape's
    # bounding box to determine the transformation
    return pya.Trans(self.shape.bbox().center())
  
  def produce_impl(self):
  
    # This is the main part of the implementation: create the layout

    # fetch the parameters convert to database units
    inner_r_dbu = self.inner_r / self.layout.dbu
    outer_r_dbu = self.outer_r / self.layout.dbu
    
    # compute the spiral
    pts = []
    da = math.pi * 2 / self.n
    dr = (self.width+self.spacing)/self.n/self.layout.dbu
    current_radius = inner_r_dbu
    current_angle = 0
    while current_radius < outer_r_dbu:
      pts.append(pya.Point.from_dpoint(pya.DPoint(current_radius * math.cos(current_angle), current_radius * math.sin(current_angle))))
      current_radius += dr
      current_angle = (current_angle+da)%(math.pi*2)
    
    # create the shape
    pt1 = pya.Point.from_dpoint(pya.DPoint(0/ self.layout.dbu,0/ self.layout.dbu))
    pt2 = pya.Point.from_dpoint(pya.DPoint(1/ self.layout.dbu,1/ self.layout.dbu))
    pt3 = pya.Point.from_dpoint(pya.DPoint(2/ self.layout.dbu,2/ self.layout.dbu))
    pt4 = pya.Point.from_dpoint(pya.DPoint(3/ self.layout.dbu,3/ self.layout.dbu))
    pts = [pt1,pt2,pt3,pt4]
    self.cell.shapes(self.l_layer).insert(pya.Path(pts,self.width/self.layout.dbu))




class Oval_shape(pya.PCellDeclarationHelper):
  """
  The PCell declaration for the circle
  """

  def __init__(self):

    # Important: initialize the super class
    super(Oval_shape, self).__init__()

    # declare the parameters
    self.param("l", self.TypeLayer, "Layer", default = pya.LayerInfo(1, 0))
    self.param("s", self.TypeShape, "", default = pya.DPoint(0, 0))
    self.param("r", self.TypeDouble, "Radius", default = 1)
    self.param("n", self.TypeInt, "Number of points", default = 64)     
    self.param("horizontal_scale", self.TypeDouble, "Vorizontal_scale", default = 1)     
    self.param("vertical_scale", self.TypeDouble, "Vorizontal_scale", default = 1)     
    # this hidden parameter is used to determine whether the radius has changed
    # or the "s" handle has been moved
    self.param("ru", self.TypeDouble, "Radius", default = 0.0, hidden = True)
    self.param("rd", self.TypeDouble, "Double radius", readonly = True)

  def display_text_impl(self):
    # Provide a descriptive text for the cell
    return "Oval_shape(L=" + str(self.l) + ",R=" + ('%.3f' % self.r) + ")"
  
  def coerce_parameters_impl(self):
  
    # We employ coerce_parameters_impl to decide whether the handle or the 
    # numeric parameter has changed (by comparing against the effective 
    # radius ru) and set ru to the effective radius. We also update the 
    # numerical value or the shape, depending on which on has not changed.
    rs = None
    if isinstance(self.s, pya.DPoint): 
      # compute distance in micron
      rs = self.s.distance(pya.DPoint(0, 0))
    if rs != None and abs(self.r-self.ru) < 1e-6:
      self.ru = rs
      self.r = rs 
    else:
      self.ru = self.r
      self.s = pya.DPoint(-self.r, 0)
    
    self.rd = 2*self.r
    
    # n must be larger or equal than 4
    if self.n <= 4:
      self.n = 4
  
  def can_create_from_shape_impl(self):
    # Implement the "Create PCell from shape" protocol: we can use any shape which 
    # has a finite bounding box
    return self.shape.is_box() or self.shape.is_polygon() or self.shape.is_path()
  
  def parameters_from_shape_impl(self):
    # Implement the "Create PCell from shape" protocol: we set r and l from the shape's 
    # bounding box width and layer
    self.r = self.shape.bbox().width() * self.layout.dbu / 2
    self.l = self.layout.get_info(self.layer)
  
  def transformation_from_shape_impl(self):
    # Implement the "Create PCell from shape" protocol: we use the center of the shape's
    # bounding box to determine the transformation
    return pya.Trans(self.shape.bbox().center())
  
  def produce_impl(self):
  
    # This is the main part of the implementation: create the layout

    # fetch the parameters
    ru_dbu = self.ru / self.layout.dbu
    
    # compute the circle
    pts = []
    da = math.pi * 2 / self.n
    for i in range(0, self.n):
      pts.append(pya.Point.from_dpoint(pya.DPoint(ru_dbu * self.horizontal_scale * math.cos(i * da), ru_dbu * self.vertical_scale* math.sin(i * da))))
    
    # create the shape
    self.cell.shapes(self.l_layer).insert(pya.Polygon(pts))


class Circle(pya.PCellDeclarationHelper):
  """
  The PCell declaration for the circle
  """

  def __init__(self):

    # Important: initialize the super class
    super(Circle, self).__init__()

    # declare the parameters
    self.param("l", self.TypeLayer, "Layer", default = pya.LayerInfo(1, 0))
    self.param("s", self.TypeShape, "", default = pya.DPoint(0, 0))
    self.param("r", self.TypeDouble, "Radius", default = 0.1)
    self.param("n", self.TypeInt, "Number of points", default = 64)     
    # this hidden parameter is used to determine whether the radius has changed
    # or the "s" handle has been moved
    self.param("ru", self.TypeDouble, "Radius", default = 0.0, hidden = True)
    self.param("rd", self.TypeDouble, "Double radius", readonly = True)

  def display_text_impl(self):
    # Provide a descriptive text for the cell
    return "Circle(L=" + str(self.l) + ",R=" + ('%.3f' % self.r) + ")"
  
  def coerce_parameters_impl(self):
  
    # We employ coerce_parameters_impl to decide whether the handle or the 
    # numeric parameter has changed (by comparing against the effective 
    # radius ru) and set ru to the effective radius. We also update the 
    # numerical value or the shape, depending on which on has not changed.
    rs = None
    if isinstance(self.s, pya.DPoint): 
      # compute distance in micron
      rs = self.s.distance(pya.DPoint(0, 0))
    if rs != None and abs(self.r-self.ru) < 1e-6:
      self.ru = rs
      self.r = rs 
    else:
      self.ru = self.r
      self.s = pya.DPoint(-self.r, 0)
    
    self.rd = 2*self.r
    
    # n must be larger or equal than 4
    if self.n <= 4:
      self.n = 4
  
  def can_create_from_shape_impl(self):
    # Implement the "Create PCell from shape" protocol: we can use any shape which 
    # has a finite bounding box
    return self.shape.is_box() or self.shape.is_polygon() or self.shape.is_path()
  
  def parameters_from_shape_impl(self):
    # Implement the "Create PCell from shape" protocol: we set r and l from the shape's 
    # bounding box width and layer
    self.r = self.shape.bbox().width() * self.layout.dbu / 2
    self.l = self.layout.get_info(self.layer)
  
  def transformation_from_shape_impl(self):
    # Implement the "Create PCell from shape" protocol: we use the center of the shape's
    # bounding box to determine the transformation
    return pya.Trans(self.shape.bbox().center())
  
  def produce_impl(self):
  
    # This is the main part of the implementation: create the layout

    # fetch the parameters
    ru_dbu = self.ru / self.layout.dbu
    
    # compute the circle
    pts = []
    da = math.pi * 2 / self.n
    for i in range(0, self.n):
      pts.append(pya.Point.from_dpoint(pya.DPoint(ru_dbu * math.cos(i * da), ru_dbu * math.sin(i * da))))
      

    # create the shape

    self.cell.shapes(self.l_layer).insert(pya.Polygon(pts))


class MyLib(pya.Library):
  """
  The library where we will put the PCell into 
  """

  def __init__(self):
  
    # Set the description
    self.description = "My First Library"
    
    # Create the PCell declarations
    self.layout().register_pcell("Circle", Circle())
    # That would be the place to put in more PCells ...
    self.layout().register_pcell("Oval_shape", Oval_shape())
    self.layout().register_pcell("Spiral", Spiral())
    self.layout().register_pcell("Euler_90_bend", Euler_90_bend())
    # Register us with the name "MyLib".
    # If a library with that name already existed, it will be replaced then.
    self.register("MyLib")

    



# Instantiate and register the library
MyLib()
