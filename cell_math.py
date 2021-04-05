# adopted from https://github.com/amccaugh/phidl/blob/master/docs/index.rst #

import numpy as np
from matplotlib import pyplot as plt
import math

def arc(radius = 10, angle = 90, num_pts = 720):
    """ Create a circular arc Path
    Parameters
    ----------
    radius : int or float
        Radius of arc
    angle : int or float
        Total angle of arc
    num_pts : int
        Number of points used per 360 degrees
    Returns
    -------
    Path
        A Path object with the specified arc
    """
    num_pts = abs(int(num_pts*angle/360))
    t = np.linspace(-90*np.pi/180, (angle-90)*np.pi/180, num_pts)
    x = radius*np.cos(t)
    y = radius*(np.sin(t)+1)
    points = np.array((x,y)).T * np.sign(angle)

    # P = Path()
    # # Manually add points & adjust start and end angles
    # P.points = points
    # P.start_angle = 0
    # P.end_angle = angle
    return points


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



def euler(radius = 10, angle = 90, p = 1.0, use_eff = True, num_pts = 720):
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





# plt.figure(figsize = (5,5))








# euler1 = euler( radius = 10,angle =  90, p = 0.1, use_eff = True, num_pts = 300)

# euler2 = euler( radius = 10,angle =  90, p = .25, use_eff = True, num_pts = 200)

# euler3 = euler( radius = 10,angle =  90, p = .5, use_eff = True, num_pts = 200)

# euler4 = euler( radius = 10,angle =  90, p = 1, use_eff = True, num_pts = 200)

# arc = arc(radius = 10, num_pts = 400)


# # print(arc)

# x1 = []
# y1 = []

# x2 = []
# y2 = []

# x3 = []
# y3 = []

# x4 = []
# y4 = []


# for i in range(len(euler1)):
#     x1.append(euler1[i][0])
#     y1.append(euler1[i][1])

# for i in range(len(euler2)):
#     x2.append(euler2[i][0])
#     y2.append(euler2[i][1])

# for i in range(len(euler3)):
#     x3.append(euler3[i][0])
#     y3.append(euler3[i][1])

# for i in range(len(euler4)):
#     x4.append(euler4[i][0])
#     y4.append(euler4[i][1])
    
       
# plt.plot(x1,y1,'--b')
# plt.plot(x2,y2,'--k')
# plt.plot(x3,y3,'--g')
# plt.plot(x4,y4,'--r')




# x2 = []
# y2 = []
# for i in range(len(arc)):
#     x2.append(arc[i][0])
#     y2.append(arc[i][1])
    
       
# plt.plot(x2,y2)


# plt.grid()
# plt.xlabel('length [um]')
# plt.ylabel('length [um]')
# plt.title('Maintaining effective radius at 10um \n sweeping P-parameter ')
# plt.legend(['euler p=0.15','euler p=0.25','euler p=0.5','euler p=1.0','arc 90° 10um radius'])








# plt.figure(figsize = (5,5))



# euler1 = euler( radius = 10,angle =  90, p = 0.1, use_eff = False, num_pts = 400)

# euler2 = euler( radius = 10,angle =  90, p = 0.25, use_eff = False, num_pts = 400)

# euler3 = euler( radius = 10,angle =  90, p = 0.5, use_eff = False, num_pts = 400)

# euler4 = euler( radius = 10,angle =  90, p = 1, use_eff = False, num_pts = 400)

# # arc2 = arc(radius = 10, num_pts = 150)

# # print(arc)

# x1 = []
# y1 = []

# x2 = []
# y2 = []

# x3 = []
# y3 = []

# x4 = []
# y4 = []


# for i in range(len(euler1)):
#     x1.append(euler1[i][0])
#     y1.append(euler1[i][1])

# for i in range(len(euler2)):
#     x2.append(euler2[i][0])
#     y2.append(euler2[i][1])

# for i in range(len(euler3)):
#     x3.append(euler3[i][0])
#     y3.append(euler3[i][1])

# for i in range(len(euler4)):
#     x4.append(euler4[i][0])
#     y4.append(euler4[i][1])
    
       
# plt.plot(x1,y1,'--b')
# plt.plot(x2,y2,'--k')
# plt.plot(x3,y3,'--g')
# plt.plot(x4,y4,'--r')




# x2 = []
# y2 = []
# for i in range(len(arc)):
#     x2.append(arc[i][0])
#     y2.append(arc[i][1])
    
       
# plt.plot(x2,y2)


# plt.grid()
# plt.xlabel('length [um]')
# plt.ylabel('length [um]')
# plt.title('Maintaining minimum radius at 10um \n sweeping P-parameter ')
# plt.legend(['euler p=0.15','euler p=0.25','euler p=0.5','euler p=1.0','arc 90° 10um radius'])





# plt.figure(figsize = (5,5))



# euler1 = euler( radius = 10,angle =  90, p = 0.3, use_eff = False, num_pts = 400)

# euler2 = euler( radius = 15,angle =  90, p = 0.3, use_eff = False, num_pts = 400)

# euler3 = euler( radius = 20,angle =  90, p = 0.3, use_eff = False, num_pts = 400)

# euler4 = euler( radius = 25,angle =  90, p = 0.3, use_eff = False, num_pts = 400)

# # arc2 = arc(radius = 10, num_pts = 150)

# # print(arc)

# x1 = []
# y1 = []

# x2 = []
# y2 = []

# x3 = []
# y3 = []

# x4 = []
# y4 = []


# for i in range(len(euler1)):
#     x1.append(euler1[i][0])
#     y1.append(euler1[i][1])

# for i in range(len(euler2)):
#     x2.append(euler2[i][0])
#     y2.append(euler2[i][1])

# for i in range(len(euler3)):
#     x3.append(euler3[i][0])
#     y3.append(euler3[i][1])

# for i in range(len(euler4)):
#     x4.append(euler4[i][0])
#     y4.append(euler4[i][1])
    
       
# plt.plot(x1,y1,'--b')
# plt.plot(x2,y2,'--k')
# plt.plot(x3,y3,'--g')
# plt.plot(x4,y4,'--r')




# x2 = []
# y2 = []
# for i in range(len(arc)):
#     x2.append(arc[i][0])
#     y2.append(arc[i][1])
    
       
# plt.plot(x2,y2)


# plt.grid()
# plt.xlabel('length [um]')
# plt.ylabel('length [um]')
# plt.title('Maintaining P-parameter at 0.3 \n sweeping minimum radius')
# plt.legend(['euler Rmin = 10','euler Rmin = 15','euler Rmin = 20','euler Rmin = 25','arc 90° 10um radius'])



# import numpy as np
# from scipy.integrate import odeint
# def clothoid_ode_rhs(state, s, kappa0, kappa1):
#     x, y, theta = state[0], state[1], state[2]
#     return np.array([np.cos(theta), np.sin(theta), kappa0 + kappa1*s])
# def eval_clothoid(x0,y0,theta0, kappa0, kappa1, s):
#     return odeint(clothoid_ode_rhs, np.array([x0,y0,theta0]), s, (kappa0, kappa1))





# x0,y0,theta0 = 0,0,0
# L = 5
# kappa0, kappa1 = 0, 0.1
# s = np.linspace(0, L, 1000)

# sol = eval_clothoid(x0, y0, theta0, kappa0, kappa1, s)

# xs, ys, thetas = sol[:,0], sol[:,1], sol[:,2] 

# plt.figure(figsize = (5,5))
# plt.plot(xs, ys, lw=3);
# plt.xlabel('x (m)');
# plt.ylabel('y (m)');
# plt.title('An awesome Clothoid');


















# plt.show()
