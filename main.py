import math
from copy import copy
import scipy.constants
from dataclasses import dataclass


# speed of light
c = scipy.constants.c
gravitational_constant = scipy.constants.gravitational_constant

@dataclass
class RectVector3:
    x : float
    y : float
    z : float
@dataclass
class PolarVector4:
    r : float
    theta : float
    phi : float
    t : float

def rect3ToPolar4(p : RectVector3, tau : float):
    r = math.sqrt(math.pow(p.x,2) + math.pow(p.y,2) + math.pow(p.z,2))
    return PolarVector4(
        r,
        math.atan(p.y / p.x),
        math.acos(p.z / r),
        tau
    )

# formal time step and bounds
# note: this is formal time (tau), t is part of the relativistic motion
tau_current = 0
delta = 1440 
tau_final = delta * 1000
## starting conditions (taken from NASA horizons)
# magnitude of velocity projected vector: 9.40096186476512 (RR)
# magnitude of velocity orthogonal vector: 45.3379365047866
# arc length converted into angle: 7.57663025389721637475430034537880544202681327076606062922278207310289 × 10^-7
rectp_start = RectVector3(
    -5.899190110070539E+07,
    -8.872862871229433E+06,
    4.685695586575562E+06
)
rectv_start = RectVector3(
    -2.883332406833908E+00,
    -4.607965057313365E+01,
    -3.501252316591184E+00
)
#p_cur = rect3ToPolar4(rectp_start, 0)
#v_cur = rect3ToPolar4(rectv_start, 1440)
p_cur = PolarVector4(
    5.983918310051619E+07,
    math.pi / 2,
    0,
    0
)
v_cur = PolarVector4(
    9.400961864765113E+00,
    0,
    7.57663025389721637475430034537880544202681327076606062922278207310289E-7,
    1
)
#v_cur.r = -9.400961864765113E+00
#v_cur.phi *= -1
#p_cur.theta = math.pi / 2
#p_cur.phi = math.pi / 2
#v_cur.theta = 1e-3
#v_cur.phi = 1e-3

print(p_cur)
print(v_cur)

entity_mass = 3.302E23 # taken from NASA Horizons
r_s = (2 * gravitational_constant * entity_mass) / (math.pow(c,2)) # schwarzchild radius
print("schwarzchild radius: " + str(r_s)) #schwarzchild radius: 0.0004904236241725966


# helper functions
def f_w(r):
    return 1 - (r_s / r)
def df_w(r):
    return r_s  / math.pow(r,2)

def f_v(r):
    return 1 / f_w(r)
def df_v(r):
    return (-1 * r_s) / math.pow(r - r_s, 2)

def clamp(value, lower, upper):
    return lower if value < lower else upper if value > upper else value

def angular_momentum(p : PolarVector4, v : PolarVector4):
    return p.r * v.r * math.sin(math.fabs(v.theta - p.theta))

print(p_cur.r / f_v(p_cur.r))
# Γ^{theta}_{phiphi} ???
def d2theta_dq2(p : PolarVector4, v : PolarVector4):
    return 0
    return (-1) * (2.0 / p.r) * v.theta * v.r + math.sin(p.theta) * math.cos(p.theta) * math.pow(v.phi,2)
def d2phi_dq2(p : PolarVector4, v : PolarVector4):
    #print(f"v r: {v.r} v theta: {v.theta} v phi: {v.phi}")
    #print((-1) * (2.0 / p.r) * v.phi * v.r)
    return (20) * (2.0 / p.r) * v.phi * v.r # - 2 * clamp(1 / math.tan(p.theta), 0, 0) * v.phi * v.theta
def d2t_dq2(p : PolarVector4, v : PolarVector4):
    return (-1) * ((4 * r_s) / ((p.r**2)*(1 - ((4 * r_s) / p.r)))) * v.t * v.r
    return (-1) *  (1.0 / f_w(p.r)) * df_w(p.r) * v.t * v.r
def d2r_dq2(p : PolarVector4, v : PolarVector4):
    f_v1 = f_v(p.r)
    print(f"p r: {p.r} p theta: {p.theta} p phi: {p.phi}")
    print(f"v r: {v.r} v theta: {v.theta} v phi: {v.phi} v t: {v.t}")
    #print((1 / f_v1) * df_v(p.r) * math.pow(v.r, 2))
    #print((p.r / f_v1) * math.pow(v.theta, 2))
    #print(((r_s * (math.sin(p.theta)**2)) / f_v1) * math.pow(v.phi, 2))
    #print((math.pow(c, 2) / (2 * f_v1)) * df_w(p.r) * math.pow(v.t, 2))
    #print((df_v(p.r) / f_v1) * (v.r**2))
    #print((p.r / f_v1) * (v.phi**2))
    #print(-((1**2) / (2*f_v1)) * (v.t **2))
    print(-((2 * r_s) / (p.r**2)) * (1 - ((4 * r_s) / p.r)) * (v.t ** 2))
    print(((2 * r_s) / ((p.r**2)*(1 - ((4 * r_s) / p.r)))) * (v.r ** 2))
    print(1 * (1 - ((4 * r_s) / p.r)) * (v.phi ** 2))
    return -((2 * r_s) / (p.r**2)) * (1 - ((4 * r_s) / p.r)) * (v.t ** 2) + \
        ((2 * r_s) / ((p.r**2)*(1 - ((4 * r_s) / p.r)))) * (v.r ** 2) + \
        1 * (1 - ((4 * r_s) / p.r)) * (v.phi ** 2)
    return (df_v(p.r) / f_v1) * (v.r**2) + (p.r / f_v1) * (v.phi**2) - \
        ((1**2) / (2*f_v1)) * (v.t **2)
    return (1 / f_v1) * df_v(p.r) * math.pow(v.r, 2) - \
            (math.pow(c, 2) / (2 * f_v1)) * df_w(p.r) * math.pow(v.t, 2)
    #return (1 / f_v1) * df_v(p.r) * math.pow(v.r, 2) + (p.r / f_v1) * math.pow(v.theta, 2) + \
    #    ((p.r * (math.sin(p.theta)**2)) / f_v1) * math.pow(v.phi, 2) - \
    #        (math.pow(c, 2) / (2 * f_v1)) * df_w(p.r) * math.pow(v.t, 2)

# file handling for output
def format_point(p : PolarVector4):
    return f"[{p.r}, {p.theta}, {p.phi}] "
f = open("output_file.py","wt")
f.write("points = [" + format_point(p_cur) + "# " + format_point(v_cur) + "\n")
while(tau_current < tau_final):
    # temp storage of state
    tau_curtemp = tau_current
    p_curtemp = copy(p_cur)
    v_curtemp = copy(v_cur)
    # simultaneous update using euler's method
    p_cur.r = p_curtemp.r + delta * v_curtemp.r
    p_cur.theta = p_curtemp.theta + delta * v_curtemp.theta
    p_cur.phi = p_curtemp.phi + delta * v_curtemp.phi
    p_cur.t = p_curtemp.t + delta * v_curtemp.t

    v_cur.r = v_curtemp.r + delta * d2r_dq2(p_curtemp, v_curtemp)
    v_cur.theta = v_curtemp.theta + delta * d2theta_dq2(p_curtemp, v_curtemp)
    v_cur.phi = v_curtemp.phi + delta * d2phi_dq2(p_curtemp, v_curtemp)
    v_cur.t = v_curtemp.t + delta * d2t_dq2(p_curtemp, v_curtemp)
    
    tau_current += delta
    # print current point
    #print(f"TIME: {tau_current}")
    #print(p_cur)
    #print(v_cur)
    print(angular_momentum(p_cur,v_cur))
    f.write(","  + format_point(p_cur) + "# " + format_point(v_cur) + "\n")

f.write("]")
f.close()