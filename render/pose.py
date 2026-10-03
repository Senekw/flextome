"""Kinematic posing of the FlexTome assembly for renders/animation (numpy, works on render.load_parts() output)."""
import numpy as np, math
from scipy.optimize import brentq
import render

PZ, R_CON, R_CRANK, L_CONROD, X_CRANK = 95.0, 60.0, 13.0, 71.2, 60.0
Z_CRANK = PZ - L_CONROD
PHI0 = math.degrees(math.asin(R_CRANK/(2*L_CONROD)))
XP0 = X_CRANK + R_CRANK*math.cos(math.radians(PHI0)); ZP0 = Z_CRANK + R_CRANK*math.sin(math.radians(PHI0))
ARM_GRP = {"arm", "cassette_upper_jaw", "cassette_thumbscrew", "histology_cassette", "paraffin_block", "pivot_clamp_plates_arm"}
CRANK_GRP = {"crank_disc", "crank_pin", "handwheel", "handle", "crankshaft"}

def rot_y(t, center, ang_deg):
    """rotate points t about an axis parallel to y through center; positive = right-hand about +y."""
    a = math.radians(ang_deg); c, s = math.cos(a), math.sin(a)
    v = t - center
    x = v[..., 0]*c + v[..., 2]*s; z = -v[..., 0]*s + v[..., 2]*c
    out = v.copy(); out[..., 0] = x; out[..., 2] = z
    return out + center

def theta_of_phi(phi_deg):
    """arm angle (deg, +up) for crank angle phi (deg, measured from the rest pose)."""
    ph = math.radians(PHI0 + phi_deg)
    pin = np.array([X_CRANK + R_CRANK*math.cos(ph), Z_CRANK + R_CRANK*math.sin(ph)])
    f = lambda th: math.hypot(R_CON*math.cos(th) - pin[0], PZ + R_CON*math.sin(th) - pin[1]) - L_CONROD
    return math.degrees(brentq(f, -0.6, 0.6)), pin

def posed(parts, phi_deg, hide=()):
    th, pin = theta_of_phi(phi_deg)
    out = []
    A = np.array([R_CON*math.cos(math.radians(th)), PZ + R_CON*math.sin(math.radians(th))])
    A0, B0 = np.array([R_CON, PZ]), np.array([XP0, ZP0])
    ang0 = math.degrees(math.atan2(B0[1]-A0[1], B0[0]-A0[0])); ang1 = math.degrees(math.atan2(pin[1]-A[1], pin[0]-A[0]))
    for p in parts:
        if p["name"] in hide: continue
        t = p["tris"]
        if p["name"] in ARM_GRP:
            t = rot_y(t, np.array([0, 0, PZ]), -th)            # +theta lifts the face (z up): rot_y with negative angle
        elif p["name"] in CRANK_GRP:
            t = rot_y(t, np.array([X_CRANK, 0, Z_CRANK]), -phi_deg)
        elif p["name"] == "conrod":
            t = t + np.array([A[0]-A0[0], 0, A[1]-A0[1]])
            t = rot_y(t, np.array([A[0], 0, A[1]]), -(ang1-ang0))
        out.append(dict(name=p["name"], color=p["color"], tris=t))
    return out, th

def check():
    import numpy as np
    for phi in (0, 90, 180, 270):
        th, pin = theta_of_phi(phi)
        A = np.array([R_CON*math.cos(math.radians(th)), PZ + R_CON*math.sin(math.radians(th))])
        print(phi, round(th, 2), round(float(np.hypot(*(A-pin))), 3))
