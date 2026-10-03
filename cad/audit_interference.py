"""Boolean interference audit of the FlexTome assembly (CadQuery).
Poses the arm group / crank group / connecting rod at 12 crank angles over one revolution (loop closure),
and the screw group at -5, 0, +5 mm of feed-screw excursion, and intersects every part pair.
Reports pair, maximum overlap volume, and the pose at which it occurs.  Intended overlaps (pins in bores, strips in
clamp pockets, nuts in pockets, pawl tooth in the drum, striker on the pawl flag) are listed in INTENDED and reported
separately from unintended ones.  Run:  python audit_interference.py   (from a directory containing flextome_cad.py)."""
import itertools, json, math, sys
import numpy as np
from scipy.optimize import brentq
import flextome_cad as fc

ARM_GRP = {"arm", "cassette_upper_jaw", "cassette_thumbscrew", "histology_cassette", "paraffin_block", "pivot_clamp_plates_arm"}
CRANK_GRP = {"crank_disc", "crank_pin"}
SCREW_GRP = {"feed_screw_M8x1_25", "feed_screw_M8x1_0", "screw_coupler", "ratchet_wheel_50T", "feed_knob", "coupler_nut_M8x1_25"}
SAME = [ARM_GRP, CRANK_GRP, SCREW_GRP]
INTENDED = {  # pairs whose overlap is a modelled fastening/engagement, with the volume cap (mm3) we accept
    ("crank_pin", "conrod"): 120, ("conrod", "crank_disc"): 5, ("pivot_strips_steel", "arm"): 80, ("pivot_strips_steel", "frame"): 80,
    ("pivot_strips_steel", "pivot_clamp_plates_arm"): 80, ("pivot_strips_steel", "pivot_clamp_plates_frame"): 80,
    ("pivot_clamp_plates_arm", "arm"): 100, ("pivot_clamp_plates_frame", "frame"): 100,
    ("ratchet_wheel_50T", "pawl_rod"): 25, ("arm", "pawl_rod"): 200,            # pawl tooth engagement; striker bar lifting the flag (rod drawn at rest)
    ("handwheel", "handle"): 200, ("blade_holder", "pressure_strip_steel"): 60, ("blade_holder", "blade_low_profile"): 5,
    ("blade_holder", "blade_clamp_screws_M4"): 250, ("pressure_strip_steel", "blade_clamp_screws_M4"): 30,
    ("feed_nuts_steel", "stage_body_flexure"): 5, ("feed_nuts_steel", "screw_coupler"): 5,
    ("stage_body_flexure", "preload_springs"): 60, ("fastener_heads_M3_M4", "thickness_selector_shield"): 5,
    ("blade_holder", "fastener_heads_M3_M4"): 5, ("stage_body_flexure", "stage_clamp_knobs_M6"): 5,
}

def theta_of(phi_deg):
    ph = math.radians(fc.PHI0 + phi_deg)
    pin = np.array([fc.X_CRANK + fc.R_CRANK*math.cos(ph), fc.Z_CRANK + fc.R_CRANK*math.sin(ph)])
    f = lambda th: np.hypot(fc.R_CON*math.cos(th) - pin[0], fc.PZ + fc.R_CON*math.sin(th) - pin[1]) - fc.L_CONROD
    return math.degrees(brentq(f, -0.6, 0.6)), pin

def vol(s):
    try: return s.val().Volume()
    except Exception: return 0.0

EXPECTED_SHELLS = {"crankshaft": 2, "crank_disc": 2, "tie_rods_M8_steel_spine": 6, "fastener_heads_M3_M4": 9}   # multi-piece parts modelled as one STL

def connectivity(meta):
    """Each exported STL must consist of exactly qty connected shells (or the documented count for grouped visual parts)."""
    import trimesh
    rows = {}
    for m in meta:
        n = len(trimesh.load(f"stl/{m['name']}.stl", force="mesh").split(only_watertight=False))
        exp = EXPECTED_SHELLS.get(m["name"], m["qty"])
        rows[m["name"]] = dict(shells=n, expected=exp, ok=(n == exp))
    return rows

def main():
    fc.PARTS.clear(); meta = fc.build()
    conn = connectivity(meta)
    print("connectivity failures:", [k for k, r in conn.items() if not r["ok"]])
    P = {p["name"]: p["solid"] for p in fc.PARTS}; names = list(P)
    A0, B0 = np.array([fc.R_CON, fc.PZ]), np.array([fc.XP0, fc.ZP0])
    ang0 = math.degrees(math.atan2(B0[1]-A0[1], B0[0]-A0[0]))
    report = {}
    for phi in np.arange(0, 360, 30):
        th, pin = theta_of(phi)
        A = np.array([fc.R_CON*math.cos(math.radians(th)), fc.PZ + fc.R_CON*math.sin(math.radians(th))])
        ang1 = math.degrees(math.atan2(pin[1]-A[1], pin[0]-A[0]))
        for dx in (-5.0, 0.0, 5.0):
            moved = {}
            for n in names:
                s = P[n]
                if n in ARM_GRP: s = s.rotate((0, 0, fc.PZ), (0, 1, fc.PZ), -th)
                elif n in CRANK_GRP: s = s.rotate((fc.X_CRANK, 0, fc.Z_CRANK), (fc.X_CRANK, 1, fc.Z_CRANK), -phi)
                elif n in SCREW_GRP: s = s.translate((dx, 0, 0))
                elif n == "conrod":
                    s = s.translate((A[0]-A0[0], 0, A[1]-A0[1])).rotate((A[0], 0, A[1]), (A[0], 1, A[1]), -(ang1-ang0))
                moved[n] = s
            for a, b in itertools.combinations(names, 2):
                if any({a, b} <= g for g in SAME): continue
                try: v = vol(moved[a].intersect(moved[b]))
                except Exception: v = -1.0
                if v > 0.5 or v < 0:
                    key = (a, b)
                    if key not in report or v > report[key]["max_mm3"]:
                        report[key] = dict(max_mm3=round(v, 1), phi=float(phi), theta=round(th, 1), screw_dx=dx)
    unintended, intended = {}, {}
    for (a, b), r in report.items():
        cap = INTENDED.get((a, b), INTENDED.get((b, a)))
        (intended if cap is not None and r["max_mm3"] <= cap else unintended)[f"{a} | {b}"] = r
    out = dict(n_parts=len(names), poses="12 crank angles x 3 screw positions", unintended=unintended, intended=intended, connectivity=conn)
    json.dump(out, open("audit_results.json", "w"), indent=1)
    print("parts:", len(names), " unintended overlaps:", len(unintended))
    for k, r in sorted(unintended.items(), key=lambda kv: -kv[1]["max_mm3"]): print("  UNINTENDED", k, r)
    for k, r in sorted(intended.items(), key=lambda kv: -kv[1]["max_mm3"]): print("  intended  ", k, r)
    return out

if __name__ == "__main__":
    main()
