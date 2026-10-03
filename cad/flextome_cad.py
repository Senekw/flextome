"""
FlexTome parametric CAD (CadQuery).  Units: mm.  Global frame: x = feed (specimen->blade), y = lateral (blade edge),
z = vertical.  Base-plate top (rear) at z=0.  Arm pivot P = (0, 0, 95).
Exports: step/flextome_assembly.step, stl/<part>.stl, parts.json (name, color, material, qty, printed?).
"""
import cadquery as cq, json, math, os
os.makedirs("stl", exist_ok=True); os.makedirs("step", exist_ok=True)

# ------------------------------------------------------------------ key parameters
PZ = 95.0            # pivot height
L_ARM = 125.0        # pivot -> block face
R_CON = 60.0         # pivot -> conrod joint on arm
R_CRANK = 13.0
L_CONROD = 71.2              # pin-to-pin conrod length
Z_CRANK = PZ - L_CONROD      # 23.8: crank centre exactly L_CONROD below the arm joint -> symmetric +/-12.5 deg swing
X_CRANK = R_CON              # 60
THETA_MAX = math.degrees(math.asin(R_CRANK / R_CON))
PHI0 = math.degrees(math.asin(R_CRANK / (2*L_CONROD)))   # crank angle (from +x) at which the arm is exactly horizontal (rest pose), 5.24 deg
XP0 = X_CRANK + R_CRANK*math.cos(math.radians(PHI0))    # crank-pin rest position
ZP0 = Z_CRANK + R_CRANK*math.sin(math.radians(PHI0))
BLADE_LEAN = 27.0            # blade back-lean from vertical [deg] (clearance ~ 23-17.5 = 5.5 deg + arc effect)
STRIP_T, STRIP_FREE, STRIP_CLAMP = 0.2, 20.0, 6.0
SX_DECK, Z_DECK = 110.0, 22.0        # front deck starts at x=110, top at z=22 (lowered 8 mm in rev. 2 to make room for the blade holder)
X_BASE0, X_BASE1 = -40.0, 195.0
Y_BASE = 75.0
Z_SCREW = 60.0               # feed-screw axis (rev. 2: lowered with the deck so the coupler clears the platform)
DXD = 3.0                    # rev. 2: drum/shield/knob shifted +3 mm so the drum clears the 10 mm upstand at -5 mm screw travel
N_TEETH = 50
R_RATCHET = 25.0


def safe_fillet(solid, selector, r):
    try:
        return solid.edges(selector).fillet(r)
    except Exception:
        return solid

def safe_chamfer(solid, selector, d):
    try:
        return solid.edges(selector).chamfer(d)
    except Exception:
        return solid

def bolt_head_z(x, y, z0, d_head=8.5, h=5.0, hex_=True):
    """socket/hex head sitting on a face at z0 (visual fastener)."""
    r = d_head/math.sqrt(3) if hex_ else d_head/2
    wp = cq.Workplane("XY", origin=(0, 0, z0)).center(x, y)
    return (wp.polygon(6, 2*r).extrude(h) if hex_ else wp.circle(r).extrude(h))

PARTS = []   # (name, solid, color(rgb), material, qty, printed)

def add(name, solid, color, material="PLA (FDM)", qty=1, printed=True):
    PARTS.append(dict(name=name, solid=solid, color=color, material=material, qty=qty, printed=printed))
    return solid

def box(xr, yr, zr):
    """axis-aligned box from ranges."""
    return (cq.Workplane("XY").box(xr[1]-xr[0], yr[1]-yr[0], zr[1]-zr[0], centered=False)
            .translate((xr[0], yr[0], zr[0])))

def cyl_y(x, z, y0, y1, d):
    return cq.Workplane("XZ", origin=(0, y0, 0)).center(x, z).circle(d/2).extrude(-(y1-y0))

def cyl_x(y, z, x0, x1, d):
    return cq.Workplane("YZ", origin=(x0, 0, 0)).center(y, z).circle(d/2).extrude(x1-x0)

def cyl_z(x, y, z0, z1, d):
    return cq.Workplane("XY", origin=(0, 0, z0)).center(x, y).circle(d/2).extrude(z1-z0)

def hex_nut(x0, x1, af, hole, axis="x", y=0, z=0):
    r = af/math.sqrt(3)
    wp = cq.Workplane("YZ", origin=(x0, 0, 0)).center(y, z).polygon(6, 2*r).extrude(x1-x0)
    return wp.cut(cyl_x(y, z, x0-1, x1+1, hole))

# ================================================================== FRAME (monolithic base + tower + deck)
def frame():
    f = safe_fillet(box((X_BASE0, X_BASE1), (-Y_BASE, Y_BASE), (-55, 0)), "|Z", 8.0)
    f = f.union(safe_fillet(box((SX_DECK, X_BASE1), (-Y_BASE, Y_BASE), (0, Z_DECK)), "|Z", 8.0))      # raised deck
    # rib pockets from below (ribbed casting style): 4 mm walls, 4 mm ribs
    pockets = None
    ys = [-Y_BASE+4, -28, 28, Y_BASE-4]
    for xs, ztop in (([X_BASE0+4, 20, 80, SX_DECK-2], -6), ([SX_DECK+2, 150, X_BASE1-4], Z_DECK-6)):
        for i in range(len(xs)-1):
            for j in range(len(ys)-1):
                x0, x1 = xs[i]+ (0 if i==0 else 2), xs[i+1] - (0 if i==len(xs)-2 else 2)
                y0, y1 = ys[j]+ (0 if j==0 else 2), ys[j+1] - (0 if j==len(ys)-2 else 2)
                p = box((x0, x1), (y0, y1), (-56, ztop))
                pockets = p if pockets is None else pockets.union(p)
    f = f.cut(pockets)
    # tower with forward gusset: 60 mm root depth tapering to 32 mm at the pivot head (FE-driven revision)
    tower = (cq.Workplane("XZ", origin=(0, 30, 0))
             .polyline([(X_BASE0, 0), (20, 0), (20, 55), (-8, 80), (-8, 118), (X_BASE0, 118)]).close().extrude(60))
    tower = safe_fillet(tower, "|Z", 6.0)
    tower = safe_chamfer(tower, ">Z", 3.0)
    f = f.union(tower)
    # pivot pockets for strip clamps (tower side): strips A (outer, u_A=(c,0,s)) and B (inner, u_B=(c,0,-s))
    f = cut_strip_pockets(f, side=-1)
    # crank pillow blocks
    for yc in (-50, 50):
        blk = box((X_CRANK-15, X_CRANK+15), (yc-7, yc+7), (0, Z_CRANK+13))
        blk = safe_fillet(blk, "|Y and >Z", 8.0)
        blk = blk.cut(cyl_y(X_CRANK, Z_CRANK, yc-8, yc+8, 22.0))     # 608 bearing seat
        f = f.union(blk)
    # deck: heat-set insert holes for stage clamp knobs (M6) and coarse-approach upstand
    for yc in (-35, 35):
        f = f.cut(cyl_z(161, yc, Z_DECK-12, Z_DECK+1, 8.0))
    f = f.union(box((X_BASE1-9, X_BASE1), (-62, -47), (Z_DECK, Z_DECK+12)).cut(cyl_x(-54, Z_DECK+4, X_BASE1-10, X_BASE1+1, 6.2)))
    # steel-spine tie rods: 2 x M8 along x at y=+-55, z=-45 (bores through the base end walls / ribs)
    for yc in (-55, 55):
        f = f.cut(cyl_x(yc, -45, X_BASE0-1, X_BASE1+1, 8.6))
    # handwheel shaft clearance cut-outs at base side for pawl? none.  Feet
    return f

def strip_frames():
    """Return list of (u, w, ycenter, kind) for strips: A outer pair width 10, B inner width 20."""
    c = math.cos(math.radians(45)); s = math.sin(math.radians(45))
    return [((c, 0, s), 10.0, -15.5, "A"), ((c, 0, s), 10.0, 15.5, "A"), ((c, 0, -s), 20.0, 0.0, "B")]

def strip_solid(u, w, yc, s0, s1, thick):
    """plate along unit vector u (in xz plane) from arc-length s0..s1 (measured from P), width w about yc, thickness thick (normal)."""
    L = s1 - s0
    ang = math.degrees(math.atan2(u[2], u[0]))
    sol = cq.Workplane("XY").box(L, w, thick).translate(((s0+s1)/2, yc, 0))
    sol = sol.rotate((0, 0, 0), (0, 1, 0), -ang)       # rotate about y so +x maps to u
    return sol.translate((0, 0, PZ))

def cut_strip_pockets(body, side):
    """side=-1: tower side (s in [-16,-10]); +1: arm side (s in [10,16]).  Pocket = strip + 3 mm clamp plate above.
    A clamp-pad boss is first added around each clamped strip end so the pad protrudes from the body face."""
    for (u, w, yc, kind) in strip_frames():
        n = (-u[2], 0, u[0])
        if side < 0: b0, b1 = -(STRIP_FREE/2+STRIP_CLAMP)-2.5, -(STRIP_FREE/2)+1.0
        else:        b0, b1 = STRIP_FREE/2-1.0, STRIP_FREE/2+STRIP_CLAMP+2.5
        boss = strip_solid(u, w+14, yc, b0, b1, 7.0).translate((n[0]*(-3.3), 0, n[2]*(-3.3)))   # boss ends at the strip seat; pad sits on top
        body = body.union(boss)
    for (u, w, yc, kind) in strip_frames():
        if side < 0: s0, s1 = -(STRIP_FREE/2+STRIP_CLAMP), -(STRIP_FREE/2)
        else:        s0, s1 = STRIP_FREE/2, STRIP_FREE/2+STRIP_CLAMP
        pocket = strip_solid(u, w+14.4, yc, s0-0.5, s1+0.5, 3.6)      # pocket takes the strip and the full-width clamp pad
        # offset pocket so strip mid-plane is at pocket bottom: shift along normal by +1.7
        n = (-u[2], 0, u[0])
        pocket = pocket.translate((n[0]*1.7, 0, n[2]*1.7))
        body = body.cut(pocket)
        # M3 heat-set insert holes (dia 4.0, 6 deep) beside the strip, through the clamp pad region
        sm = (s0+s1)/2
        for dy in (-(w/2+4.0), (w/2+4.0)):
            hole = strip_solid(u, 4.0, yc+dy, sm-2.0, sm+2.0, 14.0).translate((n[0]*(-2.0), 0, n[2]*(-2.0)))
            body = body.cut(hole)
    return body

def pivot_strips_and_clamps():
    """Returns (strips, frame-side clamp pads, arm-side clamp pads). Negative strip coordinate = frame side."""
    strips = None; clamps = [None, None]
    for (u, w, yc, kind) in strip_frames():
        st = strip_solid(u, w, yc, -(STRIP_FREE/2+STRIP_CLAMP), STRIP_FREE/2+STRIP_CLAMP, STRIP_T)
        strips = st if strips is None else strips.union(st)
        n = (-u[2], 0, u[0])
        for k, (s0, s1) in enumerate([(-(STRIP_FREE/2+STRIP_CLAMP)-0.5, -(STRIP_FREE/2)), (STRIP_FREE/2, STRIP_FREE/2+STRIP_CLAMP+0.5)]):
            cp = strip_solid(u, w+14, yc, s0, s1, 3.0).translate((n[0]*1.7, 0, n[2]*1.7))
            sm = (s0+s1)/2
            for dy in (-(w/2+4.0), (w/2+4.0)):      # M3 clearance holes
                cp = cp.cut(strip_solid(u, 3.4, yc+dy, sm-1.7, sm+1.7, 8.0).translate((n[0]*1.7, 0, n[2]*1.7)))
            clamps[k] = cp if clamps[k] is None else clamps[k].union(cp)
    return strips, clamps[0], clamps[1]

# ================================================================== ARM
def arm():
    root = safe_fillet(box((8, 32), (-20, 20), (PZ-17, PZ+20)), "|X", 4.0)
    root = cut_strip_pockets(root, side=+1)
    # tapered box beam 32 -> 110 via loft
    beam = (cq.Workplane("YZ", origin=(32, 0, PZ)).rect(40, 40)
            .workplane(offset=78).rect(30, 30).loft(combine=True))
    # lightening: hollow (3 mm wall) -- approximate by inner loft cut
    inner = (cq.Workplane("YZ", origin=(34, 0, PZ)).rect(34, 34)
             .workplane(offset=72).rect(24, 24).loft(combine=True))
    beam = beam.cut(inner)
    # vent/access window into the cavity on the underside near the root (makes the hollow an open, printable pocket)
    beam = beam.cut(box((40, 70), (-10, 10), (PZ-22, PZ-14)))
    # conrod pin boss & hole (pin along y at x=R_CON, z=PZ)
    beam = beam.union(cyl_y(R_CON, PZ, -19, 19, 12)).cut(cyl_y(R_CON, PZ, -25, 25, 5.2))
    # end plate for cassette clamp (back plate)
    back = box((110, 114), (-26, 26), (PZ-20, PZ+24))
    lower_lip = box((110, 124), (-26, 26), (PZ-20, PZ-15))         # fixed lower jaw
    posts = cyl_z(119, -23, PZ-15, PZ+34, 4).union(cyl_z(119, 23, PZ-15, PZ+34, 4))
    # striker-bar outrigger (right side): bracket + bar to x=202 at y 45..53, z 112..120
    bracket = box((86, 104), (14, 55), (PZ+13, PZ+25))
    bar = box((86, 196), (47, 55), (PZ+17, PZ+25))
    bracket = safe_fillet(bracket, "|Z", 3.0)
    bar = safe_fillet(bar, "|X", 2.0)
    a = root.union(beam).union(back).union(lower_lip).union(posts).union(bracket).union(bar)
    a = safe_fillet(a, "|X and (not <Z)", 1.5)
    return a

def upper_jaw():
    j = box((114, 124), (-26, 26), (PZ+15, PZ+22))
    j = j.cut(cyl_z(119, -23, PZ+10, PZ+30, 4.4)).cut(cyl_z(119, 23, PZ+10, PZ+30, 4.4))
    j = j.cut(cyl_z(119, 0, PZ+10, PZ+30, 5.2))      # thumbscrew clearance
    return j

def thumbscrew(x, y, z0, z1, d=5, head=14):
    return cyl_z(x, y, z0, z1, d).union(cyl_z(x, y, z1, z1+6, head))

def cassette_and_block():
    cas = box((114, 120.8), (-20, 20), (PZ-14, PZ+14))     # 40 x 28 x 6.8 standard cassette (simplified)
    blk = box((120.8, L_ARM), (-12.5, 12.5), (PZ-10, PZ+10))  # paraffin block 25 x 20 face
    return cas, blk

# ================================================================== CONROD / CRANK / HANDWHEEL
def conrod():
    # fork plates straddling the arm at pin (R_CON, PZ); lower eye at the crank pin rest position (XP0, ZP0), between the two crank webs
    plates = None
    for yc in (-23, 23):
        p = (cyl_y(R_CON, PZ, yc-3, yc+3, 18)
             .union(box((R_CON-5, R_CON+5), (yc-3, yc+3), (PZ-30, PZ)))
             .cut(cyl_y(R_CON, PZ, yc-4, yc+4, 5.2)))
        plates = p if plates is None else plates.union(p)
    cross = box((R_CON-5, R_CON+5), (-26, 26), (PZ-36, PZ-30))
    # body from cross (z=PZ-33) to lower eye, aligned toward (X_CRANK+R_CRANK, Z_CRANK)
    x1, z1 = XP0, ZP0
    x0, z0 = R_CON, PZ-33
    L = math.hypot(x1-x0, z1-z0); ang = math.degrees(math.atan2(z1-z0, x1-x0))
    body = (cq.Workplane("XY").box(L, 10, 8).translate((L/2, 0, 0))
            .rotate((0, 0, 0), (0, 1, 0), -ang).translate((x0, 0, z0)))
    eye = cyl_y(x1, z1, -5, 5, 16).cut(cyl_y(x1, z1, -6, 6, 5.2))
    return plates.union(cross).union(body).union(eye)

def crankshaft():
    # built-up crankshaft: two 8 mm shaft stubs, the conrod runs between the two webs
    return cyl_y(X_CRANK, Z_CRANK, -70, -14, 8).union(cyl_y(X_CRANK, Z_CRANK, 14, 125, 8))

def crank_disc():
    # two crank webs (printed, 44 mm dia) carrying the 5 mm crank pin; shaft stubs grub-screwed into the webs
    d = cyl_y(X_CRANK, Z_CRANK, -14, -6, 44).union(cyl_y(X_CRANK, Z_CRANK, 6, 14, 44))
    d = d.cut(cyl_y(X_CRANK, Z_CRANK, -15, -5, 8.2)).cut(cyl_y(X_CRANK, Z_CRANK, 5, 15, 8.2))
    return d

def crank_pin():
    return cyl_y(XP0, ZP0, -14, 14, 5)

def bearings():
    b = None
    for yc in (-50, 50):
        r = cyl_y(X_CRANK, Z_CRANK, yc-3.5, yc+3.5, 22).cut(cyl_y(X_CRANK, Z_CRANK, yc-4, yc+4, 8))
        b = r if b is None else b.union(r)
    return b

def handwheel():
    rim = cyl_y(X_CRANK, Z_CRANK, 95, 107, 110).cut(cyl_y(X_CRANK, Z_CRANK, 94, 108, 90))
    rim = safe_fillet(rim, "%CIRCLE", 2.5)
    hub = cyl_y(X_CRANK, Z_CRANK, 90, 118, 26)
    w = rim.union(hub)
    for k in range(5):
        a = math.radians(90 + 72*k)
        # spoke along radial direction in the x-z plane at angle a, thickness along y
        sp = (cq.Workplane("XZ", origin=(0, 101, 0)).center(0, 0).rect(8, 38).extrude(-6)
              .translate((0, 0, 31)).rotate((0, 0, 0), (0, 1, 0), -math.degrees(a)+90))
        sp = sp.translate((X_CRANK, 0, Z_CRANK))
        w = w.union(sp)
    w = w.cut(cyl_y(X_CRANK, Z_CRANK, 85, 120, 8.2))
    return w

def handle():
    return cyl_y(X_CRANK, Z_CRANK+45, 107, 135, 16).union(cyl_y(X_CRANK, Z_CRANK+45, 104, 108, 8))

def eccentric_cam():  # reserved (not used in baseline, pawl is arm-driven)
    return None

# ================================================================== BLADE MODULE (monolithic flexure stage)
def stage_body():
    """Monolithic flexure stage (rev. 2): foot clamped to the deck, two 44 mm plates, platform, hanging nut tab, and a
    gusseted upstand that carries the frame nut.  The upstand is continuous with the foot (rev. 1 had a 2 mm gap)."""
    z_foot = Z_DECK + 6                                   # 28
    z_plat0 = z_foot + 44                                 # 72: plates are 44 mm free height
    foot = safe_fillet(box((130, 196), (-45, 45), (Z_DECK, z_foot)), "|Z", 4.0)
    for yc in (-35, 35):   # clamp slots along x (coarse positioning +-10 mm)
        slot = cq.Workplane("XY", origin=(0, 0, Z_DECK-1)).center(161, yc).slot2D(26, 6.6).extrude(12)
        foot = foot.cut(slot)
    front = box((137, 138.0), (-30, 30), (z_foot, z_plat0))
    rearL = box((180, 181.0), (-30, -10), (z_foot, z_plat0))
    rearR = box((180, 181.0), (10, 30), (z_foot, z_plat0))
    platform = safe_fillet(box((128, 190), (-45, 45), (z_plat0, z_plat0+8)), "|Z", 4.0)
    platform = platform.cut(box((183, 195), (-33, 33), (z_plat0-1, z_plat0+9)))   # clearance around the fixed upstand
    for y0 in (-36, 25):                                                              # relief under the blade-clamp screw heads
        platform = platform.cut(box((127, 134.5), (y0, y0+11), (z_plat0+5, z_plat0+9)))
    # gusseted upstand (fixed) carrying the frame nut: 10 mm wall, 64 mm wide, nut pocket opens toward the coupler so the
    # cutting-force reaction pushes the nut INTO the wall
    upstand = box((184, 194), (-32, 32), (z_foot, 82))
    upstand = upstand.cut(cyl_x(0, Z_SCREW, 183, 195, 8.4))
    upstand = upstand.cut(hex_nut(183.9, 188.6, 13.4, 0.1, z=Z_SCREW))
    g_lo = (cq.Workplane("XZ", origin=(0, 16, 0)).polyline([(170, z_foot), (184, z_foot), (184, 46)]).close().extrude(32))
    g_side = None
    for y0 in (-32, 26):
        g = (cq.Workplane("XZ", origin=(0, y0+6, 0)).polyline([(172, z_foot), (184, z_foot), (184, 70), (178, 70)]).close().extrude(6))
        g_side = g if g_side is None else g_side.union(g)
    # platform hanging tab with the M8x1.0 nut pocket (opens toward -x so the +x load pushes the nut into the tab)
    tab = box((150, 160), (-24, 24), (Z_SCREW-12, z_plat0)).cut(cyl_x(0, Z_SCREW, 149, 161, 8.4))
    tab = tab.cut(hex_nut(150, 155, 13.4, 0.1, z=Z_SCREW))
    # tab gussets (rev. 2): full-width on the -x side, two side gussets on the +x side clear of the coupler (|y| > 10)
    tab = tab.union(cq.Workplane("XZ", origin=(0, 24, 0)).polyline([(140, z_plat0), (150, z_plat0), (150, 56)]).close().extrude(48))
    for y0 in (-24, 12):
        tab = tab.union(cq.Workplane("XZ", origin=(0, y0+12, 0)).polyline([(160, z_plat0), (172, z_plat0), (160, 56)]).close().extrude(12))
    # spring anchor pins (preload springs between tab and upstand)
    pins = None
    for yc in (-20, 20):
        p = cyl_x(yc, Z_SCREW-6, 158, 162, 3).union(cyl_x(yc, Z_SCREW-6, 182, 186, 3))
        pins = p if pins is None else pins.union(p)
    # pawl-rod guide column (fixed, on foot) with square bore 8.4 at (200, 66)
    col = box((194, 206), (58, 74), (z_foot, 135)).cut(box((195.8, 204.2), (61.8, 70.2), (Z_DECK, 136)))
    colbase = box((186, 206), (36, 74), (Z_DECK, z_foot))
    sidetab = box((168, 184), (-62, -45), (Z_DECK, z_foot))
    col = col.cut(box((201, 207), (57, 67), (Z_SCREW-3, Z_SCREW+30)))     # slot in the column wall for the pawl finger travel
    s = (foot.union(front).union(rearL).union(rearR).union(platform).union(upstand).union(g_lo).union(g_side)
         .union(tab).union(pins).union(col).union(colbase).union(sidetab))
    return s

LEAN_T = math.tan(math.radians(BLADE_LEAN)); LEAN_C = math.cos(math.radians(BLADE_LEAN))
Z_HOLD_TOP = 92.5          # holder / pressure-strip top: edge protrudes 2.5 mm
def inclined_plate(x_edge_offset, width, height, thick, z_top):
    """plate leaning back BLADE_LEAN from vertical; its top-front corner at (126 + x_edge_offset + (95-z_top)*tan, z_top)."""
    p = (cq.Workplane("XY").box(thick, width, height).translate((thick/2, 0, -height/2))
         .rotate((0, 0, 0), (0, 1, 0), -BLADE_LEAN))
    x0 = 126.0 + (PZ - z_top)*LEAN_T + x_edge_offset
    return p.translate((x0, 0, z_top))

N_LEAN = (LEAN_C, 0.0, math.sin(math.radians(BLADE_LEAN)))      # blade normal (front -> back, upward)
def cyl_n(s_along, y, n0, n1, d):
    """cylinder along the blade normal n; axis passes through the point s_along mm below the strip top along the blade."""
    px = 126.0 + (PZ - Z_HOLD_TOP)*LEAN_T + s_along*math.sin(math.radians(BLADE_LEAN))
    pz = Z_HOLD_TOP - s_along*LEAN_C
    cyl = cq.Workplane("YZ", origin=(n0, 0, 0)).center(y, 0).circle(d/2).extrude(n1-n0)   # along +x from n0..n1
    return cyl.rotate((0, 0, 0), (0, 1, 0), -BLADE_LEAN).translate((px, 0, pz))

def blade():
    # 80 x 8 x 0.25 blade, edge at (126, z=95), leaning back BLADE_LEAN; its rear face bears on the holder's back wall
    return inclined_plate(0.0, 80, 8, 0.25, PZ)

def pressure_strip():
    # 1.0 mm steel clamp strip in front of the blade, top 2.5 mm below the edge, 14 mm tall, screwed to the holder below the blade
    return inclined_plate(-1.0/LEAN_C, 84, 14, 1.0, Z_HOLD_TOP)

def blade_holder():
    """Rev. 2 holder: full-width back wall supports the blade's rear face over 80 mm (the cutting force pushes the blade onto
    this wall); a 0.20 mm recess for the blade so the 1.0 mm strip clamps it when the two M4 screws below the blade are tightened."""
    body = box((126, 160), (-42, 42), (80, Z_HOLD_TOP))
    # remove everything in front of the blade's rear plane
    front_cut = (cq.Workplane("XY").box(40, 90, 60).translate((-20, 0, -10))
                 .rotate((0, 0, 0), (0, 1, 0), -BLADE_LEAN)
                 .translate((126.0 + (PZ-Z_HOLD_TOP)*LEAN_T + 0.25/LEAN_C, 0, Z_HOLD_TOP)))
    body = body.cut(front_cut)
    # strip seat below the blade: 0.20 mm proud of the blade's rear plane (blade is 0.25 -> 0.05 mm clamping interference)
    seat = inclined_plate(0.05/LEAN_C, 84, 8.5, 0.20, Z_HOLD_TOP - 5.5*LEAN_C)
    body = body.union(seat)
    # M4 heat-set insert holes along the blade normal, 10.5 mm below the strip top, y = +-30
    for yc in (-30, 30):
        body = body.cut(cyl_n(10.5, yc, -0.5, 10.0, 5.6))
    # hold-down to platform: 4 x M4 clearance holes (heads on top)
    for xx in (134, 154):
        for yc in (-36, 36):
            body = body.cut(cyl_z(xx, yc, 79, 93, 4.4))
    return body

def blade_clamp_screws():
    g = None
    for yc in (-30, 30):
        head = cyl_n(10.5, yc, -1.0-4.0, -1.0, 7.0)
        shank = cyl_n(10.5, yc, -1.0, 9.0, 4.0)
        s = head.union(shank)
        g = s if g is None else g.union(s)
    return g

def blade_clamp_plate():   # retained for compatibility: now the steel clamp strip
    return pressure_strip()

def feed_screw_coarse():   # M8x1.25 threaded rod: bears on the fine-bolt head inside the coupler -> frame nut -> drum -> knob; translates +-5 mm
    return cyl_x(0, Z_SCREW, 171.3, 224+DXD, 8)

def feed_screw_fine():     # M8x1.0 hex bolt: through the stage-tab nut, hex head captured in the coupler pocket
    return cyl_x(0, Z_SCREW, 144, 165, 8).union(hex_nut(165, 170.3, 13.0, 0.1, z=Z_SCREW))

def coupler():
    """Rev. 2 coupler (13 mm long, dia 20): hex pocket for the fine-bolt head, plain bore, and a captured M8x1.25 thin nut.
    The coarse rod is screwed through the nut until its end bears on the bolt head, so the axial load path is steel-on-steel
    and the printed body carries torque only.  Length limited by the +-5 mm screw excursion between tab (x=160) and upstand (x=184)."""
    cp = cyl_x(0, Z_SCREW, 165, 178, 20)
    cp = cp.cut(hex_nut(164.9, 171, 13.4, 0.1, z=Z_SCREW)).cut(cyl_x(0, Z_SCREW, 164, 179, 8.4))
    cp = cp.cut(hex_nut(174, 178.1, 13.4, 0.1, z=Z_SCREW))
    return cp

def feed_nuts():
    """Steel nuts fixed in the stage (visual): tab nut M8x1.0 and frame nut M8x1.25."""
    return hex_nut(150.1, 154.9, 13.0, 8.1, z=Z_SCREW).union(hex_nut(184.0, 188.5, 13.0, 8.1, z=Z_SCREW))

def coupler_nut():
    """M8x1.25 thin nut captured in the coupler (travels with the screw)."""
    return hex_nut(174.1, 178.0, 13.0, 8.1, z=Z_SCREW)

def preload_springs():
    """Two extension springs (dia 5) between the tab pins and the upstand pins: keep the tab nut on one thread flank in +x
    over the whole +-1 mm travel (the symmetric flexure alone reverses its preload at mid-travel)."""
    s = None
    for yc in (-20, 20):
        sp = cyl_x(yc, Z_SCREW-6, 161, 183, 5).cut(cyl_x(yc, Z_SCREW-6, 160, 184, 3.6))
        s = sp if s is None else s.union(sp)
    return s

def ratchet_wheel():
    pts = []
    for k in range(N_TEETH):
        a0 = 2*math.pi*k/N_TEETH; a1 = 2*math.pi*(k+1)/N_TEETH
        pts.append((R_RATCHET*math.cos(a0), R_RATCHET*math.sin(a0)))
        pts.append(((R_RATCHET-2.6)*math.cos(a1), (R_RATCHET-2.6)*math.sin(a1)))
    prof = cq.Workplane("YZ", origin=(198+DXD, 0, 0)).center(0, Z_SCREW).polyline(pts).close().extrude(14)   # 14 mm drum: pawl stays engaged over +-5 mm screw travel
    w = prof.cut(cyl_x(0, Z_SCREW, 197+DXD, 213+DXD, 8.2))
    w = w.cut(hex_nut(198+DXD, 202.5+DXD, 13.4, 0.1, z=Z_SCREW)).cut(hex_nut(207.5+DXD, 212+DXD, 13.4, 0.1, z=Z_SCREW))   # jam-nut pockets both ends
    return w

def selector_shield():
    # sector 2 mm thick covering teeth; rotatable on hub (friction detent). Open sector ~ 60 deg where pawl engages
    ring = (cq.Workplane("YZ", origin=(197+DXD, 0, 0)).center(0, Z_SCREW).circle(R_RATCHET+3.0).circle(R_RATCHET+0.6).extrude(16))
    # opening sector on the +y side where the pawl finger engages (selector rotates to expose 1-6 teeth)
    cut = (cq.Workplane("YZ", origin=(196+DXD, 0, 0)).center(0, Z_SCREW)
           .polyline([(0, 0), (40*math.cos(math.radians(-35)), 40*math.sin(math.radians(-35))),
                      (40, 0), (40*math.cos(math.radians(35)), 40*math.sin(math.radians(35)))]).close().extrude(18))
    # strap leaves the ring at the top (outside the cut sector) and runs to the pawl-rod column, where an M3 screw clamps it
    strap = box((206, 212+DXD), (-6, 58), (Z_SCREW+R_RATCHET+0.8, Z_SCREW+R_RATCHET+3.0))
    strap = strap.union(box((206, 212+DXD), (52, 58), (Z_SCREW+R_RATCHET-6, Z_SCREW+R_RATCHET+3.0)))
    return ring.cut(cut).union(strap)

def feed_knob():
    k = cyl_x(0, Z_SCREW, 214+DXD, 224+DXD, 28).cut(cyl_x(0, Z_SCREW, 213+DXD, 225+DXD, 8.2))
    for i in range(12):
        a = math.radians(30*i)
        k = k.cut(cyl_x(14*math.sin(a), Z_SCREW+14*math.cos(a), 213+DXD, 225+DXD, 4))
    return k

def pawl_rod():
    rod = box((196, 204), (62, 70), (52, 147))
    finger = box((199+DXD, 203+DXD), (R_RATCHET-2.0, 66), (Z_SCREW-2, Z_SCREW+6))   # reaches wheel rim at y=+25 (engages tooth)
    foot = box((186, 204), (45, 70), (141, 147))   # flag over the striker bar (y 47..55); lifted ~22 mm over the last ~6 deg of upstroke
    return rod.union(finger).union(foot)

def tie_rods():
    r = None
    for yc in (-55, 55):
        rod = cyl_x(yc, -45, X_BASE0-8, X_BASE1+8, 8)
        nuts = hex_nut(X_BASE0-7, X_BASE0-0.5, 13, 8.1, y=yc, z=-45).union(hex_nut(X_BASE1+0.5, X_BASE1+7, 13, 8.1, y=yc, z=-45))
        part = rod.union(nuts)
        r = part if r is None else r.union(part)
    return r

def fastener_heads():
    f = None
    for yc in (-50, 50):
        for dx in (-9, 9):
            h = bolt_head_z(X_CRANK+dx, yc, Z_CRANK+13, 6.5, 3.5, hex_=False)
            f = h if f is None else f.union(h)
    for xx in (134, 154):     # blade holder hold-down screws (4 x M4)
        for yc in (-36, 36):
            f = f.union(bolt_head_z(xx, yc, Z_HOLD_TOP, 5.5, 2.5, hex_=False))
    f = f.union(cyl_x(55, Z_SCREW+R_RATCHET-1.5, 212+DXD, 215+DXD, 5.5))   # M3 head: selector strap to column
    return f

def knob_bolts():
    k = None
    for yc in (-35, 35):
        b = cyl_z(161, yc, Z_DECK-10, Z_DECK+8, 6).union(cyl_z(161, yc, Z_DECK+8, Z_DECK+16, 20))
        k = b if k is None else k.union(b)
    return k

def coarse_screw():
    return cyl_x(-54, Z_DECK+4, 184, 215, 6).union(cyl_x(-54, Z_DECK+4, 215, 223, 20))

# ================================================================== BUILD
def build():
    add("frame", frame(), (0.82, 0.82, 0.84))
    strips, clamps_f, clamps_a = pivot_strips_and_clamps()
    add("pivot_strips_steel", strips, (0.35, 0.38, 0.45), "spring steel feeler gauge 0.20 mm", 3, False)
    add("pivot_clamp_plates_frame", clamps_f, (0.25, 0.25, 0.28), "PLA (FDM) + M3 screws", 3)
    add("pivot_clamp_plates_arm", clamps_a, (0.25, 0.25, 0.28), "PLA (FDM) + M3 screws", 3)
    add("arm", arm(), (0.92, 0.55, 0.15))
    add("cassette_upper_jaw", upper_jaw(), (0.92, 0.55, 0.15))
    add("cassette_thumbscrew", thumbscrew(119, 0, PZ+15, PZ+30), (0.2, 0.2, 0.22), "M5 thumbscrew", 1, False)
    cas, blk = cassette_and_block()
    add("histology_cassette", cas, (0.95, 0.95, 0.97), "standard cassette (consumable)", 1, False)
    add("paraffin_block", blk, (0.98, 0.92, 0.75), "paraffin specimen block", 1, False)
    add("conrod", conrod(), (0.2, 0.55, 0.85))
    add("crankshaft", crankshaft(), (0.55, 0.55, 0.6), "8 mm steel shaft", 1, False)
    add("crank_disc", crank_disc(), (0.2, 0.55, 0.85))
    add("crank_pin", crank_pin(), (0.3, 0.3, 0.33), "M5 shoulder bolt", 1, False)
    add("bearings_608", bearings(), (0.45, 0.45, 0.5), "608ZZ ball bearing", 2, False)
    add("handwheel", handwheel(), (0.2, 0.55, 0.85))
    add("handle", handle(), (0.15, 0.15, 0.17))
    add("stage_body_flexure", stage_body(), (0.35, 0.75, 0.45))
    add("blade_holder", blade_holder(), (0.25, 0.6, 0.35))
    add("pressure_strip_steel", pressure_strip(), (0.6, 0.6, 0.66), "1.0 mm steel feeler-gauge strip 84x14, 2 holes", 1, False)
    add("blade_clamp_screws_M4", blade_clamp_screws(), (0.2, 0.2, 0.22), "M4 x 12 cap screw x2 + heat-set inserts", 2, False)
    add("blade_low_profile", blade(), (0.9, 0.9, 0.95), "disposable low-profile microtome blade 80x8x0.25", 1, False)
    add("feed_screw_M8x1_25", feed_screw_coarse(), (0.5, 0.5, 0.55), "M8x1.25 threaded rod 55 mm", 1, False)
    add("feed_screw_M8x1_0", feed_screw_fine(), (0.6, 0.6, 0.65), "M8x1.0 fine hex bolt 25 mm", 1, False)
    add("screw_coupler", coupler(), (0.85, 0.25, 0.25))
    add("feed_nuts_steel", feed_nuts(), (0.5, 0.5, 0.55), "M8x1.0 nut (tab) + M8x1.25 nut (frame)", 2, False)
    add("coupler_nut_M8x1_25", coupler_nut(), (0.5, 0.5, 0.55), "M8x1.25 thin nut (captured in coupler)", 1, False)
    add("preload_springs", preload_springs(), (0.55, 0.55, 0.6), "extension spring 5 x 20 mm, ~5 N", 2, False)
    add("ratchet_wheel_50T", ratchet_wheel(), (0.85, 0.25, 0.25))
    add("thickness_selector_shield", selector_shield(), (0.95, 0.75, 0.2))
    add("feed_knob", feed_knob(), (0.85, 0.25, 0.25))
    add("pawl_rod", pawl_rod(), (0.85, 0.25, 0.25))
    add("stage_clamp_knobs_M6", knob_bolts(), (0.2, 0.2, 0.22), "M6 knob bolt + heat-set insert", 2, False)
    add("fastener_heads_M3_M4", fastener_heads(), (0.25, 0.25, 0.28), "M4 x 12 socket screws x8 (visual)", 8, False)
    add("tie_rods_M8_steel_spine", tie_rods(), (0.5, 0.5, 0.55), "M8x250 threaded rod + 4 nuts (base spine)", 2, False)
    add("coarse_approach_screw_M6", coarse_screw(), (0.2, 0.2, 0.22), "M6 x 40 bolt + knob", 1, False)

    assy = cq.Assembly(name="FlexTome")
    meta = []
    for p in PARTS:
        assy.add(p["solid"], name=p["name"], color=cq.Color(*p["color"]))
        cq.exporters.export(p["solid"], f"stl/{p['name']}.stl", tolerance=0.05, angularTolerance=0.1)
        vol = p["solid"].val().Volume() if hasattr(p["solid"], "val") else 0
        meta.append(dict(name=p["name"], color=p["color"], material=p["material"], qty=p["qty"],
                         printed=p["printed"], volume_cm3=round(vol/1000, 2)))
    assy.save("step/flextome_assembly.step")
    json.dump(meta, open("parts.json", "w"), indent=1)
    return meta

if __name__ == "__main__":
    import time; t = time.time()
    meta = build()
    print(f"built {len(meta)} parts in {time.time()-t:.0f}s")
    for m in meta: print(f"{m['name']:30s} {m['volume_cm3']:8.1f} cm3  printed={m['printed']}")
