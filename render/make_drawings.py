"""Labelled assembly drawings: numbered-balloon isometric with parts legend, and annotated mechanism section.
Balloons are placed automatically on a perimeter ring in the same angular order as their targets, so leaders do not cross."""
import numpy as np, matplotlib.pyplot as plt, json
import render
parts = render.load_parts()

def ring_positions(targets, W, H, margin=0.055, cx=None, cy=None, min_sep_px=90):
    """Balloon positions on a rectangular ring: each target is projected radially from the image centre onto the ring,
    then neighbours are spread along the ring to a minimum spacing while keeping angular order (so leaders never cross)."""
    t = np.asarray(targets, float)
    cx = W/2 if cx is None else cx; cy = H/2 if cy is None else cy
    x0, x1, y0, y1 = margin*W, (1-margin)*W, margin*H, (1-margin)*H
    corners = [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
    seg_len = [x1-x0, y1-y0, x1-x0, y1-y0]; L = sum(seg_len)
    def point_at(d):
        d %= L
        for (px, py), ln, k in zip(corners, seg_len, range(4)):
            if d <= ln:
                return (px + d, py) if k == 0 else (px, py + d) if k == 1 else (px - d, py) if k == 2 else (px, py - d)
            d -= ln
        return corners[0]
    def ray_param(px, py):
        dx, dy = px - cx, py - cy
        if abs(dx) < 1e-9 and abs(dy) < 1e-9: dx = 1.0
        # intersect ray with rectangle
        cands = []
        if dx > 0: cands.append((x1 - cx)/dx)
        if dx < 0: cands.append((x0 - cx)/dx)
        if dy > 0: cands.append((y1 - cy)/dy)
        if dy < 0: cands.append((y0 - cy)/dy)
        k = min(c for c in cands if c > 0)
        ix, iy = cx + k*dx, cy + k*dy
        # perimeter distance of (ix, iy), clockwise from top-left corner
        if abs(iy - y0) < 1e-6: return ix - x0
        if abs(ix - x1) < 1e-6: return seg_len[0] + (iy - y0)
        if abs(iy - y1) < 1e-6: return seg_len[0] + seg_len[1] + (x1 - ix)
        return seg_len[0] + seg_len[1] + seg_len[2] + (y1 - iy)
    d = np.array([ray_param(px, py) for px, py in t])
    order = np.argsort(d); ds = d[order].copy()
    # spread: iterate a few passes of cyclic min-separation relaxation
    n = len(ds)
    for _ in range(50):
        moved = False
        for k in range(n):
            a, b = ds[k], ds[(k+1) % n] + (L if k == n-1 else 0)
            gap = b - a
            if gap < min_sep_px:
                shift = (min_sep_px - gap)/2
                ds[k] -= shift; 
                if k == n-1: ds[0] += shift
                else: ds[k+1] += shift
                moved = True
        if not moved: break
    out = np.empty_like(t); out[order] = np.array([point_at(v) for v in ds])
    return out

def balloon_fig(out, azim, elev, items, W=1800, H=1300, title=None, hide=(), section_y=None, legend_cols=2, fs=10, ring_margin=0.055, ring_center=(0.5, 0.5), insets=()):
    img, proj = render.render_cad(parts, "_tmp.png", azim=azim, elev=elev, W=W, H=H, hide=hide, section_y=section_y)
    tgt = np.array([proj(pt3) for (_, _, pt3) in items])
    pos = ring_positions(tgt, W, H, margin=ring_margin, cx=ring_center[0]*W, cy=ring_center[1]*H)
    rows = int(np.ceil(len(items)/legend_cols))
    LH = 40*rows + 60                          # legend strip height (px)
    IH = int(0.38*H) if insets else 0          # inset strip height (px)
    TOT = H + IH + LH
    fig = plt.figure(figsize=(W/160, TOT/160), dpi=160)
    ax = fig.add_axes([0, (IH+LH)/TOT, 1, H/TOT]); ax.imshow(img); ax.axis("off")
    for (num, label, pt3), (px, py), (tx, ty) in zip(items, tgt, pos):
        ax.annotate("", xy=(px, py), xytext=(tx, ty), arrowprops=dict(arrowstyle="-", color="0.2", lw=0.9, shrinkA=13, shrinkB=0))
        ax.plot(px, py, "o", ms=3.2, color="0.2", mec="white", mew=0.6)
        ax.text(tx, ty, str(num), ha="center", va="center", fontsize=fs+2, weight="bold",
                bbox=dict(boxstyle="circle,pad=0.3", fc="white", ec="0.2", lw=1.1))
    if title: ax.text(0.5, 0.995, title, transform=ax.transAxes, fontsize=14, weight="bold", va="top", ha="center")
    # close-up insets in their own strip below the render; source boxes marked on the render
    if insets:
        n = len(insets); gap = 0.02; wfrac = (1 - gap*(n+1))/n
        for k, (lab, (bx0, by0, bx1, by1)) in enumerate(insets):
            ax.add_patch(plt.Rectangle((bx0, by0), bx1-bx0, by1-by0, fill=False, ec="0.2", lw=1.0, ls="--"))
            ax.text(bx0+6, by0-6, lab.split(" (")[0], fontsize=fs-1, va="bottom", color="0.2", bbox=dict(fc="white", ec="none", alpha=0.8, pad=1))
            iax = fig.add_axes([gap + k*(wfrac+gap), LH/TOT + 0.10*IH/TOT, wfrac, 0.78*IH/TOT])
            iax.imshow(img[by0:by1, bx0:bx1], aspect="equal"); iax.set_xticks([]); iax.set_yticks([])
            for sp in iax.spines.values(): sp.set_edgecolor("0.2"); sp.set_linewidth(1.0)
            iax.set_title("close-up: " + lab, fontsize=fs, pad=4, loc="left")
            for (num, label, pt3), (px, py) in zip(items, tgt):
                if bx0 <= px <= bx1 and by0 <= py <= by1:
                    iax.plot(px-bx0, py-by0, "o", ms=4, color="0.2", mec="white", mew=0.7)
                    iax.annotate(str(num), xy=(px-bx0, py-by0), xytext=(14, -14), textcoords="offset points", fontsize=fs, weight="bold",
                                 bbox=dict(boxstyle="circle,pad=0.22", fc="white", ec="0.2", lw=0.8), arrowprops=dict(arrowstyle="-", color="0.2", lw=0.7))
    lax = fig.add_axes([0.03, 0.0, 0.94, LH/TOT]); lax.axis("off")
    for k, (num, label, pt3) in enumerate(items):
        col, row = k // rows, k % rows
        lax.text(col/legend_cols, 1 - (row+0.6)/(rows+0.6), f"{str(num):>3}", fontsize=fs, va="center", weight="bold")
        lax.text(col/legend_cols + 0.025, 1 - (row+0.6)/(rows+0.6), label, fontsize=fs, va="center")
    fig.savefig(out, dpi=160); plt.close(fig)

# ---------------- Figure: numbered isometric
items_iso = [
 (1, "Frame: ribbed base 235 × 150 × 55 mm, gusseted tower (one print)", (-30, -40, 60)),
 (2, "Cross-strip flexure pivot: 3 × 0.20 mm steel strips in clamp pads", (0, -20, 95)),
 (3, "Rocking arm: tapered hollow box beam, 125 mm pivot to face", (55, -22, 108)),
 (4, "Cassette clamp: lower lip, sliding jaw, M5 thumbscrew", (117, -22, 118)),
 (5, "Histology cassette and paraffin block (consumable)", (122, -14, 95)),
 (6, "Disposable low-profile blade 80 × 8 × 0.25 mm, 27° lean", (126, -28, 95)),
 (7, "Blade holder (full-width back wall), 1.0 mm clamp strip, 2 × M4", (140, -42, 90)),
 (8, "Parallelogram flexure stage (one print), plates 1.0 × 44 mm", (137, -32, 52)),
 (9, "Stage clamp knobs M6 (coarse-approach slots ±10 mm)", (161, -35, 38)),
 (10, "Coarse-approach screw M6 in side boss", (215, -54, 26)),
 (11, "Differential feed screw: M8×1.0 bolt, coupler, M8×1.25 rod", (170, 0, 60)),
 (12, "Ratchet drum 50 T × 14 mm (5 µm per tooth) with jam nuts", (208, 0, 35)),
 (13, "Thickness-selector shield (1 to 6 teeth = 5 to 30 µm)", (208, 28, 82)),
 (14, "Feed knob", (222, 0, 60)),
 (15, "Pawl rod with finger, lifted by the striker bar at top of stroke", (200, 66, 140)),
 (16, "Striker bar on arm (feed actuation)", (150, 51, 116)),
 (17, "Fork connecting rod 71.2 mm and crank webs, r = 13 mm", (60, -23, 80)),
 (18, "Built-up crankshaft, 8 mm stubs in 2 × 608 bearings", (60, -50, 25)),
 (19, "Handwheel 110 mm with handle (one turn = one section)", (60, 101, 25)),
 (20, "M8 tie rods (base pre-compression, split-base option)", (195, -55, -45)),
]
balloon_fig("fig_assembly_callouts.png", -55, 25, items_iso, title="FlexTome assembly: numbered components",
            insets=[("pivot and arm root (2, 3)", (380, 150, 760, 430)), ("blade and clamp (4 to 7)", (770, 220, 1070, 470)),
                    ("feed train (11 to 14)", (1040, 360, 1420, 720))])

# ---------------- Figure: annotated section (y = 0 cut, viewed from +y)
items_sec = [
 ("1",   "Frame tower: 60 mm root depth, FE compliance 0.21 µm/N", (-20, 0, 40)),
 ("2",   "Cross-strip pivot: rotates about the strip crossing, zero backlash", (0, 0, 95)),
 ("3",   "Arm: root clamp pads (M3) and hollow box beam", (40, 0, 95)),
 ("17",  "Connecting rod: joint at r = 60 mm; crank r = 13 mm gives ±12.5°", (60, 0, 95)),
 ("18",  "Crankshaft: two webs with pin between, axis at z = 23.8 mm", (60, 0, 25)),
 ("5",   "Cassette and block: face on 125 mm arc, 54 mm stroke", (122, 0, 95)),
 ("6",   "Blade: edge at (126, 95), lean 27°, clearance 9.5°, rake 45.5°", (126, 0, 95)),
 ("7a",  "Clamp strip (part of 7): top 2.5 mm below the edge", (127.3, 0, 92)),
 ("7b",  "Holder back wall (part of 7): carries the cutting force", (150, 0, 86)),
 ("8a",  "Stage front plate (part of 8): 1.0 × 44 mm, 60 mm wide", (137.5, 0, 50)),
 ("8b",  "Stage rear plates (part of 8): 2 × 20 mm wide", (180.5, 0, 50)),
 ("8c",  "Stage tab (part of 8): M8×1.0 nut on washer, moves ±1 mm", (155, 0, 55)),
 ("11a", "Coupler (part of 11): rod end bears on bolt head; ±5 mm", (171, 0, 60)),
 ("1a",  "Upstand (part of 1): frame nut on washer, FE 0.19 µm/N", (189, 0, 60)),
 ("12",  "Ratchet drum: 14 mm long, pawl engaged over screw travel", (208, 0, 60)),
 ("14",  "Feed knob", (222, 0, 60)),
 ("1b",  "Base rib pockets (part of 1), M8 tie rod (20) at z = −45 mm", (100, 0, -30)),
 ("1c",  "Raised deck (part of 1, z = 22 mm): clamp slots, coarse approach", (160, 0, 22)),
]
balloon_fig("fig_section_callouts.png", 90, 0, items_sec, W=1800, H=1200, section_y=0.5, hide=("handwheel", "handle"),
            title="Mechanism section at y = 0 (viewed from +y)", ring_center=(0.33, 0.42),
            insets=[("cross-strip pivot (2, 3)", (1150, 150, 1450, 400)), ("blade, strip and holder (6, 7a, 7b)", (520, 180, 760, 380)),
                    ("feed train (8c, 11a, 1a, 12, 14)", (110, 300, 620, 680))])
print("drawings done")
