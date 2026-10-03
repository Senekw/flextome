"""Labelled assembly drawings: numbered-balloon isometric with parts legend, and annotated mechanism section."""
import numpy as np, matplotlib.pyplot as plt, matplotlib as mpl, json
import render
parts = render.load_parts()
meta = {m["name"]: m for m in json.load(open("parts.json"))}

def balloon_fig(out, azim, elev, items, W=1800, H=1300, title=None, hide=(), section_y=None, legend_cols=2, fs=8.5):
    img, proj = render.render_cad(parts, "_tmp.png", azim=azim, elev=elev, W=W, H=H, hide=hide, section_y=section_y)
    fig = plt.figure(figsize=(W/160, (H+330)/160), dpi=160)
    ax = fig.add_axes([0, 330/(H+330), 1, H/(H+330)]); ax.imshow(img); ax.axis("off")
    for k, (num, label, pt3, (fx, fy)) in enumerate(items):
        px, py = proj(pt3)
        tx, ty = fx*W, (1-fy)*H
        ax.annotate("", xy=(px, py), xytext=(tx, ty), arrowprops=dict(arrowstyle="-", color="0.15", lw=0.9, shrinkA=11, shrinkB=0))
        ax.plot(px, py, "o", ms=3, color="0.15")
        ax.text(tx, ty, str(num), ha="center", va="center", fontsize=fs+1, weight="bold",
                bbox=dict(boxstyle="circle,pad=0.35", fc="white", ec="0.15", lw=1.0))
    if title: ax.text(0.01, 0.985, title, transform=ax.transAxes, fontsize=12, weight="bold", va="top")
    # legend
    lax = fig.add_axes([0.02, 0.0, 0.96, 330/(H+330)]); lax.axis("off")
    n = len(items); rows = int(np.ceil(n/legend_cols))
    for k, (num, label, pt3, _) in enumerate(items):
        col, row = k // rows, k % rows
        lax.text(col/legend_cols, 1 - (row+0.5)/(rows+0.5), f"{num:>2}  {label}", fontsize=fs, va="center", family="DejaVu Sans")
    fig.savefig(out, dpi=160); plt.close(fig)

# ---------------- Figure: numbered isometric
items_iso = [
 (1, "Frame: ribbed base 235x150x55, gusseted tower (one print)", (-30, -40, 60), (0.07, 0.62)),
 (2, "Cross-strip flexure pivot: 3 x 0.20 mm spring-steel strips + clamp pads", (0, -20, 95), (0.08, 0.88)),
 (3, "Rocking arm: tapered hollow box beam, 125 mm pivot-to-face", (55, -22, 108), (0.27, 0.96)),
 (4, "Cassette clamp: lower lip, sliding jaw, M5 thumbscrew", (117, -22, 118), (0.46, 0.99)),
 (5, "Histology cassette + paraffin block (consumable)", (122, -14, 95), (0.50, 0.83)),
 (6, "Disposable low-profile blade (80x8x0.25), 27 deg lean", (126, -28, 95), (0.60, 0.93)),
 (7, "Steel pressure strip 1.0 mm + blade holder, 2 x M4 grub screws", (140, -42, 90), (0.57, 0.66)),
 (8, "Parallelogram flexure stage (monolithic print): plates 1.0x44 mm", (137, -32, 60), (0.57, 0.45)),
 (9, "Stage clamp knobs M6 (coarse approach slots +/-10 mm)", (161, -35, 46), (0.60, 0.30)),
 (10, "Coarse-approach screw M6 (side boss)", (215, -54, 34), (0.78, 0.22)),
 (11, "Differential feed screw: M8x1.0 bolt / coupler / M8x1.25 rod", (160, 0, 70), (0.70, 0.60)),
 (12, "Ratchet drum 50 T x 14 mm (5 um per tooth) + jam nuts", (205, 0, 45), (0.80, 0.42)),
 (13, "Thickness-selector shield (1-6 teeth = 5-30 um)", (205, 28, 92), (0.86, 0.55)),
 (14, "Feed knob", (219, 0, 70), (0.90, 0.33)),
 (15, "Pawl rod with finger; lifted by striker bar at top of stroke", (200, 66, 140), (0.90, 0.86)),
 (16, "Striker bar on arm (feed actuation)", (150, 51, 116), (0.72, 0.95)),
 (17, "Connecting rod (fork) and crank disc, r = 13 mm", (60, -23, 80), (0.32, 0.50)),
 (18, "Crankshaft 8 mm in 2 x 608 bearings, pillow blocks", (60, -50, 25), (0.25, 0.28)),
 (19, "Handwheel 110 mm with handle (one turn = one section)", (60, 101, 25), (0.46, 0.14)),
 (20, "M8 tie rods (base pre-compression / split-base option)", (195, -55, -45), (0.84, 0.10)),
]
balloon_fig("fig_assembly_callouts.png", -55, 25, items_iso, title="FlexTome assembly - numbered components")

# ---------------- Figure: annotated section (y = 0 cut, viewed from +y)
items_sec = [
 (1, "Gusseted tower: 60 mm root depth, FE compliance 0.21 um/N", (-20, 0, 40), (0.93, 0.45)),
 (2, "Cross-strip pivot (X): rotation about the strip crossing; zero backlash", (0, 0, 95), (0.93, 0.82)),
 (3, "Arm root clamp pads (M3) + hollow box beam", (40, 0, 95), (0.72, 0.96)),
 (4, "Conrod joint at r = 60 mm; crank pin r = 13 mm -> theta = +/-12.5 deg", (60, 0, 95), (0.58, 0.96)),
 (5, "Crank webs (pin between them) / shaft axis z = 23.8", (60, 0, 25), (0.64, 0.16)),
 (6, "Cassette + block; face on 125 mm arc, 54 mm stroke", (122, 0, 95), (0.44, 0.96)),
 (7, "Blade edge at (126, 95): lean 27 deg, clearance 9.5 deg, rake 45.5 deg", (126, 0, 95), (0.35, 0.92)),
 (8, "Pressure strip top 2.5 mm below edge; outside the swept circle", (127.3, 0, 92), (0.19, 0.86)),
 (9, "Platform + blade holder (2 x M4 grub screws from rear)", (150, 0, 90), (0.31, 0.95)),
 (10, "Front flexure plate 1.0 x 44 mm (60 wide)", (137.5, 0, 58), (0.42, 0.42)),
 (11, "Rear flexure plates 2 x 20 wide", (180.5, 0, 58), (0.17, 0.21)),
 (12, "Stage tab with M8x1.0 nut (moves +/-1 mm)", (155, 0, 65), (0.47, 0.66)),
 (13, "Coupler (bolt head pocket + jam nut) - translates +/-5 mm", (175, 0, 70), (0.03, 0.62)),
 (14, "Frame nut M8x1.25 in fixed upstand", (189, 0, 70), (0.03, 0.72)),
 (15, "Ratchet drum 14 mm long: pawl stays engaged over screw travel", (205, 0, 70), (0.06, 0.90)),
 (16, "Feed knob", (219, 0, 70), (0.03, 0.50)),
 (17, "Rib pockets (open underside) + M8 tie rod at z = -45", (100, 0, -30), (0.50, 0.05)),
 (18, "Raised deck (z = 30): stage clamp slots, coarse approach", (160, 0, 30), (0.25, 0.09)),
]
balloon_fig("fig_section_callouts.png", 90, 0, items_sec, W=1800, H=1200, section_y=0.5, hide=("handwheel", "handle"),
            title="Mechanism section at y = 0 (viewed from +y)")
print("drawings done")
