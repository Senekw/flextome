"""Mechanics explainer figures M1-M5 (schematics with governing equations and design numbers)."""
import numpy as np, json, matplotlib.pyplot as plt, matplotlib as mpl
from matplotlib.patches import Polygon, Circle, Rectangle, FancyArrowPatch, Arc, Wedge
from scipy.optimize import brentq
apply_figure_style(sizes=(9, 8, 7))
dr = json.load(open("design_results.json")); er = json.load(open("extra_results.json"))
C = dict(frame="0.55", arm="#d98b2b", crank="#2f78b8", stage="#3aa655", steel="0.25", blade="0.1", block="#e9d8a6", acc="#c0392b")

# ---------------- geometry (mm) ----------------
PZ, L_ARM, R_CON, R_CRANK, L_CONROD = 95.0, 125.0, 60.0, 13.0, 71.2
XC, ZC = 60.0, PZ - L_CONROD   # crank centre 71.2 mm below the arm joint (matches CAD)
HALF_BLOCK = 10.0
def arm_pt(r, th): return np.array([r*np.cos(th), PZ + r*np.sin(th)])
def theta_of_phi(phi):
    pin = np.array([XC + R_CRANK*np.cos(phi), ZC + R_CRANK*np.sin(phi)])
    f = lambda th: np.linalg.norm(arm_pt(R_CON, th) - pin) - L_CONROD
    return brentq(f, -0.6, 0.6)
phi = np.linspace(0, 2*np.pi, 721)
th = np.array([theta_of_phi(p) for p in phi])
th_max = th.max(); th_min = th.min()
zf = PZ + L_ARM*np.sin(th); xf = L_ARM*np.cos(th)
stroke = zf.max() - zf.min()
# cutting window: face passes the edge height while |z_face - 95| < HALF_BLOCK
cut = np.abs(zf - PZ) < HALF_BLOCK
# split into down-stroke (cutting) and up-stroke (return)
dz = np.gradient(zf, phi)
cut_down = cut & (dz < 0)
omega = 2*np.pi            # 1 rev/s (60 rpm)
t_cut_ms = (phi[cut_down].max()-phi[cut_down].min())/omega*1000 if cut_down.any() else float("nan")
v_face = dz*omega           # mm/s at 1 rev/s
print(f"theta range {np.degrees(th_min):.2f}..{np.degrees(th_max):.2f} deg, stroke {stroke:.1f} mm, cut window {t_cut_ms:.0f} ms at 60 rpm, |v| in cut {np.abs(v_face[cut_down]).min():.0f}-{np.abs(v_face[cut_down]).max():.0f} mm/s")

# ================= M1: crank-rocker kinematics =================
fig = plt.figure(figsize=(7.2, 6.4)); gs = fig.add_gridspec(2, 2, height_ratios=[1.35, 1], hspace=0.38, wspace=0.3)
ax = fig.add_subplot(gs[0, :]); ax.set_aspect("equal")
# frame outline
ax.add_patch(Polygon([(-40, -55), (235, -55), (235, 0), (20, 0), (20, 55), (-8, 80), (-8, 118), (-40, 118)], closed=True, fc="0.92", ec="0.6", lw=0.8))
ax.add_patch(Rectangle((128, 30), 62, 58, fc="#d9f0dc", ec=C["stage"], lw=0.8))
for k, (t_, ls, a) in enumerate([(0.0, "-", 1.0), (th_max, "--", 0.55), (th_min, "--", 0.55)]):
    p0 = np.array([0, PZ]); pc = arm_pt(R_CON, t_); pf = arm_pt(L_ARM, t_)
    ax.plot([p0[0], pf[0]], [p0[1], pf[1]], color=C["arm"], lw=5 if ls == "-" else 2.5, ls=ls, alpha=a, solid_capstyle="round")
    # block on arm end
    n = np.array([-np.sin(t_), np.cos(t_)]); u = np.array([np.cos(t_), np.sin(t_)])
    blk = [pf - 12*u + HALF_BLOCK*n, pf + HALF_BLOCK*n, pf - HALF_BLOCK*n, pf - 12*u - HALF_BLOCK*n]
    ax.add_patch(Polygon(blk, closed=True, fc=C["block"], ec="0.3", lw=0.7, alpha=a))
    # crank pin for this theta: solve phi
    ph = phi[np.argmin(np.abs(th - t_) + (phi > np.pi)*0 )]
    pin = np.array([XC + R_CRANK*np.cos(ph), ZC + R_CRANK*np.sin(ph)])
    ax.plot([pin[0], pc[0]], [pin[1], pc[1]], color=C["crank"], lw=3 if ls == "-" else 1.8, ls=ls, alpha=a)
    ax.plot([XC, pin[0]], [ZC, pin[1]], color=C["crank"], lw=3 if ls == "-" else 1.8, ls=ls, alpha=a)
    ax.plot(*pin, "o", ms=4, color=C["crank"], alpha=a); ax.plot(*pc, "o", ms=4, color=C["crank"], alpha=a)
ax.add_patch(Circle((XC, ZC), R_CRANK, fill=False, ec=C["crank"], lw=0.8, ls=":"))
ax.plot(XC, ZC, "o", ms=5, color=C["crank"]); ax.plot(0, PZ, "X", ms=9, color=C["steel"])
# swept arc of face
tt = np.linspace(th_min, th_max, 50); ax.plot(L_ARM*np.cos(tt), PZ + L_ARM*np.sin(tt), color=C["acc"], lw=1, ls="-.")
# blade edge
ax.plot(126, PZ, "v", ms=7, color=C["blade"]); ax.annotate("blade edge (126, 95)", (126, PZ), (160, 118), fontsize=7, arrowprops=dict(arrowstyle="-", lw=0.6))
# dimension annotations
def dim(ax, p, q, text, off=(0, 0), color="0.2"):
    p, q = np.array(p), np.array(q)
    ax.annotate("", xy=q, xytext=p, arrowprops=dict(arrowstyle="<->", lw=0.7, color=color, shrinkA=0, shrinkB=0))
    m = (p + q)/2 + np.array(off); ax.text(*m, text, fontsize=7, color=color, ha="center", va="center", bbox=dict(fc="white", ec="none", pad=0.5))
dim(ax, (0, 132), (125, 132), "L$_{arm}$ = 125 mm", (0, 5))
dim(ax, (0, 14), (60, 14), "x$_c$ = 60", (0, 5))
dim(ax, (XC-R_CRANK, ZC-19), (XC+R_CRANK, ZC-19), "2r = 26", (0, -6))
dim(ax, (230, zf.min()), (230, zf.max()), "stroke\n54.2 mm", (-16, 0), C["acc"])
ax.text(-4, PZ-11, "flexure pivot\nP = (0, 95)", fontsize=7, ha="right", va="top")
ax.text(XC+R_CRANK+6, 58, "rod 71.2", fontsize=7, color=C["crank"], rotation=80)
ax.text(20, 108, r"$\theta=\pm$%.1f$^\circ$" % np.degrees(th_max), fontsize=8, color=C["arm"])
ax.text(132, 42, "flexure\nstage", fontsize=7, color=C["stage"], ha="center")
ax.set_xlim(-50, 245); ax.set_ylim(-60, 142); ax.set_xlabel("x (mm)"); ax.set_ylabel("z (mm)")
ax.set_title("(a) Crank-rocker linkage: crank r = 13 mm, arm joint at 60 mm, connecting rod 71.2 mm", fontsize=8, loc="left")
ax2 = fig.add_subplot(gs[1, 0])
ax2.plot(np.degrees(phi), np.degrees(th), color=C["arm"], lw=1.6)
ax2.fill_between(np.degrees(phi), -15, 15, where=cut_down, color=C["acc"], alpha=0.15, lw=0)
ax2.set_xlabel("crank angle φ (°)"); ax2.set_ylabel("arm angle θ (°)"); ax2.set_xlim(0, 360); ax2.set_ylim(-15, 15)
ax2.text(0.03, 0.9, "shaded: cutting window\n(block face passing edge)", transform=ax2.transAxes, fontsize=6.5, va="top", color=C["acc"])
ax2.set_title("(b) Rocker angle vs crank angle", fontsize=8, loc="left")
ax3 = fig.add_subplot(gs[1, 1])
ax3.plot(np.degrees(phi), v_face, color=C["crank"], lw=1.6)
ax3.fill_between(np.degrees(phi), -400, 400, where=cut_down, color=C["acc"], alpha=0.15, lw=0)
ax3.set_xlabel("crank angle φ (°)"); ax3.set_ylabel("face velocity v$_z$ (mm/s at 60 rpm)"); ax3.set_xlim(0, 360)
ax3.set_ylim(-220, 220)
ax3.text(0.03, 0.08, f"cut lasts {t_cut_ms:.0f} ms per rev\n|v| = {np.abs(v_face[cut_down]).min():.0f}–{np.abs(v_face[cut_down]).max():.0f} mm/s", transform=ax3.transAxes, fontsize=6.5)
ax3.set_title("(c) Cutting speed through the window", fontsize=8, loc="left")
fig.savefig("figM1_kinematics.png", dpi=220, bbox_inches="tight"); plt.close(fig)

# ================= M2: blade / block geometry =================
fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.6), gridspec_kw=dict(width_ratios=[1.15, 1]))
ax = axs[0]; ax.set_aspect("equal")
lean = np.radians(27.0); wedge = np.radians(35.0)
# blade: thin plate leaning 27 deg from the face plane (vertical), edge at origin, body going down-right (+x, -z)
E = np.array([0.0, 0.0]); d = np.array([np.sin(lean), -np.cos(lean)])  # along blade from edge
n = np.array([np.cos(lean), np.sin(lean)])
Lb, tb = 9.0, 0.6   # exaggerated thickness for legibility
# facet: edge to full thickness over a length tb/tan(wedge/2)... draw symmetric wedge then body
fl = tb/2/np.tan(wedge/2)
poly = [E, E + fl*d + (tb/2)*n, E + Lb*d + (tb/2)*n, E + Lb*d - (tb/2)*n, E + fl*d - (tb/2)*n]
ax.add_patch(Polygon(poly, closed=True, fc="0.75", ec="0.1", lw=0.8))
# block: face vertical at x=0 (slightly left), occupying x<0, with section ribbon above
ax.add_patch(Rectangle((-6, -8), 6, 8, fc=C["block"], ec="0.3", lw=0.7))
ax.add_patch(Rectangle((-6, 0), 6, 0.35, fc="#c9a86a", ec="none"))
ax.annotate("section ribbon (5 µm,\nexaggerated)", (-3, 0.35), (-6.2, 1.6), fontsize=6.5, arrowprops=dict(arrowstyle="-", lw=0.6))
# arrows: block motion downward along face
ax.annotate("", xy=(-3, -6.5), xytext=(-3, -3.5), arrowprops=dict(arrowstyle="->", lw=1.2, color=C["arm"]))
ax.text(-3.5, -5, "block\nmotion", fontsize=6.5, color=C["arm"], ha="right", va="center")
# pressure strip: 1 mm steel under blade, top 2.5 mm below edge measured along blade
ps0 = E + 2.5*d - (tb/2)*n
ax.add_patch(Polygon([ps0, ps0 + 7*d, ps0 + 7*d - 1.0*n, ps0 - 1.0*n], closed=True, fc="0.45", ec="0.1", lw=0.6))
ax.text(*(ps0 + 4*d - 1.9*n), "steel pressure strip 1.0 mm", fontsize=6, rotation=-63, ha="center", va="center", color="0.2")
ax.text(*(E + 6*d + 1.0*n), "blade 0.25 mm,\n35° facet", fontsize=6.5, rotation=-63, ha="center", va="center")
# angle arcs
ax.add_patch(Arc((0, 0), 6, 6, angle=0, theta1=270, theta2=270+27, color=C["acc"], lw=0.9))
ax.annotate("lean 27°", xy=(1.4, -2.7), xytext=(-5.7, -2.6), fontsize=7, color=C["acc"], arrowprops=dict(arrowstyle="-", lw=0.6, color=C["acc"]))
# clearance = angle between lower facet and face: lower facet direction = d rotated by -wedge/2 (toward +x? ) draw
ax.add_patch(Arc((0, 0), 4, 4, angle=0, theta1=270, theta2=270+np.degrees(lean - wedge/2), color=C["stage"], lw=0.9))
ax.annotate("clearance α = 9.5°\n(27° − 35°/2)", xy=(1.0, -1.9), xytext=(-5.7, -0.9), fontsize=6.5, color=C["stage"], arrowprops=dict(arrowstyle="-", lw=0.6, color=C["stage"]))
ax.add_patch(Arc((0, 0), 5, 5, angle=0, theta1=90+np.degrees(0), theta2=180 - np.degrees(lean + wedge/2) + 90 if False else 90, color="0.3", lw=0.0))
# rake: angle between upper facet and the normal to the face (horizontal, pointing +x): upper facet direction angle = -(90-27)+17.5 = -45.5 deg from +x
ax.add_patch(Arc((0, 0), 7, 7, angle=0, theta1=-45.5, theta2=0, color="#8e44ad", lw=0.9))
ax.plot([0, 4.5], [0, 0], color="0.4", lw=0.6, ls=":")
ax.text(3.2, -0.9, "rake γ = 45.5°", fontsize=6.5, color="#8e44ad")
ax.text(0.50, 0.97, "wedge β = 35°;  α + β/2 + γ = 90°\nclearance varies 4.9°–14.1° over the\n±4.6° arm angle during the cut", transform=ax.transAxes, fontsize=6.5, va="top")
ax.set_xlim(-6.5, 8.5); ax.set_ylim(-9, 4.5); ax.set_xlabel("feed axis x (mm)"); ax.set_ylabel("z (mm)")
ax.set_title("(a) Edge geometry at the cut (y = 0 section)", fontsize=8, loc="left")
ax = axs[1]
# face arc vs flat: block face points on radius 125; sagitta over 20 mm
s = np.linspace(-HALF_BLOCK, HALF_BLOCK, 100); sag = (L_ARM - np.sqrt(L_ARM**2 - s**2))*1000  # µm
ax.plot(sag, s, color=C["acc"], lw=1.6); ax.axvline(0, color="0.5", lw=0.7, ls="--")
ax.fill_betweenx(s, 0, sag, color=C["block"], alpha=0.6)
ax.set_xlabel("face deviation from flat (µm)"); ax.set_ylabel("height on block (mm)")
ax.text(120, -2, f"sagitta {sag.max():.0f} µm over 20 mm:\nthe arc is trimmed into the block\nby the arm itself (self-consistent),\nso section thickness is unaffected", fontsize=6.5, va="top")
ax.set_title("(b) Face trimmed on the 125 mm arc", fontsize=8, loc="left")
fig.savefig("figM2_blade_geometry.png", dpi=220, bbox_inches="tight"); plt.close(fig)

# ================= M3: differential screw + ratchet =================
fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.4), gridspec_kw=dict(width_ratios=[1.3, 1]))
ax = axs[0]; ax.set_aspect("equal")
# frame nut (fixed) on right, stage tab nut on left, screw through both
ax.add_patch(Rectangle((40, -14), 14, 28, fc="0.85", ec="0.3", hatch="///", lw=0.7)); ax.text(47, 17, "frame upstand\n(fixed)", ha="center", fontsize=6.5)
ax.add_patch(Rectangle((8, -14), 10, 28, fc="#d9f0dc", ec=C["stage"], lw=0.7)); ax.text(13, 17, "stage tab\n(moves)", ha="center", fontsize=6.5, color=C["stage"])
ax.add_patch(Rectangle((-10, -4), 62, 8, fc="0.6", ec="0.2", lw=0.7))
# threads: pitch 1.0 on left section, 1.25 on right
for x in np.arange(-9, 20, 1.0): ax.plot([x, x+0.5], [-4, 4], color="0.25", lw=0.5)
for x in np.arange(24, 52, 1.25): ax.plot([x, x+0.6], [-4, 4], color="0.25", lw=0.5)
ax.add_patch(Rectangle((19, -6), 5, 12, fc="0.4", ec="0.2", lw=0.6)); ax.text(21.5, -9, "coupler\n+ jam nut", ha="center", fontsize=6, va="top")
ax.add_patch(Rectangle((52, -9), 8, 18, fc=C["acc"], ec="0.2", lw=0.6)); ax.text(61, -16, "ratchet drum\n50 T × 14 mm", ha="center", fontsize=6.5, color=C["acc"], va="top")
ax.text(4, -7, "M8 × 1.0", fontsize=6.5, ha="center", va="top"); ax.text(38, -7, "M8 × 1.25", fontsize=6.5, ha="center", va="top")
ax.annotate("", xy=(50, -20), xytext=(38, -20), arrowprops=dict(arrowstyle="->", lw=1, color="0.2")); ax.text(44, -22, "screw advances 1.25 mm/rev\nrelative to the frame", fontsize=6.5, ha="center", va="top")
ax.annotate("", xy=(-6, -20), xytext=(6, -20), arrowprops=dict(arrowstyle="->", lw=1, color=C["stage"])); ax.text(0, -22, "stage nut retreats 1.0 mm/rev\nalong the screw", fontsize=6.5, ha="center", va="top", color=C["stage"])
ax.text(26, 26, "net stage feed = (1.25 − 1.00) mm/rev = 0.25 mm/rev\n÷ 50 teeth = 5.0 µm per ratchet tooth", fontsize=7.5, ha="center", weight="bold")
ax.set_xlim(-12, 72); ax.set_ylim(-34, 32); ax.axis("off")
ax.set_title("(a) Differential screw: two standard threads give 5 µm", fontsize=8, loc="left")
ax = axs[1]
teeth = np.arange(1, 7); ax.bar(teeth, 5*teeth, color=C["acc"], width=0.6)
for t_, v in zip(teeth, 5*teeth): ax.text(t_, v+0.8, f"{v} µm", ha="center", fontsize=7)
ax.set_xlabel("teeth released per stroke (selector)"); ax.set_ylabel("section thickness (µm)"); ax.set_ylim(0, 46)
ax.text(0.03, 0.97, "screw travels ±5 mm for ±1 mm stage travel\n(5:1 reduction) → 400 sections at 5 µm per reset;\npitch-error sensitivity: δt/t = (δp₁ − δp₂)/0.25 mm\n→ 0.37 µm (1 s.d.) in the budget", transform=ax.transAxes, fontsize=6.5, va="top")
ax.set_title("(b) Thickness selector", fontsize=8, loc="left")
fig.savefig("figM3_differential_screw.png", dpi=220, bbox_inches="tight"); plt.close(fig)

# ================= M4: cross-strip pivot + metrology loop =================
fig, axs = plt.subplots(1, 2, figsize=(7.6, 4.2), gridspec_kw=dict(width_ratios=[1, 1.1], wspace=0.55))
ax = axs[0]; ax.set_aspect("equal")
Lf = 20.0
for ang, w in [(45, 20), (-45, 10)]:
    a = np.radians(ang); u = np.array([np.cos(a), np.sin(a)]); n = np.array([-u[1], u[0]])
    strip = [ -Lf/2*u - 0.6*n, Lf/2*u - 0.6*n, Lf/2*u + 0.6*n, -Lf/2*u + 0.6*n]
    ax.add_patch(Polygon(strip, closed=True, fc="0.35" if w == 20 else "0.55", ec="0.1", lw=0.6))
    # clamp pads
    for sgn in (-1, 1):
        c0 = sgn*(Lf/2 + 3)*u
        ax.add_patch(Polygon([c0 - 3*u - 3*n, c0 + 3*u - 3*n, c0 + 3*u + 3*n, c0 - 3*u + 3*n], closed=True, fc=C["frame"] if sgn*np.cos(a) < 0 else C["arm"], ec="0.2", lw=0.6))
ax.plot(0, 0, "o", ms=5, color=C["acc"]); ax.text(1.0, -1.0, "instant centre\n(strip crossing)", fontsize=6.5, color=C["acc"], va="top")
ax.text(-19, 12, "frame side", fontsize=7, color="0.3"); ax.text(10, -17, "arm side", fontsize=7, color=C["arm"])
ax.text(0, -21, "3 strips (10/20/10 mm wide) × 0.20 mm,\nfree length L = 20 mm, crossing at 90°;\nmiddle strip opposed by the two outer\nstrips → symmetric, no out-of-plane twist", fontsize=6.3, ha="center", va="top")
ax.text(-20, 41, r"$k_\theta=\sum\frac{E w t^3}{12L}$ = %.2f N·m/rad" % dr["pivot"]["k_theta"] + "\n" +
        r"$\sigma_{max}=\frac{E t \theta}{2L}$ = %.0f MPa at ±12.5°" % (dr["pivot"]["sigma_max"]/1e6) + "\n" +
        r"$T$ = %.3f N·m" % dr["pivot"]["torque_max"] + "\n" +
        r"centre shift ≈ %.0f µm (repeatable)" % (dr["pivot"]["center_shift"]*1e6) + "\n" +
        r"axial $k$ = %.0f N/µm" % (dr["pivot"]["k_axial"]/1e6) + "\n" +
        r"yield SF %.1f, fatigue SF %.1f" % (dr["pivot"]["safety_yield"], dr["pivot"]["safety_fatigue"]),
        fontsize=6.3, va="top")
ax.set_xlim(-21, 21); ax.set_ylim(-36, 42); ax.axis("off")
ax.set_title("(a) Cross-strip flexure pivot", fontsize=8, loc="left")
ax = axs[1]
lc = dr["loop_compliance_um_per_N"]
order = [("blade overhang (2.5 mm)", "blade_overhang"), ("pressure strip / holder", "blade_clamp_compression"), ("stage pitch (FE)", "stage_pitch"),
         ("screw tip contact", "ball_contact"), ("screw axial", "screw_axial"), ("nut bedding in PLA", "thread_nut_in_PLA"),
         ("frame: upstand→base→tower (FE)", "frame_tower_base_FE"), ("flexure pivot axial", "pivot_axial"), ("arm axial", "arm_axial"), ("cassette clamp", "cassette_clamp")]
vals = [lc[k] for _, k in order]; labels = [l for l, _ in order]
y = np.arange(len(order))[::-1]
cols = [C["stage"]]*6 + [C["frame"], C["steel"], C["arm"], C["arm"]]
ax.barh(y, vals, color=cols, height=0.65)
for yi, v in zip(y, vals): ax.text(v + 0.004, yi, f"{v:.3f}", va="center", fontsize=6.5)
ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=6.8); ax.set_xlabel("feed-axis compliance (µm/N)"); ax.set_xlim(0, 0.36)
tot = sum(lc.values())
ax.text(0.98, 0.05, f"series sum {tot:.2f} µm/N → stiffness {1/tot:.2f} N/µm\nforce path: blade edge → stage → screw\n→ frame → pivot → arm → block face", transform=ax.transAxes, ha="right", fontsize=6.3)
ax.set_title("(b) Loop compliance by element", fontsize=8, loc="left")
fig.savefig("figM4_pivot_loop.png", dpi=220, bbox_inches="tight"); plt.close(fig)

# ================= M5: feed timing over one revolution =================
fig, axs = plt.subplots(3, 1, figsize=(7.2, 4.6), sharex=True, gridspec_kw=dict(hspace=0.12))
deg = np.degrees(phi)
axs[0].plot(deg, np.degrees(th), color=C["arm"], lw=1.6); axs[0].set_ylabel("arm angle θ (°)")
axs[0].axhline(6, color=C["acc"], lw=0.7, ls="--"); axs[0].text(2, 6.6, "striker meets pawl rod at θ ≈ 6° (block 13 mm clear of edge)", fontsize=6.3, color=C["acc"])
for a in axs: a.fill_between(deg, -1e3, 1e3, where=cut_down, color=C["acc"], alpha=0.12, lw=0)
axs[0].set_ylim(-15, 15)
# pawl lift: proportional to (theta - 6deg) * lever, capped 22.3 mm
lift = np.clip((np.degrees(th) - 6.0)/(np.degrees(th_max) - 6.0), 0, 1)*22.3
axs[1].plot(deg, lift, color=C["steel"], lw=1.6); axs[1].set_ylabel("pawl lift (mm)"); axs[1].set_ylim(-2, 26)
axs[1].text(2, 20, "pawl finger slides back over N teeth (N set by selector shield, 1–6)", fontsize=6.3)
# drum index: step at the moment the pawl drops (lift returns to 0) after the top
drop = np.where((lift[:-1] > 0) & (lift[1:] == 0))[0]
feed = np.zeros_like(deg); 
if len(drop): feed[drop[0]+1:] = 5.0
axs[2].step(deg, feed, where="post", color=C["stage"], lw=1.8); axs[2].set_ylabel("stage feed (µm)"); axs[2].set_ylim(-1, 7.5)
axs[2].text(2, 5.8, "pawl re-engages → drum turns 1 tooth (7.2°) → screw 0.25/50 mm → +5.0 µm blade advance, before the next cut", fontsize=6.3)
axs[2].set_xlabel("crank angle φ (°); one handwheel turn = one section"); axs[2].set_xlim(0, 360)
axs[0].set_title("Feed sequence per revolution (shaded: cutting window)", fontsize=8, loc="left")
fig.savefig("figM5_feed_timing.png", dpi=220, bbox_inches="tight"); plt.close(fig)
json.dump(dict(theta_max_deg=float(np.degrees(th_max)), theta_min_deg=float(np.degrees(th_min)), stroke_mm=float(stroke), t_cut_ms_60rpm=float(t_cut_ms),
               v_cut_min=float(np.abs(v_face[cut_down]).min()), v_cut_max=float(np.abs(v_face[cut_down]).max()), sagitta_um=float(sag.max()),
               cut_window_theta_deg=float(np.degrees(np.arcsin(HALF_BLOCK/L_ARM)))), open("mechanics_results.json", "w"), indent=1)
print("mechanics figs done")
