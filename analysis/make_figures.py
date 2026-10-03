import numpy as np, matplotlib as mpl, matplotlib.pyplot as plt, json
import flextome_design as fd
apply_figure_style(sizes=(9, 8, 7))
d = fd.Design()
C_FOCAL = "#c0392b"; C_A = "#2c6fbb"; C_B = "#7f8c8d"; C_C = "#e67e22"

# ---------------------------------------------------------------- Fig 2: kinematics & pivot
phi, th, z, xr = fd.crank_rocker(d)
rpm = 60.0; omega = 2*np.pi*rpm/60
dz = np.gradient(z, phi)*omega
fig, axs = plt.subplots(2, 2, figsize=(7.6, 5.6))
ax = axs[0, 0]
ax.plot(np.degrees(phi), z*1e3, color=C_A, label="vertical position z")
ax.set_xlabel("crank angle (deg)"); ax.set_ylabel("block-face position (mm)")
ax2 = ax.twinx(); ax2.plot(np.degrees(phi), -xr*1e3, color=C_C, ls="--", label="retraction from blade")
ax2.set_ylabel("retraction from blade (mm)", color=C_C); ax2.tick_params(axis="y", colors=C_C)
ax.set_title(f"Arc stroke: {(z.max()-z.min())*1e3:.0f} mm travel, {xr.max()*1e3:.1f} mm end retraction")
ax.axhspan(-10, 10, color="0.9", zorder=0); ax.text(300, -8.5, "cut zone (20 mm face)", fontsize=7, color="0.3")
ax = axs[0, 1]
ax.plot(np.degrees(phi), dz*1e3, color=C_A)
ax.set_xlabel("crank angle (deg)"); ax.set_ylabel("cutting speed (mm/s) at 60 rpm")
ax.axhline(0, color="0.6", lw=0.6)
ax.set_title(f"Cutting speed peaks at {abs(dz).max()*1e3:.0f} mm/s (60 rpm)")
ax = axs[1, 0]
Ls = [0.100, 0.125, 0.150]
zz = np.linspace(-10, 10, 200)
for L, c in zip(Ls, [C_B, C_FOCAL, C_B]):
    ang = d.clearance_nom_deg - np.degrees(np.arcsin(zz*1e-3/L))
    ax.plot(zz, ang, color=c, lw=2 if L == 0.125 else 1, label=f"arm {L*1e3:.0f} mm")
ax.axhspan(3, 15, color="#dff0d8", zorder=0)
ax.text(-9.5, 3.4, "recommended clearance band 3-15 deg", fontsize=7, color="0.3")
ax.set_xlabel("position across block face (mm)"); ax.set_ylabel("effective clearance angle (deg)")
ax.set_title("27 deg blade lean keeps clearance 4.9-14.1 deg")
ax.legend(frameon=False, loc="lower left", fontsize=7, bbox_to_anchor=(0.0, 0.08))
ax = axs[1, 1]
ts = np.linspace(0.1, 0.4, 100)*1e-3
thm = np.arcsin(d.r_crank/d.r_con)
for Lp, c in zip([0.015, 0.020, 0.025], [C_B, C_FOCAL, C_B]):
    sig = 200e9*ts*thm/(2*Lp)/1e6
    ax.plot(ts*1e3, sig, color=c, lw=2 if Lp == 0.020 else 1, label=f"free length {Lp*1e3:.0f} mm")
ax.axhline(450, color=C_C, ls="--", lw=1); ax.text(0.21, 462, "spring-steel endurance limit ~450 MPa", fontsize=7, color=C_C)
ax.plot([0.2], [200e9*0.2e-3*thm/(2*0.020)/1e6], "o", color=C_FOCAL, ms=6)
ax.set_xlabel("strip thickness (mm)"); ax.set_ylabel("peak bending stress at 12.5 deg (MPa)")
ax.set_title("0.20 mm strips: 218 MPa, below endurance limit")
ax.legend(frameon=False, loc="lower right", fontsize=7)
for a, l in zip(axs.flat, "abcd"): panel_letter(a, l)
fig.tight_layout(); fig.savefig("fig2_kinematics_pivot.png", dpi=300); plt.close(fig)

# ---------------------------------------------------------------- Fig 3: error budget
s = fd.summary(d); lc = s["loop_compliance_um_per_N"]
items = sorted(lc.items(), key=lambda kv: kv[1])
names = {"frame_tower_base_FE": "frame: gusseted tower + base (FE)", "base_rotation": "base rotation", "blade_overhang": "blade overhang 2.5 mm",
         "stage_pitch": "feed-stage pitch coupling", "thread_nut_in_PLA": "nut bedding in PLA", "tower_bending": "tower bending",
         "ball_contact": "screw-tip contact", "arm_axial": "arm axial", "pivot_axial": "pivot strips axial", "cassette_clamp": "cassette clamp",
         "screw_axial": "feed screw axial", "blade_clamp_compression": "blade holder compression"}
fig = plt.figure(figsize=(7.2, 6.0))
gs = fig.add_gridspec(2, 2, height_ratios=[1.1, 1])
ax = fig.add_subplot(gs[0, :])
ax.barh([names[k] for k, v in items], [v for k, v in items], color=[C_FOCAL if k in ("frame_tower_base_FE","blade_overhang") else C_A for k, v in items])
tot = sum(lc.values())
ax.set_xlabel("compliance along feed axis (um/N)")
ax.set_title(f"Metrology-loop compliance budget: total {tot:.2f} um/N  (stiffness {1/tot:.2f} N/um)")
for i, (k, v) in enumerate(items): ax.text(v+0.005, i, f"{v:.2f}", va="center", fontsize=7)
ax.margins(x=0.12)
ax = fig.add_subplot(gs[1, 0])
for k, col in ((1, C_A), (2, C_C)):
    mc = fd.monte_carlo(d, teeth_per_stroke=k)
    ax.hist(mc["samples"]*1e6, bins=60, color=col, alpha=0.75, label=f"{5*k} um setting: {mc['mean']*1e6:.2f} +/- {mc['std']*1e6:.2f} um (1 sd)")
ax.set_xlabel("predicted section thickness (um)"); ax.set_ylabel("count (n = 20 000 strokes)")
ax.set_title("Monte-Carlo section thickness, all error sources")
ax.legend(frameon=False, fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.28)); ax.set_ylim(0, 1300)
ax = fig.add_subplot(gs[1, 1])
mc = fd.monte_carlo(d)
contrib = mc["contrib"]; labs = {"compliance": "compliance x cut-force variation", "thread_pitch": "differential thread pitch", "pawl": "pawl seating",
                                "ratchet_tooth": "printed ratchet tooth pitch", "thermal": "thermal drift (1 s)"}
ks = sorted(contrib, key=lambda k: contrib[k])
ax.barh([labs[k] for k in ks], [contrib[k]*1e6 for k in ks], color=C_A)
for i, k in enumerate(ks): ax.text(contrib[k]*1e6+0.01, i, f"{contrib[k]*1e6:.2f}", va="center", fontsize=7)
ax.set_xlabel("1-sd contribution to thickness (um)"); ax.set_title("Compliance and thread pitch dominate"); ax.margins(x=0.2)
for a, l in zip(fig.axes, "abc"): panel_letter(a, l)
fig.tight_layout(); fig.savefig("fig3_error_budget.png", dpi=300); plt.close(fig)

# ---------------------------------------------------------------- Fig 4: stage flexure design map + stiffness sensitivity
fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.3))
ax = axs[0]
T = np.linspace(0.6, 1.6, 60)*1e-3; H = np.linspace(0.025, 0.060, 60)
TT, HH = np.meshgrid(T, H)
E = fd.MAT["PLA"]["E"]
sig = 3*E*TT*d.st_travel/HH**2/1e6
kx = 2*12*E*(d.st_w*TT**3/12)/HH**3/1e3   # N/mm
cs = ax.contourf(TT*1e3, HH*1e3, sig, levels=[0, 5, 10, 15, 20, 30, 45], cmap="YlOrRd")
cb = fig.colorbar(cs, ax=ax); cb.set_label(f"peak plate stress at +/-{d.st_travel*1e3:.0f} mm travel (MPa)")
c2 = ax.contour(TT*1e3, HH*1e3, kx, levels=[2, 5, 10, 20, 40], colors="k", linewidths=0.7)
ax.clabel(c2, fmt="%g N/mm", fontsize=6)
ax.contour(TT*1e3, HH*1e3, sig, levels=[15], colors=[C_FOCAL], linewidths=1.5)
ax.plot([d.st_t*1e3], [d.st_h*1e3], "o", color="k", ms=6); ax.text(d.st_t*1e3+0.03, d.st_h*1e3+0.8, "chosen: 1.0 x 44 mm", fontsize=7)
ax.set_xlabel("flexure plate thickness (mm)"); ax.set_ylabel("plate free height (mm)")
ax.set_title("Stage plates below 15 MPa PLA creep line")
ax = axs[1]
Fv = np.linspace(0.1, 1.0, 50)
for cf, col, lab in ((0.8315e-6, C_B, "slender 32 mm tower (FE 0.83 um/N)"), (0.2127e-6, C_FOCAL, "gusseted 60 mm-root tower (FE 0.21 um/N, chosen)")):
    dd = fd.Design(frame_fe_compliance=cf)
    Ct = sum(fd.loop_compliance(dd).values())
    std = np.sqrt((Ct*Fv)**2 + (0.371e-6)**2 + (0.1e-6)**2 + (0.066e-6)**2)
    ax.plot(Fv, std*1e6, color=col, lw=2 if cf < 0.5e-6 else 1.2, label=lab)
ax.axvline(d.F_cut_var, color="0.5", ls=":", lw=0.8); ax.text(d.F_cut_var+0.01, 0.45, "assumed 0.4 N (1 sd)", fontsize=7, color="0.4")
ax.set_xlabel("stroke-to-stroke cutting-force variation (N, 1 sd)"); ax.set_ylabel("section-thickness scatter (um, 1 sd)")
ax.set_title("Tower gusset cuts force-driven scatter by ~25%")
ax.legend(frameon=False, fontsize=7, loc="upper left")
for a, l in zip(axs, "ab"): panel_letter(a, l)
fig.tight_layout(); fig.savefig("fig4_stage_design_map.png", dpi=300); plt.close(fig)

json.dump(dict(summary={k: (v if not hasattr(v, "tolist") else v.tolist()) for k, v in s.items() if k not in ("pivot", "stage", "loop_compliance_um_per_N")},
               pivot={k: float(v) for k, v in s["pivot"].items()}, stage={k: float(v) for k, v in s["stage"].items()},
               loop_compliance_um_per_N={k: float(v) for k, v in lc.items()},
               mc5={k: float(v) for k, v in fd.monte_carlo(d).items() if k not in ("samples", "contrib")},
               mc5_contrib={k: float(v) for k, v in fd.monte_carlo(d)["contrib"].items()}),
          open("design_results.json", "w"), indent=1, default=float)
print("figures done")
