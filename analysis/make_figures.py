import numpy as np, matplotlib as mpl, matplotlib.pyplot as plt, json
import flextome_design as fd
apply_figure_style(sizes=(9, 8, 7))
d = fd.Design()
C_FOCAL = "#c0392b"; C_A = "#2c6fbb"; C_B = "#7f8c8d"; C_C = "#e67e22"

# ---------------------------------------------------------------- Fig 2: kinematics & pivot
phi, th, z, xr = fd.crank_rocker(d)
rpm = 60.0; omega = 2*np.pi*rpm/60
dz = np.gradient(z, phi)*omega
fig, axs = plt.subplots(2, 2, figsize=(8.0, 6.0))
ax = axs[0, 0]
ax.plot(np.degrees(phi), z*1e3, color=C_A, label="vertical position z")
ax.set_xlabel("crank angle φ (°)"); ax.set_ylabel("block-face position (mm)")
ax2 = ax.twinx(); ax2.plot(np.degrees(phi), -xr*1e3, color=C_C, ls="--", label="retraction from blade")
ax2.set_ylabel("retraction from blade (mm)", color=C_C); ax2.tick_params(axis="y", colors=C_C)
ax.set_title(f"Arc stroke: {(z.max()-z.min())*1e3:.0f} mm travel, {xr.max()*1e3:.1f} mm end retraction")
ax.axhspan(-10, 10, color="0.9", zorder=0); ax.text(300, -8.5, "cut zone (20 mm face)", fontsize=7, color="0.3")
ax = axs[0, 1]
ax.plot(np.degrees(phi), dz*1e3, color=C_A)
ax.set_xlabel("crank angle φ (°)"); ax.set_ylabel("cutting speed (mm/s) at 60 rpm")
ax.axhline(0, color="0.6", lw=0.6)
ax.set_title(f"Cutting speed peaks at {abs(dz).max()*1e3:.0f} mm/s (60 rpm)")
ax = axs[1, 0]
Ls = [0.100, 0.125, 0.150]
zz = np.linspace(-10, 10, 200)
for L, c in zip(Ls, [C_B, C_FOCAL, C_B]):
    ang = d.clearance_nom_deg - np.degrees(np.arcsin(zz*1e-3/L))
    ax.plot(zz, ang, color=c, lw=2 if L == 0.125 else 1, label=f"arm {L*1e3:.0f} mm")
ax.axhspan(3, 15, color="#dff0d8", zorder=0)
ax.text(-9.5, 3.4, "recommended clearance band 3–15°", fontsize=7, color="0.3")
ax.set_xlabel("position across block face (mm)"); ax.set_ylabel("effective clearance angle (°)")
ax.set_title("27° blade lean: clearance 4.9–14.1° across the face")
ax.legend(frameon=False, loc="lower left", fontsize=7, bbox_to_anchor=(0.0, 0.08))
ax = axs[1, 1]
ts = np.linspace(0.1, 0.4, 100)*1e-3
thm = np.arcsin(d.r_crank/d.r_con)
for Lp, c in zip([0.015, 0.020, 0.025], [C_B, C_FOCAL, C_B]):
    sig = 200e9*ts*thm/(2*Lp)/1e6
    ax.plot(ts*1e3, sig, color=c, lw=2 if Lp == 0.020 else 1, label=f"free length {Lp*1e3:.0f} mm")
ax.axhline(450, color=C_C, ls="--", lw=1); ax.text(0.21, 462, "spring-steel endurance limit ~450 MPa", fontsize=7, color=C_C)
ax.plot([0.2], [200e9*0.2e-3*thm/(2*0.020)/1e6], "o", color=C_FOCAL, ms=6)
ax.set_xlabel("strip thickness (mm)"); ax.set_ylabel("peak strip stress at θ = 12.5° (MPa)")
ax.set_title("0.20 mm strips: 218 MPa, below endurance limit")
ax.legend(frameon=False, loc="lower right", fontsize=7)
for a, l in zip(axs.flat, "abcd"): panel_letter(a, l)
fig.tight_layout(); fig.savefig("fig2_kinematics_pivot.png", dpi=300); plt.close(fig)

# ---------------------------------------------------------------- Fig 3: error budget
s = fd.summary(d); lc = s["loop_compliance_um_per_N"]
items = sorted(lc.items(), key=lambda kv: kv[1])
names = {"frame_tower_base_FE": "frame: tower + base on two end pads (FE)", "blade_overhang": "blade overhang 2.5 mm", "upstand_FE": "frame-nut upstand, gusseted (FE)",
         "stage_pitch": "stage pitch about the tab nut (FE)", "tab_nut_bedding": "tab nut on washer", "frame_nut_bedding": "frame nut on washer",
         "coupler_steel_contact": "rod end on bolt head", "arm_axial": "arm axial", "pivot_axial": "pivot strips axial", "cassette_clamp": "cassette clamp",
         "screw_axial": "feed screw axial", "blade_backwall_bearing": "blade on holder back wall"}
fig = plt.figure(figsize=(7.2, 6.0))
gs = fig.add_gridspec(2, 2, height_ratios=[1.1, 1])
ax = fig.add_subplot(gs[0, :])
ax.barh([names[k] for k, v in items], [v for k, v in items], color=[C_FOCAL if k in ("frame_tower_base_FE","blade_overhang","upstand_FE") else C_A for k, v in items])
tot = sum(lc.values())
ax.set_xlabel("feed-axis compliance (µm/N)")
ax.set_title(f"Loop compliance budget: {tot:.2f} µm/N total (stiffness {1/tot:.2f} N/µm)")
for i, (k, v) in enumerate(items): ax.text(v+0.005, i, f"{v:.2f}", va="center", fontsize=7)
ax.margins(x=0.12)
ax = fig.add_subplot(gs[1, 0])
for k, col in ((1, C_A), (2, C_C)):
    mc = fd.monte_carlo(d, teeth_per_stroke=k)
    ax.hist(mc["samples"]*1e6, bins=60, color=col, alpha=0.75, label=f"{5*k} µm setting: {mc['mean']*1e6:.2f} ± {mc['std']*1e6:.2f} µm (1 s.d.)")
ax.set_xlabel("predicted section thickness (µm)"); ax.set_ylabel("count (n = 20 000 strokes)")
ax.set_title("Monte-Carlo section thickness")
ax.legend(frameon=False, fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.28)); ax.set_ylim(0, 1300)
ax = fig.add_subplot(gs[1, 1])
mc = fd.monte_carlo(d)
contrib = mc["contrib"]; labs = {"compliance": "compliance x cut-force change (√2·σF)", "thread_lead": "thread lead (periodic + local)", "thread_bias": "thread lead bias (per screw pair)",
                                "pawl": "pawl seating", "ratchet_tooth": "printed ratchet tooth pitch", "thermal": "thermal drift (1 s)"}
ks = sorted(contrib, key=lambda k: contrib[k])
ax.barh([labs[k] for k in ks], [contrib[k]*1e6 for k in ks], color=C_A)
for i, k in enumerate(ks): ax.text(contrib[k]*1e6+0.01, i, f"{contrib[k]*1e6:.2f}", va="center", fontsize=7)
ax.set_xlabel("1 s.d. contribution to thickness (µm)"); ax.set_title("Force change through the loop dominates"); ax.margins(x=0.2)
for a, l in zip(fig.axes, "abc"): panel_letter(a, l)
fig.tight_layout(); fig.savefig("fig3_error_budget.png", dpi=300); plt.close(fig)

# ---------------------------------------------------------------- Fig 4: stage flexure design map + stiffness sensitivity
fig, axs = plt.subplots(1, 2, figsize=(7.6, 3.5), gridspec_kw=dict(width_ratios=[1, 1.15]))
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
ax.plot([d.st_t*1e3], [d.st_h*1e3], "o", color="k", ms=6); ax.text(d.st_t*1e3+0.03, d.st_h*1e3+0.8, "chosen: 1.0 × 44 mm", fontsize=7)
ax.set_xlabel("flexure plate thickness (mm)"); ax.set_ylabel("plate free height (mm)")
ax.set_title("Stage plates below 15 MPa PLA creep line")
ax = axs[1]
Fv = np.linspace(0.1, 1.0, 50)
def sd_curve(dd):
    Ct = sum(fd.loop_compliance(dd).values()); mc = fd.monte_carlo(dd)
    other = np.sqrt(mc["std"]**2 - (np.sqrt(2)*Ct*dd.F_cut_var)**2)
    return np.sqrt((np.sqrt(2)*Ct*Fv)**2 + other**2)*1e6, Ct*1e6
for dd, col, lw, lab in ((fd.Design(), C_FOCAL, 2.0, "rev. 2 as designed ({C:.2f} µm/N)"),
                         (fd.Design(nut_seat="insert", blade_support_overhang=2.0e-3, frame_fe_compliance=0.2107e-6), C_A, 1.2, "rev. 2 + brass inserts, 2.0 mm\noverhang, full-contact base ({C:.2f} µm/N)")):
    y, Ct = sd_curve(dd); ax.plot(Fv, y, color=col, lw=lw, label=lab.format(C=Ct))
ax.axvline(d.F_cut_var, color="0.5", ls=":", lw=0.8); ax.text(d.F_cut_var+0.015, 1.05, "assumed 0.4 N", fontsize=7, color="0.4", rotation=90, va="bottom")
ax.text(0.98, 0.04, "rev. 1 load path (7.3 µm/N upstand):\n4.8 µm at 0.4 N, off scale", transform=ax.transAxes, ha="right", fontsize=6.5, color="0.35")
ax.set_xlabel("per-stroke cutting-force variation (N, 1 s.d.)"); ax.set_ylabel("section-thickness scatter (µm, 1 s.d.)")
ax.set_title("Scatter = force variation × loop compliance"); ax.set_ylim(0, 1.8)
ax.legend(frameon=False, fontsize=6.5, loc="upper left")
for a, l in zip(axs, "ab"): panel_letter(a, l)
fig.tight_layout(); fig.savefig("fig4_stage_design_map.png", dpi=300); plt.close(fig)

json.dump(dict(summary={k: (v if not hasattr(v, "tolist") else v.tolist()) for k, v in s.items() if k not in ("pivot", "stage", "loop_compliance_um_per_N")},
               pivot={k: float(v) for k, v in s["pivot"].items()}, stage={k: float(v) for k, v in s["stage"].items()},
               loop_compliance_um_per_N={k: float(v) for k, v in lc.items()},
               mc5={k: float(v) for k, v in fd.monte_carlo(d).items() if k not in ("samples", "contrib")},
               mc5_contrib={k: float(v) for k, v in fd.monte_carlo(d)["contrib"].items()}),
          open("design_results.json", "w"), indent=1, default=float)
print("figures done")
