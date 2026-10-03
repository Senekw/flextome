"""Supplementary calculations for FlexTome (rev. 2), all previously reported without a script:
 (1) cutting-force bracket from the Atkins (2004) pressing-and-slicing wedge model,
 (2) first structural modes (lumped estimates) and the cutting excitation band,
 (3) entry transient of the loop under a ramped cutting force (2nd-order system),
 (4) improvement levers (nut seat, back-wall overhang, frame support) and the cutting-force sensitivity curve.
Run:  python extra_analysis.py  -> extra_results.json"""
import json, numpy as np
from scipy.signal import lsim, TransferFunction
import flextome_design as fd

# (1) cutting force: wedge of half-angle beta pressed into paraffin; R = fracture work (J/m2), sy = wax yield (Pa)
def F_cut(R, sy, t=5e-6, w=25e-3, beta=np.radians(17.5), mu=0.3):
    Fx_per_w = R*np.tan(beta) + sy*t*(np.tan(beta) + mu)      # feed-axis (normal to the face)
    Fz_per_w = R + sy*t*(1 + mu*np.tan(beta))                 # along the stroke
    return Fx_per_w*w, Fz_per_w*w
grid = [(R, sy) for R in (20, 50, 80) for sy in (1e6, 2e6, 3e6)]
vals = np.array([F_cut(R, sy) for R, sy in grid]); dull = np.array([F_cut(R*1.5, sy, mu=0.6) for R, sy in grid])

# (2) modes
d = fd.Design(); C0 = sum(fd.loop_compliance(d).values()); k_loop = 1/C0
E = 3.5e9; I_arm = ((35e-3)**4 - (29e-3)**4)/12; L = 0.125
k_arm_z = 3*E*I_arm/L**3; m_eff = 0.040 + 0.24*0.127
f_arm = np.sqrt(k_arm_z/m_eff)/(2*np.pi)
f_stage = np.sqrt((1/d.stage_fe_compliance)/0.08)/(2*np.pi)
f_loop = np.sqrt(k_loop/0.10)/(2*np.pi)
v_cut = 0.174; t_cut = 0.020/0.170

# (3) entry transient: force ramps 0 -> 2 N over a 0.5 mm entry chamfer at the cutting speed
m = 0.10; wn = np.sqrt(k_loop/m); zeta = 0.03
t = np.linspace(0, 0.12, 6000); F = 2.0*np.clip(v_cut*t/0.5e-3, 0, 1)
_, x, _ = lsim(TransferFunction([wn**2/k_loop], [1, 2*zeta*wn, wn**2]), F, t)

# (4) levers and sensitivity
def sd_for(dd, Fv=None):
    if Fv is not None: dd = fd.Design(**{**dd.__dict__, "F_cut_var": Fv})
    mc = fd.monte_carlo(dd); return float(mc["std"]*1e6), float(sum(fd.loop_compliance(dd).values())*1e6)
levers = {}
for name, dd in [("rev2_washer (baseline)", fd.Design()), ("rev2_bare_nuts", fd.Design(nut_seat="bare")), ("rev2_brass_inserts", fd.Design(nut_seat="insert")),
                 ("rev2_backwall_overhang_2.0", fd.Design(blade_support_overhang=2.0e-3)), ("rev2_full_contact_base", fd.Design(frame_fe_compliance=0.2107e-6)),
                 ("rev2_all_levers", fd.Design(nut_seat="insert", blade_support_overhang=2.0e-3, frame_fe_compliance=0.2107e-6)),
                 ("rev1_as_modelled (7.35 um/N upstand, bare nut, slender-stage FE)", fd.Design(upstand_fe_compliance=7.35e-6, nut_seat="bare", stage_fe_compliance=0.1086e-6, frame_fe_compliance=0.2127e-6))]:
    sd, C = sd_for(dd); levers[name] = dict(C_um_per_N=round(C, 3), sd_um=round(sd, 3))
sens = {str(Fv): round(sd_for(fd.Design(), Fv)[0], 3) for Fv in (0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0)}
out = dict(Fx_range_N=[float(vals[:, 0].min()), float(vals[:, 0].max())], Fz_range_N=[float(vals[:, 1].min()), float(vals[:, 1].max())], dull_Fx_max=float(dull[:, 0].max()),
           modes=dict(f_arm_z_Hz=float(f_arm), f_stage_x_Hz=float(f_stage), f_loop_x_Hz=float(f_loop), cut_duration_s=float(t_cut), excitation_below_Hz=float(1/t_cut)),
           transient=dict(qs_um=float(C0*2e6), peak_um=float(x.max()*1e6), overshoot_pct=float((x.max()/(C0*2)-1)*100), settle_ms=float(4/(zeta*wn)*1e3)),
           loop=dict(C_um_per_N=float(C0*1e6), k_N_per_um=float(1/(C0*1e6))), levers=levers, sd_vs_force_variation=sens)
json.dump(out, open("extra_results.json", "w"), indent=1)
print(json.dumps(out, indent=1))
