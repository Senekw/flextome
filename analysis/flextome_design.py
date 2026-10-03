"""
FlexTome design analysis: parameters, kinematics, flexure pivot, differential screw,
structural-loop compliance, thermal drift, Monte-Carlo section-thickness error budget.
All SI units unless noted.  Coordinate frame: x = feed (specimen -> blade, +x),
y = lateral (blade-edge direction), z = vertical.
"""
import numpy as np
from dataclasses import dataclass, field, asdict

# ----------------------------------------------------------------------------- materials
MAT = {
    # E [Pa], yield/UTS [Pa], density [kg/m3], CTE [1/K], fatigue endurance estimate [Pa]
    "PLA":   dict(E=3.5e9,  sy=50e6,  rho=1240, cte=68e-6,  se=12e6),
    "PETG":  dict(E=2.1e9,  sy=45e6,  rho=1270, cte=60e-6,  se=15e6),
    "steel_spring": dict(E=200e9, sy=1100e6, rho=7850, cte=11.5e-6, se=450e6),
    "steel": dict(E=200e9, sy=400e6, rho=7850, cte=11.5e-6, se=200e6),
}

@dataclass
class Design:
    # --- arm / rocking geometry
    L_arm: float = 0.125        # pivot -> block face [m]
    r_con: float = 0.060        # pivot -> connecting-rod joint on arm [m]
    r_crank: float = 0.013      # crank radius [m]
    L_conrod: float = 0.0712     # connecting rod length [m]
    z_pivot: float = 0.095      # pivot height above base-plate top [m]
    block_h: float = 0.020      # cut face height (z) [m]
    block_w: float = 0.025      # cut face width (y) [m]
    blade_overhang: float = 2.5e-3 # blade protrusion beyond pressure strip [m]
    blade_lean_deg: float = 27.0
    blade_wedge_deg: float = 35.0
    # --- cross-strip flexure pivot (spring steel feeler gauge blades)
    piv_t: float = 0.20e-3      # blade thickness [m]
    piv_L: float = 0.020        # free blade length [m]
    piv_w_total: float = 0.040  # total blade width summed over the 3 strips [m]
    piv_n: int = 3
    piv_mat: str = "steel_spring"
    # --- differential screw & ratchet
    pitch_a: float = 1.25e-3    # frame nut pitch (M8 coarse) [m]
    pitch_b: float = 1.00e-3    # stage nut pitch (M8 fine)  [m]
    n_teeth: int = 50
    # --- blade feed stage (monolithic parallelogram flexure, printed)
    st_t: float = 1.0e-3        # flexure plate thickness [m]
    st_h: float = 0.044         # flexure plate free height (bending length) [m]
    st_w: float = 0.060         # flexure plate width (y) [m]
    st_s: float = 0.043         # separation of the two plates along x [m]
    st_mat: str = "PLA"
    st_travel: float = 1.0e-3   # fine-feed range (+-) [m]; screw translates 5x this
    stage_fe_compliance: float = 0.1318e-6  # [m/N] nut-held stage at the blade, rev. 2 geometry, from fea_stage.py
    upstand_fe_compliance: float = 0.185e-6 # [m/N] gusseted frame-nut upstand under the screw reaction, from fea_upstand.py (rev. 1: 7.35e-6)
    frame_fe_compliance: float = 0.353e-6   # [m/N] pivot-to-deck, gusseted tower, base on two end pads (conservative), fea_frame.py; 0.211e-6 for full-contact base
    nut_seat: str = "washer"                # "bare" (nut face on PLA layer lines), "washer" (16 mm steel washer under each nut), "insert" (brass insert)
    blade_support_overhang: float = 2.5e-3  # blade cantilever beyond the holder back wall [m] (same as the strip side by default)
    # --- structural members (printed)
    frame_mat: str = "PLA"
    arm_sec: tuple = (0.040, 0.040, 0.003)   # outer w, h, wall thickness of box arm [m]
    tower_sec: tuple = (0.070, 0.070, 0.004) # tower box section
    base_depth: float = 0.055   # ribbed box base depth [m]
    tie_rods: bool = True       # 2 x M8 steel tie rods along the base bottom (steel spine)
    rod_offset: float = 0.015   # rod distance from base neutral axis [m]
    base_width: float = 0.150
    base_len: float = 0.230
    # --- cutting
    F_cut_nom: float = 2.0      # nominal normal (x) cutting force [N]
    F_cut_var: float = 0.4      # 1-sigma variation of the feed-axis cutting force PER STROKE [N]; thickness depends on the
                                # difference between consecutive strokes, so its s.d. is sqrt(2) x this (handled in monte_carlo)
    lead_bias_pct: float = 0.10 # cumulative lead error of each thread, uniform +-this % (ISO 6g rolled thread, ~0.025 mm over 25 mm)
    lead_periodic_um: float = 1.5  # per-revolution periodic (drunken-thread) error amplitude of each thread [um]
    lead_local_pct: float = 0.3 # random local pitch deviation of each thread, 1 s.d., % of pitch, independent per tooth step
    clearance_nom_deg: float = 9.5
    section_um: float = 5.0

    def feed_per_rev(self):
        return self.pitch_a - self.pitch_b

    def feed_per_tooth(self):
        return self.feed_per_rev() / self.n_teeth


# ----------------------------------------------------------------------------- kinematics
def crank_rocker(d: Design, n=721):
    """Position analysis of the crank-rocker driving the arm.
    Crank centre O2 at (x=d.r_con, z=z_pivot - h) below the arm joint radius; vertical conrod at mid-stroke.
    Returns arrays: crank angle phi [rad], arm angle theta [rad], tip z [m], tip x-retraction [m], clearance angle at the blade [deg]."""
    phi = np.linspace(0, 2*np.pi, n)
    # Simplified: conrod nearly vertical, arm joint radius r_con.
    # Exact loop closure: pivot P at origin; arm joint A = r_con*(cos th, sin th);
    # crank centre C = (r_con, -(L_conrod)) at mid stroke; crank pin B = C + r_crank*(cos phi, sin phi)
    C = np.array([d.r_con, -d.L_conrod])
    th = np.zeros(n)
    for i, p in enumerate(phi):
        B = C + d.r_crank*np.array([np.cos(p), np.sin(p)])
        # solve |A - B| = L_conrod for theta via Newton
        t = 0.0 if i == 0 else th[i-1]
        for _ in range(50):
            A = d.r_con*np.array([np.cos(t), np.sin(t)])
            f = np.sum((A-B)**2) - d.L_conrod**2
            dA = d.r_con*np.array([-np.sin(t), np.cos(t)])
            df = 2*np.sum((A-B)*dA)
            t -= f/df
            if abs(f) < 1e-14: break
        th[i] = t
    z_tip = d.L_arm*np.sin(th)
    x_ret = d.L_arm*(1-np.cos(th))  # retraction away from the blade (positive = back)
    return phi, th, z_tip, x_ret

def clearance_angle_variation(d: Design):
    """Effective clearance angle varies across the block face because the trajectory is an arc."""
    half = np.degrees(np.arcsin((d.block_h/2)/d.L_arm))
    return d.clearance_nom_deg - half, d.clearance_nom_deg + half, half

# ----------------------------------------------------------------------------- flexure pivot
def pivot_props(d: Design):
    m = MAT[d.piv_mat]
    I_each = (d.piv_w_total/d.piv_n)*d.piv_t**3/12
    k_theta = d.piv_n*m["E"]*I_each/d.piv_L     # [N m / rad], pure-moment loaded strips
    th_max = np.arcsin(d.r_crank/d.r_con)
    sigma_max = m["E"]*d.piv_t*th_max/(2*d.piv_L)  # bending stress at full swing
    strain_max = d.piv_t*th_max/(2*d.piv_L)
    # axial (x) stiffness: strips at +-45 deg; force along x => each strip loaded by F cos45 axially
    A_total = d.piv_w_total*d.piv_t
    k_axial = m["E"]*A_total/d.piv_L*np.cos(np.pi/4)**2*1.0
    # centre shift (Wittrick / Zelenika-De Bona, strips crossing at mid-length):  delta ~ 0.1 * L * theta^2 (upper-bound order)
    center_shift = 0.1*d.piv_L*th_max**2
    T_max = k_theta*th_max
    return dict(k_theta=k_theta, theta_max_deg=np.degrees(th_max), sigma_max=sigma_max, strain_max=strain_max,
                k_axial=k_axial, center_shift=center_shift, torque_max=T_max,
                safety_yield=m["sy"]/sigma_max, safety_fatigue=m["se"]/sigma_max)

# ----------------------------------------------------------------------------- feed stage
def stage_props(d: Design):
    m = MAT[d.st_mat]
    I = d.st_w*d.st_t**3/12
    k_x = 2*12*m["E"]*I/d.st_h**3           # two fixed-guided plates, motion along x
    F_preload_max = k_x*d.st_travel
    sigma = 3*m["E"]*d.st_t*d.st_travel/d.st_h**2   # fixed-guided beam max bending stress
    strain = sigma/m["E"]
    # rotational stiffness about y (pitch) from plate axial stiffness, plates separated by s along x
    k_ax_plate = m["E"]*(d.st_w*d.st_t)/d.st_h
    k_pitch = k_ax_plate*d.st_s**2/2         # [N m/rad]
    # stiffness along z (plates in-plane shear/axial): very high
    k_z = 2*k_ax_plate
    # stiffness along y (plates in-plane bending about z): I_y = t*w^3/12
    k_y = 2*12*m["E"]*(d.st_t*d.st_w**3/12)/d.st_h**3
    return dict(k_x=k_x, F_preload_max=F_preload_max, sigma=sigma, strain=strain, k_pitch=k_pitch,
                k_z=k_z, k_y=k_y, safety_yield=m["sy"]/sigma, safety_creep=(15e6)/sigma)   # 15 MPa sustained-stress limit: an ASSUMPTION for PLA (~0.3 x yield), not taken from a measured creep curve

# ----------------------------------------------------------------------------- structural loop compliance (x direction at the cut)
def box_I(w, h, t):
    return (w*h**3 - (w-2*t)*(h-2*t)**3)/12

def box_A(w, h, t):
    return w*h - (w-2*t)*(h-2*t)

def loop_compliance(d: Design):
    """Compliance [m/N] of each element in the metrology loop, for a unit force F_x applied at the blade/specimen
    interface (specimen pushed -x, blade pushed +x).  Returns dict element->compliance."""
    fm = MAT[d.frame_mat]
    c = {}
    # 1 arm axial (force along arm at mid-stroke)
    A_arm = box_A(*d.arm_sec)
    c["arm_axial"] = d.L_arm/(fm["E"]*A_arm)
    # 2 pivot blades axial
    c["pivot_axial"] = 1/pivot_props(d)["k_axial"]
    # 3+4 frame (tower bending + base rotation): FE-calibrated pivot-to-deck compliance when available
    if d.frame_fe_compliance is not None:
        c["frame_tower_base_FE"] = d.frame_fe_compliance
    else:
        I_t = box_I(*d.tower_sec)
        c["tower_bending"] = d.z_pivot**3/(3*fm["E"]*I_t)
        I_b = 0.4*d.base_width*d.base_depth**3/12
        EI_base = fm["E"]*I_b
        if d.tie_rods:
            A_rod = np.pi*(6.5e-3/2)**2
            EI_base += 2*MAT["steel"]["E"]*A_rod*d.rod_offset**2
        k_theta_base = EI_base/d.L_arm
        h_na = d.z_pivot + d.base_depth/2
        c["base_rotation"] = (h_na**2)/k_theta_base
    # 5 feed train (rev. 2 load path): blade -> holder -> platform -> plates (pitch) -> tab nut -> fine bolt -> bolt head on
    #    rod end (steel-steel) -> coarse rod -> frame nut -> upstand -> foot -> deck.  No ball contact exists in the design.
    A_screw = np.pi*(6.5e-3/2)**2
    c["screw_axial"] = 0.045/(200e9*A_screw)            # fine bolt + rod between the two nuts (~45 mm of M8 core)
    c["coupler_steel_contact"] = 1/5.0e7                 # rod end bearing on the bolt head, steel on steel (linearised Hertz, ~2 N)
    # nut bedding on the printed pocket floor: bare AF13 face on PLA layer lines ~0.5e7 N/m (allowance for roughness/creep);
    # a 16 mm steel washer doubles the bearing area and removes the layer-line contact -> ~1.5e7; brass insert ~2e7
    k_seat = {"bare": 0.5e7, "washer": 1.5e7, "insert": 2.0e7}[d.nut_seat]
    c["tab_nut_bedding"] = 1/k_seat
    c["frame_nut_bedding"] = 1/k_seat
    c["upstand_FE"] = d.upstand_fe_compliance
    # 6 blade overhang: blade 0.25 mm thick cantilevered beyond the full-width back wall of the holder (the cutting force
    #    pushes the blade onto this wall); width 80 mm
    E_bl = 200e9; t_bl = 0.25e-3; w_bl = 0.080; L_ov = d.blade_support_overhang
    I_bl = w_bl*t_bl**3/12
    c["blade_overhang"] = L_ov**3/(3*E_bl*I_bl)
    # 7 blade bearing on the PLA back wall (80 x 5.5 mm band, ~10 mm of PLA behind it) + holder on platform
    c["blade_backwall_bearing"] = 0.010/(fm["E"]*(0.080*0.0055))
    # 8 cassette clamp (printed jaws, compression) + paraffin block
    c["cassette_clamp"] = 0.015/(fm["E"]*(0.040*0.010))
    # stage pitch coupling: force at height h_off above stage flexure mid-plane -> rotation -> blade x motion
    # FE-calibrated (scikit-fem plane-stress, P2): screw-held stage compliance at the blade incl. platform/tab bending
    c["stage_pitch"] = d.stage_fe_compliance
    return c

# ----------------------------------------------------------------------------- thermal
def thermal_drift(d: Design, dT=1.0):
    fm = MAT[d.frame_mat]
    # loop: PLA base span (L_arm) + tower (z_pivot) ~ in series, minus steel screw in loop
    drift = fm["cte"]*(d.L_arm + 0.02)*dT + MAT["steel"]["cte"]*0.06*dT - fm["cte"]*0.0*dT
    return drift  # [m] per dT

# ----------------------------------------------------------------------------- Monte Carlo error budget
def monte_carlo(d: Design, N=20000, teeth_per_stroke=1, seed=0, n_screws=200):
    """Section-thickness error budget by simulation of consecutive strokes.
    Thickness of section i = (blade position during cut i) - (blade position during cut i-1), so every error source enters
    as a DIFFERENCE between consecutive strokes:
      * cutting force: x_i = feed_i + C*F_i  ->  C*(F_i - F_{i-1}), s.d. sqrt(2)*C*sigma_F  (the rev. 1 model applied C*F_i only)
      * thread lead error: a fixed function of position along each screw, evaluated where each nut sits at stroke i; the
        model uses a cumulative (bias) term, a per-revolution periodic term and a local random term, each drawn ONCE per
        simulated screw pair (n_screws realisations) and then walked along, so adjacent steps are correlated as on a real thread.
    Returns mean, s.d., percentiles, contributions, and the bias distribution (mean thickness of each simulated screw pair)."""
    rng = np.random.default_rng(seed)
    k = teeth_per_stroke
    step = d.feed_per_tooth()*k                       # nominal stage advance per stroke
    per_screw = N // n_screws
    thick = []; bias = []; e_lead_all = []
    for sidx in range(n_screws):
        # lead error functions for this screw pair, as a function of nut travel s (m) along each thread
        b1 = rng.uniform(-d.lead_bias_pct, d.lead_bias_pct)/100; b2 = rng.uniform(-d.lead_bias_pct, d.lead_bias_pct)/100
        A1 = d.lead_periodic_um*1e-6*rng.uniform(0.5, 1.0); A2 = d.lead_periodic_um*1e-6*rng.uniform(0.5, 1.0)
        ph1, ph2 = rng.uniform(0, 2*np.pi, 2)
        n = per_screw + 1
        # travel of each engagement point: coarse thread advances pitch_a per rev relative to the frame nut, fine thread pitch_b
        turns = np.arange(n)*k/d.n_teeth
        s1 = turns*d.pitch_a; s2 = turns*d.pitch_b
        loc1 = np.cumsum(rng.normal(0, d.lead_local_pct/100*d.pitch_a*k/d.n_teeth, n))   # local random lead walk
        loc2 = np.cumsum(rng.normal(0, d.lead_local_pct/100*d.pitch_b*k/d.n_teeth, n))
        e1 = b1*s1 + A1*np.sin(2*np.pi*s1/d.pitch_a + ph1) + loc1
        e2 = b2*s2 + A2*np.sin(2*np.pi*s2/d.pitch_b + ph2) + loc2
        x_lead = (s1 + e1) - (s2 + e2)                 # stage position from the screw train
        e_lead = np.diff(x_lead) - step                # per-stroke lead error
        e_lead_all.append(e_lead); bias.append(e_lead.mean())
    e_lead = np.concatenate(e_lead_all); N = len(e_lead)
    # ratchet tooth-pitch error (printed drum, 0.05 mm 1 s.d. on a 50 mm pitch circle) -> per-step angular error difference
    tooth = rng.normal(0, 0.05e-3/(0.050/2), N+1)
    e_tooth = np.diff(tooth)/(2*np.pi)*d.feed_per_rev()
    # cutting force through the loop compliance: difference between consecutive strokes
    C_tot = sum(loop_compliance(d).values())
    F = rng.normal(0, d.F_cut_var, N+1)
    e_comp = C_tot*np.diff(F)
    # thermal drift between consecutive strokes (1 s cadence; room drift 1 K / 30 min)
    e_th = thermal_drift(d, dT=1.0/1800)*rng.normal(1, 0.3, N)
    # pawl seating / drum settling: residual 0.1 um
    e_pawl = rng.normal(0, 0.1e-6, N)
    total = step + e_lead + e_tooth + e_comp + e_th + e_pawl
    bias = np.array(bias)
    return dict(nominal=step, mean=total.mean(), std=total.std(), p5=np.percentile(total, 5), p95=np.percentile(total, 95),
                bias_sd=bias.std(), bias_p5_pct=np.percentile(bias, 5)/step*100, bias_p95_pct=np.percentile(bias, 95)/step*100,
                contrib=dict(thread_lead=np.std(e_lead - np.repeat(bias, per_screw)), thread_bias=bias.std(), ratchet_tooth=np.std(e_tooth),
                             compliance=np.std(e_comp), thermal=np.std(e_th), pawl=np.std(e_pawl)),
                samples=total)

def summary(d: Design):
    piv = pivot_props(d); st = stage_props(d); lc = loop_compliance(d)
    C = sum(lc.values())
    cl_lo, cl_hi, half = clearance_angle_variation(d)
    phi, th, z, xr = crank_rocker(d)
    return dict(
        feed_per_rev_um=d.feed_per_rev()*1e6, feed_per_tooth_um=d.feed_per_tooth()*1e6,
        stroke_mm=(z.max()-z.min())*1e3, theta_max_deg=np.degrees(th.max()), retraction_at_ends_mm=xr.max()*1e3,
        face_sagitta_um=(d.L_arm*(1-np.cos(np.arcsin(d.block_h/2/d.L_arm))))*1e6,
        clearance_range_deg=(cl_lo, cl_hi),
        pivot=piv, stage=st, loop_compliance_um_per_N={k: v*1e6 for k, v in lc.items()},
        loop_stiffness_N_per_um=1/(C*1e6), thermal_um_per_K=thermal_drift(d)*1e6,
    )

if __name__ == "__main__":
    d = Design()
    import pprint; pprint.pprint(summary(d))
    mc = monte_carlo(d); mc.pop("samples"); pprint.pprint(mc)
