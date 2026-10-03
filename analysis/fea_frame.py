"""2D plane-stress FE (scikit-fem, P2) of the FlexTome frame side section (x-z), unit thickness with E scaled by
effective out-of-plane width / 150 mm.  Ribbed base: web fraction 0.4 (walls+ribs) except 6 mm solid top skin.
Tie rods: 2 x M8 steel (A = 2*33.2 mm^2) at z=-45 -> represented as a 6.5 mm tall steel band with width scaling 2*6.5/150.
Load: self-equilibrated couple: +F at pivot (x=0,z=95) along x, -F at the stage-foot clamp zone (x 150..172, z=30).
Output: relative x displacement pivot vs deck zone per N -> base+tower compliance."""
import numpy as np, json
from skfem import MeshTri, Basis, ElementVector, ElementTriP2, ElementTriP0, BilinearForm, asm, condense, solve
from skfem.helpers import ddot, sym_grad, trace

E_PLA = 3500.0; E_ST = 200000.0; NU = 0.36; REF_W = 150.0
GUSSET = True
def props(cx, cz):
    """return (E_eff) per element centroid; 0 outside."""
    E = np.zeros_like(cx)
    # base region x -40..195, z -55..0 ; deck x 110..195 z 0..30 ; tower x -40..-8 z 0..118 (width 60)
    base = (cx >= -40) & (cx <= 195) & (cz >= -55) & (cz <= 0)
    deck = (cx >= 110) & (cx <= 195) & (cz >= 0) & (cz <= 30)
    if GUSSET:
        # polygon (-40,0),(20,0),(20,55),(-8,80),(-8,118),(-40,118)
        tower = (cx >= -40) & (cz >= 0) & (cz <= 118) & ((cx <= -8) | ((cz <= 55) & (cx <= 20)) | ((cz > 55) & (cz < 80) & (cx <= 20 - 28*(cz-55)/25)))
    else:
        tower = (cx >= -40) & (cx <= -8) & (cz >= 0) & (cz <= 118)
    skin = base & (cz >= -6)
    deckskin = deck & (cz >= 24)
    E[base] = E_PLA*0.4            # ribbed web fraction
    E[skin] = E_PLA*1.0
    E[deck] = E_PLA*0.4
    E[deckskin] = E_PLA*1.0
    E[tower] = E_PLA*(60/REF_W)
    rods = (cx >= -40) & (cx <= 195) & (cz >= -48.25) & (cz <= -41.75)
    E[rods] = E_ST*(2*6.5/REF_W)
    return E

def run(h=1.0, tie_rods=True, F=1.0):
    xs = np.arange(-40, 195.001, h); zs = np.arange(-55, 118.001, h)
    m = MeshTri.init_tensor(xs, zs)
    c = m.p[:, m.t].mean(axis=1); E = props(c[0], c[1])
    if not tie_rods:
        rods = (c[1] >= -48.25) & (c[1] <= -41.75) & (E > 0); E[rods] = E_PLA*0.4
    keep = E > 0
    m = MeshTri(m.p, m.t[:, keep]).remove_unused_nodes()
    c = m.p[:, m.t].mean(axis=1); E = props(c[0], c[1])
    if not tie_rods:
        rods = (c[1] >= -48.25) & (c[1] <= -41.75); E[rods] = E_PLA*0.4
    basis = Basis(m, ElementVector(ElementTriP2()))
    Ef = basis.with_element(ElementTriP0()).interpolate(E)
    @BilinearForm
    def stiff(u, v, w):
        Eq = w["E"]; lam = Eq*NU/((1+NU)*(1-NU)); mu = Eq/(2*(1+NU))
        return 2*mu*ddot(sym_grad(u), sym_grad(v)) + lam*trace(sym_grad(u))*trace(sym_grad(v))
    K = asm(stiff, basis, E=Ef)
    # supports: bottom z=-55 fixed in z (table); one node fixed in x (rigid body) at x=80
    bot = basis.get_dofs(lambda x: x[1] < -55+1e-6)
    D = np.concatenate([bot.nodal["u^2"], basis.get_dofs(lambda x: (x[1] < -55+1e-6) & (np.abs(x[0]-80) < 0.6)).nodal["u^1"]])
    piv = basis.get_dofs(lambda x: (x[0] < -8+1e-6) & (x[0] > -9) & (np.abs(x[1]-95) < 3)).nodal["u^1"]
    deckz = basis.get_dofs(lambda x: (x[0] > 150) & (x[0] < 172) & (x[1] > 30-1e-6)).nodal["u^1"]
    f = np.zeros(basis.N); f[piv] = F/len(piv); f[deckz] = -F/len(deckz)
    u = solve(*condense(K, f, D=D))
    rel = (u[piv].mean() - u[deckz].mean())/REF_W   # mm per N (unit thickness -> /150 mm)
    return m, basis, u, rel

if __name__ == "__main__":
    out = {}
    import sys
    for g in (False, True):
        GUSSET = g
        for tr in (True, False):
            m, basis, u, rel = run(tie_rods=tr)
            key = ("gusset" if g else "slender") + ("_tierods" if tr else "_norods")
            out[key] = rel*1e3
            print(key, "frame compliance pivot-deck: %.3f um/N" % (rel*1e3), "elements", m.t.shape[1])
    GUSSET = True
    json.dump(out, open("fea_frame_results.json", "w"))
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    from matplotlib.tri import Triangulation
    m, basis, u, rel = run(tie_rods=True)
    ux = u[basis.nodal_dofs[0]]/REF_W*1e3; uz = u[basis.nodal_dofs[1]]/REF_W*1e3
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    sc = 25/np.abs(ux).max()
    tc = ax.tripcolor(Triangulation(m.p[0]+sc*ux, m.p[1]+sc*uz, m.t.T), ux, shading="gouraud", cmap="coolwarm")
    ax.triplot(Triangulation(m.p[0], m.p[1], m.t.T), color="0.75", lw=0.1)
    ax.set_aspect("equal"); ax.set_xlabel("x (mm)"); ax.set_ylabel("z (mm)")
    ax.set_title(f"Frame under 1 N pivot-to-deck couple: relative compliance {rel*1e3:.2f} um/N (deformation x{sc:.0f})", fontsize=9)
    cb = fig.colorbar(tc, ax=ax, shrink=0.8); cb.set_label("x-displacement (um per N)", fontsize=8)
    fig.tight_layout(); fig.savefig("fig6_fea_frame.png", dpi=300)
