"""Plane-stress FE (scikit-fem, P2) of the stage foot + frame-nut upstand in the x-z plane.
Load: +x force at the frame-nut pocket floor (the screw reaction).  Support: foot bottom fixed (clamped to the deck).
Out-of-plane widths differ per region, so E is scaled per element by width/REF_W.  Units: mm, N, MPa.
Two geometries: rev. 1 (6 mm wall, 32 wide, 2 mm gap closed for the model, as the critic analysed) and rev. 2 (10 mm wall,
64 wide, lower gusset 32 wide and two 6 mm side gussets)."""
import numpy as np, json, sys
from skfem import MeshTri, Basis, ElementVector, ElementTriP2, ElementTriP0, BilinearForm, asm, condense, solve
from skfem.helpers import ddot, sym_grad, trace

E0 = 3500.0; NU = 0.36; REF_W = 60.0
Z_FOOT0, Z_FOOT1 = 22.0, 28.0

def regions(rev):
    if rev == 1:   # rev. 1 geometry: deck at 30, foot 30-36, 6 mm wall; the 2 mm gap of the rev. 1 CAD is closed here (best case)
        return [((130, 192, 30.0, 36.0), 90.0), ((186, 192, 36.0, 92), 32.0)], (188.5, 193.0), 70.0
    reg = [((130, 196, Z_FOOT0, Z_FOOT1), 90.0), ((184, 194, Z_FOOT1, 82), 64.0)]
    return reg, (184.0, 188.6), 60.0

def width_field(cx, cz, rev):
    reg, _, _ = regions(rev)
    w = np.zeros_like(cx)
    for (x0, x1, z0, z1), wd in reg:
        m = (cx >= x0) & (cx <= x1) & (cz >= z0) & (cz <= z1); w[m] = wd
    if rev == 2:
        # lower central gusset: triangle (170,28)-(184,28)-(184,46), width 32
        m = (cx >= 170) & (cx <= 184) & (cz >= Z_FOOT1) & (cz <= Z_FOOT1 + (cx-170)/14*18) & (w == 0); w[m] = 32.0
        # side gussets: trapezoid (172,28)-(184,28)-(184,70)-(178,70), 2 x 6 wide
        m = (cx >= 172) & (cx <= 184) & (cz >= Z_FOOT1) & (cz <= 70) & (cx >= 172 + (cz-Z_FOOT1)/42*6) & (w == 0); w[m] = 12.0
        # where the lower gusset and side gussets coexist with the wall they are already counted in the wall width
    return w

def run(rev, h=0.5, F=1.0):
    reg, nut_x, nut_z = regions(rev)
    z_lo = 30.0 if rev == 1 else Z_FOOT0
    xs = np.arange(128, 198.0001, h); zs = np.arange(z_lo, 94.0001, h)
    m = MeshTri.init_tensor(xs, zs)
    c = m.p[:, m.t].mean(axis=1); w = width_field(c[0], c[1], rev)
    m = MeshTri(m.p, m.t[:, w > 0]).remove_unused_nodes() if hasattr(MeshTri, "remove_unused_nodes") else MeshTri(m.p, m.t[:, w > 0])
    c = m.p[:, m.t].mean(axis=1); w = width_field(c[0], c[1], rev)
    basis = Basis(m, ElementVector(ElementTriP2())); b0 = basis.with_element(ElementTriP0())
    Ef = b0.interpolate(E0*w/REF_W)
    lam = lambda E: E*NU/((1+NU)*(1-NU)); mu = lambda E: E/(2*(1+NU))
    @BilinearForm
    def stiff(u, v, wq):
        E = wq["E"]; return 2*mu(E)*ddot(sym_grad(u), sym_grad(v)) + lam(E)*trace(sym_grad(u))*trace(sym_grad(v))
    K = asm(stiff, basis, E=Ef)
    D = basis.get_dofs(lambda x: x[1] < z_lo + 1e-6).all()
    # load on the nut-pocket floor: nodes at x = nut_x[1] (pocket floor plane), |z - nut_z| <= 6.5 (nut AF/2)
    x_floor = nut_x[0] if rev == 1 else nut_x[1]
    ln = basis.get_dofs(lambda x: (np.abs(x[0]-x_floor) < h/2+1e-6) & (np.abs(x[1]-nut_z) <= 6.5)).nodal["u^1"]
    f = np.zeros(basis.N); f[ln] = F/len(ln)
    u = solve(*condense(K, f, D=D))
    ux = u[ln].mean()/REF_W          # mm per N (unit-thickness model scaled to REF_W)
    return dict(rev=rev, compliance_um_per_N=float(ux*1e3), n_elements=int(m.t.shape[1])), (m, basis, u)

if __name__ == "__main__":
    out = {}
    for rev in (1, 2):
        r, _ = run(rev); out[f"rev{rev}"] = r
        print(f"rev {rev}: upstand compliance at the frame nut = {r['compliance_um_per_N']:.3f} um/N  ({r['n_elements']} elements)")
    # analytic check for rev 1: cantilever 6 x 32 mm, L = 70-28 = 42 mm above the foot top
    I = 32*6**3/12; L = 70-36; print("rev 1 cantilever estimate (34 mm above the foot): %.2f um/N" % (L**3/(3*E0*I)*1e3))
    json.dump(out, open("fea_upstand_results.json", "w"), indent=1)
