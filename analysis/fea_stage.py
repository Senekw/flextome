"""2D plane-stress FE model (scikit-fem) of the FlexTome feed stage in the x-z plane.
Out-of-plane widths differ per region, so E is scaled per element by width/60 mm (reference width)."""
import numpy as np, json
import skfem
from skfem import MeshTri, Basis, ElementVector, ElementTriP1, ElementTriP2, ElementTriP0, BilinearForm, asm, condense, solve
from skfem.helpers import ddot, sym_grad, trace, eye

E0 = 3.5e9; NU = 0.36            # PLA
REF_W = 60.0
# regions: (x0,x1,z0,z1,width_mm)
REG = [
    (130, 192, 30, 36, 90),      # foot (plus column base)
    (137, 138, 36, 80, 60),      # front flexure plate
    (180, 181, 36, 80, 40),      # rear plates (2 x 20)
    (128, 184, 80, 88, 90),      # platform (clearance to upstand at 186)
    (127.5, 160, 88, 92.5, 84),  # blade holder body
    (166, 176, 58, 80, 24),      # nut tab
    (186, 192, 36, 92, 32),      # rear upstand (fixed, carries frame nut)
]
def inside(cx, cz):
    w = np.zeros_like(cx)
    for (x0, x1, z0, z1, wd) in REG:
        m = (cx >= x0) & (cx <= x1) & (cz >= z0) & (cz <= z1)
        w[m] = wd
    return w

def run(h=0.5, F=1.0, load_z=92.5, nut_fixed=True):
    xs = np.arange(126, 192.0001, h); zs = np.arange(30, 93.0001, h)
    m = MeshTri.init_tensor(xs, zs)
    c = m.p[:, m.t].mean(axis=1)
    w = inside(c[0], c[1])
    keep = w > 0
    m = MeshTri(m.p, m.t[:, keep]).remove_unused_nodes() if hasattr(MeshTri, "remove_unused_nodes") else MeshTri(m.p, m.t[:, keep])
    c = m.p[:, m.t].mean(axis=1); w = inside(c[0], c[1])
    basis = Basis(m, ElementVector(ElementTriP2()))
    b0 = basis.with_element(ElementTriP0())
    Efield = b0.interpolate(E0*w/REF_W)
    lam_f = lambda E: E*NU/((1+NU)*(1-NU))     # plane stress Lame
    mu_f = lambda E: E/(2*(1+NU))
    @BilinearForm
    def stiff(u, v, wq):
        E = wq["E"]
        return 2*mu_f(E)*ddot(sym_grad(u), sym_grad(v)) + lam_f(E)*trace(sym_grad(u))*trace(sym_grad(v))
    K = asm(stiff, basis, E=Efield)
    # BCs: foot bottom z=30 fixed (all dofs); nut tab x-dofs fixed (screw holds it) if nut_fixed
    dofs_fix = basis.get_dofs(lambda x: x[1] < 30 + 1e-6).all()
    if nut_fixed:
        tab = basis.get_dofs(lambda x: (x[0] > 166-1e-6) & (x[0] < 176+1e-6) & (x[1] > 64) & (x[1] < 76))
        dofs_fix = np.unique(np.concatenate([dofs_fix, tab.nodal["u^1"]]))
    # load: F (N) in +x spread over nodes at top-front of holder (x<=129, z >= load_z-0.5)
    ln = basis.get_dofs(lambda x: (x[0] < 129.01) & (x[1] > load_z-0.51)).nodal["u^1"]
    f = np.zeros(basis.N); f[ln] = F/len(ln)
    u = solve(*condense(K, f, D=dofs_fix))
    # displacement at load nodes (x), in metres -> mm input means u in mm (E in Pa, lengths mm -> u in mm * 1e-3? ) handle units:
    # K assembled with lengths in mm and E in Pa=N/m^2 -> stiffness in N/mm * 1e-6 ... convert: K_true = K*1e-3 (N/mm). so u_true[mm] = u/1e-3*? simpler: redo with consistent units below.
    return m, basis, u, ln, w

if __name__ == "__main__":
    # Use consistent units: mm, N, MPa.  E0 given in Pa -> MPa.
    E0 = E0/1e6
    m, basis, u, ln, w = run()
    ux = u[ln].mean()/REF_W      # mm per 1 N, divided by the 60 mm reference depth (unit-thickness model)
    print("stage compliance at blade (x) = %.4f um/N" % (ux*1e3))
    # free stage (no nut) -> k_x of the parallelogram
    m2, b2, u2, ln2, _ = run(nut_fixed=False)
    kx = REF_W/(u2[ln2].mean())
    print("parallelogram k_x = %.2f N/mm (analytic 4.1-4.9)" % kx)
    json.dump(dict(stage_compliance_um_per_N=ux*1e3, k_x_N_per_mm=kx, n_elements=int(m.t.shape[1]), element="P2 plane stress"), open("fea_results.json", "w"))
    # plot deformed shape
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    from matplotlib.tri import Triangulation
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.4))
    for ax, (mm, uu, title) in zip(axs, [(m2, u2, "free parallelogram, 1 N at blade"), (m, u, "screw-held stage, 1 N at blade")]):
        tri = Triangulation(mm.p[0], mm.p[1], mm.t.T)
        U = uu.reshape(-1, 2) if uu.size == 2*mm.p.shape[1] else None
        ux_n = uu[basis.nodal_dofs[0]] if mm is m else uu[b2.nodal_dofs[0]]
        uz_n = uu[basis.nodal_dofs[1]] if mm is m else uu[b2.nodal_dofs[1]]
        scale = 0.25*60/max(abs(ux_n).max(), 1e-12)
        tc = ax.tripcolor(Triangulation(mm.p[0]+scale*ux_n, mm.p[1]+scale*uz_n, mm.t.T), ux_n*1e3/REF_W, shading="gouraud", cmap="viridis")
        ax.triplot(tri, color="0.7", lw=0.15)
        ax.set_aspect("equal"); ax.set_title(title, fontsize=9); ax.set_xlabel("x (mm)"); ax.set_ylabel("z (mm)")
        cb = fig.colorbar(tc, ax=ax, shrink=0.8); cb.set_label("x-displacement (um per N)", fontsize=8)
    fig.suptitle("Plane-stress FE of the feed stage (deformation exaggerated)", fontsize=10)
    fig.tight_layout(); fig.savefig("fig5_fea_stage.png", dpi=300)
