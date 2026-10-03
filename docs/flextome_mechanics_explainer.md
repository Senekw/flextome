# FlexTome, how the machine works (mechanics walkthrough)

Companion to the paper (Section "Mechanism walkthrough"). Coordinates: *x* = feed axis (specimen → blade), *y* = along the blade edge, *z* = up; pivot axis is the *y*-parallel line through P = (0, 0, 95) mm.

![Numbered assembly](fig_assembly_callouts.png)
![Mechanism section](fig_section_callouts.png)
![Exploded view](fig_render_exploded.png)

## 1. Drive train: handwheel → crankshaft → connecting rod → arm
- **Handwheel (19)** Ø110 mm, five spokes, offset handle. One turn = one section.
- **Built-up crankshaft (18)**: two 8 mm steel stubs in two 608 bearings (pillow blocks at y = ±50) grub-screwed into two printed crank webs (Ø44) carrying a 5 mm crank pin at r = 13 mm. The **fork connecting rod (17)** runs *between* the webs, so it never crosses the shaft when the pin is below the axis.
- Crank centre C = (60, 23.8) is exactly one rod length ℓ = 71.2 mm below the arm joint A (60 mm from the pivot), so the swing is symmetric:
  `sin θmax = r / r_con = 13/60 → θ = ±12.5°`, stroke `s = 2·125·sin θmax = 54.2 mm`.
- Loop closure `|A(θ) − B(φ)| = ℓ` gives θ(φ) (fig M1b). Arm horizontal at φ₀ = asin(r/2ℓ) = 5.2° (CAD rest pose).
- At 60 rpm the face crosses the edge at 153–174 mm/s; a 20 mm face is cut in 118 ms (excitation < 9 Hz, first structural mode 357 Hz).
- Handwheel effort is dominated by the cut itself: pivot torque 0.058 N·m → ~0.013 N·m at the shaft (0.23 N at the handle).

![Kinematics](figM1_kinematics.png)

## 2. The cut: blade and block geometry
- Disposable low-profile blade 80 × 8 × 0.25 mm, 35° facet, clamped behind a 1.0 mm steel pressure strip, edge 2.5 mm proud, holder leaning 27° back from the face.
- `clearance α = 27° − 17.5° = 9.5°`, `rake γ = 90° − 27° − 17.5° = 45.5°`, `α + β/2 + γ = 90°`.
- Face rides a 125 mm arc → clearance sweeps 9.5° ± 4.6° (4.9°–14.1°) during the cut; face is cylindrical, sagitta `h = 125 − √(125² − 10²) = 401 µm` over 20 mm. Every stroke follows the same arc about the same flexure-defined centre, so the arc is trimmed in once and successive sections are uniform.
- The pressure-strip top sits 2.5 mm below the edge along the blade, just outside the block's swept circle (this is why a 2.0 mm overhang is not available).

![Blade geometry](figM2_blade_geometry.png)

## 3. The arm pivot: a cross-strip flexure
- Three 0.20 mm feeler-gauge strips: two 10 mm wide at y = ±15 in one plane, one 20 mm wide at y = 0 in the perpendicular plane, crossing at the pivot axis, 20 mm free length, clamped under printed pads with M3 screws.
- `kθ = Σ E w t³ / 12L = 0.267 N·m/rad`, `σmax = E t θmax / 2L = 218 MPa` (yield SF 5.0, fatigue SF 2.1), `Tmax = 0.058 N·m`.
- Axially 40 N/µm (0.025 µm/N in the loop). Centre shift ≈ 0.1·L·θ² ≈ 95 µm at full swing; this is deterministic in θ, so it reshapes the arc but not the thickness.
- Middle strip twice the width of the outer two, placed between them → both strip planes balanced about y = 0, no twist under off-centre load.

## 4. Feed train: differential screw, drum ratchet, selector
- The **blade** moves, not the specimen. Stage (8) = monolithic parallelogram: platform on 1.0 × 44 mm PLA plates (front 60 wide, rear 2 × 20). Hanging tab carries an M8×1.0 nut; fixed upstand carries an M8×1.25 nut. One screw assembly runs through both (M8×1.0 bolt → printed coupler + jam nut → M8×1.25 rod → 50-tooth drum (12) → jam nuts → knob (14)).
- `Δx/rev = p₁ − p₂ = 1.25 − 1.00 = 0.25 mm`, `Δx/tooth = 0.25/50 = 5.0 µm`.
- Flexure preload 4.9 N at full travel keeps the tab on one thread flank → no backlash.
- Consequences: (i) the screw body translates 5× the stage → 14 mm **drum**, not a wheel, keeps the pawl engaged over ±5 mm (±1 mm stage travel = 400 sections per reset); (ii) only the pitch *difference* matters: `δt/t = (δp₁ − δp₂)/0.25 mm` → 0.37 µm (1 s.d.) with 1 % pitch errors; (iii) both threads are ordinary rolled hardware.
- Selector shield (13) on the stage column masks all but k = 1…6 teeth → 5k µm per stroke. Coarse approach: loosen two M6 knobs (9), push with the M6 side screw (10) through ±10 mm slots.

![Differential screw](figM3_differential_screw.png)

## 5. Feed actuation: one turn = one section
- Pawl rod (15) in a guide on the stage column, finger on the drum under gravity + light spring. Striker bar (16) on the arm reaches over the rod's flag.
- Up-stroke: past θ ≈ 6° (block 13 mm clear of the edge) the bar lifts the rod ≤ 22 mm; the finger rides back over k teeth.
- Down-stroke: past θ ≈ 6° the bar releases; the finger drops behind the k-th tooth and rod weight + spring index the drum k × 7.2° → blade advances 5k µm **before** the face reaches the edge. The drum is at rest during the cut.
- No return-stroke retraction (3.0 mm natural retraction at the stroke ends only); this is a documented limitation shared with the classical Cambridge rocker.

![Feed timing](figM5_feed_timing.png)

## 6. The force path (metrology loop)
Only parts between the blade edge and the block face can change thickness: edge → overhang → strip/holder → platform → flexure plates (pitch about the nut) → tab nut → screw → frame nut in PLA → upstand/deck/base → tower → strip pivot → arm → cassette clamp → block.

`C_loop = Σ Cᵢ = 1.00 µm/N → k_loop = 1.0 N/µm`

| element | µm/N |
|---|---|
| blade overhang (2.5 mm) | 0.250 |
| frame: upstand → base → tower (FE) | 0.213 |
| nut bedding in PLA | 0.200 |
| stage pitch (FE) | 0.109 |
| screw-tip contact | 0.100 |
| arm axial | 0.080 |
| flexure pivot axial | 0.025 |
| cassette clamp | 0.011 |
| screw axial | 0.009 |
| pressure strip / holder | 0.003 |

A steady 2 N feed-axis force opens the loop 2 µm, but a steady deflection is trimmed into the face; only stroke-to-stroke *variation* matters (0.4 N s.d. → 0.40 µm). The cheapest upgrade is a brass insert for the frame nut (1.00 → 0.85 µm/N; scatter 0.56 → 0.52 µm).

![Pivot and loop](figM4_pivot_loop.png)

## 7. Operating sequence
1. Clamp the cassette under the sliding jaw (M5 thumbscrew), block face proud of the jaw.
2. Loosen the two M6 knobs, push the stage to ~0.5 mm from the face with the coarse screw, re-clamp.
3. Selector to 6 (30 µm); turn until the whole face is trimmed to the arc.
4. Selector to 1 (5 µm); each turn yields one section, carried off as a ribbon.
5. After ~400 sections wind the feed knob back 20 turns and repeat step 2.

## 8. What the interference audit changed
Blade holder vs swept circle (→ 27° lean, 2.5 mm overhang); striker bar through pawl rod (→ outboard bar + flag); coarse boss into stage foot; conrod eye position; screw body translating 5× (→ drum, coupler gap, shield on column); conrod sweeping through a continuous crankshaft at bottom dead centre (→ two-web built-up crankshaft); 1.2 mm crank-centre offset giving −11.4°/+13.7° (→ crank centre at z = 23.8, ±12.5°). Residual intersections are intended (pins, strips, pawl tooth 16 mm³, striker on flag between 6° and 12.5°).
