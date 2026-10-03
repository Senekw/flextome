# FlexTome, how the machine works (mechanics walkthrough, revision 2)

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
- Disposable low-profile blade 80 × 8 × 0.25 mm, 35° facet. Its rear face bears on a full-width printed back wall (this is what carries the cutting force, which pushes the blade away from the specimen); a 1.0 mm steel clamp strip in front is pulled onto a 0.20 mm blade recess by two M4 screws below the blade. Edge 2.5 mm proud, holder leaning 27° back from the face.
- `clearance α = 27° − 17.5° = 9.5°`, `rake γ = 90° − 27° − 17.5° = 45.5°`, `α + β/2 + γ = 90°`.
- Face rides a 125 mm arc → clearance sweeps 9.5° ± 4.6° (4.9°–14.1°) during the cut; face is cylindrical, sagitta `h = 125 − √(125² − 10²) = 401 µm` over 20 mm. Every stroke follows the same arc about the same flexure-defined centre, so the arc is trimmed in once and successive sections are uniform.
- The clamp-strip top sits 2.5 mm below the edge along the blade, just outside the block's swept circle; the back wall behind the blade is not bound by the circle and could rise to 2.0 mm (a lever).

![Blade geometry](figM2_blade_geometry.png)

## 3. The arm pivot: a cross-strip flexure
- Three 0.20 mm feeler-gauge strips: two 10 mm wide at y = ±15 in one plane, one 20 mm wide at y = 0 in the perpendicular plane, crossing at the pivot axis, 20 mm free length, clamped under printed pads with M3 screws.
- `kθ = Σ E w t³ / 12L = 0.267 N·m/rad`, `σmax = E t θmax / 2L = 218 MPa` (yield SF 5.0, fatigue SF 2.1), `Tmax = 0.058 N·m`.
- Axially 40 N/µm (0.025 µm/N in the loop). Centre shift ≈ 0.1·L·θ² ≈ 95 µm at full swing; this is deterministic in θ, so it reshapes the arc but not the thickness.
- Middle strip twice the width of the outer two, placed between them → both strip planes balanced about y = 0, no twist under off-centre load.

## 4. Feed train: differential screw, drum ratchet, selector
- The **blade** moves, not the specimen. Stage (8) = monolithic parallelogram: platform on 1.0 × 44 mm PLA plates (front 60 wide, rear 2 × 20). A gusseted hanging tab carries an M8×1.0 nut on a steel washer; a 10 mm thick, 64 mm wide gusseted upstand (continuous with the foot) carries an M8×1.25 nut on a washer. One screw assembly runs through both: M8×1.0 hex bolt → 13 mm printed coupler (bolt head in a hex pocket, captured M8×1.25 thin nut) → M8×1.25 rod screwed in until its end bears on the bolt head → frame nut → 50-tooth drum (12) → jam nuts → knob (14). The axial path is steel throughout; the printed coupler carries only torque.
- `Δx/rev = p₁ − p₂ = 1.25 − 1.00 = 0.25 mm`, `Δx/tooth = 0.25/50 = 5.0 µm`.
- Two extension springs (≈10 N total) between tab and upstand pull the tab in +x, the direction of the cutting force, so the nut stays on one flank over the whole ±1 mm travel (the symmetric flexure's own 4.9 N restoring force would reverse at mid-travel).
- Consequences: (i) the screw body translates 5× the stage → 14 mm **drum**, not a wheel, keeps the pawl engaged over ±5 mm (±1 mm stage travel = 400 sections per reset); (ii) only the pitch *difference* matters, but it is a fifth of either pitch, so relative lead errors are amplified five-fold: `δt/t = (δp₁ − δp₂)/0.25 mm` → a per-machine bias within ±0.7 % (calibrate once) plus 0.18 µm (1 s.d.) of scatter; (iii) both threads are ordinary rolled hardware.
- Selector shield (13) on the stage column masks all but k = 1…6 teeth → 5k µm per stroke. Coarse approach: loosen two M6 knobs (9), push with the M6 side screw (10) through ±10 mm slots.

![Differential screw](figM3_differential_screw.png)

## 5. Feed actuation: one turn = one section
- Pawl rod (15) in a guide on the stage column, finger on the drum under gravity + light spring. Striker bar (16) on the arm reaches over the rod's flag.
- Up-stroke: past θ ≈ 6° (block 13 mm clear of the edge) the bar lifts the rod ≤ 22 mm; the finger rides back over k teeth.
- Down-stroke: past θ ≈ 6° the bar releases; the finger drops behind the k-th tooth and rod weight + spring index the drum k × 7.2° → blade advances 5k µm **before** the face reaches the edge. The drum is at rest during the cut.
- No return-stroke retraction (3.0 mm natural retraction at the stroke ends only); this is a documented limitation shared with the classical Cambridge rocker.

![Feed timing](figM5_feed_timing.png)

## 6. The force path (metrology loop)
Only parts between the blade edge and the block face can change thickness: edge → overhang → back wall and holder → platform → flexure plates (pitch about the tab nut) → tab nut on washer → bolt, rod end on bolt head, rod → frame nut on washer → gusseted upstand → foot, deck, base → tower → strip pivot → arm → cassette clamp → block.

`C_loop = Σ Cᵢ = 1.20 µm/N → k_loop = 0.83 N/µm`

| element | µm/N |
|---|---|
| frame: tower + base on two end pads (FE; 0.21 with full contact) | 0.353 |
| blade overhang (2.5 mm) | 0.250 |
| gusseted frame-nut upstand (FE; rev. 1 was 7.3) | 0.185 |
| stage pitch about the tab nut (FE) | 0.132 |
| arm axial | 0.080 |
| tab nut on washer | 0.067 |
| frame nut on washer | 0.067 |
| flexure pivot axial | 0.025 |
| rod end on bolt head | 0.020 |
| cassette clamp | 0.011 |
| fine bolt + rod axial | 0.007 |
| blade on back wall | 0.006 |

A steady 2 N feed-axis force opens the loop 2.4 µm, but a steady deflection is trimmed into the face; only the stroke-to-stroke *change* matters (0.4 N per-stroke s.d. → √2 × 0.4 N of change → 0.68 µm). Predicted section at the 5 µm setting: 5.01 ± 0.72 µm (90 % interval 3.84 to 6.19 µm). Levers: back-wall overhang 2.0 mm, base shimmed to full contact, brass inserts → 0.90 µm/N and 0.56 µm together.

![Pivot and loop](figM4_pivot_loop.png)

## 7. Operating sequence
1. Clamp the cassette under the sliding jaw (M5 thumbscrew), block face proud of the jaw.
2. Loosen the two M6 knobs, push the stage to ~0.5 mm from the face with the coarse screw, re-clamp.
3. Selector to 6 (30 µm); turn until the whole face is trimmed to the arc.
4. Selector to 1 (5 µm); each turn yields one section, carried off as a ribbon.
5. After ~400 sections wind the feed knob back 20 turns and repeat step 2.

All numbers above are model predictions; the machine has not been built or measured.

## 8. What the audits and the external review changed
Interference audit (`audit_interference.py`, 12 crank angles × 3 screw positions, plus a mesh-connectivity check): blade holder vs swept circle (→ 27° lean, 2.5 mm overhang); striker bar through pawl rod (→ outboard bar + flag); coarse boss into stage foot; conrod eye position; screw body translating 5× (→ drum, coupler gap, shield on column); conrod sweeping through a continuous crankshaft (→ two-web crankshaft); 1.2 mm crank-centre offset (→ ±12.5°); in revision 2: drum vs thicker upstand (→ +3 mm), coupler vs lowered platform (→ screw axis z = 60), knob heads vs side gussets, clamp-screw heads vs platform (→ relief pockets), enclosed pivot pockets and sealed arm cavity.

External review of revision 1 (the decisive one): the frame-nut upstand was disconnected from the foot and only 6 mm thick (7.3 µm/N on its own, so the real loop was ≈8.5 µm/N and the predicted scatter ≈4.8 µm); the blade was pressed onto two grub-screw tips instead of a wall; the coupler was 0.7 mm thick with no rod engagement; the budget had a ball contact that did not exist and one nut too few; the flexure preload reversed at mid-travel; the Monte-Carlo applied C·F instead of C·ΔF and redrew thread pitch every stroke. Every item is addressed in revision 2 (Table 2 of the report).
