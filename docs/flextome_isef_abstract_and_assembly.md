# FlexTome, ISEF project abstract (≤250 words) and assembly instructions

## Abstract (Engineering Mechanics / Biomedical Engineering)

Cancer diagnosis depends on 3–5 µm paraffin sections, but the rotary microtomes that cut them cost US$9,000–22,700, and much of sub-Saharan Africa has one pathologist per million people. The thinnest published open-source sectioning device cuts at 50 µm. I designed FlexTome, a rocking-arm microtome printed on a desktop FDM printer and finished with hardware-store parts for US$60. Three decisions make 5 µm sectioning plausible in printed plastic: the arm pivot is a cross-strip flexure of spring-steel feeler-gauge strips, so the knife-to-specimen loop has no sliding joint and no backlash; the 5 µm feed comes from a differential screw (M8 × 1.25 against M8 × 1.0, 0.25 mm/rev) indexed by a 50-tooth ratchet that the arm advances at top-dead-centre, so the micrometre chain is entirely steel; and the disposable blade rides on a monolithic printed flexure stage. I built a parametric CAD model (32 parts), an analytic stiffness and error budget, plane-stress finite-element models of the stage and frame, and a 20,000-stroke Monte-Carlo thickness simulation. The frame FE revealed a tower compliance four times larger than beam theory predicted and drove a gusset redesign; an interference audit over the full stroke caught and fixed three collisions. The design predicts a loop stiffness of 1.0 N/µm and sections of 5.0 ± 0.6 µm (1 s.d.), comparable to the >±10 % scatter measured on commercial instruments, with all structural resonances ≥357 Hz against <10 Hz cutting excitation. A seven-step validation protocol accompanies the open-source design files.

---

## Assembly instructions

### A. Print (≈27 h total on a 0.4 mm nozzle)
| Part | Material / settings | Orientation |
|---|---|---|
| frame (235 × 150 × 173 mm) | PLA, 0.2 mm layers, 4 perimeters, 100 % infill in tower/gusset and deck, 30 % elsewhere | Upside-down (rib pockets open upward), or split at x = 110 mm for small beds |
| arm | PLA, 0.2 mm, 4 perimeters, 40 % gyroid | Beam horizontal, clamp end up |
| stage_body_flexure | PETG or ASA preferred (PLA acceptable), 0.1 mm layers, 100 % infill | Flexure plates vertical |
| blade_holder, cassette_upper_jaw, conrod, crank_disc, handwheel, handle, feed_knob, ratchet_wheel_50T (drum), thickness_selector_shield, pawl_rod, screw_coupler, pivot_clamp_plates ×6 | PLA, 0.2 mm, 100 % infill for the ratchet drum and coupler | Flat faces down |

### B. Hardware to buy (see `flextome_bom.csv`)
Feeler-gauge set (0.20 mm and 1.0 mm blades used), 2 × 608ZZ bearings, 8 mm shaft 200 mm, 2 × M8 × 250 threaded rod + 4 nuts, M8 × 1.25 rod 50 mm + 1 nut, M8 × 1.0 bolt 30 mm + 2 fine nuts, 12 × M3 screws + heat-set inserts, 2 × M4 grub screws, M5 thumbscrew, M5 shoulder bolt, 2 × M4 set screws, 2 × M6 knob bolts + inserts, M6 × 40 bolt, one tension spring (~2 N/mm), one light compression spring, one box of disposable low-profile blades.

### C. Build sequence
1. **Frame prep.** Press heat-set inserts into the deck (2 × M6) and pivot clamp pads (6 × M3). Pass the two M8 × 250 tie rods through the base bores and tighten the end nuts hand-tight + ¼ turn (pre-compression against creep; also joins a split base).
2. **Pivot strips.** From the 0.20 mm feeler blade, shear three strips: 10 × 32 mm, 20 × 32 mm, 10 × 32 mm. Deburr edges with 600-grit paper. Lay each in its angled pocket on the tower head (6 mm clamped), fit the printed clamp plate, torque the M3 screws evenly. The 10 mm strips lie in the outer pockets (one orientation), the 20 mm strip in the centre pocket (opposite orientation), forming an X when viewed from the side.
3. **Arm.** Offer the arm root to the free strip ends; clamp with the three remaining plates. Check that the arm swings ±12° freely with no rubbing on the tower or gusset.
4. **Crank drive.** Press the 608 bearings into the pillow blocks, slide the 8 mm shaft through, fit the crank disc (M4 set screw) centred under the arm, fit the M5 shoulder bolt as crank pin, fit the handwheel hub (M4 set screw) and the handle. Connect the conrod: fork over the arm pin at 60 mm from the pivot, lower eye on the crank pin. Turn the handwheel: the arm must rock smoothly through one full cycle.
5. **Blade stage.** Press the M8 × 1.25 nut into the rear upstand pocket and the M8 × 1.0 nut into the hanging tab pocket. Thread the M8 × 1.0 bolt through the tab nut from the front; seat its head in the coupler's hex pocket. Thread the M8 × 1.25 rod through the frame nut into the coupler and lock with a jam nut. Slide the ratchet drum onto the rod behind the upstand, lock with jam nuts at both ends; fit the knob. Fit the pawl rod in its column (finger toward the drum, flag over the striker bar's path), with the light compression spring under the flag. Fit the selector shield onto the column strap.
6. **Blade.** Slide the 1.0 mm feeler strip (84 × 11 mm) into the front slot, then a blade behind it (edge up, protruding 2.5 mm above the strip). Tighten the two M4 grub screws from the rear until the blade cannot be moved by finger pressure. **Always handle blades with forceps.**
7. **Mount the stage.** Place the stage on the deck slots, fit the two M6 knob bolts loosely, fit the M6 coarse-approach screw in its side boss against the stage tab.
8. **Counterbalance.** Hook the tension spring between the arm's underside hook and the base so the arm rests near mid-stroke.
9. **Alignment.** Clamp a blank paraffin block in the cassette clamp. Back the fine feed to mid-travel (coupler centred in its gap). Wind the coarse screw until the block face is ~0.5 mm from the edge, lock the knob bolts. Set the selector to 6 teeth (30 µm) and trim until the full face cuts; then select 1 tooth (5 µm) and cut a ribbon.

### D. Safety
The blade is surgical-grade; fit and remove it with forceps, cover it with its guard when not sectioning, and keep fingers out of the arm's swept volume while the handwheel turns.
