"""Bill of materials for FlexTome (rev. 2).  Printed-part lines come from parts.json (solid-modelled volume, PLA 1.24 g/cm3,
US$20/kg, an upper bound since infill is <100 %); hardware lines are typed here with indicative retail unit prices.
Run after flextome_cad.py:  python make_bom.py  -> flextome_bom.csv, bom_results.json"""
import json, csv
RHO, PRICE = 1.24, 20.0/1000      # g/cm3, USD/g
meta = json.load(open("parts.json"))
rows = []
for m in meta:
    if not m["printed"]: continue
    g = m["volume_cm3"]*RHO
    rows.append(dict(item=m["name"], category="3D-printed (PLA)", qty=m["qty"], spec=f"{m['volume_cm3']:.2f} cm3 modelled solid",
                     unit_cost_usd=round(g*PRICE, 2), source="filament @ $20/kg", mass_g=round(g, 1)))
HW = [  # item, category, qty, spec, unit USD, source
 ("spring-steel feeler gauge strips 0.20 mm (from a 32-blade set)", "hardware", 1, "set; 3 strips cut to 10/20/10 x 32 mm", 6.00, "hardware store / online"),
 ("steel clamp strip 1.0 mm (feeler gauge blade, 84 x 14 mm, 2 holes)", "hardware", 1, "from the same set", 0.00, "included above"),
 ("608ZZ ball bearing", "hardware", 2, "8 x 22 x 7 mm", 1.00, "skate / 3D-printer supplier"),
 ("8 mm steel shaft stubs (2 x 60 mm)", "hardware", 1, "smooth rod, cut", 3.00, "3D-printer supplier"),
 ("M8 x 1.25 threaded rod, 250 mm (base tie rods)", "hardware", 2, "zinc steel", 1.50, "hardware store"),
 ("M8 x 1.25 threaded rod, 55 mm (feed screw, rear)", "hardware", 1, "cut from rod", 0.30, "hardware store"),
 ("M8 x 1.0 fine-pitch hex bolt, 25 mm (feed screw, front)", "hardware", 1, "class 8.8", 1.20, "automotive / online"),
 ("M8 x 1.25 nut", "hardware", 5, "frame nut + 4 tie-rod nuts", 0.10, "hardware store"),
 ("M8 x 1.25 thin nut (coupler, captured)", "hardware", 1, "4 mm", 0.15, "hardware store"),
 ("M8 x 1.0 fine nut (stage tab)", "hardware", 1, "", 0.40, "automotive / online"),
 ("M8 steel washer, 16 mm OD (under tab nut and frame nut)", "hardware", 2, "form A", 0.05, "hardware store"),
 ("disposable low-profile microtome blade", "consumable", 1, "80 x 8 x 0.25 mm (any brand)", 0.60, "lab supplier ($30/50)"),
 ("M3 x 8 screws + heat-set inserts (pivot clamp pads, selector strap)", "hardware", 13, "", 0.10, "online"),
 ("M4 x 12 cap screws + heat-set inserts (blade clamp strip)", "hardware", 2, "", 0.15, "online"),
 ("M4 x 12 cap screws + inserts (blade holder to platform)", "hardware", 4, "", 0.15, "online"),
 ("M4 set screws (crank webs on shaft stubs, handwheel hub)", "hardware", 3, "", 0.15, "online"),
 ("M5 x 30 thumbscrew (cassette clamp)", "hardware", 1, "", 0.80, "online"),
 ("M5 shoulder bolt / crank pin", "hardware", 1, "5 x 30 mm", 0.80, "online"),
 ("M6 knob bolts + heat-set inserts (stage clamp)", "hardware", 2, "", 1.20, "online"),
 ("M6 x 40 bolt + printed knob (coarse approach)", "hardware", 1, "", 0.40, "hardware store"),
 ("extension springs 5 x 20 mm, ~5 N (stage preload)", "hardware", 2, "", 0.40, "hardware store"),
 ("compression spring (pawl rod return, light)", "hardware", 1, "", 0.30, "hardware store"),
 ("cyanoacrylate / epoxy, PTFE grease, thread locker", "consumable", 1, "", 2.00, "hardware store"),
]
for it, cat, q, spec, c, src in HW:
    rows.append(dict(item=it, category=cat, qty=q, spec=spec, unit_cost_usd=c, source=src, mass_g=""))
for r in rows: r["line_cost_usd"] = round(r["qty"]*r["unit_cost_usd"], 2)
with open("flextome_bom.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
by_cat = {}
for r in rows: by_cat[r["category"]] = round(by_cat.get(r["category"], 0) + r["line_cost_usd"], 2)
printed = [r for r in rows if r["category"].startswith("3D")]
out = dict(total_usd=round(sum(r["line_cost_usd"] for r in rows), 2), by_category=by_cat,
           printed_part_types=len(printed), printed_pieces=sum(r["qty"] for r in printed),
           printed_mass_g=round(sum(r["mass_g"]*r["qty"] for r in printed), 1),
           printed_volume_L=round(sum(m["volume_cm3"]*m["qty"] for m in meta if m["printed"])/1000, 2),
           n_components=len(rows))
json.dump(out, open("bom_results.json", "w"), indent=1)
print(out)
