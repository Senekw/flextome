"""Render the full CAD view set with render_cad (flat shading + feature edges) and the stroke animation."""
import numpy as np, render, pose
from PIL import Image
parts = render.load_parts()
render.render_cad(parts, "fig_render_iso_right.png", azim=-55, elev=25)
render.render_cad(parts, "fig_render_iso_left.png", azim=125, elev=22)
render.render_cad(parts, "fig_render_top.png", azim=-90, elev=89)
render.render_cad(parts, "fig_render_front.png", azim=-90, elev=0)
render.render_cad(parts, "fig_render_section.png", azim=90, elev=0, section_y=0.5, hide=("handwheel", "handle"))
# exploded view: lift groups in distinct directions
ex = {"arm": (0, 0, 70), "cassette_upper_jaw": (0, 0, 110), "cassette_thumbscrew": (0, 0, 130), "histology_cassette": (0, 0, 95), "paraffin_block": (0, 0, 95),
      "pivot_clamp_plates_arm": (0, 0, 70), "pivot_strips_steel": (0, 0, 40), "pivot_clamp_plates_frame": (0, 0, 20),
      "conrod": (0, -60, 30), "crank_disc": (0, -60, 0), "crank_pin": (0, -60, 0), "crankshaft": (0, -60, 0), "bearings_608": (0, -60, 0),
      "handwheel": (0, 70, 0), "handle": (0, 70, 0),
      "stage_body_flexure": (60, 0, 0), "blade_holder": (60, 0, 45), "pressure_strip_steel": (60, 0, 60), "blade_grub_screws_M4": (60, 0, 45), "blade_low_profile": (60, 0, 75),
      "feed_screw_M8x1_25": (120, 0, 0), "feed_screw_M8x1_0": (120, 0, 0), "screw_coupler": (120, 0, 0), "ratchet_wheel_50T": (120, 0, 0),
      "thickness_selector_shield": (120, 50, 0), "feed_knob": (120, 0, 0), "pawl_rod": (120, 50, 40), "stage_clamp_knobs_M6": (60, 0, 40),
      "fastener_heads_M3_M4": (0, 0, 30), "tie_rods_M8_steel_spine": (0, 0, -50), "coarse_approach_screw_M6": (60, -50, 0)}
render.render_cad(parts, "fig_render_exploded.png", azim=-50, elev=22, explode=ex, W=2200, H=1500)
# poses
for phi, name in ((90, "pose_top"), (0, "pose_mid"), (270, "pose_bottom")):
    pp, th = pose.posed(parts, phi)
    render.render_cad(pp, f"{name}.png", azim=-55, elev=22, title=f"crank {phi}°, arm θ = {th:+.1f}°", W=1500, H=1050)
# animation (24 frames, cutting stroke = decreasing theta)
frames = []
for phi in np.linspace(0, 360, 24, endpoint=False):
    pp, th = pose.posed(parts, phi)
    img, _ = render.render_cad(pp, "_f.png", azim=-55, elev=22, W=900, H=640, title=f"θ = {th:+.1f}°")
    frames.append(Image.fromarray((np.clip(img, 0, 1)*255).astype(np.uint8) if img.dtype != np.uint8 else img).convert("P", palette=Image.ADAPTIVE, colors=128))
frames[0].save("flextome_stroke_animation.gif", save_all=True, append_images=frames[1:], duration=70, loop=0)
print("renders done")
