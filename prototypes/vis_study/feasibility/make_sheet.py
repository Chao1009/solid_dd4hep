"""Overview sheet of the best images of the SoLID visualization study (3x3 panels)."""
import os
import sys
from PIL import Image, ImageDraw, ImageFont

V = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(V, "feasibility", "feasibility_overview.png")
PANELS = [
    ("polyhedra/img/g4_PVDIS/cutaway.png", None,
     "Q1 (a) no DAWNFILE: Geant4 polyhedra dump -> solid_lineview, PVDIS (0 px diff vs DAWNFILE)"),
    ("polyhedra/img/g4_JPsi_lite/section_yz.png", None,
     "Q1 (b) same route, SoLID-Lite J/psi x=0 section (0.064 % px diff, vertex-on-plane)"),
    ("geant4_native/samples/02_vector_hlr_vs_lineview.png", (900, 0, 1800, 560),
     "Q1 (c) Geant4-native TSG_OFFSCREEN gl2ps + recolour: close, 4.3 MB SVG"),
    ("events_animation/out/jpsi1_side.png", None,
     "Q2 (d) event overlay on solid_lineview: tracks, hits, module energy (log)"),
    ("events_animation/out/jpsi1_cutaway.png", None,
     "Q2 (e) same event, hidden-line cutaway at the event plane"),
    ("events_animation/out/jpsi0_beam.png", None,
     "Q2 (f) beam view through the LAEC, modules coloured by energy"),
    ("events_animation/out/frames_jpsi1_side/f0070.png", None,
     "Q3 (g) animation frame t = 17 ns (0.04-0.2 s/frame; GIF/WebP/WebM/SVG)"),
    ("geant4_native/samples/03_tsg_zb_cutaway_tracks.png", (60, 0, 1600, 1000),
     "Q2/Q3 (h) Geant4-native TSG colour render + trajectories"),
    ("web_interactive/shots2/jsroot_ev2_side_ortho.png", None,
     "Q2 (i) JSROOT in the browser: ortho view, module energy colours"),
]
CW, CH, CAP, GAP = 780, 470, 30, 12
try:
    font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 15)
except OSError:
    font = ImageFont.load_default()
W = 3 * CW + 4 * GAP
H = 3 * (CH + CAP) + 4 * GAP
sheet = Image.new("RGB", (W, H), "white")
draw = ImageDraw.Draw(sheet)
for i, (rel, crop, cap) in enumerate(PANELS):
    r, c = divmod(i, 3)
    x0 = GAP + c * (CW + GAP)
    y0 = GAP + r * (CH + CAP + GAP)
    im = Image.open(os.path.join(V, rel))
    if im.mode in ("RGBA", "LA", "P"):
        bg = Image.new("RGB", im.size, "white")
        im = im.convert("RGBA")
        bg.paste(im, mask=im.split()[3])
        im = bg
    else:
        im = im.convert("RGB")
    if crop:
        im = im.crop(crop)
    im.thumbnail((CW - 8, CH - 8), Image.LANCZOS)
    sheet.paste(im, (x0 + (CW - im.width) // 2, y0 + CAP + (CH - im.height) // 2))
    draw.rectangle([x0, y0 + CAP, x0 + CW - 1, y0 + CAP + CH - 1], outline=(190, 190, 190))
    draw.text((x0 + 4, y0 + 7), cap, fill=(20, 20, 20), font=font)
rgb = sheet
sheet = rgb.quantize(colors=256, method=Image.Quantize.FASTOCTREE, dither=Image.Dither.FLOYDSTEINBERG)
sheet.save(OUT, optimize=True)
print(OUT, os.path.getsize(OUT), sheet.size)
