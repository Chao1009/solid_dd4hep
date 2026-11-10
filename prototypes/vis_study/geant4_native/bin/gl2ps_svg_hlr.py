"""Turn a TSG/gl2ps SVG of an opaque surface+edges (hlhsr) scene into a vector hidden-line drawing.

gl2ps (BSP_SORT) writes all primitives back to front, so painting the surface triangles white and the
geometry edge polylines black hides the lines that lie behind surfaces (painter's algorithm).
Geometry edges are recognised by their stroke colour (the uniform colour given to all volumes with
/vis/geometry/set/colour all -1 R G B 1); other polylines (trajectories) keep their colours.

usage: gl2ps_svg_hlr.py IN.svg OUT.svg [--shaded] [--lw 0.6] [--edge-colour ededed] [--traj-on-top]
  --traj-on-top: move the non-geometry polylines (trajectories) to the end, i.e. draw them over everything
"""
import collections, re, sys
src, dst = sys.argv[1], sys.argv[2]
args = sys.argv[3:]
shaded = '--shaded' in args
lw = args[args.index('--lw') + 1] if '--lw' in args else '0.6'
s = open(src).read()
strokes = collections.Counter(re.findall(r'<polyline[^>]*?stroke="#([0-9a-fA-F]{6})"', s))
edge = args[args.index('--edge-colour') + 1] if '--edge-colour' in args else strokes.most_common(1)[0][0]
first_poly = s.find('<polygon')            # first polygon = gl2ps background rectangle
head, body = s[:first_poly], s[first_poly:]
bg_end = body.find('/>') + 2
bg, body = body[:bg_end], body[bg_end:]
if shaded:   # keep the shading, hairline stroke in the fill colour hides triangle seams
    body = re.sub(r'<polygon fill="#([0-9a-fA-F]{6})"', r'<polygon fill="#\1" stroke="#\1" stroke-width="0.25"', body)
else:
    body = re.sub(r'<polygon fill="#[0-9a-fA-F]{6}"', '<polygon fill="#ffffff" stroke="#ffffff" stroke-width="0.25"', body)
body = re.sub(r'(<polyline[^>]*?)stroke="#' + edge + r'" stroke-width="[0-9.]+"',
              lambda m: m.group(1) + f'stroke="#000000" stroke-width="{lw}"', body)
body = body.replace(' shape-rendering="crispEdges"', '')
if '--traj-on-top' in args:
    end = body.rfind('</g>')
    pat = re.compile(r'<polyline(?![^>]*stroke="#000000")[^>]*/>\n?')
    trajs = pat.findall(body[:end])
    body = pat.sub('', body[:end]) + ''.join(trajs) + body[end:]
open(dst, 'w').write(head + bg + body)
print(f'{dst}: {len(head + bg + body) // 1024} KiB, edge colour #{edge}, '
      f'{s.count("<polygon")} polygons, {s.count("<polyline")} polylines')
