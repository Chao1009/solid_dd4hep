import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tgeopoly
ROOT = tgeopoly.cling()
mgr = ROOT.TGeoManager('t', 't')
shapes = {
  'Tube(360)': ROOT.TGeoTube(10, 20, 5),
  'TubeSeg(6deg)': ROOT.TGeoTubeSeg(10, 20, 5, 0, 6),
  'TubeSeg(360)': ROOT.TGeoTubeSeg(10, 20, 5, 0, 360),
  'Sphere(th0-90)': ROOT.TGeoSphere(280, 283, 0, 90, 0, 360),
  'Torus': ROOT.TGeoTorus(100, 0, 10, 0, 360),
}
x = np.array([0, 10, 10, 6, 6, 0], float); y = np.array([0, 0, 10, 10, 4, 4], float)   # concave L
xt = ROOT.TGeoXtru(2); xt.DefinePolygon(6, x, y); xt.DefineSection(0, -1); xt.DefineSection(1, 1)
shapes['Xtru(L, concave)'] = xt
for nseg in (20, 72):
    for k, s in shapes.items():
        m = ROOT.solidtgeo.tessellate(s, nseg)
        pol = np.frombuffer(m.pol.data(), dtype=np.int32, count=m.pol.size())
        sizes = []
        i = 0
        while i < len(pol):
            sizes.append(int(pol[i])); i += 1 + pol[i]
        print('nseg %2d  %-18s nv %5d  npol %5d  max polygon size %3d' % (nseg, k, m.nv, m.npol, max(sizes) if sizes else 0))
