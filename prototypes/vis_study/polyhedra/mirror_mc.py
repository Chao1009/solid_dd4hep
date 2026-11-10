import sys, numpy as np, ROOT, dd4hep
desc = ROOT.dd4hep.Detector.getInstance(); desc.fromXML(sys.argv[1])
shape = None
for v in ROOT.gGeoManager.GetListOfVolumes():
    if v.GetName().startswith(sys.argv[2]):
        shape = v.GetShape(); break
print(shape.GetName(), shape.ClassName())
lo = np.array([85.0, -24.0, 12.0]); hi = np.array([222.0, 24.0, 50.0])   # cm, generous box around both meshes
rng = np.random.default_rng(1)
found_lo = np.full(3, 1e9); found_hi = np.full(3, -1e9); nin = 0
N = 2_000_000
for it in range(10):
    P = (lo + (hi - lo) * rng.random((N, 3))).astype(np.float64)
    flat = np.ascontiguousarray(P.reshape(-1))
    inside = np.zeros(N, dtype=np.bool_)
    shape.Contains_v(flat, inside, N)
    Q = P[inside]
    nin += len(Q)
    if len(Q):
        found_lo = np.minimum(found_lo, Q.min(0)); found_hi = np.maximum(found_hi, Q.max(0))
vol = nin / (10 * N) * np.prod(hi - lo)
print('MC inside %d, volume %.1f cm3; extent min %s max %s [mm]' % (nin, vol, (found_lo * 10).round(1), (found_hi * 10).round(1)))
