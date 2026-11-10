"""imgdiff.py A.png B.png [diff.png]: pixel statistics + diff image (red = only in A, blue = only in B)."""
import sys, re
import numpy as np
from PIL import Image
a = np.asarray(Image.open(sys.argv[1]).convert('RGB')).astype(int)
b = np.asarray(Image.open(sys.argv[2]).convert('RGB')).astype(int)
if a.shape != b.shape:
    print('shape differs', a.shape, b.shape)
    h = min(a.shape[0], b.shape[0]); w = min(a.shape[1], b.shape[1]); a = a[:h, :w]; b = b[:h, :w]
d = np.abs(a - b).max(2)
ink_a = a.min(2) < 200; ink_b = b.min(2) < 200
n = d.size
print('differing pixels (>8/255): %d of %d (%.4f%%); ink A %d, ink B %d, ink only A %d, only B %d'
      % ((d > 8).sum(), n, 100.0 * (d > 8).sum() / n, ink_a.sum(), ink_b.sum(), (ink_a & ~ink_b).sum(), (ink_b & ~ink_a).sum()))
if len(sys.argv) > 3:
    out = np.full(a.shape, 255, np.uint8)
    both = ink_a & ink_b
    out[both] = (170, 170, 170)
    out[ink_a & ~ink_b] = (230, 0, 0)
    out[ink_b & ~ink_a] = (0, 0, 230)
    im = Image.fromarray(out)
    im.thumbnail((1500, 1500))
    im.save(sys.argv[3], optimize=True)
