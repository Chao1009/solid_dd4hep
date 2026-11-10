#!/usr/bin/env python3
"""Encode a directory of PNG frames (f0000.png ...) into GIF / WebP (Pillow) and MP4 / WebM (GStreamer).

usage: encode.py FRAMEDIR OUTBASE [--fps 12] [--hold 2.0] [--formats gif,webp,mp4,webm]
ffmpeg and imageio are not available (host or eic_xl 26.10); GStreamer 1.20 on the host has vp8/vp9 encoders
and an MP4 muxer, but no H.264 encoder.
"""
import argparse
import glob
import os
import shutil
import subprocess
import tempfile
import time

from PIL import Image

ap = argparse.ArgumentParser()
ap.add_argument("framedir")
ap.add_argument("outbase")
ap.add_argument("--fps", type=float, default=12)
ap.add_argument("--hold", type=float, default=2.0, help="seconds to hold the last frame")
ap.add_argument("--formats", default="gif,webp,mp4,webm")
ap.add_argument("--gif-colors", type=int, default=256)
ap.add_argument("--gif-scale", type=float, default=1.0)
ap.add_argument("--webp-lossless", action="store_true")
a = ap.parse_args()

files = sorted(glob.glob(os.path.join(a.framedir, "f*.png")))
dt = 1000.0 / a.fps
durs = [int(round(dt))] * len(files)
durs[-1] += int(a.hold * 1000)
res = {}

if "gif" in a.formats:
    t0 = time.perf_counter()
    import numpy as np
    # one global palette (no flicker): the most frequent exact colours of every 8th frame + the last frames.
    # The drawings are flat fills + anti-aliased lines, so the frequent colours are the fills, the line
    # colours and the main anti-aliasing shades; the rest map to the nearest entry.
    sample = sorted(set(files[::8] + files[-3:]))
    px = np.concatenate([np.asarray(Image.open(f).convert("RGB")).reshape(-1, 3) for f in sample])
    key = (px[:, 0].astype(np.int64) << 16) | (px[:, 1].astype(np.int64) << 8) | px[:, 2]
    u, c = np.unique(key, return_counts=True)
    # reserved entries: the particle colours and their blends with white (thin anti-aliased lines are rare
    # colours and would otherwise snap to a wrong hue)
    keyc = []
    for h in ("#2a78d6", "#eb6834", "#1baf7a", "#52514e", "#898781", "#0b0b0b"):
        rgb = np.array([int(h[i:i + 2], 16) for i in (1, 3, 5)])
        for w in (0.0, 0.35, 0.6, 0.8):
            keyc.append(np.round(rgb * (1 - w) + 255 * w).astype(np.int64))
    # the energy colour ramp (ColorBrewer Reds steps 1..8 as in evdisp.energy_color), 40 levels
    reds = np.array([[int(h[i:i + 2], 16) for i in (1, 3, 5)] for h in
                     ("#fff5f0", "#fee0d2", "#fcbba1", "#fc9272", "#fb6a4a", "#ef3b2c", "#cb181d", "#a50f15",
                      "#67000d")], float)
    for uu in np.linspace(0, 1, 40):
        x = 1 + uu * (len(reds) - 2)
        i = min(int(x), len(reds) - 2)
        keyc.append(np.round(reds[i] * (1 - (x - i)) + reds[i + 1] * (x - i)).astype(np.int64))
    keyc = np.unique(np.array(keyc), axis=0)
    kk = (keyc[:, 0] << 16) | (keyc[:, 1] << 8) | keyc[:, 2]
    order = u[np.argsort(-c, kind="stable")]
    order = order[~np.isin(order, kk)]
    top = np.concatenate([kk, order[:min(a.gif_colors, 255) - len(kk)]])     # index 255: transparent
    pal_rgb = np.stack([(top >> 16) & 255, (top >> 8) & 255, top & 255], 1).astype(np.uint8)
    flat = list(pal_rgb.reshape(-1)) + [0] * (3 * (256 - len(pal_rgb)))
    P = pal_rgb.astype(np.int64)
    ims = []
    prev = None
    for f in files:          # exact nearest-palette mapping (Pillow's quantize uses a coarse colour cache)
        im = Image.open(f).convert("RGB")
        if a.gif_scale != 1.0:
            im = im.resize((int(im.width * a.gif_scale), int(im.height * a.gif_scale)), Image.LANCZOS)
        arr = np.asarray(im).reshape(-1, 3).astype(np.int64)
        k = (arr[:, 0] << 16) | (arr[:, 1] << 8) | arr[:, 2]
        uk, inv = np.unique(k, return_inverse=True)
        uc = np.stack([(uk >> 16) & 255, (uk >> 8) & 255, uk & 255], 1)
        d = ((uc[:, None, :] - P[None, :, :]) ** 2).sum(2)
        idx = d.argmin(1).astype(np.uint8)[inv.ravel()].reshape(im.height, im.width)
        out = idx.copy()
        if prev is not None:             # unchanged pixels -> transparent: LZW packs the static scene away
            out[idx == prev] = 255
        prev = idx
        pim = Image.fromarray(out, mode="P")
        pim.putpalette(flat)
        ims.append(pim)
    p = a.outbase + ".gif"
    ims[0].save(p, save_all=True, append_images=ims[1:], duration=durs, loop=0, optimize=False, disposal=1,
                transparency=255)
    res["gif"] = (p, time.perf_counter() - t0)

if "webp" in a.formats:
    t0 = time.perf_counter()
    ims = [Image.open(f).convert("RGB") for f in files]
    p = a.outbase + ".webp"
    ims[0].save(p, save_all=True, append_images=ims[1:], duration=durs, loop=0, quality=90, method=4,
                lossless=a.webp_lossless)
    res["webp"] = (p, time.perf_counter() - t0)


def gst(muxer, ext, enc):
    # hold the last frame: duplicate it in a temp copy of the sequence
    t0 = time.perf_counter()
    tmp = tempfile.mkdtemp(prefix="seq_", dir=os.path.dirname(os.path.abspath(a.outbase)))
    n = len(files)
    extra = int(round(a.hold * a.fps))
    for i in range(n + extra):
        os.symlink(os.path.abspath(files[min(i, n - 1)]), os.path.join(tmp, "f%04d.png" % i))
    p = a.outbase + ext
    fr = "%d/1" % int(a.fps)
    cmd = ["gst-launch-1.0", "-q", "multifilesrc", "location=%s/f%%04d.png" % tmp, "index=0",
           "caps=image/png,framerate=%s" % fr, "!", "pngdec", "!", "videoconvert", "!",
           "video/x-raw,format=I420", "!"] + enc + ["!", muxer, "!", "filesink", "location=%s" % p]
    r = subprocess.run(cmd, capture_output=True, text=True)
    shutil.rmtree(tmp)
    if r.returncode != 0 or not os.path.exists(p) or os.path.getsize(p) == 0:
        print("gstreamer %s failed: %s" % (ext, (r.stderr or r.stdout)[-400:]))
        return None
    return p, time.perf_counter() - t0


VP9 = ["vp9enc", "deadline=0", "cpu-used=4", "end-usage=cq", "cq-level=30", "threads=2",
       "keyframe-max-dist=60"]
if "mp4" in a.formats:
    # MP4 needs H.264 / H.265 / AV1 / MPEG-4 part 2 (mp4mux 1.20 does not take VP8/VP9); try what may exist
    MP4_ENCODERS = [("x264enc", ["x264enc", "speed-preset=medium", "tune=animation", "pass=quant", "quantizer=20"]),
                    ("openh264enc", ["openh264enc", "bitrate=4000000"]),
                    ("avenc_h264_omx", ["avenc_h264_omx"]),
                    ("avenc_mpeg4", ["avenc_mpeg4", "bitrate=4000000"]),
                    ("svtav1enc", ["svtav1enc"]), ("rav1enc", ["rav1enc"]), ("av1enc", ["av1enc"])]
    found = [(n, e) for n, e in MP4_ENCODERS
             if subprocess.run(["gst-inspect-1.0", n], capture_output=True).returncode == 0]
    if found:
        r = gst("mp4mux", ".mp4", found[0][1])
        if r:
            res["mp4"] = r
    else:
        print("mp4: no encoder element among %s; ffmpeg not installed -> MP4 skipped "
              "(e.g. ffmpeg -framerate %g -i %s/f%%04d.png -c:v libx264 -pix_fmt yuv420p -crf 20 %s.mp4)"
              % (", ".join(n for n, _ in MP4_ENCODERS), a.fps, a.framedir, a.outbase))
        r = gst("mp4mux", ".mp4", VP9)       # documents the failure: mp4mux 1.20 has no video/x-vp9 sink caps
        if r:
            res["mp4"] = r
if "webm" in a.formats:
    r = gst("webmmux", ".webm", VP9)
    if r:
        res["webm"] = r

w, h = Image.open(files[0]).size
print("%d frames %dx%d at %g fps (+%.1f s hold)" % (len(files), w, h, a.fps, a.hold))
for k, (p, t) in res.items():
    print("%-5s %8.1f kB  %6.2f s  %s" % (k, os.path.getsize(p) / 1024, t, p))
