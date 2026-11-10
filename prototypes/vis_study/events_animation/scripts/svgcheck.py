#!/usr/bin/env python3
"""Render an animated (SMIL) SVG from evdisp.py at fixed event times with headless Firefox.

usage: svgcheck.py ANIM.svg OUTDIR T1_ns T2_ns ...

The SVG is inlined in an HTML page whose script pauses the SMIL timeline and seeks it to
t_ns * duration / t_end (data-t-end-ns / data-duration-s on the <svg> root); Firefox writes a screenshot.
Firefox is the Ubuntu snap here: it cannot see the host /tmp, so page, profile and screenshot live in a
temporary directory under ~/snap/firefox/common (removed afterwards) and the PNG is moved to OUTDIR.
"""
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time

svg, outdir = sys.argv[1], sys.argv[2]
times = [float(x) for x in sys.argv[3:]]
os.makedirs(outdir, exist_ok=True)
txt = open(svg, encoding="utf-8").read()
t_end = float(re.search(r'data-t-end-ns="([\d.]+)"', txt).group(1))
dur = float(re.search(r'data-duration-s="([\d.]+)"', txt).group(1))
W, H = [int(v) for v in re.search(r'viewBox="0 0 (\d+) (\d+)"', txt).groups()]
width = 1600
height = int(round(H * width / W))
body = txt.split("?>", 1)[1] if txt.startswith("<?xml") else txt
body = re.sub(r'width="[\d.]+pt" height="[\d.]+pt"', 'width="%d" height="%d"' % (width, height), body, count=1)
base = os.path.splitext(os.path.basename(svg))[0]
snapdir = os.path.expanduser("~/snap/firefox/common")
tmp = tempfile.mkdtemp(prefix="svgcheck_", dir=snapdir if os.path.isdir(snapdir) else None)
prof = os.path.join(tmp, "profile")
os.makedirs(prof)
try:
    for t in times:
        ts = t * dur / t_end
        html = os.path.join(tmp, "page.html")
        png_tmp = os.path.join(tmp, "shot.png")
        png = os.path.join(outdir, "%s_t%.1f.png" % (base, t))
        with open(html, "w", encoding="utf-8") as fh:
            fh.write('<!doctype html><html><head><meta charset="utf-8"><title>svg check</title></head>'
                     '<body style="margin:0;background:#fff">' + body +
                     '<script>const s=document.querySelector("svg");s.pauseAnimations();s.setCurrentTime(%.4f);'
                     '</script></body></html>' % ts)
        t0 = time.time()
        for attempt in range(3):     # headless Firefox occasionally hangs (another Firefox instance running)
            try:
                r = subprocess.run(["firefox", "--headless", "--no-remote", "--new-instance", "--profile", prof,
                                    "--screenshot", png_tmp, "--window-size=%d,%d" % (width, height),
                                    "file://" + html], capture_output=True, text=True, timeout=60)
                break
            except subprocess.TimeoutExpired as exc:
                print("  attempt %d timed out" % (attempt + 1), flush=True)
                r = exc
                shutil.rmtree(prof, ignore_errors=True)
                os.makedirs(prof)
        ok = os.path.exists(png_tmp)
        if ok:
            shutil.move(png_tmp, png)
        print("t = %5.1f ns -> SMIL clock %.3f s: %s (%.1f s)%s" % (
            t, ts, png if ok else "FAILED", time.time() - t0, "" if ok else " " + str(getattr(r, "stderr", ""))[-300:]), flush=True)
finally:
    shutil.rmtree(tmp, ignore_errors=True)
