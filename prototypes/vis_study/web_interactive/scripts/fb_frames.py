#!/usr/bin/env python3
"""Drive the Firebird frontend in headless Firefox: load geometry + DEX event, pick a camera preset,
then press the time 'Step' button N times and grab a PNG per step (time animation frames).

usage: fb_frames.py URL OUTDIR [--frames 30] [--camera 'Camera to right'] [--w 1400 --h 860]
"""
import argparse, json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from wd import start_driver  # noqa: E402

CLICK = """
const [label, text] = arguments;
const b = [...document.querySelectorAll('button')].find(e =>
   (label === null || (e.getAttribute('aria-label') || e.title || e.getAttribute('mattooltip') || '') === label) &&
   (text === null || (e.innerText || '').trim() === text));
if (!b) return 'not found: ' + label + '/' + text;
b.click(); return 'ok';
"""
TIME = "const i=[...document.querySelectorAll('input[type=text]')].find(e=>/^-?[0-9.]+$/.test(e.value)); return i?i.value:null;"

ap = argparse.ArgumentParser()
ap.add_argument('url'); ap.add_argument('outdir')
ap.add_argument('--frames', type=int, default=30)
ap.add_argument('--camera', default=None, help="menu text of a camera preset, e.g. 'Camera to right'")
ap.add_argument('--w', type=int, default=1400); ap.add_argument('--h', type=int, default=860)
ap.add_argument('--step-sleep', type=float, default=0.4)
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)
t0 = time.time()
p, wd = start_driver()
log = []
try:
    wd.new_session(a.w, a.h); wd.set_rect(a.w, a.h); wd.get(a.url)
    while time.time() - t0 < 180 and not wd.js('return !!(window.firebird && window.firebird.ready)'):
        time.sleep(0.5)
    print(f'ready after {time.time() - t0:.1f}s', wd.js('return window.firebird'))
    time.sleep(2)
    if a.camera:
        print('camera menu', wd.js(CLICK, ['Toggle camera controls', 'videocam'])); time.sleep(0.5)
        print('camera', wd.js(CLICK, [None, a.camera])); time.sleep(1.5)
        wd.js("document.dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape'}));"); time.sleep(0.5)
    wd.screenshot(os.path.join(a.outdir, 'start.png'))
    print('rewind', wd.js(CLICK, ['Rewind', 'first_page'])); time.sleep(1)
    tf = time.time()
    for i in range(a.frames):
        r = wd.js(CLICK, ['Step', 'chevron_right'])
        time.sleep(a.step_sleep)
        fn = os.path.join(a.outdir, f'f{i:03d}.png')
        wd.screenshot(fn)
        log.append((i, wd.js(TIME), r))
    print(f'{a.frames} frames in {time.time() - tf:.1f}s; total {time.time() - t0:.1f}s')
    print('time values:', [x[1] for x in log])
finally:
    wd.quit(); p.terminate()
