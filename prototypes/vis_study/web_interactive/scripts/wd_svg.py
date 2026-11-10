import sys, json, time
sys.path.insert(0, sys.argv[0].rsplit('/', 1)[0])
from wd import start_driver
url, out = sys.argv[1], sys.argv[2]
t0 = time.time()
p, wd = start_driver()
try:
    wd.new_session(1200, 750)
    wd.set_rect(1200, 750)
    wd.get(url)
    for _ in range(240):
        if wd.js('return window.__done === true'):
            break
        time.sleep(0.5)
    time.sleep(2)
    st = wd.js('return window.__status')
    svg = wd.js('const s=document.querySelector("#d svg svg") || document.querySelector("#d svg"); return s ? s.outerHTML : null')
    print('status', st, 'svg chars', len(svg) if svg else None, f'{time.time()-t0:.1f}s')
    if svg:
        open(out, 'w').write(svg)
    wd.screenshot(out.replace('.svg', '_screen.png'))
finally:
    wd.quit(); p.terminate()
