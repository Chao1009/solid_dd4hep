"""usage: probe_fetch.py PAGE_URL TARGET_URL [prefs-json]  - can a page at PAGE_URL fetch TARGET_URL? (15 s JS timeout)"""
import json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from wd import start_driver
page, target = sys.argv[1], sys.argv[2]
prefs = json.loads(sys.argv[3]) if len(sys.argv) > 3 else {}
p, wd = start_driver()
try:
    wd.new_session(800, 600, prefs); wd.get(page); time.sleep(3)
    r = wd.js_async("""const done = arguments[arguments.length-1];
      const t = new Promise(res => setTimeout(() => res('timeout 15 s (no answer)'), 15000));
      const f = fetch(arguments[0]).then(r => r.text().then(x => 'status ' + r.status + ' len ' + x.length)).catch(e => 'error: ' + e);
      Promise.race([f, t]).then(done);""", [target])
    print('prefs', prefs, 'fetch from', page, '->', r)
finally:
    wd.quit(); p.terminate()
