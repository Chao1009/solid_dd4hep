#!/usr/bin/env python3
"""Minimal W3C WebDriver client for headless Firefox (snap firefox.geckodriver).

Usage examples:
  wd.py shot URL out.png [--w 1600 --h 1000] [--wait-js 'return !!window.x'] [--timeout 120]
           [--sleep 5] [--js 'script run after ready'] [--log console.txt]
"""
import argparse, base64, json, os, socket, subprocess, sys, time, urllib.request, urllib.error


def free_port():
    s = socket.socket(); s.bind(("127.0.0.1", 0)); p = s.getsockname()[1]; s.close(); return p


class WD:
    def __init__(self, port):
        self.base = f"http://127.0.0.1:{port}"
        self.sid = None

    def req(self, method, path, body=None, timeout=300):
        data = None if body is None else json.dumps(body).encode()
        r = urllib.request.Request(self.base + path, data=data, method=method,
                                   headers={"Content-Type": "application/json"})
        try:
            with urllib.request.urlopen(r, timeout=timeout) as f:
                return json.loads(f.read().decode())
        except urllib.error.HTTPError as e:
            msg = e.read().decode()
            raise RuntimeError(f"{method} {path}: {e.code} {msg[:2000]}")

    def new_session(self, w, h, prefs=None):
        caps = {"capabilities": {"alwaysMatch": {
            "browserName": "firefox",
            "moz:firefoxOptions": {"args": ["-headless", f"--width={w}", f"--height={h}"],
                                   "prefs": prefs or {}},
            "pageLoadStrategy": "normal",
            "timeouts": {"script": 600000, "pageLoad": 300000}}}}
        r = self.req("POST", "/session", caps)
        self.sid = r["value"]["sessionId"]
        return r["value"]["capabilities"]

    def s(self, path):
        return f"/session/{self.sid}{path}"

    def set_rect(self, w, h):
        return self.req("POST", self.s("/window/rect"), {"width": w, "height": h})

    def get(self, url):
        return self.req("POST", self.s("/url"), {"url": url})

    def js(self, script, args=None):
        return self.req("POST", self.s("/execute/sync"), {"script": script, "args": args or []})["value"]

    def js_async(self, script, args=None):
        return self.req("POST", self.s("/execute/async"), {"script": script, "args": args or []})["value"]

    def screenshot(self, out):
        b64 = self.req("GET", self.s("/screenshot"))["value"]
        with open(out, "wb") as f:
            f.write(base64.b64decode(b64))

    def quit(self):
        try:
            self.req("DELETE", f"/session/{self.sid}", timeout=30)
        except Exception as e:
            print("quit:", e, file=sys.stderr)


def start_driver():
    port = free_port()
    log = open(os.environ.get("WD_LOG", os.devnull), "w")
    p = subprocess.Popen(["firefox.geckodriver", "--port", str(port), "--log", "warn"],
                         stdout=log, stderr=log)
    for _ in range(100):
        try:
            urllib.request.urlopen(f"http://127.0.0.1:{port}/status", timeout=1)
            break
        except Exception:
            time.sleep(0.2)
    return p, WD(port)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["shot"])
    ap.add_argument("url")
    ap.add_argument("out")
    ap.add_argument("--w", type=int, default=1600)
    ap.add_argument("--h", type=int, default=1000)
    ap.add_argument("--wait-js", default=None, help="JS returning truthy when ready")
    ap.add_argument("--timeout", type=float, default=180)
    ap.add_argument("--sleep", type=float, default=2.0, help="extra settle time after ready")
    ap.add_argument("--js", action="append", default=[], help="JS to run after ready (repeatable)")
    ap.add_argument("--js-sleep", type=float, default=1.0)
    ap.add_argument("--report-js", default=None, help="JS whose return value is printed at the end")
    ap.add_argument("--prefs", default=None, help="JSON dict of firefox prefs")
    ap.add_argument("--frames", type=int, default=0, help="also grab N frames, running --frame-js before each")
    ap.add_argument("--frame-js", default=None)
    ap.add_argument("--frame-pattern", default=None, help="printf pattern for frame files")
    ap.add_argument("--frame-sleep", type=float, default=0.3)
    a = ap.parse_args()
    prefs = json.loads(a.prefs) if a.prefs else {}
    t0 = time.time()
    p, wd = start_driver()
    try:
        caps = wd.new_session(a.w, a.h, prefs)
        print("browser", caps.get("browserName"), caps.get("browserVersion"), f"{time.time()-t0:.1f}s")
        wd.set_rect(a.w, a.h)
        wd.get(a.url)
        t1 = time.time()
        ok = True
        if a.wait_js:
            ok = False
            while time.time() - t1 < a.timeout:
                try:
                    if wd.js(a.wait_js):
                        ok = True
                        break
                except Exception as e:
                    pass
                time.sleep(0.5)
        print(f"ready={ok} after {time.time()-t1:.1f}s")
        time.sleep(a.sleep)
        for s in a.js:
            r = wd.js(s)
            print("js ->", json.dumps(r)[:500] if r is not None else None)
            time.sleep(a.js_sleep)
        wd.screenshot(a.out)
        print("saved", a.out, f"total {time.time()-t0:.1f}s")
        for i in range(a.frames):
            if a.frame_js:
                wd.js(a.frame_js, [i])
            time.sleep(a.frame_sleep)
            wd.screenshot(a.frame_pattern % i)
        if a.frames:
            print("frames", a.frames, f"total {time.time()-t0:.1f}s")
        if a.report_js:
            r = wd.js(a.report_js)
            print("report:", json.dumps(r, indent=1)[:6000])
    finally:
        wd.quit()
        p.terminate()
        try:
            p.wait(10)
        except Exception:
            p.kill()


if __name__ == "__main__":
    main()
