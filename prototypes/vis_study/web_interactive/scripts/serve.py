#!/usr/bin/env python3
"""Tiny static web server for the display tests: URL-prefix mounts, HTTP Range (needed by JSROOT and
Firebird's in-browser ROOT reader), CORS, and an SPA fallback to index.html (Angular routes).

usage: serve.py PORT /=DIR [/solid/=DIR ...] [--spa]
"""
import http.server, mimetypes, os, re, sys, urllib.parse

PORT = int(sys.argv[1])
MOUNTS = []
SPA = '--spa' in sys.argv
for arg in sys.argv[2:]:
    if '=' in arg:
        pre, d = arg.split('=', 1)
        MOUNTS.append((pre, os.path.realpath(d)))
MOUNTS.sort(key=lambda m: -len(m[0]))
mimetypes.add_type('application/javascript', '.mjs')
mimetypes.add_type('application/json', '.jsonc')
mimetypes.add_type('application/octet-stream', '.root')


class H(http.server.BaseHTTPRequestHandler):
    def log_message(self, fmt, *args):
        if os.environ.get('SERVE_LOG'):
            sys.stderr.write('%s\n' % (fmt % args))

    def resolve(self):
        path = urllib.parse.unquote(urllib.parse.urlsplit(self.path).path)
        for pre, d in MOUNTS:
            if path.startswith(pre):
                rel = path[len(pre):].lstrip('/')
                full = os.path.realpath(os.path.join(d, rel))
                if not (full == d or full.startswith(d + os.sep)):
                    return None
                if os.path.isdir(full):
                    full = os.path.join(full, 'index.html') if os.path.exists(os.path.join(full, 'index.html')) \
                        else os.path.join(full, 'index.htm')
                if not os.path.isfile(full) and SPA and pre == '/' and '.' not in os.path.basename(rel):
                    full = os.path.join(d, 'index.html')
                return full if os.path.isfile(full) else None
        return None

    def cors(self):
        self.send_header('Access-Control-Allow-Origin', '*')
        self.send_header('Access-Control-Allow-Headers', 'Range, Content-Type')
        self.send_header('Access-Control-Expose-Headers', 'Content-Range, Content-Length, Accept-Ranges')

    def do_OPTIONS(self):
        self.send_response(204); self.cors(); self.end_headers()

    def do_HEAD(self):
        self.do_GET(head=True)

    def do_GET(self, head=False):
        full = self.resolve()
        if not full:
            self.send_response(404); self.cors(); self.end_headers(); return
        size = os.path.getsize(full)
        ctype = mimetypes.guess_type(full)[0] or 'application/octet-stream'
        rng = self.headers.get('Range')
        start, end = 0, size - 1
        status = 200
        if rng:
            m = re.match(r'bytes=(\d*)-(\d*)', rng.split(',')[0].strip())
            if m:
                if m.group(1):
                    start = int(m.group(1)); end = int(m.group(2)) if m.group(2) else size - 1
                else:
                    start = size - int(m.group(2)); end = size - 1
                end = min(end, size - 1); status = 206
        self.send_response(status)
        self.cors()
        self.send_header('Content-Type', ctype)
        self.send_header('Accept-Ranges', 'bytes')
        self.send_header('Cache-Control', 'no-store')
        if status == 206:
            self.send_header('Content-Range', f'bytes {start}-{end}/{size}')
        self.send_header('Content-Length', str(end - start + 1))
        self.end_headers()
        if head:
            return
        with open(full, 'rb') as f:
            f.seek(start)
            left = end - start + 1
            while left > 0:
                chunk = f.read(min(left, 1 << 20))
                if not chunk:
                    break
                self.wfile.write(chunk); left -= len(chunk)


http.server.ThreadingHTTPServer.allow_reuse_address = True
http.server.ThreadingHTTPServer(('127.0.0.1', PORT), H).serve_forever()
