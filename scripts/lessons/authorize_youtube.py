"""Run on the computer with your browser: python3 authorize_youtube.py CLIENT.json"""
import base64
import hashlib
import json
import os
from pathlib import Path
import secrets
import sys
from http.server import BaseHTTPRequestHandler, HTTPServer
from urllib.parse import parse_qs, urlencode, urlsplit
from urllib.request import Request, urlopen
import webbrowser


def main(client_file):
    client = json.loads(Path(client_file).read_text())["installed"]
    state = secrets.token_urlsafe(32)
    verifier = secrets.token_urlsafe(48)
    challenge = base64.urlsafe_b64encode(hashlib.sha256(verifier.encode()).digest()).decode().rstrip("=")

    class Callback(BaseHTTPRequestHandler):
        def log_message(self, *args):
            pass  # The URL contains a one-time credential.

        def do_GET(self):
            parts = urlsplit(self.path)
            query = parse_qs(parts.query)
            if parts.path != "/" or query.get("state") != [state]:
                self.send_error(400)
                return
            self.server.result = query
            self.send_response(200)
            self.end_headers()
            self.wfile.write(b"Return to your terminal.")

    with HTTPServer(("127.0.0.1", 0), Callback) as server:
        redirect = f"http://127.0.0.1:{server.server_port}/"
        params = dict(client_id=client["client_id"], redirect_uri=redirect,
                      response_type="code", scope="https://www.googleapis.com/auth/youtube.force-ssl",
                      state=state, code_challenge=challenge, code_challenge_method="S256",
                      access_type="offline", prompt="select_account consent")
        url = "https://accounts.google.com/o/oauth2/v2/auth?" + urlencode(params)
        print("Select the Google account and channel for @PhDQtami.\n" + url)
        webbrowser.open(url)
        while not hasattr(server, "result"):
            server.handle_request()
        if "code" not in server.result:
            raise SystemExit("Authorization was declined.")
        code = server.result["code"][0]

    fields = dict(client_id=client["client_id"], client_secret=client["client_secret"],
                  grant_type="authorization_code", code=code, redirect_uri=redirect,
                  code_verifier=verifier)
    with urlopen(Request("https://oauth2.googleapis.com/token", data=urlencode(fields).encode())) as response:
        token = json.load(response)
    bundle = {key: client[key] for key in ("client_id", "client_secret")}
    bundle["refresh_token"] = token["refresh_token"]
    dest = Path.home() / ".config" / "linearknot-youtube" / "authorization.json"
    dest.parent.mkdir(mode=0o700, parents=True, exist_ok=True)
    import tempfile
    with tempfile.NamedTemporaryFile(mode="w", dir=dest.parent, delete=False) as f:
        json.dump(bundle, f)
    os.replace(f.name, dest)
    print("Authorization saved to", dest)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    main(sys.argv[1])
