"""Check unique external PDF URLs listed by verify_release.py."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import json
from pathlib import Path
from urllib.request import Request, urlopen
from urllib.error import HTTPError


def check(url):
    def request(method):
        headers = {"User-Agent": "Manuscript-reference-check/1.0"}
        req = Request(url, method=method, headers=headers)
        with urlopen(req, timeout=20) as response:
            if method == "GET":
                response.read(1024)
            return {"url": url, "status": response.status, "final_url": response.url}
    try:
        return request("HEAD")
    except Exception:
        try:
            return request("GET")
        except HTTPError as exc:
            return {"url": url, "status": exc.code, "error": str(exc)}
        except Exception as exc:
            return {"url": url, "status": None, "error": str(exc)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("report", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    report = json.loads(args.report.read_text(encoding="utf-8"))
    urls = report["external_urls"]
    with ThreadPoolExecutor(max_workers=8) as pool:
        results = list(pool.map(check, urls))
    out = {"checked_utc": datetime.now(timezone.utc).isoformat(), "links": results}
    args.output.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    failures = [r for r in results if r.get("status") != 200]
    print(json.dumps({"total": len(results), "not_200": failures}, indent=2))


if __name__ == "__main__":
    main()
