"""Check compiled references, local dependencies, and repository PDF links.

Requires PyMuPDF (pip install pymupdf). HTTP checking is intentionally separate:
repository links are validated against the actual checkout, including URL decoding.
"""
from pathlib import Path
from urllib.parse import urlparse, unquote
import argparse
import json
import re
import sys
import fitz


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    exact_files = {p.relative_to(root).as_posix() for p in root.rglob("*") if p.is_file()}
    exact_dirs = {p.relative_to(root).as_posix() for p in root.rglob("*") if p.is_dir()}
    report = {"papers": [], "errors": [], "external_urls": [], "repository_urls": []}
    external, repository = set(), set()
    for markdown in [root / "README.md", *sorted((root / "manuscripts").rglob("*.md"))]:
        for target in re.findall(r"\[[^\]]*\]\(([^)]+)\)", markdown.read_text(encoding="utf-8")):
            parsed = urlparse(target)
            if parsed.scheme or not parsed.path:
                continue
            candidate = (markdown.parent / unquote(parsed.path)).resolve()
            try:
                relative = candidate.relative_to(root).as_posix()
            except ValueError:
                report["errors"].append(f"Markdown target outside repository: {markdown.name}: {target}")
                continue
            if relative not in exact_files and relative not in exact_dirs:
                report["errors"].append(f"Missing Markdown target: {markdown.relative_to(root)}: {target}")
    for source in sorted((root / "manuscripts").glob("0*/*.tex")):
        text = source.read_text(encoding="utf-8")
        clean = re.sub(r"(?<!\\)%[^\n]*", "", text)
        if re.search(r"(?:[CF]:[/\\]|file://|sandbox:|Desktop[/\\]|Sandbox[/\\])", clean, re.I):
            report["errors"].append(f"Private filesystem reference: {source.relative_to(root)}")
        for match in re.findall(r"\\(?:input|include)\{([^}]+)\}", clean):
            path = source.parent / match
            if not path.suffix:
                path = path.with_suffix(".tex")
            if not path.is_file():
                report["errors"].append(f"Missing LaTeX dependency: {source.name}: {match}")
        pdf = source.with_suffix(".pdf")
        if not pdf.is_file():
            report["errors"].append(f"Missing PDF: {pdf.relative_to(root)}")
            continue
        doc = fitz.open(pdf)
        item = {"source": source.relative_to(root).as_posix(), "pdf": pdf.relative_to(root).as_posix(),
                "pages": len(doc), "links": 0, "overfull_warnings": []}
        for page in doc:
            for link in page.get_links():
                uri = link.get("uri")
                if not uri:
                    if link.get("kind") == fitz.LINK_GOTO and not (0 <= link.get("page", -1) < len(doc)):
                        report["errors"].append(f"Invalid internal PDF destination: {pdf.name}")
                    continue
                item["links"] += 1
                parsed = urlparse(uri)
                if parsed.scheme not in ("http", "https", "mailto"):
                    report["errors"].append(f"Non-public PDF URL: {pdf.name}: {uri}")
                    continue
                prefixes = ("/DiarHaidary/Results/blob/main/", "/DiarHaidary/Results/tree/main/")
                prefix = next((p for p in prefixes if parsed.path.startswith(p)), None)
                if parsed.netloc == "github.com" and prefix:
                    repository.add(uri)
                    target = unquote(parsed.path[len(prefix):])
                    path = root / target
                    valid = (path.is_file() and target in exact_files) if "/blob/" in prefix else (path.is_dir() and target in exact_dirs)
                    if not valid:
                        report["errors"].append(f"Missing repository target: {pdf.name}: {uri}")
                else:
                    external.add(uri)
        logpath = source.with_suffix(".log")
        if logpath.exists():
            log = logpath.read_text(encoding="utf-8", errors="replace")
            for warning in re.findall(r"LaTeX Warning:[^\n]*(?:\n[^\n]*)?", log):
                if "undefined" in warning.lower() or "multiply defined" in warning.lower():
                    report["errors"].append(f"{source.name}: {warning}")
            item["overfull_warnings"] = re.findall(r"Overfull \\[hv]box[^\n]*", log)
            if item["overfull_warnings"]:
                report["errors"].append(f"Overfull box in {source.name}")
        report["papers"].append(item)
    report["external_urls"] = sorted(external)
    report["repository_urls"] = sorted(repository)
    result = json.dumps(report, indent=2, ensure_ascii=False)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(result + "\n", encoding="utf-8")
    print(result)
    sys.exit(bool(report["errors"]))


if __name__ == "__main__":
    main()
