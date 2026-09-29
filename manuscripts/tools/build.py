"""Build the seven manuscript PDFs from their committed LaTeX sources."""
from pathlib import Path
import argparse
import os
import shutil
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdflatex", default=os.environ.get("PDFLATEX", "pdflatex"))
    args = parser.parse_args()
    compiler = shutil.which(args.pdflatex)
    if not compiler:
        sys.exit("pdflatex was not found. Install a TeX distribution or pass --pdflatex PATH.")
    root = Path(__file__).resolve().parents[1]
    names = (
        "00-Overview/Program_Overview.tex",
        "01-Random-Frame-RPCholesky/Random_Frame_RPCholesky.tex",
        "02-Product-Probe-Drift/Product_Probe_Drift_v2.tex",
        "03-Exact-Gaussian-Range-Laws/Exact_Gaussian_Range_Laws.tex",
        "04-SparseStack-Companion/SparseStack_Constants_Nystrom_Adaptive.tex",
        "05-SRHT-Companion/SRHT_Prescriptions_Nystrom_Adaptive.tex",
        "06-Robust-Sketched-RPCholesky/RPCholesky_Approximate_Sketched_Pivots.tex",
    )
    for name in names:
        source = root / name
        for _ in range(3):
            result = subprocess.run(
                [compiler, "-interaction=nonstopmode", "-halt-on-error", source.name],
                cwd=source.parent, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                text=True, encoding="utf-8", errors="replace",
            )
            if result.returncode:
                sys.exit(f"Build failed: {source}\n{result.stdout[-6000:]}")
        log = source.with_suffix(".log").read_text(encoding="utf-8", errors="replace")
        unresolved = [line for line in log.splitlines() if
                      ("undefined" in line.lower() and "warning" in line.lower())
                      or "Rerun to get cross-references right" in line]
        if unresolved:
            sys.exit(f"Unresolved references in {source.name}:\n" + "\n".join(unresolved))
        print(f"Built {source.parent.name}/{source.with_suffix('.pdf').name}")


if __name__ == "__main__":
    main()
