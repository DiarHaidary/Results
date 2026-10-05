"""Run the manuscript's deterministic checks using this Python interpreter."""
from pathlib import Path
import subprocess
import sys

root = Path(__file__).resolve().parent
for name in ["check_srht.py", "check_srht_terminal.py", "check_kr_spectrum.py",
             "check_kr_fourth.py", "check_spherical.py"]:
    subprocess.run([sys.executable, str(root / name)], cwd=root, check=True)
print("All five manuscript check programs passed.")
