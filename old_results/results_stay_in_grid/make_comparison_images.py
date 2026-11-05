"""Side-by-side classification maps: the paper's run (with the sink state) above, the
stay-in-grid run below, and a third row that shows only the cells whose class changed (in
their new colour; unchanged pixels are faded).

    python make_comparison_images.py
"""
from pathlib import Path

import numpy as np
from PIL import Image

HERE = Path(__file__).resolve().parent
PAPER = next((d / "paper_results" for d in HERE.parents if (d / "paper_results").is_dir()),
             HERE.parent / "paper_results") / "results_hscc_2027_newbellmans_12workers_current_code"
OUT = HERE / "comparison_images"
OUT.mkdir(exist_ok=True)
CASES = [("fig4a_dubins120_std", "ra_nodiscount", "dubins_res120_"),
         ("fig4b_dubins120_discount", "ra_discount", "dubins_res120_"),
         ("fig5_avoid80", "avoid_nodiscount", "dubins_res80_"),
         ("fig2_dubins20_discount", "ra_discount", "dubins_res20_")]


def run_dir(root, mode, prefix, restricted):
    for d in sorted((root / mode).glob(prefix + "*")):
        if d.name.endswith("_restrictgrid") == restricted:
            return d
    return None


for name, mode, prefix in CASES:
    mine, ref = run_dir(HERE, mode, prefix, True), run_dir(PAPER, mode, prefix, False)
    if mine is None:
        continue
    phases = sorted(int(p.stem.split("_")[3]) for p in mine.glob("value_function_phase_*_complete.png"))
    k = phases[-1]
    a = np.asarray(Image.open(ref / f"value_function_phase_{k}_complete.png").convert("RGB"))
    b = np.asarray(Image.open(mine / f"value_function_phase_{k}_complete.png").convert("RGB"))
    h = a.shape[0]
    top = int(h * 0.662)                      # the third row: cell classification
    a, b = a[top:], b[top:]
    changed = np.any(a != b, axis=2)
    faded = (0.25 * b + 0.75 * 255).astype(np.uint8)
    diff = np.where(changed[..., None], b, faded)
    img = Image.fromarray(np.vstack([a, b, diff]))
    img = img.resize((img.width // 3, img.height // 3), Image.LANCZOS)
    out = OUT / f"{name}_phase{k}.png"
    img.save(out)
    print(f"{out.name}: {changed.mean() * 100:.2f}% of the pixels changed")
