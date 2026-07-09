"""
TOPOLOGICA ablation sweep -- Step 1: enumerate the 18-configuration sweep
grid used by run_sweep.py.

Design (per universal-ablation-engine Part 2 -- this is primarily a
hyperparameter SWEEP, not a component ablation, except for the imputation
RULE which is a genuine discrete-component ablation):

  1. Gene-count sweep:      HVG in {500, 1000, 2000, 4000}, PC=50, zero-fill, n_perm=500
  2. PC-count sweep:        PC in {10, 25, 50, 100}, HVG=2000, zero-fill, n_perm=500
  3. Imputation ablation:   rule in {zero, mean, median, drop}, HVG=2000, PC=50, n_perm=500
  4. Permutation convergence: n_perm in {100, 200, 500, 1000, 2000}, HVG=2000, PC=50, zero
  5. Interaction corners:   HVG in {500, 4000} x PC in {10, 100}, zero-fill, n_perm=500

The pre-registered default (HVG2000_PC50_zero_np500) is the shared center
point and appears exactly once (configs are de-duplicated by tag), which is
why the union of the five items above yields 18 configs, not 4+4+4+5+4=21.

Output: configs.json (list of dicts consumed by run_sweep.py)
Run: python 01_make_configs.py
"""
import json
import os

_ROOT = os.path.dirname(os.path.abspath(__file__))

configs = []
seen_tags = set()


def add(tag, imputation_rule, n_hvg, n_pcs, n_perm_pipeline, n_draws_gaussian):
    if tag in seen_tags:
        return
    seen_tags.add(tag)
    configs.append({
        'tag': tag, 'imputation_rule': imputation_rule, 'n_hvg': n_hvg, 'n_pcs': n_pcs,
        'n_perm_pipeline': n_perm_pipeline, 'n_draws_gaussian': n_draws_gaussian
    })


def build_configs():
    # Item 1: gene count sweep (PC=50, zero, np=500)
    for hvg in [500, 1000, 2000, 4000]:
        add(f"HVG{hvg}_PC50_zero_np500", 'zero', hvg, 50, 500, 500)

    # Item 2: PC count sweep (HVG=2000, zero, np=500)
    for pc in [10, 25, 50, 100]:
        add(f"HVG2000_PC{pc}_zero_np500", 'zero', 2000, pc, 500, 500)

    # Item 3: imputation ablation (HVG=2000, PC=50, np=500)
    for rule in ['zero', 'mean', 'median', 'drop']:
        add(f"HVG2000_PC50_{rule}_np500", rule, 2000, 50, 500, 500)

    # Item 4: permutation count convergence (HVG=2000, PC=50, zero, gauss fixed at 500)
    for np_ in [100, 200, 500, 1000, 2000]:
        add(f"HVG2000_PC50_zero_np{np_}", 'zero', 2000, 50, np_, 500)

    # Item 5: interaction check corners (HVG in {500,4000} x PC in {10,100}, zero, np=500)
    for hvg in [500, 4000]:
        for pc in [10, 100]:
            add(f"HVG{hvg}_PC{pc}_zero_np500", 'zero', hvg, pc, 500, 500)

    return configs


if __name__ == "__main__":
    cfgs = build_configs()
    print(f"Total unique configs: {len(cfgs)}")
    for c in cfgs:
        print(" ", c['tag'])
    assert len(cfgs) == 18, f"expected 18 configs, got {len(cfgs)}"
    out_path = os.path.join(_ROOT, "configs.json")
    with open(out_path, "w") as f:
        json.dump(cfgs, f, indent=2)
    print(f"Saved {out_path}")
