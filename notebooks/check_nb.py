"""Report executed/errored cells of notebooks (helper for the build)."""
import json, sys, glob
for f in sys.argv[1:] or sorted(glob.glob('*/*.ipynb')):
    nb = json.load(open(f))
    cells = [c for c in nb['cells'] if c['cell_type'] == 'code']
    run = sum(1 for c in cells if c.get('execution_count'))
    errs = [(i, o) for i, c in enumerate(cells) for o in c.get('outputs', []) if o.get('output_type') == 'error']
    print(f'{f}: {run}/{len(cells)} code cells executed', '' if not errs else f"| ERROR cell {errs[0][0]}: {errs[0][1]['ename']}: {errs[0][1]['evalue'][:300]}")
