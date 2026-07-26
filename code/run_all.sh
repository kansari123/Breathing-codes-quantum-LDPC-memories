#!/bin/bash
# Reproduces every result in the breathing-codes v0+v1 study.
# Total runtime ~10-15 min on a laptop. pip install -r requirements.txt first.
set -e
echo "=== 1. Code verification (n,k via GF2 rank; d via ILP; ~1 min) ==="
python3 cellulation.py
echo "=== 2. Circuit sanity (stim rejects any non-deterministic detector) ==="
python3 circuits.py
echo "=== 3. Kill-test 1: exact fault distances (hernia check; ~1 min) ==="
python3 hernia.py
echo "=== 4. Kill-test 2: Monte Carlo economics (50k shots/pt; ~6-8 min) ==="
python3 economics.py
echo "=== 5. Kill-test 3: movement accounting ==="
python3 kt3_moves.py
echo "=== 6. v1-A: deep-breath seam locality + KT3 rematch ==="
python3 v1_deep.py
echo "=== 7. v1-B: partial breathing (column dial + Bring control; ~2 min) ==="
python3 v1_partial.py
echo "=== 8. v1-B addendum: targeted-logical distance ==="
python3 v1_targeted_check.py
echo "All results reproduced."

# v1-C: exhale-to-operate (structural only, no stim needed)
python3 exhale_ops.py
python3 exhale_ops2.py
