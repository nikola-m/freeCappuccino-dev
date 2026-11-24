#!/bin/sh
# run_validation.sh - generate and run the CLSVOF validation cases, print the numbers quoted in README.md.
#
#   sh run_validation.sh [quick|full]      (default: quick, about 10 minutes on one core; full adds the finest grids)
#
# Needs: bin/cappuccino built with 'make' in the repository root, python3 with numpy. Work directory: ./run (created).
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
EXE=${EXE:-$HERE/../../bin/cappuccino}
MODE=${1:-quick}
W=${WORK:-$PWD/run}
mkdir -p "$W"; cd "$W"

runcase () {   # runcase <dir>  : run the solver in the case directory
  ( cd "$1" && "$EXE" input.nml monitor restart > run.log 2>&1 )
}

echo "=== 1. Solid body rotation of a disk, one revolution (volume of fluid advection, prescribed flow) ==="
for n in 32 64 $( [ "$MODE" = full ] && echo 128 ); do
  python3 "$HERE/mkcase.py" advect-rotate rot$n --n $n --steps $((10*n)) > /dev/null; runcase rot$n
  echo "n = $n"; python3 "$HERE/clsvof_post.py" shape rot$n --center 0.5 0.75 0 --radius 0.15 | sed -n 2,3p
done

echo "=== 2. Reversing vortex (LeVeque), T = 8 ==="
for n in 32 64 $( [ "$MODE" = full ] && echo 128 ); do
  python3 "$HERE/mkcase.py" advect-vortex vor$n --n $n --steps $((8*n)) > /dev/null; runcase vor$n
  echo "n = $n"; python3 "$HERE/clsvof_post.py" shape vor$n --center 0.5 0.75 0 --radius 0.15 | sed -n 2,3p
done

echo "=== 3. Static drop, Laplace pressure jump (sigma = 1, R = 0.25) ==="
for n in 16 32 64; do
  python3 "$HERE/mkcase.py" static-drop sd$n --n $n --steps 200 > /dev/null; runcase sd$n
  echo "2D n = $n"; python3 "$HERE/clsvof_post.py" laplace sd$n --ic cylinder --radius 0.25 --sigma 1.0
done
python3 "$HERE/mkcase.py" static-drop sd3d --n 24 --ndim 3 --steps 100 > /dev/null; runcase sd3d
echo "3D n = 24"; python3 "$HERE/clsvof_post.py" laplace sd3d --ic sphere --radius 0.25 --sigma 1.0

echo "=== 4. Rising bubble, Hysing et al. (2009) test case 1 ==="
for n in 20 40 $( [ "$MODE" = full ] && echo 80 ); do
  python3 "$HERE/mkcase.py" rising-bubble rb$n --n $n > /dev/null; runcase rb$n
  echo "h = 1/$n"; python3 "$HERE/clsvof_post.py" hysing rb$n
done

echo "=== 5. Falling drop in a closed box ==="
python3 "$HERE/mkcase.py" falling-drop fall2d --n 32 --ndim 2 --steps 300 > /dev/null; runcase fall2d
python3 "$HERE/clsvof_post.py" fall fall2d --gravity 1.0 --density-ratio 1000 --added-mass 1.0
