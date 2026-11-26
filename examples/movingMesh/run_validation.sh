#!/bin/sh
# run_validation.sh - generate and run the moving mesh validation cases, print the numbers quoted in README.md.
#
#   sh run_validation.sh        (about 2 minutes on one core)
#
# Needs: bin/cappuccino built with 'make' in the repository root, python3 with numpy. Work directory: ./run (created).
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
EXE=${EXE:-$HERE/../../bin/cappuccino}
W=${WORK:-$PWD/run}
mkdir -p "$W"; cd "$W"

runcase () {   # runcase <dir> [restart file name] : run the solver in the case directory
  ( cd "$1" && rm -rf vtk/[0-9]* restart monitor && "$EXE" input.nml monitor restart > run.log 2>&1 )
}
gcl () {       # gcl <dir> : largest GCL error and smallest cell volume over the run
  grep "GCL error" "$1/monitor" | sed 's/.*min vol = *\([^ ]*\) .*GCL error = *\(.*\)/\1 \2/' \
    | awk 'NR==1{mv=$1;g=$2}{if($1<mv)mv=$1; if($2>g)g=$2}END{printf "  GCL error (max over the run) = %.2e   smallest cell volume (min over the run) = %.3e\n", g, mv}'
}
mass () {      # mass <dir> : largest global continuity error over the run (closed domains only: no outflow)
  grep "continuity errors" "$1/monitor" | sed 's/.*global = *\([^,]*\),.*/\1/' \
    | awk 'BEGIN{m=0}{x=$1+0; if(x<0)x=-x; if(x>m)m=x}END{printf "  global continuity error (max over the run) = %.2e\n", m}'
}

echo "=== 1. Piston, constant velocity: exact plug flow, u = U, p = 0 (2D 40x10, SIMPLE) ==="
python3 "$HERE/mkcase.py" piston pis2d --motion translate --U 0.5 --steps 40 --dt 0.01 > /dev/null; runcase pis2d
python3 "$HERE/piston_check.py" pis2d || true; gcl pis2d

echo "=== 2. The same with PISO ==="
python3 "$HERE/mkcase.py" piston pispiso --motion translate --U 0.5 --steps 40 --dt 0.01 --solver piso > /dev/null; runcase pispiso
python3 "$HERE/piston_check.py" pispiso || true

echo "=== 3. The same in 3D (20 x 6 x 6 cells, one cell column is the 2D case) ==="
python3 "$HERE/mkcase.py" piston pis3d --ndim 3 --nx 20 --ny 6 --motion translate --U 0.5 --steps 20 --dt 0.01 > /dev/null; runcase pis3d
python3 "$HERE/piston_check.py" pis3d || true; gcl pis3d

echo "=== 4. Oscillating piston, x = 0.1 sin(2 pi t), one period, two meshes ==="
for n in 40 80; do
  python3 "$HERE/mkcase.py" piston osc$n --nx $n --ny $((n/4)) --motion oscillate --amp 0.1 --freq 1.0 --steps 100 --dt 0.01 > /dev/null; runcase osc$n
  echo "nx = $n"; python3 "$HERE/piston_check.py" osc$n || true; gcl osc$n
done

echo "=== 5. Square oscillating in a closed box (RBF, 8 cells per unit length, 100 steps) ==="
python3 "$HERE/mkcase.py" body bodyosc --motion oscillate --steps 100 > /dev/null; runcase bodyosc; gcl bodyosc; mass bodyosc
echo "=== 6. Square rotating in a closed box (1 rad) ==="
python3 "$HERE/mkcase.py" body bodyrot --motion rotate --steps 50 > /dev/null; runcase bodyrot; gcl bodyrot; mass bodyrot
echo "=== done ==="
