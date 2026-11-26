#!/bin/bash
# run_validation.sh - runs the all-speeds SIMPLE test set and compares with the existing solvers.
#
#   ./run_validation.sh [path/to/cappuccino] [workdir]
#
# 1. cavity 32x32 Re 100: stock SIMPLE vs all-speeds SIMPLE (incompressible mode, AMG for p'): iterations to convergence
# 2. nozzle (subsonic, isentropic, steady with pseudo time): Mach number, p, T against the quasi-1D solution (nozzle_check.py)
# 3. NACA0012 M = 0.5, Re = 5000 (coarse C-grid): all-speeds SIMPLE (dual time) vs block-coupled solver: drag
# Python 3 with numpy is needed for the case generators and the post-processing.
set -e
BIN=${1:-$(cd "$(dirname "$0")/../.." && pwd)/bin/cappuccino}
WORK=${2:-/tmp/asimple_validation}
HERE=$(cd "$(dirname "$0")" && pwd)
CPL=$HERE/../coupledPU
mkdir -p "$WORK"; cd "$WORK"

echo "== 1. cavity"
python3 "$HERE/mkcase.py" cavity cav --nx 32 --ny 32 > /dev/null
( cd cav
  $BIN input-simple.nml mon1 rst1 > /dev/null 2>&1; echo "stock SIMPLE         : $(grep -c Timestep mon1) iterations"
  $BIN input-asimple.nml mon2 rst2 > /dev/null 2>&1; echo "all-speeds SIMPLE+AMG: $(grep -c Timestep mon2) iterations" )

echo "== 2. nozzle"
python3 "$HERE/mkcase.py" nozzle noz > /dev/null
( cd noz; $BIN input-asimple.nml mon rst > /dev/null 2>&1; python3 "$CPL/nozzle_check.py" . rst | tail -3 )

echo "== 3. NACA0012"
python3 "$HERE/mkcase.py" naca naca > /dev/null
( cd naca
  $BIN input-asimple.nml mon1 rst1 > /dev/null 2>&1; echo "all-speeds SIMPLE:"; python3 "$CPL/postairfoil.py" . rst1 | sed -n 1,2p
  $BIN input-coupled.nml mon2 rst2 > /dev/null 2>&1; echo "block-coupled:";     python3 "$CPL/postairfoil.py" . rst2 | sed -n 1,2p )

echo "== 4. nozzle, compressible BC kinds of the namelist (M_in = 0.3, asimple_urfP = 0.2)"
for k in stagnation totalpressure massflow; do
  python3 "$HERE/mkcase.py" nozzle noz-$k --inlet $k --Min 0.3 > /dev/null
  ( cd noz-$k; $BIN input-asimple.nml mon rst > /dev/null 2>&1; echo "inlet $k: $(python3 "$CPL/nozzle_check.py" . rst | sed -n 1p | cut -c1-100)" )
done
for o in compoutlet open riemann; do
  python3 "$HERE/mkcase.py" nozzle noz-$o --inlet stagnation --outlet $o --Min 0.3 > /dev/null
  ( cd noz-$o; $BIN input-asimple.nml mon rst > /dev/null 2>&1; echo "outlet $o: $(python3 "$CPL/nozzle_check.py" . rst | sed -n 1p | cut -c1-100)" )
done

echo "== 5. circular-arc bump channel (GAMM type), M = 0.5 / 0.85 / 1.4"
for m in 0.5 0.85 1.4; do
  python3 "$HERE/mkcase.py" bump bump-$m --M $m > /dev/null
  ( cd bump-$m; $BIN input-asimple.nml mon rst > /dev/null 2>&1; echo "M = $m: $(grep -c Timestep mon) iterations"; python3 "$HERE/bump_check.py" . rst )
done

echo "== 6. channel Re 20000, wall functions (incompressible)"
python3 "$HERE/mkcase.py" channel chw --nx 40 --ny 10 --Re 20000 > /dev/null
for w in reichardt spl unified; do
  ( cd chw; sed "s/ASIMPLE = T/ASIMPLE = T\n asimple_wallFunc = '$w'\n ltest = T/" input-asimple.nml > in-$w.nml
    $BIN in-$w.nml mon-$w rst-$w > /dev/null 2>&1; echo "$w: $(grep -c Timestep mon-$w) iterations, $(grep 'wall function:' mon-$w | tail -1)" )
done
