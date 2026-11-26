#!/usr/bin/env python3
"""
mkcase.py - test cases for the all-speeds SIMPLE solver (module simple_allspeeds, switch ASIMPLE = T).

  python3 mkcase.py naca    <dir> [--nt 101 --nw 40 --nr 40 --M 0.5 --Re 5000 --alpha 0 --dt 1e-4 --inner 6 --steps 1500]
  python3 mkcase.py nozzle  <dir> [--nx 60 --ny 6 --Min 0.2 --inlet compinlet|stagnation|totalpressure|massflow
                                   --outlet pressure|compoutlet|open|riemann]
  python3 mkcase.py bump    <dir> [--M 0.5|0.85|1.4 --nx 180 --ny 60]   (circular-arc bump channel, see mkbump.py)
  python3 mkcase.py cavity  <dir> [--nx 32 --ny 32 --Re 100]
  python3 mkcase.py channel <dir> [--nx 40 --ny 10 --Re 50 --outlet pressure|outlet]

The meshes and the flow set-up come from the generators of ../coupledPU (mkairfoil.py, mknozzle.py, mkcase.py); here the
input files for the all-speeds SIMPLE are derived from theirs (input files written by this script have short names: the code
reads the input file name into a 20-character string).

naca    input-asimple.nml    dual-time stepping: BDF1, dt = --dt, --inner outer iterations per time step (robust, recommended)
        input-asimple-st.nml steady with the pseudo-time term (asimple_pseudoCFL), converges slowly on the C-grid
        input-coupled.nml    the existing block-coupled all-speeds solver, for comparison
nozzle  input-asimple.nml    steady, pseudo-time CFL 10 with ramp; check with ../coupledPU/nozzle_check.py
cavity, channel  input-simple.nml (stock SIMPLE), input-asimple.nml (all-speeds SIMPLE, incompressible mode, AMG for p')
"""
import argparse, os, re, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
CPL = os.path.join(HERE, '..', 'coupledPU')


def set_params(text, **kw):
    """Set (or add) 'key = value' entries of the namelist text. Values are written as given."""
    for k, v in kw.items():
        pat = re.compile(r'^(\s*)' + re.escape(k) + r'\s*=.*$', re.I | re.M)
        if pat.search(text):
            text = pat.sub(lambda m: '%s%s = %s' % (m.group(1), k, v), text, count=1)
        else:
            text = text.replace('\n/', '\n%24s = %s\n/' % (k, v), 1) if '\n/' in text else text + '\n%s = %s\n' % (k, v)
    return text


def run(script, args):
    r = subprocess.run([sys.executable, os.path.join(CPL, script)] + args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if r.returncode != 0:
        print(r.stdout); sys.exit(1)
    return r.stdout


def write(path, text):
    with open(path, 'w') as f:
        f.write(text)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('case', choices=['naca', 'nozzle', 'bump', 'cavity', 'channel'])
    ap.add_argument('dir')
    ap.add_argument('--nt', default='101'); ap.add_argument('--nw', default='40'); ap.add_argument('--nr', default='40')
    ap.add_argument('--M', default='0.5'); ap.add_argument('--Re', default=None); ap.add_argument('--alpha', default='0')
    ap.add_argument('--dt', default='1e-4'); ap.add_argument('--inner', default='6'); ap.add_argument('--steps', default='1500')
    ap.add_argument('--nx', default=None); ap.add_argument('--ny', default=None); ap.add_argument('--Min', default='0.2')
    ap.add_argument('--outlet', default='pressure'); ap.add_argument('--inlet', default='compinlet')
    a = ap.parse_args()
    d = a.dir

    if a.case == 'naca':
        out = run('mkairfoil.py', [d, '--nt', a.nt, '--nw', a.nw, '--nr', a.nr, '--M', a.M, '--Re', a.Re or '5000', '--alpha', a.alpha])
        print(out.strip().splitlines()[-1])
        base = open(os.path.join(d, 'input-allspeeds.nml')).read()
        write(os.path.join(d, 'input-coupled.nml'), base)
        # Segregated solver: SIMPLE-type, energy equation on, no block-coupled solver
        common = dict(SIMPLE='T', COUPLED='F', calcEn='T', npcor='1', ASIMPLE='T', asimple_urfP='0.7', asimple_tolRelP='0.05')
        dual = set_params(base, ltransient='T', bdf='T', TIMESTEP=a.dt, NUMSTEP=a.steps, NZAPIS=a.steps, MAXIT=a.inner, **common)
        dual = set_params(dual, asimple_urfP='0.3', title="'NACA M=%s all-speeds SIMPLE, dual time stepping (BDF1)'" % a.M)
        write(os.path.join(d, 'input-asimple.nml'), dual)
        st = set_params(base, NUMSTEP='4000', NZAPIS='4000', asimple_pseudoCFL='3.0', asimple_cflRamp='200', **common)
        st = set_params(st, asimple_urfP='0.3', title="'NACA M=%s all-speeds SIMPLE, steady, pseudo-time'" % a.M)
        write(os.path.join(d, 'input-asimple-st.nml'), st)

    elif a.case == 'nozzle':
        out = run('mknozzle.py', [d, '--inlet', a.inlet, '--outlet', a.outlet, '--Min', a.Min, '--nx', a.nx or '60', '--ny', a.ny or '6'])
        base = open(os.path.join(d, 'input-allspeeds.nml')).read()
        write(os.path.join(d, 'input-coupled.nml'), base)
        st = set_params(base, SIMPLE='T', COUPLED='F', calcEn='T', npcor='1', ASIMPLE='T', asimple_urfP='0.7' if a.inlet == 'compinlet' else '0.2', asimple_pseudoCFL='10.0',
                        asimple_cflRamp='100', NUMSTEP='3000', NZAPIS='3000', MAXIT='1')
        write(os.path.join(d, 'input-asimple.nml'), st)

    elif a.case == 'bump':
        args = [d, '--M', a.M]
        if a.nx: args += ['--nx', a.nx]
        if a.ny: args += ['--ny', a.ny]
        r = subprocess.run([sys.executable, os.path.join(HERE, 'mkbump.py')] + args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        print(r.stdout.strip())
        if r.returncode != 0: sys.exit(1)

    else:
        args = [a.case, d]
        if a.nx: args += ['--nx', a.nx]
        if a.ny: args += ['--ny', a.ny]
        if a.Re: args += ['--Re', a.Re]
        if a.case == 'channel': args += ['--outlet', a.outlet]
        run('mkcase.py', args)
        base = open(os.path.join(d, 'input-simple.nml')).read()
        write(os.path.join(d, 'input-asimple.nml'), set_params(base, ASIMPLE='T', npcor='1'))

    print('Case written to', d)


if __name__ == '__main__':
    main()
