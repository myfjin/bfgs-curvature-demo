# BFGS curvature-condition bug demo

A tiny, self-contained demo of a real defect in two BFGS/L-BFGS implementations
from a pattern library: the curvature condition is checked with `Math.abs(s·y)`,
which lets **negative** `s·y` through.  On a convex least-squares objective this
never matters; on a non-convex objective it can make the inverse-Hessian
approximation indefinite and send the solver uphill.

## Files

- `bfgs.js` — buggy implementation, `Math.abs(ys) < 1e-18` intact.
- `bfgs_fixed.js` — same structure as `bfgs.js`, but the curvature check is
  `if (ys < 1e-18) break;` and the log line says `[FIXED]` instead of `[BUGGY]`.

The only functional difference between the two files is the curvature condition:
`Math.abs(ys) < 1e-18` becomes `ys < 1e-18`.  Everything else that differs is a
comment or a log line.

- `probe.js` — runs both on a convex and a non-convex objective and prints divergence.
- `run.sh` — one command to run the demo.

## Run it

```bash
git clone <repo-url> bfgs-curvature-demo
cd bfgs-curvature-demo
./run.sh
```

Requires Node.js (v14+).  No `npm install` — the code uses only built-ins.

## What you should see

On the convex quadratic, buggy and fixed implementations agree:

```
=== Convex quadratic (bug invisible) ===
buggy final:   1.734747, 2.602121
fixed final:   1.734747, 2.602121
||buggy - fixed||: 0.0000e+0
```

On the non-convex objective they diverge, and the buggy version prints the
negative `ys` values it is accepting:

```
[BUGGY]  step 2: ys = -2.265e-5, rho = -4.415e+4  (negative curvature accepted)
...
buggy final:   -7.706312, 28.549681
fixed final:   0.115412, 0.607006
buggy f(x):    6.6686e+4
fixed f(x):    -4.9820e-2
||buggy - fixed||: 2.9017e+1
```

The summary at the bottom states the conclusion:

```
Convex: BUGGY and FIXED agree — the defect is hidden on convex problems.
Non-convex: BUGGY and FIXED DIVERGE — the Math.abs() check accepts negative curvature.
```

## Why this matters

BFGS requires `s·y > 0` to preserve positive definiteness of the inverse-Hessian
approximation.  A check on `|s·y|` only guards against division by zero; it does
not guard against a negative inner product.

The non-convex objective in `probe.js` was chosen specifically because `s·y`
becomes negative on it.  That is the point: the defect is invisible on the convex
least-squares problem the original pattern was tested against, and only appears
on an objective that exercises the condition.  Constructing such a case is how
you test a guard, not a trick.

You do not have to trust this repo: the BFGS curvature condition `s·y > 0` is
standard and appears in Nocedal & Wright, *Numerical Optimization* (2nd ed.),
§6.1 and §7.3.  Compare it to the line in `bfgs.js`.

## License

MIT — this is a standalone educational demo.
