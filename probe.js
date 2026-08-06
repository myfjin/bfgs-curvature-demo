// Demonstrates the BFGS curvature-condition bug.
//
// BFGS requires s·y > 0 to keep the inverse-Hessian approximation positive
// definite.  The buggy implementation tests |s·y| < 1e-18, which lets a
// negative curvature through.  On a convex quadratic this never matters; on a
// non-convex objective it does.

const { bfgs: bfgsBuggy } = require('./bfgs.js');
const { bfgs: bfgsFixed } = require('./bfgs_fixed.js');

function dot(a, b) {
    let s = 0.0;
    for (let i = 0; i < a.length; i++) s += a[i] * b[i];
    return s;
}

function norm2(v) {
    return Math.sqrt(dot(v, v));
}

// Convex quadratic: f(x) = (x1 - 2)^2 + (x2 - 3)^2
function gradConvex(x) {
    return [2 * (x[0] - 2), 2 * (x[1] - 3)];
}

// Non-convex objective used to expose the bug.
// f(x) = 0.5*x1^2 - x1*x2 + 0.1*x2^4  has regions of negative curvature.
// Gradient below is exact for this f.
function gradNonConvex(x) {
    return [x[0] - x[1], -x[0] + 0.4 * x[1] * x[1] * x[1]];
}

function objectiveNonConvex(x) {
    return 0.5 * x[0] * x[0] - x[0] * x[1] + 0.1 * x[1] * x[1] * x[1] * x[1];
}

function runCase(name, x0, gradient, objective) {
    console.log(`\n=== ${name} ===`);
    const buggy = bfgsBuggy(x0.slice(), gradient, 200);
    const fixed = bfgsFixed(x0.slice(), gradient, 200);

    console.log('buggy final:  ', buggy.map(v => v.toFixed(6)).join(', '));
    console.log('fixed final:  ', fixed.map(v => v.toFixed(6)).join(', '));

    if (objective) {
        console.log('buggy f(x):   ', objective(buggy).toExponential(4));
        console.log('fixed f(x):   ', objective(fixed).toExponential(4));
    }

    const dist = norm2([buggy[0] - fixed[0], buggy[1] - fixed[1]]);
    console.log('||buggy - fixed||:', dist.toExponential(4));
    return dist;
}

console.log('BFGS curvature-condition demo');
console.log('Bug:  if (Math.abs(s·y) < 1e-18) break;   lets negative s·y through');
console.log('Fix:  if (s·y < 1e-18) break;              rejects it\n');

const convexDiff = runCase('Convex quadratic (bug invisible)', [0, 0], gradConvex, null);
const nonConvexDiff = runCase('Non-convex objective (bug visible)', [0.1, 0.6], gradNonConvex, objectiveNonConvex);

console.log('\n=== SUMMARY ===');
if (convexDiff < 1e-6) {
    console.log('Convex: BUGGY and FIXED agree — the defect is hidden on convex problems.');
} else {
    console.log('Unexpected: convex results diverged.');
}

if (nonConvexDiff > 1e-3) {
    console.log('Non-convex: BUGGY and FIXED DIVERGE — the Math.abs() check accepts negative curvature.');
} else {
    console.log('Unexpected: non-convex results did not diverge; try a different starting point.');
}
