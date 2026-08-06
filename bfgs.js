// Standalone BFGS demo — buggy curvature check.
//
// The defect: BFGS requires s·y > 0 to keep the inverse-Hessian approximation
// positive definite.  The line below tests |s·y| < 1e-18, which lets a NEGATIVE
// curvature through (e.g. s·y = -5 is non-zero, so rho = -0.2).  On a convex
// least-squares objective this never fires; the demo uses a non-convex objective
// where s·y can become negative, so the bug becomes visible.

function dot(a, b) {
    let s = 0.0;
    for (let i = 0; i < a.length; i++) s += a[i] * b[i];
    return s;
}

function matVec(M, v) {
    const out = [];
    for (let i = 0; i < M.length; i++) {
        let s = 0.0;
        for (let j = 0; j < M[0].length; j++) s += M[i][j] * v[j];
        out.push(s);
    }
    return out;
}

function bfgs(theta0, gradient, maxIter) {
    const m = theta0.length;
    let theta = theta0.slice();
    let g = gradient(theta);

    const H = [];
    for (let i = 0; i < m; i++) {
        const row = [];
        for (let j = 0; j < m; j++) row.push(i === j ? 1.0 : 0.0);
        H.push(row);
    }

    for (let step = 0; step < maxIter; step++) {
        if (dot(g, g) < 1e-24) break;

        const p = matVec(H, g).map(v => -v);
        const alpha = 0.01; // fixed small step so the non-convex landscape is traversed faithfully

        const sVec = [];
        for (let i = 0; i < m; i++) {
            sVec.push(alpha * p[i]);
            theta[i] += sVec[i];
        }
        const gNew = gradient(theta);
        const yVec = [];
        for (let i = 0; i < m; i++) yVec.push(gNew[i] - g[i]);
        g = gNew;

        const ys = dot(yVec, sVec);
        // === LINE UNDER TEST ===
        // Bug: Math.abs(ys) allows a negative s·y through.  The fixed version uses (ys < 1e-18).
        if (Math.abs(ys) < 1e-18) break;
        // === END LINE UNDER TEST ===
        const rho = 1.0 / ys;

        if (rho < 0) {
            console.log(`[BUGGY]  step ${step}: ys = ${ys.toExponential(3)}, rho = ${rho.toExponential(3)}  (negative curvature accepted)`);
        }

        const Hy = matVec(H, yVec);
        const syH = [];
        for (let j = 0; j < m; j++) {
            let s = 0.0;
            for (let i = 0; i < m; i++) s += sVec[i] * H[i][j];
            syH.push(s);
        }
        for (let i = 0; i < m; i++) {
            for (let j = 0; j < m; j++) {
                H[i][j] += rho * sVec[i] * sVec[j]
                    - rho * sVec[i] * Hy[j]
                    - rho * yVec[i] * syH[j]
                    + rho * rho * yVec[i] * dot(sVec, Hy) * sVec[j];
            }
        }
    }
    return theta;
}

module.exports = { bfgs, dot };
