#!/usr/bin/env bash
# Run the BFGS curvature-condition demo.
# Requires Node.js (v14+). No npm install needed — std-lib only.
set -e
cd "$(dirname "$0")"
node probe.js
