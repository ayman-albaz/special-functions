import std/[math, monotimes, strutils, times]
import special_functions

template bench(name: string, body: untyped) =
  let start = getMonoTime()
  body
  let elapsed = getMonoTime() - start
  let ms = elapsed.inMicroseconds.float / 1000.0
  echo name, ": ", ms.formatFloat(ffDecimal, 3), " ms"

const N = 100_000

const sep = "--------------------------------------------------"

echo "Special Functions — Nim benchmark (", N, " iterations each)"
echo sep

# ----- Gamma -----
bench "lowerIncompleteGamma(2, 0.5)":
  for i in 0 ..< N:
    discard lowerIncompleteGamma(2.0, 0.5)

bench "regularizedLowerIncompleteGamma(2, 0.5)":
  for i in 0 ..< N:
    discard regularizedLowerIncompleteGamma(2.0, 0.5)

bench "inverseLowerIncompleteGamma(2, 0.5)":
  for i in 0 ..< N:
    discard inverseLowerIncompleteGamma(2.0, 0.5)

bench "lgamma(5)":
  for i in 0 ..< N:
    discard lgamma(5.0)

# ----- Beta -----
bench "beta(2, 2)":
  for i in 0 ..< N:
    discard beta(2.0, 2.0)

bench "regularizedLowerIncompleteBeta(2, 2, 0.5)":
  for i in 0 ..< N:
    discard regularizedLowerIncompleteBeta(2.0, 2.0, 0.5)

bench "inverseRegularizedLowerIncompleteBeta(2, 2, 0.5)":
  for i in 0 ..< N:
    discard inverseRegularizedLowerIncompleteBeta(2.0, 2.0, 0.5)

# ----- Erf -----
bench "erf(1)":
  for i in 0 ..< N:
    discard erf(1.0)

bench "erfcx(3)":
  for i in 0 ..< N:
    discard erfcx(3.0)

bench "inverseErf(0.5)":
  for i in 0 ..< N:
    discard inverseErf(0.5)

# ----- Polygamma -----
bench "digamma(1)":
  for i in 0 ..< N:
    discard digamma(1.0)

bench "trigamma(1)":
  for i in 0 ..< N:
    discard trigamma(1.0)

# ----- Precomputed G/B (batch) -----
echo sep

let G = gamma(2.0)
bench "lowerIncompleteGamma(2, 0.5, G) [precomputed]":
  for i in 0 ..< N:
    discard lowerIncompleteGamma(2.0, 0.5, G)

let B = beta(2.0, 2.0)
bench "lowerIncompleteBeta(2, 2, 0.5, B) [precomputed]":
  for i in 0 ..< N:
    discard lowerIncompleteBeta(2.0, 2.0, 0.5, B)

# ----- float32 -----
echo ""
echo "float32 (", N, " iterations each)"
echo sep

bench "lowerIncompleteGamma(2f32, 0.5f32)":
  for i in 0 ..< N:
    discard lowerIncompleteGamma(2.0'f32, 0.5'f32)

bench "inverseLowerIncompleteGamma(2f32, 0.5f32)":
  for i in 0 ..< N:
    discard inverseLowerIncompleteGamma(2.0'f32, 0.5'f32)

bench "lgamma(5f32)":
  for i in 0 ..< N:
    discard lgamma(5.0'f32)

bench "beta(2f32, 2f32)":
  for i in 0 ..< N:
    discard beta(2.0'f32, 2.0'f32)

bench "regularizedLowerIncompleteBeta(2f32, 2f32, 0.5f32)":
  for i in 0 ..< N:
    discard regularizedLowerIncompleteBeta(2.0'f32, 2.0'f32, 0.5'f32)

bench "erf(1f32)":
  for i in 0 ..< N:
    discard erf(1.0'f32)

bench "erfcx(3f32)":
  for i in 0 ..< N:
    discard erfcx(3.0'f32)

bench "inverseErf(0.5f32)":
  for i in 0 ..< N:
    discard inverseErf(0.5'f32)

bench "digamma(1f32)":
  for i in 0 ..< N:
    discard digamma(1.0'f32)

bench "trigamma(1f32)":
  for i in 0 ..< N:
    discard trigamma(1.0'f32)
