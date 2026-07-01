# Statistical Correctness Gate

## Motivation

The main motivation for this gate is to eliminate grinding that came from the
Fiat-Shamir heuristic. A deterministic challenge derived from the submitted
artifact is cheap and reproducible, but it gives contestants a fixed target:
they can keep changing `ops.bin` until the derived sample set misses the
circuit's bad inputs.

We do not have the same proof-generation cost constraint here, so the trusted
evaluator chooses fresh OS randomness after `ops.bin` exists and mixes it with
`hash(ops.bin)`. The printed sample seed still makes a finished run
reproducible, but the sample stream is not available while the untrusted circuit
is being generated.

The previous local gate accepted only circuits with zero classical failures on
9024 sampled inputs. That was a compact way to make a 99%-correct circuit
unlikely to pass if the samples were not grindable:

```text
Pr[Binomial(9024, 0.01) = 0] = 0.99^9024 = 2^-130.844...
```

The new pre-promotion gate uses 50048 trusted samples and accepts a bounded
number of sampled output/phase failures. `50048 = 782 * 64`, so every sample
batch is full and the 64-lane simulator's gate counters do not include inactive
lanes from a partial final batch. Ancilla cleanup remains exact: any leftover
non-register qubit still fails the run.

## Sample Source

`eval_circuit` hashes the raw `ops.bin` bytes after the untrusted
`build_circuit` stage has emitted the artifact. If no seed is supplied, it then
chooses 256 bits from the OS CSPRNG and derives the printed sample seed as:

```text
SHAKE256(osrng_256_bits || SHAKE256(ops.bin))
```

The test-input stream is then derived from the sample seed and `hash(ops.bin)`.
The OS randomness is chosen after the artifact exists, so contestant code cannot
know the final sample stream while building `ops.bin`.

For reproducibility, the evaluator prints the seed and accepts:

```bash
./target/release/eval_circuit --sample-seed <printed-seed>
```

Use a fixed seed only after `ops.bin` already exists. Do not pass a promotion
seed through the full benchmark wrapper in an adversarial setting: even though
`benchmark.sh` forwards arguments only to the trusted evaluator, untrusted build
code may still be able to inspect parent process metadata on some systems and
grind `ops.bin` against a known seed.

## Bound

Let `X ~ Binomial(n = 50048, p = 0.01)` be the number of sampled output/phase
failures from a circuit with a 1% per-input sampled failure rate. Then:

```text
mean      = n p              = 500.48
std dev   = sqrt(n p (1-p))  = sqrt(495.4752) = 22.259272...
```

The largest acceptance cutoff with `Pr[X <= cutoff] <= 2^-128` is 238:

```text
Pr[X <= 238] = 2^-128.811...
Pr[X <= 239] = 2^-127.732...
```

So the gate accepts at most 238 sampled output/phase failures out of 50048
shots, counting a shot with both an output mismatch and phase garbage only
once. This matches the old anti-grinding target while giving approximate
circuits a large sampled test before the promotion test. Ancilla errors are
not part of this allowance; any leftover non-register qubit remains a hard
failure.

For `Y ~ Binomial(50048, p)`, the false-rejection probability for a circuit
with true sampled output/phase failure rate `p` is `Pr[Y > 238]`:

| True failure rate `p` | Correctness | Mean failures | `Pr[Y > 238]` |
|---:|---:|---:|---:|
| 0.100% | 99.900% | 50.0 | 2^-271.948 |
| 0.200% | 99.800% | 100.1 | 2^-104.478 |
| 0.300% | 99.700% | 150.1 | 1.37032e-11 |
| 0.400% | 99.600% | 200.2 | 0.004093 |
| 0.476% | 99.524% | 238.0 | 0.4828 |
| 0.500% | 99.500% | 250.2 | 0.7702 |
| 0.600% | 99.400% | 300.3 | 1 - 2^-13.181 |
| 0.617% | 99.383% | 308.8 | 1 - 2^-16.007 |
| 0.750% | 99.250% | 375.4 | 1 - 2^-45.846 |
| 0.900% | 99.100% | 450.4 | 1 - 2^-92.287 |
| 1.000% | 99.000% | 500.5 | 1 - 2^-128.811 |

So a circuit with 99.9% sampled output/phase correctness has negligible
false-rejection probability under this gate.
