# Statistical Correctness Gate

The previous local gate accepted only circuits with zero classical failures on
9024 sampled inputs. That was a compact way to make a 99%-correct circuit very
unlikely to pass:

```text
Pr[Binomial(9024, 0.01) = 0] = 0.99^9024 = 2^-130.844...
```

The new pre-promotion gate uses 100000 trusted samples and accepts a bounded
number of sampled output/phase failures. Ancilla cleanup remains exact:
any leftover non-register qubit still fails the run.

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

Let `X ~ Binomial(n = 100000, p = 0.01)` be the number of sampled output/phase
failures from a circuit with a 1% per-input sampled failure rate. Then:

```text
mean      = n p              = 1000
std dev   = sqrt(n p (1-p))  = sqrt(990) = 31.464265...
```

The largest acceptance cutoff with `Pr[X <= cutoff] <= 2^-128` is 617:

```text
Pr[X <= 617] = 2^-128.368...
Pr[X <= 618] = 2^-127.665...
```

So the gate accepts at most 617 sampled output/phase failures out of 100000
shots, counting a shot with both an output mismatch and phase garbage only
once. This matches the old anti-grinding target while giving approximate
circuits a large sampled test before the promotion test. Ancilla errors are
not part of this allowance; any leftover non-register qubit remains a hard
failure.

For `Y ~ Binomial(100000, p)`, the false-rejection probability for a circuit
with true sampled output/phase failure rate `p` is `Pr[Y > 617]`:

| True failure rate `p` | Correctness | Mean failures | `Pr[Y > 617]` |
|---:|---:|---:|---:|
| 0.100% | 99.900% | 100.0 | 2^-884.190 |
| 0.200% | 99.800% | 200.0 | 2^-409.477 |
| 0.300% | 99.700% | 300.0 | 2^-191.317 |
| 0.400% | 99.600% | 400.0 | 2^-78.169 |
| 0.500% | 99.500% | 500.0 | 1.81282e-7 |
| 0.600% | 99.400% | 600.0 | 0.2358 |
| 0.617% | 99.383% | 617.0 | 0.4893 |
| 0.750% | 99.250% | 750.0 | 1 - 2^-21.763 |
| 0.900% | 99.100% | 900.0 | 1 - 2^-77.120 |
| 1.000% | 99.000% | 1000.0 | 1 - 2^-128.368 |

So a circuit with 99.9% sampled output/phase correctness has negligible
false-rejection probability under this gate.
