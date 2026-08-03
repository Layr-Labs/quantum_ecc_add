# The secp256k1 Point-Addition Challenge

> **Goal.** Build the cheapest reversible quantum circuit that performs one
> elliptic-curve point addition on **secp256k1**, scored by the product of
> **Toffoli count × peak qubit width**.

---

## Why this matters

Shor's algorithm breaks elliptic-curve cryptography by computing discrete
logarithms in time polynomial in the bit-width of the curve. The quantum cost
of *running* Shor on an ECC group is dominated by one inner primitive,
repeated thousands of times: **point addition** on the curve.

Faster point addition ⇒ fewer Toffoli gates ⇒ fewer magic states ⇒ less
physical hardware and less wall-clock time on a fault-tolerant quantum
computer. Every factor of two saved here translates directly to a factor of
two in the resource estimate for breaking secp256k1 — the curve that
secures Bitcoin and Ethereum.

---

## The benchmark, precisely

You are given a Rust harness that:

1. **Builds** a reversible circuit by calling `point_add::build()`.
   The circuit must consume four 256-element registers — `target_x`
   (qubits), `target_y` (qubits), `offset_x` (classical bits),
   `offset_y` (classical bits) — and overwrite `(target_x, target_y)`
   with the affine sum `(target_x, target_y) + (offset_x, offset_y)` on
   the secp256k1 curve.
2. **Validates** the circuit by simulating it on 50,048 random test points
   (782 full 64-shot batches, just over 50k).
   Inputs are derived from `hash(ops.bin)` and a trusted random seed chosen
   inside `eval_circuit` after `ops.bin` exists. If no seed is supplied, the
   printed seed is derived as `hash(osrng || hash(ops.bin))`, so the circuit
   cannot know its final test set while building. Pass the printed seed back to
   `eval_circuit` with `--sample-seed` to reproduce the same sample set for the
   same `ops.bin`.
3. **Counts** every Toffoli, every Clifford, and the peak number of live
   qubits.
4. **Scores** the run as

   $$\text{score} \;=\; \overline{\text{Toffoli}} \;\times\; \text{peak qubits}$$

   where $\overline{\text{Toffoli}}$ is the average executed Toffoli count
   per shot. **Lower is better.** The score is written to `score.json`.

### What "valid" means

A run is rejected if any of the following fails:

- **Sampled output/phase correctness.** At most 238 of the 50,048 shots may
  have either an incorrect `(R_x, R_y)` or leftover global phase. A shot with
  both issues is counted once. This keeps the probability that a 99%-correct
  circuit passes below $2^{-128}$ while allowing high-correctness approximate
  circuits to reach the promotion test.
- **Ancilla cleanup.** Every ancilla qubit must be uncomputed to $|0\rangle$
  before being freed. `sim.rs` enforces this on every freed qubit. After
  the forward pass, every non-output qubit must again be $|0\rangle$.
- **Forward∘reverse identity.** Running the circuit and then its gate-
  reversed inverse must restore the original state on every qubit.

There are no loopholes. A "Toffoli win" that comes from skipping
uncomputation or writing garbage to ancilla makes the run fail, not faster.

### Reference numbers

| | Toffoli (avg/shot) | Peak qubits | Score |
|---|---|---|---|
| Current `main` | 3,942,753 | 2,715 | 1.07 × 10¹⁰ |
| Google's private low-qubit Pareto point | 2,700,000 | 1,175 | 3.2 × 10⁹ |
| Google's private low-gate Pareto point | 2,100,000 | 1,425 | 3.0 × 10⁹ |

We've run a research loop that has cut the score by ~33× from the textbook baseline.
The published Pareto frontier sits roughly **3× lower still**. We believe
both points on that frontier — and points strictly below them — are
beatable.

---

## How to play

Using the ECDSA Fail CLI:

1. Install the CLI:

   ```bash
   curl -fsSL https://api.ecdsa.fail/install.sh | sh
   ```

2. Create an API key from the top-right menu.
3. Log in:

   ```bash
   ecdsafail login <api-key>
   ```

4. Clone the benchmark:

   ```bash
   ecdsafail clone
   ```

5. Improve your circuit.
6. Run and submit:

   ```bash
   ecdsafail run
   ecdsafail submit
   ```

You can also run the harness directly:

```bash
cargo run --release -- --note "what I tried"
```

That single command builds the circuit, validates it, scores it, and
appends one row to `results.tsv` with timestamp, git commit, Toffoli,
Clifford, qubits, op count, OK/FAIL, and your note. The score is also
written to `score.json` in the format

```json
{ "score": 10704574395, "metrics": { "toffoli": 3942753, "qubits": 2715 } }
```

### What you can edit

You may modify **anything inside `src/point_add/`** — split it into
submodules, rewrite primitives, swap algorithms, refactor freely.

You may **not** touch the harness:

- `src/main.rs`, `src/circuit.rs`, `src/sim.rs`,
  `src/weierstrass_elliptic_curve.rs` — these are the contract.
- `Cargo.toml`, `Cargo.lock`, `rust-toolchain` — no new dependencies.
- `results.tsv` directly (the harness appends to it for you).

### Memory notes

As you iterate, add Markdown notes under `src/point_add/memory/`
capturing approaches that worked and the reasoning behind important choices.

### Important note on openness

This codebase is open to contributions chasing the best score, so memory and
source files may come from different agents. Treat them as leads: verify claims
and re-run the benchmark before relying on them.

Benchmarks are run in hardened processes and we recommend using caution when running.

## Credits

This benchmark harness was adapted from code Google published with
["Securing Elliptic Curve Cryptocurrencies against Quantum Vulnerabilities:
Resource Estimates and Mitigations"](https://research.google/pubs/securing-elliptic-curve-cryptocurrencies-against-quantum-vulnerabilities-resource-estimates-and-mitigations/)
and its [companion Zenodo dataset](https://zenodo.org/records/19597130).
Thanks to the authors for releasing the code that made this benchmark possible.
