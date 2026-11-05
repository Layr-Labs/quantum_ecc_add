# Constant arithmetic and workspace experiments, 2026-11-05

Promoted results:

- `7f9be43e-057a-43dc-ace8-0488da946c61`: 3,943,928 Toffolis,
  2,712 qubits, product 10,695,932,736. Source `d6e4800`.
- `629b3e92-6f74-43a6-a04b-b8242c978b2e`: same Toffolis,
  2,711 qubits, product 10,691,988,808. Source `c0a5567`.

Both passed the full 9,024-shot local validator and remote scoring. See their
public submission notes for the detailed experiment history and limitations.
The original starting frontier was `d9f34597-f6d9-4b97-ae9a-2fc507ad3b80`,
source `c4d7d58`: 3,942,753 Toffolis, 2,714 qubits, product 10,700,631,642.

## Useful changes

The direct controlled-constant add and subtract routines originally used three
Toffolis to calculate a majority at a set constant bit. The identity
`maj(a,d,c) = ((a XOR c) AND (d XOR c)) XOR c` reduces this to one Toffoli,
with CNOTs to transform and restore inputs. The existing measurement phase
correction still describes the same Boolean majority function.

Using these direct routines avoids loaded constant registers. The meaningful
peak changes required updating the zero-accumulator modular-add specialization,
not just the ordinary modular adder. Direct modular negation removes another
large constant register and exposes inversion workspace as the next peak.

Basis-constant propagation over the final stream removes hundreds of thousands
of identity gates. It treats quantum inputs as unknown, does not predict
measurement outcomes, keeps all phase gates and measurements, and conservatively
invalidates target values under conditional operations. Reset checks stay in
the output. This pass operates only on symbolic constants, never test data.

Kaliski history can be allocated one qubit per forward step and freed one per
backward step. History lifetime changes alone did not lower the measured peak.
Combining them with step-4 workspace borrowing did: the existing denominator
bound leaves `v_w[width..]` at zero, where `width=n` before iteration n and
`width=2*n-iteration` afterward. Step 4 accesses only the low v_w slice, so
those high wires can supply part of its clean temporary register. Uncompute
the full logical temporary, but free only its separately allocated owned
subset. The borrowed wires remain part of v_w and must be zero again before
later denominator operations.

## Negative results and correctness cautions

Disjoint-wire self-inverse cancellation removed 49,900 Clifford gates and no
Toffolis. A bounded commuting-X-family extension found no additional Toffolis.
Neither pass was retained.

Direct halving without majority optimization added 8,964 Toffolis. The
one-Toffoli majority rewrite eliminated that overhead. Many direct-arithmetic
integration steps reduced Clifford count without moving the product peak.

Some variants with the inherited small-r shortcut ending at iteration 262
failed phase validation in one batch despite zero coordinate and ancilla
failures. The promoted variant ends that shortcut at 256, consistent with its
stated n-step bound. It retains the inherited pair iteration and bulk-prefix
settings; finite-shot validation is not a proof that those settings cover all
inputs. Do not infer universal phase correctness from one passing circuit
hash. No failing variant was submitted as an improvement.

A conservative reference with 512 inversion steps, a 256-step bulk prefix,
and small-r threshold 256 passed all 9,024 shots but scored 13,570,521,082
(4,808,831 Toffolis, 2,822 qubits). Enabling dirty-workspace arithmetic on that
reference also passed but worsened the score to 14,800,856,642. These are useful
correctness comparisons, not promoted implementations.

The comparator carry-in copy can be removed: the first MAJ carry and its phase
correction can use the restored first input bit directly. This passed the full
validator at the same 2,711-qubit score, before later overflow-control work.

Do not edit the harness or score data to pursue these experiments. Keep scratch
backups outside the checkout. Reproduce with the pinned Rust toolchain and the
normal sandboxed `yukon run` pipeline; remote scoring remains authoritative.
