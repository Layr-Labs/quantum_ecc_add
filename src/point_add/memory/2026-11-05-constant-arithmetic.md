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

## Subsequent promoted frontier and measured controlled arithmetic

Further promotions in this session:

- `9f3aeba3-a6c6-40c4-93a8-7969f054d94c`: 3,943,905 Toffolis,
  2,711 qubits, score 10,691,926,455 (overflow bit used as reduction control).
- `9d0e5b43-42ce-44ea-8816-6d7810210eaf`: 3,943,927 Toffolis,
  2,710 qubits, score 10,688,042,170 (zero-carry and constant-only adders).
- `c0194a44-6c4d-4ec2-af62-d8ebab8e0e6f`: 3,940,907 Toffolis,
  2,588 qubits, score 10,199,067,316 (zero-extended Karatsuba merge).
- `44efa26f-fc51-49cb-96d3-27c80fce1015`: 3,941,951 Toffolis,
  2,565 qubits, score 10,111,104,315 (borrowed low-half sums, segmented merge).
- `dabbec86-75aa-4a53-9998-dc713f92c373`: 3,942,299 Toffolis,
  2,483 qubits, score 9,788,728,417 (measured controlled inversion arithmetic).
- `451b38a0-52e4-47f1-804e-ecaf79763071`: 3,938,207 Toffolis,
  2,480 qubits, score 9,766,753,360 (fast post-shift modular addition).

The controlled-ripple structural source is Algorithm 4 of Remaud & Vandaele,
https://arxiv.org/html/2501.16802v2 . The local measurement adaptation saves
original nonlinear ladder targets, measures modified targets on cleanup,
corrects the saved-value Z and control-product CZ phases, and restores the
input wires from their saved values. It uses n-1 temporary qubits and 2n-1
Toffolis per n-bit controlled modular-power-of-two addition. Subtraction uses
complement/add/complement. Its integration replaces both controlled operations
of Kaliski step 4 and their backward counterparts. Independent exhaustive and
large-width tests, including zero-before-reset assertions, pass.

A fully unitary ancilla-free version also passed the complete validator but
scored 11,396,376,563 (4,589,761 Toffolis, 2,483 qubits), so it was not retained.

The shorter-addend merge can use a segmented path: add into the low slice plus
an overflow, conditionally increment the high tail, then compare low output
against the unchanged addend to clear the overflow. This trades a comparator
for fewer simultaneous carry wires. In-place middle-product half-sums borrow
the low halves and restore them with subtraction after the middle product.
The complete inherited multiplication kernel requires distinct input wires;
it is not a squaring primitive when x and y are literally the same register.

The shift spill additions can likewise use shorter-addend arithmetic instead
of full zero padding. Generic modular shifts may contain measurements; retain
the separate phase-clean shift implementations for callers requiring pure
unitary gate reversal. The full shift roundtrip tests check multiplication by
2^k modulo p, exact reverse identity, phase, and pre-reset zero values.

## Concrete inherited termination witness

A scratch diagnostic evaluator (not a harness edit) examined the failed direct
constant-wrapper stream. It found nonzero qubit 1021 immediately before reset
at batch 113, operation 19,975,038, lane mask 0x10000000. Builder metadata maps
that qubit to v_w[0] at the end of the inherited 400-step second inversion.
Thus that stream had an incomplete inversion cleanup. This explains a class
of rare phase failures with correct coordinates and final zero ancillas.
Successful 9,024-shot runs do not prove those inherited iteration limits safe
for all inputs. No failing variant was submitted; the public notes explicitly
record this limitation. A stronger completion bound remains future work.

## Adaptive workspace and a larger completion margin

`0d343168-1229-4a10-9734-6307707664c9` promoted the shorter spill shift at
3,918,719 Toffolis, 2,479 qubits, score 9,714,504,401 (source `bed83c8`).

A subsequent local candidate passes at 4,067,987 Toffolis, 2,350 qubits,
score 9,559,769,450 with both inversion limits increased to 408. It chooses
pure or measured adders/comparators from static builder occupancy, using a
2,350-qubit soft geometry budget. This is an algorithm-selection budget,
not an input-dependent branch or mathematical truncation threshold.

Sparse constants can be segmented into a low addition and controlled tail
increment. For a controlled increment, the low-block carry-out can be erased
using control AND (low result == 0). For decrement use control AND (low result
== all ones). The existing measurement-clean zero test supplies that predicate
without loading a full comparison constant. Recurse only through the intended
high correction wrappers; keep core ripple functions available to avoid
unbounded decomposition. Dense constants fall back to their core arithmetic.

The sparse negation identity is `p-v = ~(v + (2^256-p)-1)`. This avoids treating
p+1 as a dense addend. It preserves the predecessor's noncanonical p output
when v=0, which its callers already expect.

Lambda is zero after the normal second-inversion multiplication and cleanup.
Free it there, before the backward inversion, rather than retaining 256 zero
wires to the end. Guard the final free to avoid freeing its recycled IDs twice.
Other experimental architecture branches retain their original final cleanup.

Noncompetitive but passing local comparisons:

- Sparse unconditional constants alone: score 9,717,694,874 at 2,479 qubits.
- 2,450 soft budget, 408 steps: score 9,805,741,294 at 2,461 qubits.
- 2,350 budget with blocked increment but no adaptive controlled ladder or
  early lambda release: score 9,892,136,511 at 2,459 qubits.

The 408-step setting is a larger fixed finite margin than 404/400. It is not
the theoretical 512-step bound and does not establish universal completion.
Independent tests now exercise both low-occupancy measured and high-occupancy
unitary routes, sparse controlled constants, high-occupancy field operations,
modular shifts, and strict zero-before-reset assertions.

## Direct classical offsets and clean scaling

Subsequent promoted steps include direct classical-offset field arithmetic,
classical preprocessing of 2*Qx and 3*Qx, and removal of copied pad/carry wires
in the multiplication kernel. Their public notes contain the full results.
The carry-one and biased-kernel tests preserve the multiplier's actual
control-zero increment 2^n-x, not an uncorrected signed subtraction.

The clean scaling implementation in clean_scale.rs removes the Bennett spill
and reduction flags from an in-place power-of-two field shift. For
p=2^256-c, output R, quotient q, and original captured spill s:

- R = 2^k*v - q*p, with 0<=q<2^k.
- R mod 2^k = c*q mod 2^k, so q follows from the inverse of odd c.
- s = q - delta, with delta = [R < c*q].

The forward helper turns spill into q, erases delta by comparing R with c*q,
and subtracts inverse(c)*R_low from q modulo 2^k. The reverse helper computes
q and delta from its input, reconstructs spill/flag, and invokes the existing
explicit reverse shift. Product temporaries are computed and uncomputed
with the sparse c = 1+2^4-2^6+2^10+2^32 decomposition. No target values are
read by the builder; q is computed by quantum-controlled constant arithmetic.

The tested implementation uses blocks up to 64 bits and replaces only lambda
normalization walks. It does not replace coefficient doubling inside Kaliski.
Independent tests include zero, p-1, values near p/2, and inputs mapping to
small outputs such as one, so the reduction flag's rare branch is exercised.
The 408-step normalization and its reverse are checked with exact classical
field references, zero phase, all ancillas zero, and zero-before-reset checks.
