# Partial measured Cuccaro carries

Starting frontier: 59a85cc6-bb2a-41ac-b90d-49fe31a248bf, source c1bb169.
Baseline reproduced at 3,845,263 Toffolis, 2,335 qubits, product 8,978,689,105;
17 release primitive tests passed.

The generic carry-in fast Cuccaro pair now allocates only the measured carry
prefix that fits WORKSPACE_SOFT_LIMIT. Remaining steps use in-place unitary
MAJ/UMA. Both forms leave identical transformed inputs and carry values at
step boundaries; every omitted measured carry adds one cleanup Toffoli.
Constant-one add/sub specializations route overflowing calls through this
hybrid pair with a prepared and restored carry-in qubit.

Lowering the original budget alone increased Toffolis and left the peak
unchanged. Changing the generic pair alone left the constant-one peak.
After both changes, a construction scan over budgets 2321..2360 selected
2333: 3,846,793 Toffolis, 2,333 qubits, predicted product 8,974,568,069.
2334 gives 8,977,224,522; 2332 gives 8,982,600,484. More aggressive budgets
lose on gates and can expose a different shift peak. Budget decisions use
static builder occupancy, never quantum data or test outcomes.

The biased controlled kernel test now forces every carry capacity from zero
to full for widths 1..5 and exhaustively checks arithmetic, preserved inputs,
zero phase, inverse identity and strict zero-before-reset behavior. The full
17-test suite passes. Full point-add validation and remote scoring remain
necessary before treating the predicted product as a promoted result.

Next direction: the add_nbit_qq_fast/sub_nbit_qq_fast wrappers still switch
entirely to unitary arithmetic when a full carry register does not fit.
They could instead prepare a zero carry-in and use the same hybrid pair.
Check width and cleanup before retaining that change. Existing inversion
limits remain 408 rather than a universal 512-step completion bound; this
carry-storage improvement does not resolve that inherited limitation.

Full local validation passed all 9,024 shots at the predicted 8,974,568,069,
with zero coordinate, phase or ancilla failures. Average executed Cliffords
are 20,245,223.775; emitted operations 28,134,869, classical bits 1,868,860.

The first hybrid candidate promoted as 92758ece-6022-4b20-a851-22673be9a355,
source 0cf07bf8945264dad21dcaba1c549277a76eee10. Normal frontier sync completed.
The ordinary add/sub wrapper follow-up passes 9,024 shots at 3,842,901
Toffolis, 2,333 qubits, product 8,965,488,033. It prepares a zero carry-in
and calls the hybrid pair instead of the all-unitary fallback; short addends
need no physical top padding. The seventeen-test suite passes with added
high-occupancy equal-width/short-addend coverage through 256 bits. A geometry
scan over 2240..2340 still selects 2333 among tested values.

The wrapper follow-up promoted as 49b0fd0f-5fc2-4363-bd12-8c900180ca29,
source 7e588e2d70b5eb6a798af0049a35750d01be9ea3, and normal sync completed.
Partial carry storage in cmp_lt_into_fast then passed 18 primitive tests and
all 9,024 full shots at 3,840,383 Toffolis, 2,333 qubits, 8,959,613,539.
It must retain the first scratch (the specialized first step controls u[0]);
only later steps can target u[i] directly. With no room, keep the old pure
unitary comparator. Exhaustive tests check every capacity and both initial
flag values at widths 1..6, input preservation, predicate XOR, double-call
identity, zero phase, zero workspace and strict pre-reset checks.

The partial comparator promoted as 5ff96417-3814-4097-8f89-7965d1c65f2e,
source 81c9fe8c7f19869f8af0bbf4d1e38a9402048b98; normal sync completed.
Removing the scoped comparator's copied carry input passes exhaustive tests
but changes neither Toffolis nor peak. That source experiment is cached,
not retained in the scoring candidate.

Budget 2331 alone predicts 3,843,131 Toffolis, product 8,958,338,361.
Blocked controlled increment/decrement can split off only enough low bits to
let the high suffix fit a full ripple: with room=L-active, require split>=n-room.
Cap by the existing available-prefix allowance and half-width progress bound.
This avoids unnecessary low-block overflow cleanup. Combined with budget2331,
18 primitive tests and all 9,024 full shots pass at 3,840,467 Toffolis,
2,331 qubits, product 8,952,128,577; zero coordinate/phase/ancilla failures.
The sparse controlled-constant test now varies room around full, half and
tight capacities, including cases that actually select the changed split.

Increment splitting promoted as cdedc44b-c527-4ded-9009-b40192d190b5,
source 49b0256aee3a33bfb952e905ccfa8993b395814f; normal sync completed.
A new geometry scan selects 2332 among 2270..2340 step-two plus 2331/2333.
The 2332 follow-up passes 9,024 shots at 3,838,706 Toffolis, 2,332 qubits,
product 8,951,862,392. It saves 1,761 Toffolis for one added qubit.
The proposed always-first-measured comparator fallback passes tests but
has no full-circuit resource effect here; excluded and cached separately.

Budget2332 promoted as 0d431d02-3ccd-484d-a989-5cf48bf72fb6, source
094d5af84ac3deb49a91aa0d8f8246cc95e2e91d, and normal sync completed.
Small clean shifts 2/4 are gate-expensive; larger shifts can lower gates but
add too much width at the first product reduction (2,203 live qubits).
New product_fold.rs instead folds raw 512-bit words twice with c=2^32+977.
First T=L+cH fits289 bits, q=T>>256<=c. Then S=T_low+cq<2p. Retain q33 and
the reduction bit across a callback, explicitly undo reduction and both folds.
It reuses clean_scale short_term and solinas_product, with compact66-bit cq.
At the first known-zero output multiplication, copy the canonical residue
then uncompute the fold and raw Karatsuba product. Independent arbitrary-limb
boundary/random tests at four occupancies pass, with exact raw-word restore,
phase, inverse XOR-copy and strict reset checks. Nineteen tests pass.
Full 9,024-shot validation passes at 3,824,279 Toffolis, 2,333 qubits,
product 8,922,042,907, zero coordinate/phase/ancilla failures. Next: modular
add/subtract callbacks at the other reduction sites with more spare width.
