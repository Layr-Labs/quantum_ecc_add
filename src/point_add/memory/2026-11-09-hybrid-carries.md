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
