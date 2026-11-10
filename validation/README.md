# Validation record

The regression fixture `golden_synthetic.json` contains eight q values
at each of three synthetic parameter anchors from the six-coordinate
study used in the manuscript. Those curves were generated and frozen
before the package was extracted. It contains no measured observations.
The original complete-shell cutoff was 256 and quadrature used
3072 transverse and 768 axial nodes up to 4000 A. The public test
evaluates the same finite representation, not a reoptimized curve.

The private workspace release audit H111b also compares the chemical
profile coefficients, the six retained reference-parameter-set curves,
and the axially constant and zero-displacement control curves. At matched
settings, all tested calculated curves agreed to a maximum relative
difference below 1.1e-14. This is implementation regression, not a claim
of physical accuracy or uniform numerical convergence. No measured-data
fit or PCA was rerun for this check.

Two earlier audit attempts are retained in the private analysis record.
H111 incorrectly compared the reference fit with the first stored control
row. H111a corrected the row selection but did not yet decode the archive's
null axial length as the constant-along-axis limit. H111b corrected those
validation bookkeeping issues. The package evaluator was unchanged.

The public test suite separately checks independent Gaussian resolution
integration and radial SLD quadrature, component amount conservation,
direct/series covariance agreement, zero-displacement invariance,
cache-free evaluation and recovery of a shared synthetic structure with
two contrast amplitudes and separate backgrounds. These do not prove
global identifiability, convergence for every possible input, or the
correctness of a physical interpretation for another sample.

The runnable example includes its seed and known truth. On the initial
validated environment (Python with NumPy 1.26.4 and SciPy 1.13.1), the
seeded example fit converged. Increasing cutoffs and nodes by 1.5 and
doubling the shell cutoff changed that curve by less than 0.00061 times
the supplied synthetic dI. This is one selected example and a combined
refinement check. It is not a universal default-setting error bound.
