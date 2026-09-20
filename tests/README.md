# Testing guide

This directory holds the maintained, automated test suite. Name test modules `test_<behavior>.py` and run them from the repository root with:

```sh
python -m unittest discover -s tests -p "test_*.py"
```

Tests must be self-contained and must not depend on ignored local data. Prefer deterministic inputs, state physical units explicitly, and use tolerant floating-point comparisons such as `numpy.testing.assert_allclose`.

When a test preserves knowledge from a `scratch/` experiment, state the claim and its context in the test docstring, then use a small fixture and assertions that demonstrate it. Keep longer explanations in `doc/` and link them from the test when useful. The test must remain part of the normal automated suite.

## Test data and assertions

Large production VR/ATOM files are not required for most checks. Test file handling in layers:

1. Use tiny synthetic fixtures that preserve the real header, grid, atom, and numeric-field structure, including boundary cases and malformed input.
2. Assert parsed metadata, shapes, dtypes, atom counts/order, and representative numerical values rather than only checking that a command completed.
3. Add read-write-read round-trip tests where supported and compare the resulting in-memory structures with documented tolerances.
4. Test domain invariants relevant to the change, such as grid dimensions, coordinate transforms, periodicity, atom identity/order, and conservation or normalization properties.
5. Keep any tracked golden fixture small. Compare normalized numerical content or a documented summary/checksum so irrelevant whitespace or formatting does not hide a numerical regression.

## Full-size integration checks

For changes that can only be validated with a full-size calculation, first cover the core behavior with a small automated test. Then run a representative integration calculation and report the input dataset, command, selected numerical results, and tolerances in the pull request. Never make the normal test suite depend on an untracked local dataset.
