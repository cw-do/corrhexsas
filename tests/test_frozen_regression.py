"""Independent frozen synthetic curves, generated before this package existed."""
import json
from pathlib import Path
import numpy as np
from threadpoolctl import threadpool_limits
from corrhexsas import Evaluator, Numerics, Structure, Disk


def test_frozen_manuscript_synthetic_anchors():
    fixture = json.loads((Path(__file__).parents[1]/"validation/golden_synthetic.json").read_text())
    with threadpool_limits(limits=1):
        engine = Evaluator(fixture["q"], fixture["sigma_q"], Numerics(**fixture["numerics"]))
        for case in fixture["cases"]:
            prediction = engine.evaluate(Structure(**case["structure"]), Disk(case["radius"]))[0]
            np.testing.assert_allclose(prediction, case["intensity"], rtol=2e-10, atol=1e-12)
