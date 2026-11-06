import json
from pathlib import Path
import tempfile
import unittest

from summarize_frozen_transfers import verify_observer


class ObserverParity(unittest.TestCase):
    def test_accept_identical_reject_changed_output(self):
        with tempfile.TemporaryDirectory() as temporary:
            a, b = [Path(temporary)/name for name in ('control', 'capture')]
            for p in (a, b):
                p.mkdir()
                (p/'metrics.json').write_text(json.dumps(
                    dict(status='PASS', calls=3, bytes=1234, raw_bytes=5678)))
                (p/'solution.csv').write_text('global_dof,value\n0,0.125\n1,0.25\n')
                (p/'convergence.csv').write_text(
                    'solver_calls,relative_residual,interface_jump\n3,1e-8,0\n')
            self.assertTrue(verify_observer(a,b)['numerically_identical'])
            (b/'solution.csv').write_text('global_dof,value\n0,0.125\n1,0.25001\n')
            with self.assertRaisesRegex(ValueError,'parity'):
                verify_observer(a,b)
            (b/'metrics.json').write_text(json.dumps(
                dict(status='PASS', calls=4, bytes=1234, raw_bytes=5678)))
            with self.assertRaisesRegex(ValueError,'calls'):
                verify_observer(a,b)


if __name__ == '__main__':
    unittest.main()
