import os
import pytest
import numpy as np
from uegccdutils import Dataset

dir = os.path.basename(__file__)[5:-3] # infer directory from test name

def test_output(data_regression):
	data = Dataset(dir).test_summary()
	data_regression.check(data)

def test_fort80_consistency():
	data = Dataset(dir)
	for fctr in data.CCD_structure_factors:
		sum = fctr["exchange S(G)"] \
			+ fctr["non-exchange S(G)"]
		assert np.allclose(sum, fctr["S(G)"])
