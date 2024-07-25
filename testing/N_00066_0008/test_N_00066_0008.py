import os
import pytest
from uegccdutils import Dataset

dir = os.path.basename(__file__)[5:-3] # infer directory from test name

def test_output(data_regression):
	data = Dataset(dir).test_summary()
	data_regression.check(data)