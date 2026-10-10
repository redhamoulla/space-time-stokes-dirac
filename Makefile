PYTHON ?= python

.PHONY: verify verify-original verify-extension

verify: verify-original verify-extension

verify-original:
	$(PYTHON) experiments/verify_transmission_line_v12.py

verify-extension:
	$(PYTHON) experiments/ph_gauges.py
	$(PYTHON) experiments/verify_stage_geometry.py
	$(PYTHON) experiments/validate_ph_extension.py


