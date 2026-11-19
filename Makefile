install:
	pip install -r requirements.txt
	pip install -e .

test:
	pytest

synthetic:
	python scripts/run_synthetic_demo.py --config configs/base.yaml

clean:
	find . -type d -name __pycache__ -prune -exec rm -rf {} +
	rm -rf .pytest_cache build dist *.egg-info src/*.egg-info
