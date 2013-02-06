chars:
	@echo "Chars out of 2366:"
	@tail -n +8 README.md | wc -m

spell-check:
	aspell check README.md

strip:
	python strip.py

.PHONY: chars strip spell-check
