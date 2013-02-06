chars:
	@echo "Chars out of 2366:"
	@tail -n +8 README.md | wc -m

.PHONY: chars
