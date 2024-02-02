OUTDIR=out/

setup:
	mkdir $(OUTDIR)
.PHONY: setup

clean:
	rm -rf $(SRCDIR)__pycache__
.PHONY: clean


