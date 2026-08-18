.PHONY: verify paper clean

verify:
	python experiments/verify_transmission_line_v12.py

paper:
	mkdir -p build
	latexmk -pdf -interaction=nonstopmode -halt-on-error \
		-outdir=build manuscript/moulla_2026_space_time_stokes_dirac.tex

clean:
	latexmk -C -outdir=build manuscript/moulla_2026_space_time_stokes_dirac.tex

