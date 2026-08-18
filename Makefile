# Pipeline driver. Each target is one pipeline step; targets are ordered
# by dependency, so `make all` runs the whole analysis from an empty clone.
#
# Override the interpreters if they are not on PATH:
#   make all PYTHON=/path/to/python RSCRIPT=/path/to/Rscript

PYTHON  ?= python
RSCRIPT ?= Rscript
DRAWS   ?= 1000
JOBS    ?= 4

.PHONY: all download harmonize signatures gate score null axis ceiling \
        meta headline figures forest test verify check smoke clean \
        clean-results help

all: download harmonize signatures gate score null axis ceiling meta \
     headline figures check

download:                       ## fetch ICB cohorts from ORCESTRA (~0.9 GB)
	$(RSCRIPT) R/00_download_cohorts.R

harmonize:                      ## MAE -> flat expression + clinical tables
	$(RSCRIPT) R/01_harmonize.R

signatures:                     ## assemble the 102-signature library
	$(PYTHON) scripts/01_build_signatures.py

gate:                           ## apply the pre-specified cohort gate
	$(PYTHON) scripts/02_cohort_gate.py

score:                          ## score every signature in every cohort
	$(PYTHON) scripts/03_score_signatures.py

null:                           ## the random-signature null (slowest step)
	$(PYTHON) scripts/04_null_calibration.py --draws $(DRAWS) --jobs $(JOBS)

axis:                           ## single-axis collapse + ESTIMATE confound
	$(PYTHON) scripts/05_axis_and_confound.py

ceiling:                        ## leave-one-cohort-out transferable ceiling
	$(PYTHON) scripts/06_transferable_ceiling.py

meta:                           ## random-effects meta-analysis per signature
	$(RSCRIPT) R/07_meta_analysis.R

forest:                         ## forest plots for the canonical predictors
	$(RSCRIPT) R/07_meta_analysis.R --forest \
	  T_cell_inflamed_GEP_Ayers_et_al,IFNG_6gene_Ayers,Cytolytic_Activity_Rooney_et_al

headline:                       ## collect every manuscript number to JSON
	$(PYTHON) scripts/08_headline_numbers.py

figures: forest                 ## regenerate manuscript figures
	$(PYTHON) scripts/09_figures.py

check:                          ## reconcile rerun numbers vs the manuscript
	$(PYTHON) scripts/10_check_manuscript.py

smoke:                          ## fast end-to-end check (100 draws)
	$(MAKE) harmonize signatures gate score axis ceiling meta headline \
	  DRAWS=100
	$(PYTHON) scripts/04_null_calibration.py --draws 100 --jobs $(JOBS)
	@echo "smoke run complete -- numbers are NOT publication values"

test:                           ## unit tests (no data download needed)
	$(PYTHON) -m pytest tests/ -q

verify:                         ## check the IOBR roster against upstream
	$(PYTHON) scripts/verify_roster.py

clean-results:                  ## remove generated results and figures
	rm -rf results/* figures/*.png
	touch results/.gitkeep figures/.gitkeep

clean: clean-results            ## also remove downloaded and harmonized data
	rm -rf data/raw/* data/harmonized/*
	touch data/raw/.gitkeep data/harmonized/.gitkeep

help:
	@grep -E '^[a-z-]+:.*?## .*$$' $(MAKEFILE_LIST) \
	  | awk 'BEGIN {FS = ":.*?## "}; {printf "  %-12s %s\n", $$1, $$2}'
