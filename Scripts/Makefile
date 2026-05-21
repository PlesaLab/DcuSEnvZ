#===============================================================================
SHELL := /bin/bash

# Define Directories
PROJ_DIR := /gpfs/projects/plesalab/llippert/DcuS_38_nextseq
INPUT_DIR := $(PROJ_DIR)/GC3F-LL-11236
REF_DIR := $(PROJ_DIR)/refs
CUT1_DIR := $(PROJ_DIR)/cut1
FINAL_BC_DIR := $(CUT1_DIR)/Final_BC

# Create Directories
$(shell mkdir -p cut1/Final_BC)

# Define Variables
SORT := sort -T ./ --parallel=40

# Define Software
HTS_PRIMERS := ./hts_Primers
STARCODE := ./starcode
BBMERGE := bbmerge.sh

# Input Files
PAIRS := PreSort_1:S37 PreSort_2:S38 Undetermined:S0

#===============================================================================

# Primary targets
all: starcode

# Declare phony targets
.PHONY: all starcode

# Use secondary expansion for pattern rules
.SECONDEXPANSION:

# Generate lists of output files for all samples
MERGED_READS := $(foreach pair, $(PAIRS), $(PROJ_DIR)/$(subst :,_,${pair})_mrg.fastq.gz)
TRIMMED_READS := $(foreach pair, $(PAIRS), $(PROJ_DIR)/$(subst :,_,${pair})_mrg_flip_trim_SE.fastq.gz)
CUT_SEQS := $(foreach pair, $(PAIRS), $(CUT1_DIR)/$(subst :,_,${pair})_mrg_flip_trim_SE.cut.txt)
COUNT_FILES := $(foreach pair, $(PAIRS), $(CUT1_DIR)/$(subst :,_,${pair}).count.txt)
FIX_COUNT_FILES := $(foreach pair, $(PAIRS), $(CUT1_DIR)/$(subst :,_,${pair}).fix.count.txt)

# Ensure intermediate files are not deleted
.SECONDARY: $(MERGED_READS) $(TRIMMED_READS) $(CUT_SEQS) $(COUNT_FILES) $(FIX_COUNT_FILES)

#===============================================================================
# Pattern Rules

# Step 1: Merge reads (using bbmerge)
$(PROJ_DIR)/%_mrg.fastq.gz:
	@lib=$$(echo $* | cut -d_ -f1); \
	sample=$$(echo $* | cut -d_ -f2); \
	$(BBMERGE) \
	    in1=$(INPUT_DIR)/$${lib}_$${sample}_R1_001.fastq.gz \
	    in2=$(INPUT_DIR)/$${lib}_$${sample}_R2_001.fastq.gz \
	    out=$@ \
	    outu1=$(PROJ_DIR)/$${lib}_$${sample}_R1_unmerged.fastq.gz \
	    outu2=$(PROJ_DIR)/$${lib}_$${sample}_R2_unmerged.fastq.gz

# Step 2: Primer trimming (using hts_Primers)
$(PROJ_DIR)/%_mrg_flip_trim_SE.fastq.gz: $(PROJ_DIR)/%_mrg.fastq.gz
	@lib=$$(echo $* | cut -d_ -f1); \
	sample=$$(echo $* | cut -d_ -f2); \
	$(HTS_PRIMERS) \
	    -U $< \
	    -f $(PROJ_DIR)/$${lib}_$${sample}_mrg_flip_trim \
	    -P $(REF_DIR)/HKABCstub_FWD.fasta \
	    -Q $(REF_DIR)/HKABCstub_REV.fasta \
	    -l 5 -x -e 6 -d 6 \
	    -L $(PROJ_DIR)/$${lib}_$${sample}_flip_stats.txt -F

# Step 3: Extract sequences (24 bp barcode)
$(CUT1_DIR)/%_mrg_flip_trim_SE.cut.txt: $(PROJ_DIR)/%_mrg_flip_trim_SE.fastq.gz
	gzip -cd $< | awk 'NR%4 ==2' | cut -c1-24 | \
	    awk 'BEGIN{FIELDWIDTHS="24"} {print $$1}' > $@

# Step 4: Count extracted sequences
$(CUT1_DIR)/%.count.txt: $(CUT1_DIR)/%_mrg_flip_trim_SE.cut.txt
	cat $< | $(SORT) | uniq -c | \
	    awk '{print $$2 "\t" $$1}' > $@

# Step 5: Prepare starcode input
$(CUT1_DIR)/%.fix.count.txt: $(CUT1_DIR)/%.count.txt
	tail -n +2 $< > $@

# Step 6: Run starcode (using Levenshtein distance = 1)
starcode: $(FIX_COUNT_FILES)
	@for pair in $(PAIRS); do \
	    lib=$$(echo $$pair | cut -d: -f1); \
	    sample=$$(echo $$pair | cut -d: -f2); \
	    $(STARCODE) -d1 --sphere \
	        -i $(CUT1_DIR)/$${lib}_$${sample}.fix.count.txt \
	        --output $(FINAL_BC_DIR)/$${lib}_$${sample}_collapse_d1.tsv \
	        --print-clusters & \
	done; wait