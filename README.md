# bbgoat-sv-scripts

Analysis scripts for the Data Note "Structural variants from two Black Bengal goat PacBio HiFi genomes".

- Reference: ARS1.2 (GCF_001704415.2). Animals: 803M (male, run D01) and D863F (female, run C01).
- Raw HiFi reads: SRR41146742 (803M) and SRR41146743 (D863F), BioProject PRJNA1502507. The reads are on hold and will be released at acceptance of the manuscript. The SV call sets will be released under CC0 in GigaDB at acceptance.

## Contents
- `scripts/`: the numbered scripts, committed exactly as run (first commit). The second commit only quotes `rm` arguments. Each script's header comment states its purpose and rules.
  Files: 01_align_hifi.sh, 05_control.sh, 05_novelty.sh, 28_hifi_vs_shortread.sh, 29_hifi_events_depth.sh, 30_crosstab.sh, 31_novelty_strata.sh, 32_asm_concordance.sh
- `COMMANDS.md`: calling and merging steps that were run by hand and have no script.
- `tests/smoke.sh`: syntax and lint check only. It does not test the analysis.

## Limits
Paths are as run on one server (project root /data/bbg_pangenome) and are not portable without editing. The scripts need the inputs listed in the Data Note and the conda environments named there.

Licence: MIT.
