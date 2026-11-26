# VirMarkDB

<p align="center">
  <img src="website/assets/logo_virmark.png" alt="VirMarkDB logo" width="600">
</p>



**VirMarkDB** is a reference database providing protein sequences of viral marker genes covering multiple viral groups, including *Nucleocytoviricota*, *Orthornavirae* and *Pantevenvirales*.

It is designed for viral marker detection, marker-based phylogeny, and taxonomic assignment workflows using manually curated viral marker-gene references and profile HMMs.

## Nextflow workflows

`main.nf` orchestrates three named subworkflows under `subworkflows/`:

| Subworkflow | Role | Default |
|---|---|---|
| `DATABASE_GENERATION` | Manifest, references, genomes, Prodigal, HMM searches and final database | Enabled |
| `ASSESSMENT` | NR pools, self-alignment, filtering, database comparison and missing-protein clustering | Disabled |
| `STATISTICS` | Database summaries, marker recovery, lengths, scores and taxonomic coverage | Enabled |

```bash
# Generate the database and statistics.
nextflow run main.nf -profile conda -resume

# Include assessment against an existing local BLAST NR database prefix.
nextflow run main.nf -profile conda --run_assessment \
  --nr_database /shared/bank/nr/blast/nr -resume

# Generate the database without statistics.
nextflow run main.nf -profile conda --run_statistics false -resume
```

`main.nf` is the pipeline entry point. Scheduler and software environments are
controlled by Nextflow profiles, not nested `sbatch` or `module load` calls.

Subworkflows exchange task outputs directly rather than reading published
folders. Generation emits the database, manifest, marker map and genome-marker
map. Assessment emits analysis and cluster results. Statistics emits its
`analysis/` directory. Results are published beneath `output/database_generation`,
`output/assessment` and `output/statistics`, respectively.

Assessment requires access to the NR database on the execution host and downloads
NCBI BLAST taxonomy data. Extraction rules come from `input/assessment_rules.tsv`,
joined to the cleaned mapfile by `marker_group_id`. The supplied rules preserve
the original title patterns and length ranges for the ten curated markers.
Only markers present in the mapfile are extracted; unused rules are ignored.
Missing rules, duplicate selected rules, inconsistent groups and invalid numeric
limits stop the task before NR access.

| Rule column | Meaning |
|---|---|
| `marker_group_id` | Exact identifier from the mapfile |
| `group_id` | Must match the mapfile group |
| `taxid` | NCBI taxon to extract from NR |
| `include_pattern` | Case-insensitive AWK regular expression for protein titles |
| `exclude_pattern` | Case-insensitive AWK regular expression for unwanted titles |
| `min_length` | Minimum protein length, inclusive |
| `max_length` | Maximum protein length, inclusive |

To assess a new marker, add it to the mapfile and add a matching row to the rule
TSV; no script changes are needed. All seven fields are required. Column order
may vary, and multiple mapfile rows for the same marker/group share one rule.
Use `--assessment_rules /path/to/rules.tsv` to supply another rule file. The
selected rules are copied to `pool/metadata/selected_rules.tsv` for traceability.

 Self-alignment
retains proteins with a non-self match at identity >= 0.3, alignment length >= 100
and E-value <= 1e-5. MMseqs explicitly reports fractional identity. Missing
proteins are clustered with CD-HIT at 95% identity. Empty pools and searches
produce empty outputs without invoking alignment or clustering tools.

Statistics derive marker definitions and eligible genomes from the generated
maps; markers with no retained sequences receive zero recovery. Denominators
cover all manifest-eligible genomes, including those outside the current
500-accession download limit. Generic family coverage plots replace the old
plots hardcoded for three viral groups. Website source files are not overwritten;
tables and PDF/PNG figures are published in `output/statistics/analysis/`.

## Nextflow genome processing

```bash
nextflow run main.nf -profile conda --genome_batch_size 20 \
  --private_genomes /path/to/private_genomes -resume
```

The manifest provides each genome's `virus_id` and `Source`. NCBI genomes are
fetched with Entrez Direct in batches; `genome_download_parallelism` controls
the number of concurrent download tasks (1 to 3, default 3). Downloads wait for
reference collection to finish, so Entrez Direct and rentrez do not compete
for the same NCBI request budget. Each download worker pauses between genomes;
reference fetching retries failures up to three times. Other users sharing
the same IP can still cause NCBI rate limits. The current `main.nf` limits this
branch to the first 500 unique accessions. Private genomes are read from
`<private_genomes>/<virus_id>.fna`.

Private genome processing is optional and disabled by default. Without
`--private_genomes`, `PRIVATE_GENOMES` is not invoked and only the NCBI branch
is processed, even if the manifest lists private genomes. With the option,
only manifest rows with `Source = PRIVATE` are processed. If there are none,
no private task runs; Nextflow may display `-` for that process. A listed
private genome with a missing FASTA raises an error.

Both branches emit genome IDs alongside their files. `main.nf` combines these
`(virus_id, genome_file)` pairs and passes batches to Prodigal as
`(batch_id, genome_ids, genome_files)`. IDs come from the manifest throughout;
Prodigal reads the staged `<virus_id>.fna` files directly.

All results are published to `<outdir_generation>/viruses/<virus_id>/`:

- `<virus_id>.fna`: input genome;
- `<virus_id>.faa`: predicted proteins;
- `<virus_id>.ffn`: predicted coding sequences;
- `<virus_id>.gff`: annotations.

`outdir_generation` defaults to `output/database_generation`. There are no batch
or intermediate directory names in the published paths. Prodigal output
channels retain the batch ID and genome IDs with their corresponding files.
Use `-resume` to reuse completed tasks.

`COLLECT_REF` uses the cleaned marker map emitted by `GENERATE_MANIFEST`.
The columns have distinct roles:

| Column | Example | Role |
|---|---|---|
| `marker` | `dnapol` | Protein/marker name |
| `group_id` | `nucleocytoviricota` | Viral group and reference source |
| `marker_group_id` | `dnapol_nucleocytoviricota` | Reference identifier, used for the FASTA filename |
| `manual_file` | `dnapol_nucleocytoviricota.tsv` | Optional accession file under `input/manual_reference_protein/` |

One loop in `collect_ref.R` reads each marker, loads its published reference
sequences, adds the accessions from `manual_file`, removes duplicate sequences,
and writes `references_protein/<marker_group_id>.fasta`. Only markers in the
map are prepared. Missing or empty manual files add no sequences. RdRp is a
special case: existing RdRp-scan profiles are downloaded instead of a FASTA.

Reference FASTAs, downloaded sources and alignments stay in Nextflow's `work/`
directories. `HMM_PROFILE` receives `(marker_group_id, reference_fasta)` pairs,
aligns each set with MAFFT (`--ep 0 --genafpair --maxiterate 1000`), then runs
`hmmbuild`. RdRp-scan bypasses this step because its profiles are already built.

Only final profiles are published under `outdir_generation`:

```text
hmms/
├── dnapol_nucleocytoviricota.hmm
├── g23_pantevenvirales.hmm
├── ...
└── rdrpscan/
    └── RdRp_HMM_profile_CLUSTALO.db.h3*  (four companion files)
```

`HMM_PROFILE.out.profiles` emits `(marker_group_id, hmm_file)` pairs;
RdRp-scan is available through `COLLECT_REF.out.rdrpscan`.

`GENERATE_MANIFEST` also writes `genome_marker_map.tsv`, using `taxo_rank` and
`taxo_value` to associate genomes with markers. After Prodigal, proteins from
all processed batches are pooled by marker. Empty protein files are excluded.
`HMM_SEARCH` runs `hmmsearch` with the generated profiles, or `hmmscan` with
the pre-pressed RdRp database. The RdRp reference ID comes from the mapfile.
The default threshold is `--hmm_evalue 1e-5`, matching the original scripts.

Searches publish only `search/<marker_group_id>.tbl` (`--domtblout`) and
`search/<marker_group_id>.log`. Protein pools remain in `work/`. A marker with
no proteins to search launches no task; a search with no hits still produces
a valid table containing HMMER comments. The HMMER target/query columns are
reversed between `hmmsearch` and `hmmscan`: downstream parsing must account
for the search tool. Results are available as `HMM_SEARCH.out.hits` and `.logs`.
Existing FASTAs or sources published by earlier runs are not automatically
removed.

`FINAL_RESULTS` gathers HMMER tables and Prodigal protein/CDS files after all
searches finish. It reuses `manifest_genomes.tsv` for genome metadata,
`marker_taxo_map.tsv` for marker definitions, and `genome_marker_map.tsv` to
validate genome/marker associations. These are the files emitted by
`GENERATE_MANIFEST`; taxonomy rules and marker IDs are not reconstructed. It keeps the best profile per ORF and the best ORF per genome
and marker, plus additional copies reaching 80% of the best score and length.
RdRp keeps one copy and exports the domain envelope. Selected proteins shorter
than 25% of their marker's median length are removed. Manual exclusions come
from `input/remove_specific_ORF_sequences.txt`.

Final outputs are published under `<outdir_generation>/VirMarkDB/`:

- `markers/<group_id>/<marker>/`: protein and nucleotide FASTAs and ORF tables;
- `virus_informations/`: marker composition, taxonomy and genome metadata;
- `logs/`: filtering counts, shared ORFs, family summaries and plots.

Coding sequences come from Prodigal `.ffn` files. ORFs are linked to manifest
IDs through their source FASTA filenames, including private genomes. ORF IDs
must be unique across genomes. Runs without hits still write empty metadata
and filtering tables; plots are generated only when markers are selected.
`FINAL_RESULTS.out.database` emits the complete database directory.

Reference collection requires network access to Zenodo, NCBI and GitHub.

## Documentation

Full documentation is available at:

**[https://virmarkdb.github.io/VirMarkDB/](https://virmarkdb.github.io/VirMarkDB/)**

The website provides information about:

- database content and taxonomic scope;
- marker genes and marker-group strategy;
- workflow and software versions;
- file formats and release structure;
- coverage statistics and analyses;
- contribution guidelines;
- citation, license and known limitations.

## Releases

Public release archives will be deposited on Zenodo.

| Item | Status |
|---|---|
| Release archive | v0.1.0 |
| Zenodo DOI | https://doi.org/10.5281/zenodo.23712240 |

## Citation

Citation metadata are provided in [`CITATION.cff`](CITATION.cff).

Estoup P, Baudoux A-C, Corre E, Blanc G & Henry N. 2026. *VirMarkDB: a reference database of viral marker-gene sequences*. Version 0.1.0. Zenodo. https://doi.org/10.5281/zenodo.23712240 

## License

VirMarkDB is distributed under the CC-BY-4.0 license. See [`LICENSE`](LICENSE).

## Contact

For questions, corrections, missing taxa, suspicious marker assignments or reproducibility issues, please use the GitHub issue tracker:

[https://github.com/VirMarkDB/VirMarkDB/issues](https://github.com/VirMarkDB/VirMarkDB/issues)
