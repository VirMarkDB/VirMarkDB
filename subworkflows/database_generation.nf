#!/usr/bin/env nextflow

include { GENERATE_MANIFEST }    from '../modules/database_generation/generate_manifest.nf'
include { DOWNLOAD_GENOMES }     from '../modules/database_generation/download_genomes.nf'
include { PRIVATE_GENOMES }      from '../modules/database_generation/private_genomes.nf'
include { PRODIGAL }             from '../modules/database_generation/prodigal.nf'
include { COLLECT_REF }          from '../modules/database_generation/collect_ref.nf'
include { HMM_PROFILE }          from '../modules/database_generation/hmm_profile.nf'
include { FINAL_RESULTS }        from '../modules/database_generation/final_results.nf'
include { HMM_SEARCH }           from '../modules/database_generation/hmm_search.nf'


workflow DATABASE_GENERATION {

    main:


    inputs = Channel.fromPath(params.general_input, checkIfExists: true)
    script_r = file("${params.scripts_generation}/generate_manifest.R")

    GENERATE_MANIFEST(inputs, script_r, params.ictv_version)

    // The cleaned marker map defines which references to collect.
    COLLECT_REF(
        GENERATE_MANIFEST.out.marker_map,
        file("${params.general_input}/manual_reference_protein", checkIfExists: true),
        file("${params.general_input}/reference_data", checkIfExists: true),
        file("${params.scripts_generation}/collect_ref.R", checkIfExists: true)
    )

    // Keep marker_group_id attached to each reference FASTA.
    reference_proteins = COLLECT_REF.out.proteins
        .flatten()
        .map { fasta -> tuple(fasta.baseName, fasta) }

    // Build profiles from reference FASTAs; RdRp-scan profiles are already built.
    HMM_PROFILE(reference_proteins)

    // Handle NCBI downloads and private genomes separately.
    manifest = GENERATE_MANIFEST.out.manifest
        .splitCsv(header: true, sep: '\t')

    genome_batches_ncbi = manifest
        .filter { it.Source != 'PRIVATE' }
        .map { it.virus_id }
        .unique()
        .take(500)
        .collate(params.genome_batch_size.toInteger())
        .toList()
        .flatMap { batches -> batches.withIndex().collect { ids, i -> tuple(i + 1, ids) } }

    // Wait for reference collection to avoid concurrent rentrez and efetch requests.
    DOWNLOAD_GENOMES(genome_batches_ncbi, COLLECT_REF.out.sources.first())

    // Private genomes are optional; do not invoke the process unless enabled.
    private_files = Channel.empty()
    if (params.private_genomes) {
        private_genomes = manifest
            .filter { it.Source == 'PRIVATE' }
            .map { row -> tuple(row.virus_id, file("${params.private_genomes}/${row.virus_id}.fna", checkIfExists: true)) }
            .unique()

        PRIVATE_GENOMES(private_genomes)
        private_files = PRIVATE_GENOMES.out.genomes
    }

    // Both sources emit (virus_id, genome_file) pairs.
    genomes = DOWNLOAD_GENOMES.out.genomes
        .flatMap { ids, files ->
            ids.collect { id -> tuple(id, files.find { it.name == "${id}.fna" }) }
        }
        .mix(private_files)

    // Pass genome IDs and their files to Prodigal in batches.
    genome_batches = genomes
        .collate(params.genome_batch_size.toInteger())
        .toList()
        .flatMap { batches ->
            batches.withIndex().collect { batch, i ->
                tuple(i + 1, batch.collect { it[0] }, batch.collect { it[1] })
            }
        }

    PRODIGAL(genome_batches)

    // The manifest defines which markers to search for each genome.
    genome_markers = GENERATE_MANIFEST.out.genome_markers
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(row.virus_id, row.marker_group_id) }

    proteins = PRODIGAL.out.proteins
        .flatMap { batch_id, ids, files ->
            ids.collect { id -> tuple(id, files.find { it.name == "${id}.faa" }) }
        }
        .filter { id, fasta -> fasta.size() > 0 }

    // Pool proteins from eligible genomes for each marker.
    protein_pools = proteins
        .combine(genome_markers, by: 0)
        .map { id, fasta, ref_id -> tuple(ref_id, fasta) }
        .groupTuple()

    // Choose the search tool from the reference type.
    rdrp_profiles = GENERATE_MANIFEST.out.marker_map
        .splitCsv(header: true, sep: '\t')
        .filter { row -> row.group_id == 'orthornavirae' && row.marker == 'rdrp' }
        .map { row -> row.marker_group_id }
        .unique()
        .combine(COLLECT_REF.out.rdrpscan)
        .map { ref_id, database -> tuple(ref_id, 'hmmscan', database) }

    search_profiles = HMM_PROFILE.out.profiles
        .map { ref_id, hmm -> tuple(ref_id, 'hmmsearch', hmm) }
        .mix(rdrp_profiles)

    HMM_SEARCH(protein_pools.combine(search_profiles, by: 0))

    // Wait for all searches and predictions, including runs with no HMMER hits.
    FINAL_RESULTS(
        GENERATE_MANIFEST.out.manifest,
        GENERATE_MANIFEST.out.marker_map,
        GENERATE_MANIFEST.out.genome_markers,
        HMM_SEARCH.out.hits.map { ref_id, table -> table }.toList(),
        PRODIGAL.out.proteins.flatMap { batch_id, ids, files -> files }.toList(),
        PRODIGAL.out.cds.flatMap { batch_id, ids, files -> files }.toList(),
        file("${params.general_input}/remove_specific_ORF_sequences.txt", checkIfExists: true),
        file("${params.scripts_generation}/final_results.R", checkIfExists: true)
    )

    emit:
    database = FINAL_RESULTS.out.database
    manifest = GENERATE_MANIFEST.out.manifest
    marker_map = GENERATE_MANIFEST.out.marker_map
    genome_markers = GENERATE_MANIFEST.out.genome_markers
}
