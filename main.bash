#!/usr/bin/env bash

# Generate VirMarkDB database 
    # Generate a manifest list of virus
    module load r/4.5.2
    Rscript scripts/database_generation/01_generate_manifest.R
    # Search ORF in viral genomes
    bash scripts/database_generation/02_prodigal.bash
    
    # Generate HMM profile from references protein
        # Make an alignment of references sequences
        bash scripts/database_generation/03_align_markers.bash
        # Build hmm profile
        bash scripts/database_generation/04_hmmbuild.bash
        # Search protein with HMM
        bash scripts/database_generation/05_hmmsearch.bash
    
    # Generate output files
    Rscript scripts/database_generation/06_final_results.R



# VirMarkDB assessment   
    # Generate a pool of potentential protein from ncbi
        # Load bamford protein from the local nr database of bamfordvirae
        sbatch scripts/assessment/07_scrap_viral_protein_ncbi.bash
        #  Filter out mis annotated protein with an auto alignment
        sbatch scripts/assessment/08_align_pool_vs_pool.bash
        module load r/4.5.2
        Rscript scripts/assessment/09_results_align_pool_vs_pool.R

    # Test VirMarkDB with algniment of a pool of viral protein
        # Blast with mmseqs each marker from VirMarkDB to each pool prot
        sbatch scripts/assessment/10_align_pool_VirMarkDB.bash
        # Analyze the results
        module load r/4.5.2
        Rscript scripts/assessment/11_analyse_results_align_pool_VirMarkDB.R
        # Clusterise missing protein 
        bash scripts/assessment/12_clusterise_missing_prot.bash
    
    # General statistics
    Rscript scripts/assessment/13_summary_satistics.R
