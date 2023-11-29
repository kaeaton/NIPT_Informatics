#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
from pathlib import Path

# constants

# Each number in the list provides the count past which all data will be trimmed off (inclusive).
# (in other words, any read longer than the number of nucleotides listed will be trimmed out of the results.)
# This should enhance the fetal fraction.

# This list allows you to add as many sizes (in nucleotides) as you like.
# SIZE_SELECTION_LIST = [130, 135, 136, 140]
SIZE_SELECTION_LIST = [135, 140]

# The number of nucleotides in the UMI. This shouldn't change often, but in case it needs to...
UMI_LENGTH = 8

# The number of threads for your Intel CPU. An average core has two channels, so two threads from this
# process can be run per core. You want to save one core for running the system while you run this process.
# Example: if you have 4 cores and 2 channels per core, you could run 8 threads,
# but it will process faster if you use <= 6 threads and save one (or more) cores for running the system.
# If your CPU is composed of higher and lower power cores (ARM architecture and some newer Intel CPUs)
# it gets more complicated. Apple runs threads in groups of four, routing first to the high-powered cores,
# so your number of threads should be a multiple of four.
NUMBER_OF_THREADS = 8

# The reference genome used for alignment and deduplication. This has only been tested with hg19 and hg38.
REFERENCE_GENOME = Path("./hg38/hg38.fa")

# round up the read 1 files
read_1_data_files = Path("./samples").glob('**/*_R1*.fastq.gz')

for i in read_1_data_files:

    # define the data file pairs and sample ID
    read_1 = i
    read_2 = Path(str(i).replace('R1', 'R2'))
    sample_id = str(os.path.basename(i)).split('_S')[0]

    # create results folder for this sample
    results_folder = Path("./results")
    if not os.path.isdir(results_folder):
        os.makedirs(results_folder)

    print(f"Sample {sample_id}\t")
    print(str(os.path.basename(read_1)), str(os.path.basename(read_2)))

    for j in SIZE_SELECTION_LIST:

        print(f"Sample {sample_id}, selecting for {j}bp length reads or shorter.")

        ## Fastp trimming ##
        read_1_trimmed = Path(results_folder, f'{sample_id}_R1_size_limit_{j}bp.fastq.gz')
        read_2_trimmed = Path(results_folder, f'{sample_id}_R2_size_limit_{j}bp.fastq.gz')

        print('Trimming files')

        trim = f"fastp -i {read_1} -o {read_1_trimmed} -I {read_2} -O {read_2_trimmed} \
                               --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
                               --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
                               -U --umi_prefix=UMI --umi_loc=per_read --umi_len={UMI_LENGTH} -c  \
                               --length_limit {j}"


        os.system(trim)

        print(f'\n** fastp trim of sample {sample_id}_{j}bp completed. **\n')

        ## BWA mem alignment ##
        reads_combined_sam = Path(results_folder, f'{sample_id}_{j}bp_aligned.sam')

        bwa_mem_header_string = f'"@RG\\tID:1\\tDS:SEQCAP_EZ\\tPL:ILLUMINA\\tLB:{sample_id}\\tSM:{sample_id}"'
        align = f"bwa mem -M -t {NUMBER_OF_THREADS} -R {bwa_mem_header_string} {REFERENCE_GENOME} " \
                          f"{read_1_trimmed} {read_2_trimmed} > {reads_combined_sam}"

        os.system(align)
        print(f'\n** bwa mem alignment of sample {sample_id}_{j}bp completed. **\n')

        ## GATK tools ##

        # GATK - convert from sam to bam
        reads_combined_bam = Path(results_folder, f'{sample_id}_{j}bp_aligned.bam')

        gatk_sam_to_bam_file = f"gatk SamFormatConverter -I {reads_combined_sam} -O {reads_combined_bam}"

        os.system(gatk_sam_to_bam_file)
        print(f'\n** gatk conversion of sample {sample_id}_{j}bp from sam file to bam file completed. **\n')

        # GATK - Fix Mate Information
        reads_combined_and_fixed_bam = Path(results_folder, f'{sample_id}_{j}bp_fixed.bam')

        gatk_fix = f"gatk FixMateInformation -I {reads_combined_bam} -O {reads_combined_and_fixed_bam}"

        os.system(gatk_fix)
        print(f'\n** gatk fixed read mates of sample {sample_id}_{j}bp in bam file. **\n')

        # GATK - Sort Results
        reads_sorted_bam = Path(results_folder, f'{sample_id}_{j}bp_sorted.bam')

        gatk_sort = f"gatk SortSam -I {reads_combined_and_fixed_bam} -O {reads_sorted_bam} -SO coordinate"

        os.system(gatk_sort)
        print(f'\n** gatk sorted results of sample {sample_id}_{j}bp. **\n')

        ## GenCore ##
        deduplicated_results = Path(results_folder, f'{sample_id}_{j}bp_deduplicated.bam')

        gencore_umi = f"gencore -i {reads_sorted_bam} -o {deduplicated_results} -r {REFERENCE_GENOME} -s 1 -u UMI"
        os.system(gencore_umi)
        print(f'\n** gencore deduplication for sample {sample_id}_{j}bp is complete. **\n')

        # File cleanup
        os.remove(read_1_trimmed)
        os.remove(read_2_trimmed)
        os.remove(reads_combined_sam)
        os.remove(reads_combined_bam)
        os.remove(reads_combined_and_fixed_bam)
        os.remove(reads_sorted_bam)

print("All samples complete.")