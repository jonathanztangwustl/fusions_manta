version 1.0

import "mutect2.wdl" as m2

# WDL tasks and workflow for investigating CH fusions with Terra
# Also runs mutect2 and deletion detection to search for mutations in subset of genes

# Samtools task for fusions.py
task samtools_fusions {
    input {
        File fusions_bed
        File full_cram
        File full_cram_crai
        File ref_fa
        File ref_fa_fai
    }
    Int cores = 1
    Float cram_size = size([full_cram], "GB")
    Int runtime_size = 4 + round(cram_size)
    runtime {
        memory: "4GB"
        cpu: cores
        preemptible: 1
        docker: "quay.io/biocontainers/samtools:1.2--0"
        disks: "local-disk ~{runtime_size} SSD"
        bootDiskSizeGb: runtime_size
    }
    command <<<
        set -o pipefail
        set -o errexit
        ln -s ~{full_cram} full.cram
        ln -s ~{full_cram_crai} full.cram.crai
        samtools view -b -L ~{fusions_bed} -M -T ~{ref_fa} -t ~{ref_fa_fai} -F 1036 -f 1 -q 30 -o fusions.bam full.cram
    >>>
    output {
        File fusions_bam = "fusions.bam"
    }
}

# Samtools task for mutect
task samtools_mutect {
    input {
        File subset_bed
        File full_cram
        File full_cram_crai
        File ref_fa
        File ref_fa_fai
    }
    Int cores = 1
    Float bam_size = size([full_cram], "GB")
    Int runtime_size = 4 + round(bam_size)
    runtime {
        memory: "4GB"
        cpu: cores
        preemptible: 1
        docker: "quay.io/biocontainers/samtools:1.2--0"
        disks: "local-disk ~{runtime_size} SSD"
        bootDiskSizeGb: runtime_size
    }
    command <<<
        set -o pipefail
        set -o errexit
        ln -s ~{full_cram} full.cram
        ln -s ~{full_cram_crai} full.cram.crai
        samtools view -b -L ~{subset_bed} -M -T ~{ref_fa} -t ~{ref_fa_fai} -F 1036 -f 1 -o subset.bam full.cram
        samtools index -b subset.bam
    >>>
    output {
        File subset_bam = "subset.bam"
        File subset_bai = "subset.bam.bai"
    }
}

# Fusions.py task
task fusions {
    input {
        File fusions_py
        File fusions_bam
        File ROIs
        File gene_ref_bed
    }
    Int cores = 4
    Float bam_size = size([fusions_bam], "GB")
    Int runtime_size = 4 + round(bam_size)
    runtime {
        memory: "4GB"
        cpu: cores
        preemptible: 1
        docker: "chrisamiller/genomic-analysis:0.2"
        disks: "local-disk ~{runtime_size} SSD"
        bootDiskSizeGb: runtime_size
    }
    command <<<
        python \
        ~{fusions_py} \
        ~{fusions_bam} \
        ~{ROIs} \
        ~{gene_ref_bed} \
        fusions_out.txt
    >>>
    output {
        File fusions_out = "fusions_out.txt"
    }
}

# CNV - Indexcov
task indexcov {
    input {
        File full_cram_crai
        File ref_fa_fai
    }
    Int cores = 1
    Float crai_size = size([full_cram_crai], "GB")
    Int runtime_size = 4 + round(crai_size)
    runtime {
        memory: "4GB"
        cpu: cores
        preemptible: 1
        docker: "quay.io/biocontainers/goleft:0.2.4--0"
        disks: "local-disk ~{runtime_size} SSD"
        bootDiskSizeGb: runtime_size
    }
    command <<<
        goleft indexcov --extranormalize -d indexcov_out --fai ~{ref_fa_fai} ~{full_cram_crai}
        tar -cvf indexcov.tar indexcov_out/
        gzip indexcov.tar
    >>>
    output {
        File indexcov_out = "indexcov.tar.gz"
    }
}

# CNV - mosdepth
task mosdepth {
    input {
        File ref_fa
        File ref_fa_fai
        File full_cram
        File full_cram_crai
        File subset_mosdepth
    }
    Int cores = 1
    Float bam_size = size([full_cram], "GB")
    Int runtime_size = 4 + round(bam_size)
    runtime {
        memory: "8GB"
        cpu: cores
        preemptible: 1
        docker: "quay.io/biocontainers/mosdepth:0.2.5--hb763d49_0"
        disks: "local-disk ~{runtime_size} SSD"
        bootDiskSizeGb: runtime_size
    }
    command <<<
        ln -s ~{full_cram} full.cram
        ln -s ~{full_cram_crai} full.cram.crai
        ln -s ~{ref_fa} ref.fa
        ln -s ~{ref_fa_fai} ref.fa.fai
        mkdir mosdepth
        mosdepth -b ~{subset_mosdepth} -n -f ref.fa mosdepth/cnv full.cram
        tar -cvf mosdepth.tar mosdepth/
        gzip mosdepth.tar
    >>>
    output {
        File mosdepth_out = "mosdepth.tar.gz"
    }
}

# Workflow to call samtools, fusions, mutations, and deletions
workflow fusions_mutations {
    input {
        File wf_fusions_bed
        File wf_full_cram
        File wf_full_cram_crai
        File wf_fusions_py
        File wf_ROIs
        File wf_subset_bed
        File wf_subset_mosdepth
        File wf_gene_ref_bed
        File wf_ref_fa
        File wf_ref_fa_fai
        File wf_ref_dict
        File Mutect2_intervals
        Int Mutect2_scatter_count
        String Mutect2_gatk_docker
        File Mutect2_gatk_override
        File wf_subset_mosdepth
    }
    call indexcov {
        input:
        full_cram_crai=wf_full_cram_crai,
        ref_fa_fai=wf_ref_fa_fai
    }
    call mosdepth {
        input:
        ref_fa=wf_ref_fa,
        ref_fa_fai=wf_ref_fa_fai,
        full_cram=wf_full_cram,
        full_cram_crai=wf_full_cram_crai,
        subset_mosdepth=wf_subset_mosdepth
    }
    call samtools_fusions {
        input:
        fusions_bed=wf_fusions_bed,
        full_cram=wf_full_cram,
        full_cram_crai=wf_full_cram_crai,
        ref_fa=wf_ref_fa,
        ref_fa_fai=wf_ref_fa_fai
    }
    call fusions {
        input:
        fusions_py=wf_fusions_py,
        fusions_bam=samtools_fusions.fusions_bam,
        ROIs=wf_ROIs,
        gene_ref_bed=wf_gene_ref_bed
    }
    call samtools_mutect {
        input:
        subset_bed=wf_subset_bed,
        full_cram=wf_full_cram,
        full_cram_crai=wf_full_cram_crai,
        ref_fa=wf_ref_fa,
        ref_fa_fai=wf_ref_fa_fai
    }
    call m2.Mutect2{
        input:
        intervals=Mutect2_intervals,
        ref_fasta=wf_ref_fa,
        ref_fai=wf_ref_fa_fai,
        ref_dict=wf_ref_dict,
        tumor_reads=samtools_mutect.subset_bam,
        tumor_reads_index=samtools_mutect.subset_bai,
        scatter_count=Mutect2_scatter_count,
        gatk_docker=Mutect2_gatk_docker,
        gatk_override=Mutect2_gatk_override
    }
    output {
        File fusions_bam = samtools_fusions.fusions_bam
        File fusions_out = fusions.fusions_out
        File subset_bam = samtools_mutect.subset_bam
        File subset_bai = samtools_mutect.subset_bai
        File filtered_vcf = Mutect2.filtered_vcf
        File filtered_vcf_idx = Mutect2.filtered_vcf_idx
        File filtering_stats = Mutect2.filtering_stats
        File mutect_stats = Mutect2.mutect_stats
        File? contamination_table = Mutect2.contamination_table
        File? funcotated_file = Mutect2.funcotated_file
        File? funcotated_file_index = Mutect2.funcotated_file_index
        File? bamout = Mutect2.bamout
        File? bamout_index = Mutect2.bamout_index
        File? maf_segments = Mutect2.maf_segments
        File? read_orientation_model_params = Mutect2.read_orientation_model_params
        File indexcov_out = indexcov.indexcov_out
        File mosdepth_out = mosdepth.mosdepth_out
    }
}
