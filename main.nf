#!/usr/bin/env/ nextflow

nextflow.enable.dsl = 2

process trim_adapter {
    input:
        tuple val(sample_id), file(reads)
    output:
        tuple val(sample_id), file("*_val_{1,2}.fq.gz")
    """
    trim_galore -q 20 --phred33 --paired --illumina --trim-n --gzip ${reads[0]} ${reads[1]}
    """
}

process trim_linker {
    input:
        tuple val(sample_id), file(reads)
    output:
        tuple val(sample_id), file("${sample_id}_{1,2}.valid.fastq")
    
    """
    trimLinker -e 2 -t 12 -m 1 -k 2 -l 16 -o . -n ${sample_id} -A ${params.linker_A} -B ${params.linker_B} ${reads[0]} ${reads[1]}
    """
}

process bwa_mapping {
    input:
        tuple val(sample_id), file(reads)
    output:
        tuple val(sample_id), file("${sample_id}.bam")
    """
    bwa mem -SP5M -t 20 ${params.bwa_index_prefix} ${reads[0]} ${reads[1]} | samtools view -b > ${sample_id}.bam
    """
}

process bam_to_pairs {
    input:
        tuple val(sample_id), file(bam)
    output:
        tuple val(sample_id), file("${sample_id}.pairs.gz")
    """
    pairtools parse -c ${params.genome_size} -o ${sample_id}.pairs.gz --min-mapq ${params.min_mapq} --drop-sam ${bam}
    """
}

process sort_pairs {
    input:
        tuple val(sample_id), file(pairs)
    output:
        tuple val(sample_id), file("${sample_id}.sorted.pairs.gz")
    """
    pairtools sort -o ${sample_id}.sorted.pairs.gz ${pairs}
    """
}

process dedup_pairs {
    input:
        tuple val(sample_id), file(pairs)
    output:
        tuple val(sample_id), file("${sample_id}.dedup.pairs.gz"), file("${sample_id}.dedup.stats.txt")
    """
    pairtools dedup --mark-dups --output-stats ${sample_id}.dedup.stats.txt -o ${sample_id}.dedup.pairs.gz ${pairs}
    """
}

process select_UU_pairs {
    input:
        tuple val(sample_id), file(pairs)
    output:
        tuple val(sample_id), file("${sample_id}.UU.pairs.gz")

    """
    pairtools select '(pair_type=="UU")' -o ${sample_id}.UU.pairs.gz ${pairs}
    """
}

process stat_UU_pairs {
    input:
        tuple val(sample_id), file(pairs)
    output:
        tuple val(sample_id), file("${sample_id}.UU.stats.txt")
    """
    pairtools stats -o ${sample_id}.UU.stats.txt ${pairs}
    """
}

process remove_noise {
    input:
        tuple val(sample_id), file(pairs)
    output:
        tuple val(sample_id), file("${sample_id}.valid.pairs"), file("${sample_id}.noise_reduce.log")

    publishDir "results/remove_noise", mode: 'symlink'
    """
    python ${params.script_dir}/noise_reduce.py ${pairs} ${sample_id} -r ${params.digest} --log-file ${sample_id}.noise_reduce.log
    """
}

process pairs_to_validPairs {
    input:
        tuple val(sample_id), file(pairs)
    output:
        tuple val(sample_id), file("${sample_id}.validPairs")
    
    """
    cat ${pairs} | grep -v "^#" | awk 'BEGIN{FS="\\t";OFS="\\t"}{print \$1,\$2,\$3,\$6,\$4,\$5,\$7}' > ${sample_id}.validPairs
    """
}

process cooler_csort {
    input:
        tuple val(sample_id), file(pairs)
    output:
        tuple val(sample_id), file("${sample_id}.pairs.gz"), file("${sample_id}.pairs.gz.px2")
    
    publishDir "results/indexed_pairs", mode: 'symlink'
    """
    cooler csort -c1 2 -p1 3 -c2 4 -p2 5 -o ${sample_id}.pairs.gz ${pairs} ${params.genome_size}
    """
}

process cooler_cload {
    input:
        tuple val(sample_id), file(pairs), file(pairs_ix)
    output:
        tuple val(sample_id), file("${sample_id}.cool")

    """
    cooler cload pairix ${params.genome_size}:${params.resolutions[0]} ${pairs} ${sample_id}.cool
    """
}

process cooler_zoomify {
    input:
        tuple val(sample_id), file(cool)
    output:
        tuple val(sample_id), file("${sample_id}.balance.mcool")

    publishDir "results/mcool", mode: 'symlink'
    """
    cooler zoomify -r ${params.resolutions.join(",")} -o ${sample_id}.balance.mcool --balance ${cool}
    """
}

process fithichip {
    input:
        tuple val(sample_id), file(pairs)
    output:
        tuple val(sample_id), path("${sample_id}.fithichip.out")
    
    publishDir "results/fithichip", mode: 'symlink'

    script:
        config = "${sample_id}_fithichip.config"
        result = "${sample_id}.fithichip.out"
        """
        python ${params.script_dir}/format_fithichip_config.py ${pairs} ${params.chipdata} ${result} ${params.resolution_fithichip} ${params.genome_size} ${params.use_P2P_background} > ${config}
        FitHiChIP_HiCPro.sh -C ${config}
        """
}

process collect_stats {
    input:
        tuple val(sample_id), file(stats), val(type)
    output:
        tuple val(sample_id), file("${sample_id}.stat")

    publishDir "results/stats", mode: 'symlink'

    script:
        out = "${sample_id}.stat"
        cmd = stats.indexed().collect() { index, f -> 
            """
            echo "--------------------------" >> ${out}
            echo ${type[index]} >> ${out}
            echo "--------------------------" >> ${out}
            cat ${f} >> ${out}
            echo >> ${out}
            """
        }
        cmd.join("\n")
}

process count_fastq {
    input:
        tuple val(sample_id), file(fastq), val(type)
    output:
        tuple val(sample_id), file("${sample_id}.${type}.count"), val(type)
    
    script:
        out = "${sample_id}.${type}.count"
        if (fastq ==~ /.*\.gz$/) { cat = "zcat" } else { cat = "cat" }
        """
        l=\$(${cat} ${fastq} | wc -l | awk '{print \$1}')
        echo \$((\$l / 4)) > ${out}
        """
}

process count_pairs {
    input:
        tuple val(sample_id), file(pairs), val(type)
    output:
        tuple val(sample_id), file("${sample_id}.${type}.count"), val(type)

    script:
        out = "${sample_id}.${type}.count"
        if (pairs ==~ /.*\.gz$/) { cat = "zcat" } else { cat = "cat" }
        """
        ${cat} ${pairs} | grep -v "^#" | awk 'BEGIN{c=0;t=0;OFS="\t"}{if(\$2==\$4){c+=1}else{t+=1}}END{print NR,c,t}' > ${out}
        """
}

process collect_counts {
    input:
        tuple val(sample_id), file(count), val(type)
    output:
        tuple val(sample_id), file("${sample_id}.count")

    publishDir "results/count", mode: 'symlink'

    script:
        out = "${sample_id}.count"
        cmd = count.indexed().collect() { index, f -> 
            """
            echo -n "${type[index]}\t" >> ${out}
            cat ${f} >> ${out}
            """
        }
        cmd.join("\n")
}


workflow {
    // main flow
    reads = Channel.fromFilePairs(params.input_reads) 
    tr_adapter = reads | trim_adapter
    pets = tr_adapter | trim_linker 
    bam = pets | bwa_mapping
    pairs = bam | bam_to_pairs
    deduped_outputs = pairs | sort_pairs | dedup_pairs
    deduped_pairs = deduped_outputs.map { t -> t[0..1] } 
    uu_pairs = deduped_pairs | select_UU_pairs
    rn_outputs = uu_pairs | remove_noise
    rn_pairs = rn_outputs.map { t -> [t[0], t[1]] }
    rn_pairs | cooler_csort | cooler_cload | cooler_zoomify
    rn_pairs | pairs_to_validPairs | fithichip

    // collect stats
    stats_dedup = deduped_outputs.map{ t -> [t[0], t[2], "dedup_stat"] }
    stats_uu = (uu_pairs | stat_UU_pairs).map{ t -> [t[0], t[1], "UU_stats"] }
    stats_rn = rn_outputs.map{ t -> [t[0], t[2], "remove_noise"] }
    stats_dedup.mix(stats_uu, stats_rn).groupTuple() | collect_stats

    // collect counts
    cnt_reads = reads.map{ t -> [t[0], t[1][0], "raw_reads"] }
    cnt_adapter = tr_adapter.map{ t -> [t[0], t[1][0], "trim_adapter"] }
    cnt_pets = pets.map{ t -> [t[0], t[1][0], "trim_linker"] }
    cnt_fastq = cnt_reads.mix(cnt_adapter, cnt_pets) | count_fastq
    cnt_b2p = pairs.map{ t -> [t[0], t[1], "bam2pairs"] }
    cnt_dedup = deduped_pairs.map{ t -> [t[0], t[1], "dedup"]}
    cnt_uu = uu_pairs.map{ t -> [t[0], t[1], "UU"]}
    cnt_rn = rn_pairs.map{ t -> [t[0], t[1], "remove_noise"]}
    cnt_pairs = cnt_b2p.mix(cnt_dedup, cnt_uu, cnt_rn) | count_pairs
    cnt_fastq.mix(cnt_pairs).groupTuple() | collect_counts
}