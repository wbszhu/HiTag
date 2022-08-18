#!/usr/bin/env/ nextflow

nextflow.enable.dsl = 2

process trim_adapter {
    input:
        tuple val(sample_id), file(reads)
    output:
        tuple val(sample_id), file("*_val_{1,2}.fq.gz")
    """
    trim_galore -q 20 --phred33 --paired --Nextera --trim-n --gzip ${reads[0]} ${reads[1]}
    """
}


process trim_linker {
    input:
        tuple val(sample_id), file(reads)
    output:
        tuple val(sample_id), file("${sample_id}_{1,2}.valid.fastq.gz")
    
    """
    python ${params.script_dir}/trim_linker_V2.py ${reads[0]} ${reads[1]} ${sample_id}_1.valid.fastq.gz ${sample_id}_2.valid.fastq.gz ${params.linker_A} ${params.linker_B} 
    """
}


process bwa_mem {
    input:
        tuple val(sample_id), file(reads)
    output:
        tuple val(sample_id), file("${sample_id}_mem.sam")
    """
    bwa mem -SP5M -t 20 ${params.bwa_index_prefix} ${reads[0]} ${reads[1]} > ${sample_id}_mem.sam
    """
}

process mem_sam_to_pairs {
    input:
        tuple val(sample_id), file(sam)
    output:
        tuple val(sample_id), file("${sample_id}_{1,2}.fq"), file("${sample_id}.mem.pairs")
    """
    python ${params.script_dir}/mem_pairedSam.py ${sam} ${sample_id}_1.fq ${sample_id}_2.fq ${sample_id}.mem.pairs
    """
}

process bwa_aln {
    input:
        tuple val(sample_id), val(reads_type), file(reads)
    output:
        tuple val(sample_id), val(reads_type), file("${sample_id}.${reads_type}.sai"), file(reads)
    """
    bwa aln -t 20 ${params.bwa_index_prefix} ${reads} > ${sample_id}.${reads_type}.sai
    """
}


process bwa_samse {
    input:
        tuple val(sample_id), val(reads_type), file(sai), file(reads)
    output:
        tuple val(sample_id), val(reads_type), file("${sample_id}.${reads_type}.sam")
    """
    bwa samse  ${params.bwa_index_prefix} ${sai} ${reads} > ${sample_id}.${reads_type}.sam
    """
}


process aln_sam_to_pairs{
    input:
        tuple val(sample_id), file(sam_1), file(sam_2)
    output:
        tuple val(sample_id), file("${sample_id}.aln.pairs")
    """
    python ${params.script_dir}/aln_pairedSam.py  ${sam_1} ${sam_2} ${sample_id}.aln.pairs
    """
}


process cat_allpairs{
    input:
        tuple val(sample_id), file(pairs_1), file(pairs_2)
    output:
        tuple val(sample_id), file("${sample_id}.pairs")
    """
    echo "## pairs format" > ${sample_id}.pairs
    cat ${pairs_1} >> ${sample_id}.pairs
    cat ${pairs_2} >> ${sample_id}.pairs
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
        bash FitHiChIP_HiCPro.sh -C ${config}
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
    mem_results = pets | bwa_mem | mem_sam_to_pairs
    mem_results.multiMap { r ->
        pairs: [r[0], r[2]]
        unmaped_fq_1: [r[0], "R1", r[1][0]]
        unmaped_fq_2: [r[0], "R2", r[1][1]]
 
    }.set{ mem_res }
    aln_sams = mem_res.unmaped_fq_1.concat(
        mem_res.unmaped_fq_2) | bwa_aln | bwa_samse
    aln_sams.branch {
        r1: it[1] == "R1"
        r2: it[1] == "R2"
    }.set{ aln_sam_sep }
    sams = aln_sam_sep.r1.map{
            t -> [t[0], t[2]]
        }.join(aln_sam_sep.r2.map{
            t -> [t[0], t[2]]
        })
    aln_pairs = sams | aln_sam_to_pairs
    pairs = aln_pairs.join(mem_res.pairs) | cat_allpairs
    deduped_outputs = pairs | sort_pairs | dedup_pairs
    deduped_pairs = deduped_outputs.map { t -> t[0..1] } 
    rn_outputs = deduped_pairs | remove_noise
    rn_pairs = rn_outputs.map { t -> [t[0], t[1]] }
    rn_pairs | cooler_csort | cooler_cload | cooler_zoomify
    rn_pairs | pairs_to_validPairs | fithichip

    // collect stats
    stats_dedup = deduped_outputs.map{ t -> [t[0], t[2], "dedup_stat"] }
    stats_uu = (deduped_pairs | stat_UU_pairs).map{ t -> [t[0], t[1], "UU_stats"] }
    stats_rn = rn_outputs.map{ t -> [t[0], t[2], "remove_noise"] }
    stats_dedup.mix(stats_uu, stats_rn).groupTuple() | collect_stats

    // collect counts
    cnt_reads = reads.map{ t -> [t[0], t[1][0], "raw_reads"] }
    cnt_adapter = tr_adapter.map{ t -> [t[0], t[1][0], "trim_adapter"] }
    cnt_pets = pets.map{ t -> [t[0], t[1][0], "trim_linker"] }
    cnt_fastq = cnt_reads.mix(cnt_adapter, cnt_pets) | count_fastq
    cnt_b2p = pairs.map{ t -> [t[0], t[1], "bam2pairs"] }
    cnt_dedup = deduped_pairs.map{ t -> [t[0], t[1], "dedup"]}
    cnt_rn = rn_pairs.map{ t -> [t[0], t[1], "remove_noise"]}
    cnt_pairs = cnt_b2p.mix(cnt_dedup, cnt_rn) | count_pairs
    cnt_fastq.mix(cnt_pairs).groupTuple() | collect_counts
}
