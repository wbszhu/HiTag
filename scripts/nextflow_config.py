import sys

Sample=sys.argv[1]
linker_A=sys.argv[2]
linker_B=sys.argv[3]
min_mapq=sys.argv[4]
genome_size=sys.argv[5]
bwa_index_prefix=sys.argv[6]
digest=sys.argv[7]

content_to_addvalue = {
    "Sample":Sample,
    "linker_A":linker_A,
    "linker_B":linker_B,
    "min_mapq":min_mapq,
    "genome_size":genome_size,
    "bwa_index_prefix":bwa_index_prefix,
    "digest":digest,
}


temp = """
params {{
String dir = new File(".").getAbsolutePath()

    input_reads = "${{dir}}/../00.data/{Sample}_{{1,2}}.fq.gz"

    linker_A = "{linker_A}"
    linker_B = "{linker_B}"
    min_mapq = "{min_mapq}"

    genome_size = "{genome_size}"
    bwa_index_prefix = "{bwa_index_prefix}"

    resolutions = [5000, 10000, 25000, 40000, 100000]

    chipdata = 
    use_P2P_background = 
    digest = "{digest}"

    script_dir = "${{dir}}/scripts"
}}
"""

new_file_content = temp.format(**content_to_addvalue)

print(new_file_content)
