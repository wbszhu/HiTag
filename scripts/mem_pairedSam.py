import click
import pysam
from pysam import Samfile

@click.command("extractUU to bed")
@click.argument("samfile")
@click.argument("fq1")
@click.argument("fq2")
@click.argument("bedpe")
@click.option('--mapq_thresh', default=20, show_default=True)
def parse_sam(samfile, mapq_thresh, fq1, fq2, bedpe):
    with open(fq1, 'w') as f1, open(fq2, 'w') as f2, open(bedpe, 'w') as f3:
        samf = pysam.AlignmentFile(samfile, mode='rb')
        while True:
            try:
                rec = read_rec(samf)
            except StopIteration:
                break
            R1 = rec["R1"]
            R2 = rec["R2"]
            if (unique(R1, mapq_thresh) == "unique") and (unique(R2, mapq_thresh) == "unique"):
                bedpe = bed2pe(R1, R2)
                if bedpe != "other":
                    f3.write(bedpe + '\n')
            else:
                seqid1, seq1, quality1 = sam_to_fq(R1)
                seqid2, seq2, quality2 = sam_to_fq(R2)
                write_fastq(seqid1, seq1, quality1, f1)
                write_fastq(seqid2, seq2, quality2, f2)            


def read_rec(samf):
    lines = []
    for _, line in zip(range(2), samf): # two line once
        lines.append(line)
    if len(lines) < 2 :
        raise StopIteration
    rec = {'R1': lines[0],'R2': lines[1]}
    return rec


def unique(read, mapq_thresh):
    if (read.mapq >= mapq_thresh) & (read.flag & 4 ==0):
        return 'unique'
    else:
        return 'other'


def sam_to_fq(read):
    seqid = read.qname
    seq = read.seq
    quality = read.qual
    return seqid, seq, quality


def sam_to_bed(read):
    name = read.qname
    chr_, start= read.reference_name, read.reference_start
    strand = '-' if read.is_reverse else '+'
    return name, chr_, start, strand

def bed2pe(read1, read2):
    name1, chr_1, start_1, strand_1 = sam_to_bed(read1)
    name2, chr_2, start_2, strand_2 = sam_to_bed(read2)
    bedpe = str()
    if (name1 == name2) and (start_1 < start_2):
        bedpe = "\t".join([name1, chr_1, str(start_1), chr_2, str(start_2), strand_1, strand_2, "UU"])
    elif (name1 == name2) and (start_1 > start_2):
        bedpe = "\t".join([name1, chr_2, str(start_2), chr_1, str(start_1), strand_2, strand_1, "UU"])
    else:
        bedpe = "other"
    return bedpe


def write_fastq(seqid, seq, quality, fh):
    """
    Write Fastq ord to file.
    """
    fh.write("@"+seqid+"\n")
    fh.write(seq+"\n")
    fh.write("+\n")
    fh.write(quality+"\n")


if __name__ == "__main__":
    parse_sam()