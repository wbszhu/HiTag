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
        for rec in read_rec(samf):
            R1 = rec["R1"]
            R2 = rec["R2"]
            if is_unique(R1, mapq_thresh) and is_unique(R2, mapq_thresh):
                bedpe = bed2pe(R1[0], R2[0])
                f3.write(bedpe + '\n')
            else:
                seqid1, seq1, quality1 = sam_to_fq(R1[0])
                seqid2, seq2, quality2 = sam_to_fq(R2[0])
                write_fastq(seqid1, seq1, quality1, f1)
                write_fastq(seqid2, seq2, quality2, f2)            


def new_rec():
    return {
        'R1': [],
        'R2': [],
    }


def rec_append_read(rec, read):
    if read.is_read1:
        rec['R1'].append(read)
    else:
        rec['R2'].append(read)


def read_rec(samf):
    prev_qname = "_init_"
    rec = new_rec()
    for read in samf:
        if (read.qname == prev_qname) or (prev_qname == "_init_"):
            rec_append_read(rec, read)
        else:
            yield rec
            rec = new_rec()
            rec_append_read(rec, read)
        prev_qname = read.qname
    yield rec


def is_unique(rec, mapq_thresh):
    if len(rec) > 1:
        return False
    read = rec[0]
    if (read.mapq >= mapq_thresh) & (read.flag & 4 ==0):
        return True
    else:
        return False


def sam_to_fq(read):
    seqid = read.qname
    seq = read.query_sequence
    quality = read.qual
    return seqid, seq, quality


def sam_to_bed(read):
    name = read.qname
    chr_, start = read.reference_name, read.reference_start
    strand = '-' if read.is_reverse else '+'
    return name, chr_, start, strand


def bed2pe(read1, read2):
    read1, read2 = sorted([read1, read2], key=lambda r: r.reference_start)
    name1, chr_1, start_1, strand_1 = sam_to_bed(read1)
    name2, chr_2, start_2, strand_2 = sam_to_bed(read2)
    assert name1 == name2
    bedpe = "\t".join([
        name1, chr_1, str(start_1),
        chr_2, str(start_2), strand_1, strand_2, "UU"])
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
