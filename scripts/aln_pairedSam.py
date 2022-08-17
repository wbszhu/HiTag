import click
import pysam
from pysam import Samfile

@click.command("")
@click.argument("sam1")
@click.argument("sam2")
@click.argument("bedpe")
@click.option('--mapq_thresh', default=20, show_default=True)
def aln_UUpe(sam1, sam2, mapq_thresh, bedpe):
    with open(bedpe, 'w') as out:
        samf1 = pysam.AlignmentFile(sam1, mode='rb')
        samf2 = pysam.AlignmentFile(sam2, mode='rb')
        while True:
            try:
                rec1 = read_rec(samf1)[0]
                rec2 = read_rec(samf2)[0]
            except StopIteration:
                break
            type1, name1, chr1, start1, strand1 = sam_to_bed(rec1, mapq_thresh)
            type2, name2, chr2, start2, strand2 = sam_to_bed(rec2, mapq_thresh)
            if (type1 == "unique" and type2 == "unique"):
                if (rec1.qname == rec2.qname):
                    if (start1 < start2):
                        bedpe = "\t".join([name1, chr1, str(start1), chr2, str(start2), strand1, strand2, "UU"])
                        out.write(bedpe + '\n')
                    elif (start1 > start2):
                        bedpe = "\t".join([name1, chr2, str(start2), chr1, str(start1), strand2, strand1, "UU"])
                        out.write(bedpe + '\n')


def read_rec(samf):
    rec = []
    for _, line in zip(range(1), samf): # one line once
        rec.append(line)
    if len(rec) < 1 :
        raise StopIteration
    return rec


def sam_to_bed(read, mapq_thresh):
    if (read.mapq >= mapq_thresh) & (read.flag & 4 ==0):
        name = read.qname
        chr_, start= read.reference_name, read.reference_start
        strand = '-' if read.is_reverse else '+'
        type = "unique"
        return type, name, chr_, start, strand
    else:
        name = read.qname
        chr_, start= read.reference_name, read.reference_start
        strand = '-' if read.is_reverse else '+'
        type = "other"
        return type, name, chr_, start, strand


#def bed2pe(read1, read2):
#    name1, chr_1, start_1, strand_1 = sam_to_bed(read1)
#    name2, chr_2, start_2, strand_2 = sam_to_bed(read2)
#    bedpe = str()
#    if (name1 == name2) and (start_1 < start_2):
#        bedpe = "\t".join([name1, chr_1, str(start_1), chr_2, str(start_2), strand_1, strand_2, "UU"])
#    elif (name1 == name2) and (start_1 > start_2):
#        bedpe = "\t".join([name1, chr_2, str(start_2), chr_1, str(start_1), strand_2, strand_1, "UU"])
#    return bedpe


if __name__ == "__main__":
    aln_UUpe()
