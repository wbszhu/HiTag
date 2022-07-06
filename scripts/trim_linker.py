from cutadapt._align import Aligner
import click
import gzip
import io


@click.command("")
@click.argument("fq_1")
@click.argument("fq_2")
@click.argument("out_f1")
@click.argument("out_f2")
@click.argument("linkera")
@click.argument("linkerb")
def trim_linker(fq_1, fq_2, linkera, linkerb, out_f1, out_f2):
    with open_fastq(fq_1) as f1, open_fastq(fq_2) as f2:
        o1 = write_fastq(out_f1)
        o2 = write_fastq(out_f2)
        while True:
            rec1 = read_rec(f1)
            record1, _ = linker_detect(rec1, linkera, linkerb)
            out_rec1 = dict2string(record1)
            rec2 = read_rec(f2)
            record2, _ = linker_detect(rec2, linkera, linkerb)
            out_rec2 = dict2string(record2)
            if (len(rec1["name"]) == 0) or (len(rec2["name"]) == 0):
                break
            o1.write(out_rec1 + "\n")
            o2.write(out_rec2 + "\n")


def read_rec(f):
    lines = []
    for i in range(4):
        line = f.readline().strip()
        lines.append(line)
    rec = {
        'name': lines[0],
        'seq': lines[1],
        'opt': lines[2],
        'qual': lines[3]
    }
    return rec


def dict2string(rec):
    nam = rec["name"]
    seq = rec["seq"]
    opt = rec["opt"]
    qua = rec["qual"].strip("\n")
    out_str = '\n'.join([nam, seq, opt, qua])
    return out_str


def linker_detect(record, linkera, linkerb):
    sequ = record.get("seq")
    aligner = Aligner(sequ, 0.13)
    resulta = aligner.locate(linkera)
    resultb = aligner.locate(linkerb)
    id_ = record["name"]
    if resulta is not None:
        st = resulta[0]
        if len(sequ[:st]) > 16:
            record["seq"]  = sequ[:st]
            record["qual"] = record["qual"][:st] + "\n"
    elif resultb is not None:
        st = resultb[0]
        if len(sequ[:st]) > 16:
            record["seq"]  = sequ[:st]
            record["qual"] = record["qual"][:st] + "\n"
    return record, id_


def open_fastq(fname, mode='r'):
    if fname.endswith('.gz'):
        fh = gzip.open(fname, mode=mode+'b')
        fh = io.TextIOWrapper(fh)
    else:
        fh = open(fname, mode=mode)
    return fh


def write_fastq(fname, mode='wt'):
    fh = gzip.open(fname, mode)
    return fh


if __name__ == "__main__":
    trim_linker()
