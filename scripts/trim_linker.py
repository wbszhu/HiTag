import sys
import logging
import gzip
import io
from queue import Empty
import multiprocessing as mp
import subprocess as subp

from cutadapt._align import Aligner
import click


LOGGING_FMT = "%(name)-20s %(levelname)-7s @ %(asctime)s: %(message)s"
LOGGING_DATE_FMT = "%m/%d/%y %H:%M:%S"
log = logging.getLogger(__name__)
handler = logging.StreamHandler(sys.stderr)
handler.setFormatter(logging.Formatter(fmt=LOGGING_FMT, datefmt=LOGGING_DATE_FMT))
log.addHandler(handler)
log.setLevel(logging.DEBUG)


@click.command("")
@click.argument("fq_1")
@click.argument("fq_2")
@click.argument("out_f1")
@click.argument("out_f2")
@click.argument("linkera")
@click.argument("linkerb")
@click.option( "--processes", "-p",
    default=1,
    help="Use how many processes to run")
@click.option( "--chunk_size", "-c",
    default=100000,
    help="Use how many processes to run")
def trim_linker(fq_1, fq_2, linkera, linkerb, out_f1, out_f2, processes, chunk_size):
    task_queue = mp.Queue()
    workers = [
        mp.Process(target=worker, args=(task_queue, linkera, linkerb, out_f1, out_f2))
        for _ in range(processes)
    ]
    log.info("%d workers spawned for trim linker"%len(workers))

    for w in workers:
        w.start()
    with open_fastq(fq_1) as f1, open_fastq(fq_2) as f2:
        while True:
            chunk = read_chunk(f1, f2, chunk_size)
            if len(chunk) == 0:
                break
            task_queue.put(chunk)
    for _ in workers:
        task_queue.put(None)  # End flag
    for w in workers:
        w.join()

    tmp_files_1 = [f"{out_f1}.tmp.{w.pid}" for w in workers]
    tmp_files_2 = [f"{out_f2}.tmp.{w.pid}" for w in workers]
    log.info(f"Merge tmp files read1: {tmp_files_1}")
    log.info(f"Merge tmp files read2: {tmp_files_2}")
    p1 = subp.Popen(f"cat {' '.join(tmp_files_1)} > {out_f1}", shell=True)
    p2 = subp.Popen(f"cat {' '.join(tmp_files_2)} > {out_f2}", shell=True)
    p1.wait()
    p2.wait()
    log.info("Remove tmp files.")
    subp.check_call(f"rm {' '.join(tmp_files_1 + tmp_files_2)}", shell=True)



def read_rec(f):
    lines = []
    for _ in range(4):
        line = f.readline().strip()
        lines.append(line)
    rec = {
        'name': lines[0],
        'seq': lines[1],
        'opt': lines[2],
        'qual': lines[3]
    }
    return rec


def read_chunk(f1, f2, chunk_size):
    chunk = []
    for _ in range(chunk_size):
        rec1 = read_rec(f1)
        rec2 = read_rec(f2)
        if (len(rec1['name']) == 0) or (len(rec2['name']) == 0):
            return chunk
        chunk.append((rec1, rec2))
    return chunk


def worker(task_queue, linker_a, linker_b, out_path_1, out_path_2):
    pid = mp.current_process().pid
    log.info(f"Worker {pid} start")
    tmp1 = f"{out_path_1}.tmp.{pid}"
    tmp2 = f"{out_path_2}.tmp.{pid}"
    o1 = open_fastq(tmp1, 'w', force_gz=out_path_1.endswith('.gz'))
    o2 = open_fastq(tmp2, 'w', force_gz=out_path_2.endswith('.gz'))
    while True:
        try:
            chunk = task_queue.get()
            if chunk is None:
                raise Empty
        except Empty:
            o1.close()
            o2.close()
            log.info(f"Worker {pid} stoped.")
            break
        for rec1, rec2 in chunk:
            record1, rec1_has_linker = linker_detect(rec1, linker_a, linker_b)
            record2, rec2_has_linker = linker_detect(rec2, linker_a, linker_b)
            if (len(record1["seq"])>=16) and (len(record2["seq"])>=16):
                out_rec1 = dict2string(record1)
                out_rec2 = dict2string(record2)
                o1.write(out_rec1 + "\n")
                o2.write(out_rec2 + "\n")


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
    has_linker = False
    if resulta is not None:
        st = resulta[0]
        if (resulta[1]-resulta[0]) >= 13:
            record["seq"]  = sequ[:st]
            record["qual"] = record["qual"][:st] + "\n"
            has_linker = True    
    if (not has_linker) and (resultb is not None):
        st = resultb[0]
        if (resultb[1]-resultb[0]) >= 13:
            record["seq"]  = sequ[:st]
            record["qual"] = record["qual"][:st] + "\n"
            has_linker = True
    return record, has_linker


def open_fastq(fname, mode="r", force_gz=False):
    if fname.endswith(".gz") or force_gz:
        fh = gzip.open(fname, mode)
        fh = io.TextIOWrapper(fh)
    else:
        fh = open(fname, mode)
    return fh


if __name__ == "__main__":
    trim_linker()
