import sys
import io
import os
import gzip
import logging
import subprocess
from queue import Empty
import multiprocessing as mp

import click
import numpy as np
import h5py


LOGGING_FMT = "%(name)-20s %(levelname)-7s @ %(asctime)s: %(message)s"
LOGGING_DATE_FMT = "%m/%d/%y %H:%M:%S"
log = logging.getLogger(__name__)
handler = logging.StreamHandler(sys.stderr)
handler.setFormatter(logging.Formatter(fmt=LOGGING_FMT, datefmt=LOGGING_DATE_FMT))
log.addHandler(handler)
log.setLevel(logging.DEBUG)


types = ("valid", "self-circle", 'dangling-end', "dump", "re-ligation")


def open_file(fname, mode='r'):
    if fname.endswith('.gz'):
        fh = gzip.open(fname, mode=mode+'b')
        fh = io.TextIOWrapper(fh)
    else:
        fh = open(fname, mode=mode)
    return fh


def merge_tmp_files(tmp_files, merged_file, header=""):
    if os.path.exists(merged_file):
        os.remove(merged_file)
    if header:
        with open(merged_file, 'w') as f:
            f.write(header)
    # merge tmp files
    cmd = "cat " + " ".join(tmp_files) + " >> " + merged_file
    subprocess.check_call(cmd, shell=True)
    # remove tmp files
    cmd = "rm " + " ".join(tmp_files)
    subprocess.check_call(cmd, shell=True)


def parse_line_pairs(line):
    """ parse pairs file format file line
    return read_id, chr1, pos1, chr2, pos2, strand1, strand2 """
    items = line.strip().split()
    read_id, chr1, pos1, chr2, pos2, strand1, strand2 = items[:7]
    pos1, pos2 = int(pos1), int(pos2)
    return read_id, chr1, pos1, chr2, pos2, strand1, strand2


def parse_line_bed6(line):
    """ parse bed6 line
    return chr, start, end, name, score, strand
    """
    items = line.strip().split()
    chr_, start, end, name, score, strand = items[:6]
    start, end, score = [int(i) for i in (start, end, score)]
    return chr_, start, end, name, score, strand


def read_and_skip_headers(file, comment="#"):
    headers = []
    for line in file:
        if line.startswith(comment):
            headers.append(line)
        else:
            break
    file.seek(len("".join(headers)))
    return headers


TIME_OUT = 1
CHUNK_SIZE = 10000


def load_rest_sites(frag_file):
    """
    load restriction cutting sites(positions) into memory.
    """
    rest_sites = {}
    if frag_file.endswith(".gz") or frag_file.endswith(".bed"):
        with open_file(frag_file) as f:
            for idx, line in enumerate(f):
                line = line.strip()
                rest_seq = None
                if line.startswith("# rest_seq"):
                    rest_seq = line.split()[-1]
                    continue
                chr_, start, end, _, _, _ = parse_line_bed6(line)
                if start == 0:  # remove first
                    continue
                rest_sites.setdefault(chr_, [])
                rest_sites[chr_].append(start)
        for k in rest_sites.keys():
            rest_list = rest_sites[k]
            rest_list.pop()  # remove last
            rest_sites[k] = np.asarray(rest_list)
    else:  # hdf5 file
        with h5py.File(frag_file, 'r') as f:
            rest_seq = f.attrs['rest_seq']
            chromosomes = list(f['chromosomes'].keys())
            for chr_ in chromosomes:
                fragments = f['chromosomes'][chr_][()]
                rest_sites[chr_] = fragments[1:-1]
    return rest_sites, rest_seq


def find_frag(rest_sites, pos):
    """
    Assign genome interval to fragment.
    Parameters
    ----------
    rest_sites : `numpy.ndarray`
        restriction sites(positions) of a chromosome.
    pos : int
        PET mapped position
    Return
    ------
    n : int
        Index number of fragment.
    """
    search = pos + 1  # search point
    sites = rest_sites
    frag_idx = sites.searchsorted(search, side='right')
    return frag_idx


def read_chunk(file, chunk_size=CHUNK_SIZE):
    """ read and parse bedpe file,
    return a list of bedpe items
    """
    chunk = []
    for i in range(chunk_size):
        line = file.readline()
        if not line:
            if i == 0:
                raise StopIteration
            else:
                return chunk
        items = parse_line_pairs(line)
        chunk.append(items)
    return chunk


def pairs_type(rest_sites, pairs_items, threshold_span):
    """
    judge interaction type,
    Parameters
    ----------
    rest_sites : dict
        Each chromosome's restriction start positions.
    pairs_items : list
        pairs fields.
    threshold_span : int
        If span large than this, will not check.
        Use -1 to force assign PET to fragment.
    Return
    ------
    types: {"valid", "self-circle", "re-ligation", "dangling-end", "dump"}
    frag1: int
    frag2: int
    """
    read_id, chr1, pos1, chr2, pos2, strand1, strand2 = pairs_items[:7]

    # inter chromosome interaction
    if chr1 != chr2:
        if threshold_span != -1:
            return "valid", None, None
        else:
            frag1 = find_frag(rest_sites[chr1], pos1)
            frag2 = find_frag(rest_sites[chr2], pos2)
            return "valid", frag1, frag2

    # intra chromosome interaction
    span =  abs(pos1 - pos2)
    if span > threshold_span and threshold_span != -1:
        # safe span: interaction distance enough long
        if threshold_span != -1:
            return "valid", None, None
    else:
        # unsafe span, need to check
        sites = rest_sites[chr1]
        frag1 = find_frag(sites, pos1)
        frag2 = find_frag(sites, pos2)
        tp = None

        if frag1 == frag2:
            if (strand1 == "-") and (strand2 == "+"):
                tp = "self-circle"
            elif (strand1 == "+") and (strand2 == "-"):
                tp = "dangling-end"
            else:
                tp =  "dump"
        elif (frag2 - frag1 == 1):
            tp = "re-ligation"
        else:
            tp = "valid"
        return tp, frag1, frag2


def worker(task_queue,
           rest_sites,
           files,
           threshold_span, counts, lock):
    current = mp.current_process().pid

    handles = {tp: open(f, 'w') for tp, f in files.items()}

    l_t_counts = { k: 0 for k in types }
    l_p_counts = { k: 0 for k in [ oc for oc in [i+j for i in "+-" for j in "+-"]] }

    while 1:
        try:
            chunk = task_queue.get()
            if chunk is None:
                raise Empty
        except Empty:
            for f in handles.values():
                f.close()
            lock.acquire()
            for key, l_counts in zip(['type', 'position'], [l_t_counts, l_p_counts]):
                counts_ = counts[key]
                for k, v in l_counts.items():  # update 
                    counts_[k] += v
                counts[key] = counts_  # save to manager dict
            lock.release()
            log.debug("Process-%d done."%current)
            break

        for items in chunk:
            try:
                type_, frag1, frag2 = pairs_type(rest_sites, items, threshold_span)
            except KeyError as e:
                log.warning(str(e))
                continue

            out_line = "\t".join([str(i) for i in items])
            if frag1:
                out_line += "\t{}\t{}".format(frag1, frag2)
            out_line += "\n"

            for tp in types:
                if tp == type_:
                    handles[tp].write(out_line)
                    l_t_counts[tp] += 1

            if frag1:
                position = items[-2] + items[-1]
                l_p_counts[position] += 1


def log_counts(counts, log_file):
    log.info("Noise reduce result count:")
    t_counts = counts['type']
    total = sum(list(t_counts.values()))
    for tp in types:
        n = t_counts[tp]
        r = 0 if total == 0 else n / total
        msg = f"{tp}\t{n}\tpercent\t{r:.2%}"
        log.info("\t"+msg)
    msg = "total\t{}".format(total)
    log.info("\t" + msg)

    with open(log_file, 'w') as f:
        f.write("# type count\n")
        for tp in types:
            n = t_counts[tp]
            f.write(f"{tp}\t{n}\n")
        f.write("total\t{}\n".format(total))
        f.write("\n")
        f.write("# position count\n")
        for k, v in counts['position'].items():
            f.write("{}\t{}\n".format(k, v))
        f.write("\n")


@click.command(name="noise_reduce")
@click.argument("pairs")
@click.argument("output")
@click.option("--restriction", "-r",
    required=True,
    help="File which recorded the position of all restriction sites(DNA fragments)"
        "in the genome. BED6 format or hdf5 file(default) is supported."
        "Which can be produced by 'dlohic extract_fragments' command")
@click.option("--processes", "-p", 
    default=1,
    help="Use how many processes to run.")
@click.option("--threshold-span", "-s",
    default=-1,
    show_default=True,
    help="Threshold of pair span. Use -1 to force check.")
@click.option("--log-file",
    default="noise_reduce.log",
    help="Sperate log file for storage some count information.")
def _main(pairs, output,
         restriction, processes,
         threshold_span,
         log_file):
    """
    Remove noise in HiC data.
    \b
    types:
        1. self-circle
        2. dangling-end
        3. dump
        4. re-ligation
    """

    log.info("noise reduce on file %s"%pairs)

    log.info("loading restriction sites from file: {}".format(restriction))
    rest_sites, _ = load_rest_sites(restriction)

    # init counts
    lock = mp.Lock()
    counts = mp.Manager().dict()
    counts['type'] = { k: 0 for k in types }
    counts['position'] = {  # permutation of "+-"
        k: 0 for k in [
            oc for oc in [i+j for i in "+-" for j in "+-"]
        ]
    }

    task_queue = mp.Queue()
    workers = [mp.Process(target=worker,
                          args=(task_queue,
                                rest_sites,
                                {tp: f"{output}.{tp}.tmp.{i}" for tp in types},
                                threshold_span,
                                counts, lock))
               for i in range(processes)]

    log.info("%d workers spawned for noise refuce"%len(workers))

    for w in workers:
        w.start()

    with open_file(pairs) as f:
        headers = read_and_skip_headers(f)
        while 1:
            try:
                task_queue.put(read_chunk(f))
            except StopIteration:
                break

    for w in workers:
        task_queue.put(None)

    for w in workers:
        w.join()

    log.info("merging tmporary files.")
    # merge tmp files
    for tp in types:
        tmp_files = [f"{output}.{tp}.tmp.{i}" for i in range(processes)]
        if tp == "valid":
            header = "".join(headers)
        else:
            header = ""
        merge_tmp_files(tmp_files, f"{output}.{tp}.pairs", header)

    log_counts(counts, log_file)
    log.info("Noise reduce done.")


main = _main.callback


if __name__ == "__main__":
    eval("_main()")
