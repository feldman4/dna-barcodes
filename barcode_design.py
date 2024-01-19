# MIT License

# Copyright (c) 2020 David Jacob Feldman

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


import argparse
from collections import Counter, defaultdict
from itertools import product
import logging
import multiprocessing.pool as mpp
import regex as re
import time

import Levenshtein
import numpy as np
import pandas as pd
import scipy.sparse
from tqdm.auto import tqdm


logger = logging.getLogger(__name__)


def generate_all_barcodes(n):
    if n > 11:
        raise ValueError
        
    return [''.join(x) for x in product(*(n*[list('ACTG')]))]


def calculate_gc(s):
    s = s.upper()
    return (s.count('G') + s.count('C')) / len(s)


def create_barcode_set(n, k, homopolymer, gc_min, gc_max, exclude=None,
    limit=None, progress=tqdm, cores=1):
    if 1 < max(gc_min, gc_max):
        print('Reminder: gc_min and gc_max parameters are on a 0..1 scale')
    
    df_bcs = (pd.DataFrame({'barcode': generate_all_barcodes(n)})
     .assign(gc=lambda x: x['barcode'].apply(calculate_gc))
    )

    logger.info(f'Generated {len(df_bcs)} barcodes of length {n}')

    barcodes = (df_bcs
     .query('@gc_min <= gc <= @gc_max')
     .loc[lambda x: ~(x['barcode'].apply(lambda y:
        has_homopolymer(y, homopolymer)))]
     ['barcode'].pipe(list))

    logger.info(f'Retained {len(barcodes)} barcodes after '
        'filtering for GC content and homopolymers')

    if exclude is not None:
        pat = re.compile(exclude)
        barcodes = [bc for bc in barcodes if not pat.findall(bc)]
        logger.info(f'Retained {len(barcodes)} barcodes after removing matches '
            f'to "{exclude}"')

    if limit:
        rs = np.random.RandomState(0)
        barcodes = rs.choice(barcodes, limit, replace=False)
    logger.info('Assigning barcodes to hash buckets...')
    hash_buckets = build_khash(barcodes, k)
    logger.info('Calculating distances within buckets...')

    if cores > 1:
        from multiprocessing import Pool
        with Pool(processes=cores) as p:
            work = [([x], k, None) for x in hash_buckets]
            D = {}
            with progress(total=len(hash_buckets)) as pbar:
                for d in p.istarmap(sparse_dist, work):
                    pbar.update()
                    D.update(d)
    else:
        D = sparse_dist(hash_buckets, k, progress=progress)
    cm = sparse_view(barcodes, D)

    logger.info(f'Selecting barcodes with minimum edit distance {k}...')
    group_ids = [0] * len(barcodes)
    selected = maxy_clique_groups(cm, group_ids)

    selected_barcodes = [barcodes[x] for x in selected]
    logger.info(f'Selected {len(selected_barcodes)} barcodes')
    
    return (df_bcs
            .query('barcode == @selected_barcodes')
            .assign(n=n, k=k, homopolymer=homopolymer, 
                gc_min=gc_min, gc_max=gc_max))


def prune(barcodes, k, cores=1, progress=tqdm):
    logger.info('Assigning barcodes to hash buckets...')
    hash_buckets = build_khash(barcodes, k)
    logger.info('Calculating distances within buckets...')

    if cores > 1:
        from multiprocessing import Pool
        with Pool(processes=cores) as p:
            work = [([x], k, None) for x in hash_buckets]
            D = {}
            with progress(total=len(hash_buckets)) as pbar:
                for d in p.istarmap(sparse_dist, work):
                    pbar.update()
                    D.update(d)
    else:
        D = sparse_dist(hash_buckets, k, progress=progress)
    cm = sparse_view(barcodes, D)

    logger.info(f'Selecting barcodes with minimum edit distance {k}...')
    group_ids = [0] * len(barcodes)
    selected = maxy_clique_groups(cm, group_ids)

    selected_barcodes = [barcodes[x] for x in selected]
    logger.info(f'Selected {len(selected_barcodes)} barcodes')

    return selected_barcodes



def check_barcode_set(barcodes, k, distance='Levenshtein', max_to_check=1e6):
    """Returns list of barcode pairs that fail to satisfy distance k.
    Default distance is 
    """
    if distance == 'Levenshtein':
        distance = Levenshtein.distance
    failures = []
    for a in barcodes:
        for b in barcodes:
            d = distance(a, b)
            max_to_check -= 1
            if a != b and d < k:
                failures += [(a, b, d)]
            if max_to_check == 0:
                return failures
    return failures


def parse_args():
    description = 'Design DNA barcodes with guaranteed Levenshtein edit distance'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('--distance', type=int, default=2,
        help='minimum pairwise edit distance')
    parser.add_argument('--length', type=int, default=7, 
        help='barcode length')
    parser.add_argument('--gc_min', type=int, default=30,
        help='minimum GC percentage')
    parser.add_argument('--gc_max', type=int, default=70,
        help='maximum GC percentage')
    parser.add_argument('--homopolymer', type=int, default=4,
        help='maximum homopolymer')
    parser.add_argument('--limit', type=int, default=0,
        help='limit number of sequences input to design process, 0=no limit')
    parser.add_argument('--max_to_check', type=int, default=int(1e6),
        help='maximum number of barcode pairs to verify edit distance')
    parser.add_argument('--verbosity', type=int, default=2,
        help='logging level: <=2 logs info, <=3 logs warnings')
    parser.add_argument('--exclude', type=str, default='',
        help=('pattern to exclude in barcodes, e.g., --exclude="GTTC|ATTC" '
            'or --exclude="(?:GTCC){s<2}" (see https://pypi.org/project/regex/ '
            'for fuzzy matching syntax)'))
    parser.add_argument('--cores', type=int, default=1,
        help='number of parallel processes for distance calculation')

    return parser.parse_args()


def handle_failures(barcodes, k, max_to_check=1e6, num_failures_to_print=10):
    max_to_check = int(max_to_check)
    
    logger.info(f'Validating barcodes...')
    failures = check_barcode_set(barcodes, k, max_to_check=max_to_check)
    if failures:
        logger.warning('!! Failures detected !!')
        
        for failure in failures[:num_failures_to_print]:
            logger.warning('{} {} distance={}'.format(*failure))
        if len(failures) > num_failures_to_print:
            logger.warning('...')

    else:
        if max_to_check > (len(barcodes)**2):
            logger.info('All barcodes passed!')
        else:
            logger.info(f'Checked {max_to_check:,} barcode pairs, all passed!')

    return failures


def istarmap(self, func, iterable, chunksize=1):
    """https://stackoverflow.com/questions/57354700/starmap-combined-with-tqdm
    """
    self._check_running()
    if chunksize < 1:
        raise ValueError(
            "Chunksize must be 1+, not {0:n}".format(
                chunksize))

    task_batches = mpp.Pool._get_tasks(func, iterable, chunksize)
    result = mpp.IMapIterator(self)
    self._taskqueue.put(
        (
            self._guarded_task_generation(result._job,
                                          mpp.starmapstar,
                                          task_batches),
            result._set_length
        ))
    return (item for chunk in result for item in chunk)


mpp.Pool.istarmap = istarmap


def main():
    args = parse_args()

    logging.basicConfig(format='%(asctime)s -- %(message)s',
                   datefmt='%Y-%m-%d %H:%M:%S')
    logging_level = args.verbosity * 10
    logging.root.setLevel(logging_level)
    logging.root.handlers[0].addFilter(lambda x: 'NumExpr' not in x.msg)
    progress = tqdm if logging_level <= 20 else None

    limit = None if args.limit == 0 else args.limit
    exclude = None if args.exclude == '' else args.exclude
    df_bcs = create_barcode_set(args.length, args.distance, args.homopolymer, 
        args.gc_min/100, args.gc_max/100, 
        cores=args.cores, exclude=exclude, limit=limit, progress=progress)

    failures = handle_failures(
        df_bcs['barcode'], args.distance, args.max_to_check)

    failure_tag = '.failure' if failures else ''

    filename = f'barcodes_n{args.length}_k{args.distance}{failure_tag}.csv'
    filename = timestamp(filename)
    df_bcs.to_csv(filename, index=None)

    logger.info(f'Output written to {filename}')


### FUNCTIONS BELOW FROM OpticalPooledScreens repository ###
# https://github.com/feldman4/OpticalPooledScreens/blob/master/ops/pool_design.py ###

def has_homopolymer(x, n):
    a = 'A'*n in x
    t = 'T'*n in x
    g = 'G'*n in x
    c = 'C'*n in x
    return a | t | g | c


def build_khash(xs, k, return_dict=False):
    D = defaultdict(list)
    for x in xs:
        for h in khash(x, k):
             D[h].append(x)

    D = {k: sorted(set(v)) for k,v in D.items()}
    if return_dict:
        return D
    else:
        hash_buckets = list(D.values())
        return hash_buckets


def khash(s, k):
    """Divide a string into substrings suitable for checking edit distance of 
    `k`. Two strings of the same length with Levenshtein edit distance less 
    than `k` will share at least one substring. 
    """
    n = len(s)
    window = int(np.ceil((n - k) / float(k)))
    s = s + s
    arr = []
    for i in range(n):
        for j in (0, 1):
            arr += [((i + j) % n, s[i:i+window])]
    return arr


def sparse_dist(hash_buckets, k, distance=None, progress=None):
    """Entries less than k only.
    """
    if distance is None:
        distance = Levenshtein.distance
    if progress is None:
        progress = lambda x: x
    D = {}
    for xs in progress(hash_buckets):
        for i, a in enumerate(xs):
            for b in xs[i+1:]:
                d = distance(a,b)
                if d < k:
                    key = tuple(sorted((a,b)))
                    D[key] = d
    return D


def sparse_view(xs, D, symmetric=True):
    """string barcodes
    """
    assert len(xs) == len(set(xs))
    mapper = {x: i for i, x in enumerate(xs)}
    f = lambda x: mapper[x]
    if len(D) == 0:
        i, j, data = [], [], []
    else:
        i, j, data = zip(*[(f(a), f(b), v) for (a, b), v in D.items()])
        # sparse matrix uses zero for missing values
        data = np.array(data) >= 0
        i = np.array(i)
        j = np.array(j)

    n = len(xs)
    cm = scipy.sparse.coo_matrix((data, (i, j)), shape=(n, n))

    if symmetric:
        cm = (cm + cm.T).tocsr()
        
    return cm


def maxy_clique_groups(cm, group_ids, verbose=False):
    """Prioritizes groups with the fewest selected barcodes.
    Prioritizing groups with the fewest remaining barcodes could give
    better results.
    """

    # counts => group_id
    d1 = defaultdict(set)
    for id_, counts in Counter(group_ids).items():
        d1[counts] |= {id_}

    # group_id => indices
    d2 = defaultdict(list)
    for i, id_ in enumerate(group_ids):
        d2[id_] += [i]
    # .pop() takes from the end of the list
    d2 = {k: v[::-1] for k,v in d2.items()}

    # group_id => # selected
    d3 = Counter()

    selected = []
    available = np.array(range(len(group_ids)))

    while d1:
        if verbose and (len(selected) % 1000) == 0:
            print(len(selected))
    #     assert cm[selected, :][:, selected].sum() == 0

        # pick a group_id from the lowest bin
        count = min(d1.keys())
        id_ = d1[count].pop()

        # remove bin if empty
        if len(d1[count]) == 0:
            d1.pop(count)

        # discard indices until we find a new one
        index = None
        while d2[id_]:
            index = d2[id_].pop()
            # approach 1: check for conflict every time
            # cm[index, selected].sum() == 0
            # approach 2: keep an array of available indices
            if index in available:
                break
        else:
            index = None

        # keep index
        if index:
            selected.append(index)
            d3[id_] += 1
            available = available[available != index]
            # get rid of incompatible barcodes
            remove = cm[index, available].indices
            mask = np.ones(len(available), dtype=bool)
            mask[remove] = False
            available = available[mask]


        # move group_id to another bin
        n = len(d2[id_])
        if n > 0:
            d1[n] |= {id_}

    return selected


def timestamp(filename='', fmt='%Y%m%d_%H%M%S', sep='.'):
    stamp = time.strftime(fmt)
    pat= r'(.*)\.(.*)'
    match = re.findall(pat, filename)
    if match:
        return sep.join([match[0][0], stamp, match[0][1]])
    elif filename:
        return sep.join([filename, stamp])
    else:
        return stamp


if __name__ == '__main__':
    main()
