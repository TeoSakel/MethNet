#!/usr/bin/env python3
"""
Clone Xenahubs Repositories.

This script creates a local (heterozygote) clone of the xenahubs.net browser.
The difference between the 2 is that files are ugziped where possible.

Usage:
    clone_xena.py <hub> [options]
    clone_xena.py (--list-hubs)
    clone_xena.py (-h | --help)

Options:
    -h --help            Show this message.
    --list-hubs          List known hubs
    -l --list            List datasets in <hub>. Formated as <type>: <name>
    -t TYPE --type=TYPE  Filter datasets based on TYPE
    -r RE --regex RE     Filter datasets based on regular expression RE
    -o DIR --output DIR  Directory to write results as DIR/TYPE/NAME.tsv [default: .]
"""
import sys
import zlib
import time
import re
import warnings
from pathlib import Path
from requests import get
from docopt import docopt
with warnings.catch_warnings():
    warnings.filterwarnings("ignore", category=FutureWarning)
    import xenaPython as xena

# TODO: allow a dry-run option to just download urls?
#       problem is that GET always downloads the page

HUBS = ('gdc', 'icgc', 'pancanatlas', 'tcga', 'toil', 'ucscpublic')


def hub_url(hub):
    """Returns hub url given hub name"""
    return 'https://{}.xenahubs.net'.format(hub)


def add_suffix(path, suffix):
    return path.parent / (path.name + suffix)


def try_get(hub, dat):
    """Tries to download dataset from hub. If it fails it retries with the .gz suffix"""
    url = '/'.join((hub_url(hub), 'download', dat))
    response = get(url)
    meta = get(url + '.json')
    if response.status_code == 200:
        return response, meta
    return get(url + '.gz'), meta


def get_path(outdir, fpath):
    fpath = Path(fpath)
    parts = fpath.parts
    if parts[0] == "https:":  # it's url
        fpath = Path('').joinpath(*parts[3:])  # drop https://<hub>/download
    try:
        fpath = fpath.relative_to('probeMap')
    except ValueError:
        pass
    return outdir / fpath


def write_response(response, meta, outdir):
    """Write response as a .tsv file and metadata as .json in outdir"""
    # Decide on path
    fpath = get_path(outdir, response.url)
    fpath.parent.mkdir(parents=True, exist_ok=True)
    # uncompress if needed
    content = response.content
    if fpath.suffix == '.gz':
        try:
            content = zlib.decompress(content, 16+zlib.MAX_WBITS)
            fpath = fpath.with_suffix('.tsv')
        except zlib.error:
            pass
    elif not fpath.suffix == '.tsv':
        fpath = add_suffix(fpath, '.tsv')
    # write
    fpath.write_bytes(content)
    fpath.with_suffix('.json').write_bytes(meta.content)


def clone(hubdir, datasets):
    responsive_urls = {}
    for dtype, url in datasets:
        print(dtype, ":", url)
        # prepare dir and file names
        typedir = hubdir / dtype
        fpath = get_path(typedir, url)
        if add_suffix(fpath, '.json').exists() or fpath.with_suffix('.json').exists():
            print('already exists')
            continue
        # get response and write
        response, meta = try_get(hub, url)
        if not response.status_code == 200:
            print('server not responding')
            continue
        time.sleep(2)  # to avoid ban from the server
        write_response(response, meta, typedir)
        # record successes
        try:
            responsive_urls[dtype].extend([response.url, meta.url])
        except KeyError:
            responsive_urls[dtype] = [response.url, meta.url]

    # write in each directory the urls of all the files for wget if needed
    for dtype, urls in responsive_urls.items():
        with open(typedir / 'urls.txt', 'w') as fid:
            fid.write('\n'.join(urls))

if __name__ == "__main__":
    args = docopt(__doc__)
    # case 0: show help
    if args['--help']:
        print(__doc__)
        sys.exit(0)
    # case 1: show hubs
    if args['--list-hubs']:
        print('\n'.join(HUBS))
        sys.exit(0)

    hub = args['<hub>']
    dtype = args['--type']
    regex = args['--regex']

    # filter
    host = hub_url(hub)
    cohorts = xena.all_cohorts(host, [])
    datasets = [(dset['type'], dset['name']) for dset in xena.dataset_list(host, cohorts)
                if type(dset['type']) is str] # some datasets have None instead of strings
    if dtype:
        datasets = [dset for dset in datasets if re.match(dtype, dset[0], flags=re.IGNORECASE)]
    if regex:
        datasets = [dset for dset in datasets if re.search(regex, dset[1])]

    # case 2: show datasets
    if args['--list']:
        datasets.sort()
        for dset in datasets:
            print('{}:\t{}'.format(*dset))
        sys.exit(0)

    # case 3: main!
    hubdir = Path(args['--output']) / Path(hub)
    clone(hubdir, datasets)

