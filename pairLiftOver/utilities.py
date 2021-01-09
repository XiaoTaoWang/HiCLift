import subprocess, sys, os, io, logging, cooler, pairLiftOver
from pairtools import _headerops
from pyliftover import LiftOver
from pairLiftOver.io import open_pairs, _pixel_to_reads, _pairs_write

log = logging.getLogger(__name__)

def extract_chrom_sizes(fil):

    chromsizes = []
    with open(fil, 'r') as source:
        for line in source:
            tmp = line.rstrip().split()
            chromsizes.append((tmp[0], int(tmp[1])))
    
    return chromsizes

def make_mapping_table(chroms_path, lo, resolution=200):

    chromsizes = extract_chrom_sizes(chroms_path)
    D = {}
    for c, Len in chromsizes:
        for s in range(0, Len, resolution):
            key = (c, s)
            tmp = lo.convert_coordinate(c, s)
            if tmp is None:
                continue
            if len(tmp)==1:
                value = (tmp[0][0], tmp[0][1])
                D[key] = value
    
    return D

def liftover(in_path, out_pre, in_format, out_format, in_chroms, out_chroms, in_assembly, out_assembly,
    chain_file, resolution=500, nproc_in=8, nproc_out=8, tmpdir='/tmp', memory='4G'):
    
    tmpdir = os.path.abspath(os.path.expanduser(tmpdir))
    if not os.path.exists(tmpdir):
        os.makedirs(tmpdir)
    
    outfolder, out_pre = os.path.split(out_pre)
    out_path = os.path.join(tmpdir, '{0}.pairs.gz'.format(out_pre))

    instream = open_pairs(in_path, mode='r', data_format=in_format, nproc=nproc_in)
    outstream = open_pairs(out_path, mode='w', data_format='pairs', nproc=nproc_out)
    
    # write header
    log.info('Writing headers ...')
    chromsizes = extract_chrom_sizes(out_chroms)
    header = _headerops.make_standard_pairsheader(
        assembly=out_assembly, chromsizes=chromsizes,
        columns=['readID', 'chrom1', 'pos1', 'chrom2', 'pos2', 'strand1', 'strand2'],
        shape='upper triangle'
    )
    header.append('#pairLiftOver: coordinates transformed from {0}'.format(in_assembly))

    outstream.writelines((l+'\n' for l in header))
    outstream.flush()
    
    # build the mapping table at the given resolution
    if not chain_file is None:
        lo = LiftOver(chain_file)
    else:
        lo = LiftOver(in_assembly, out_assembly)
    if not resolution is None:
        log.info('Building the mapping table at the resolution: {0}'.format(resolution))
        mapping_table = make_mapping_table(in_chroms, lo, resolution)
    else:
        mapping_table = None
    
    chrom_index = dict(_headerops.get_chrom_order(out_chroms))

    if in_format in ['cooler', 'juicer']:
        body_stream = instream
    else:
        _, body_stream = _headerops.get_header(instream)
    
    # sort command
    command = r'''/bin/bash -c 'export LC_COLLATE=C; export LANG=C; sort -k 2,2 -k 4,4 -k 3,3n -k 5,5n --stable {0} {1} -S {2} {3}'''.format(
        '--parallel={0}'.format(nproc_out), '--temporary-directory={0}'.format(tmpdir), memory,'--compress-program=lz4c'
    )
    command += "'"
    
    log.info('Converting, sorting, and compressing ...')
    total_count = 0
    mapped_count = 0
    with subprocess.Popen(command, stdin=subprocess.PIPE, bufsize=-1, shell=True, stdout=outstream) as process:
        stdin_wrapper = io.TextIOWrapper(process.stdin, 'utf-8')
        for line in body_stream:
            if in_format in ['cooler', 'juicer']:
                total_count, mapped_count = _pixel_to_reads(stdin_wrapper,
                                                            line,
                                                            chrom_index,
                                                            mapping_table,
                                                            lo, resolution,
                                                            in_format,
                                                            total_count,
                                                            mapped_count)
            elif in_format in ['pairs', 'hic-pro']:
                total_count, mapped_count = _pairs_write(stdin_wrapper,
                                                         line,
                                                         chrom_index,
                                                         mapping_table,
                                                         lo, resolution,
                                                         in_format,
                                                         total_count,
                                                         mapped_count)
        stdin_wrapper.flush()
        process.communicate()
    
    log.info('{0:,} / {1:,} pairs were uniquely mapped to {2}'.format(mapped_count, total_count, out_assembly))
    if instream != sys.stdin:
        instream.close()
    if outstream != sys.stdout:
        outstream.close()
    
    # handle with different output formats
    if out_format == 'pairs':
        dest = os.path.join(outfolder, os.path.split(out_path)[1])
        command = ['mv', out_path, dest]
        subprocess.check_call(' '.join(command), shell=True)
    else:
        command = ['pairix', out_path]
        subprocess.check_call(' '.join(command), shell=True)
        if out_format == 'cooler':
            clr = cooler.Cooler(in_path)
            binsize = clr.binsize
            bin_label = ':'.join([out_chroms, str(binsize)])
            unit, denom = ('kb', 1000) if (binsize / 1000 < 1000) else ('mb', 1000000)
            reslabel = str(binsize // denom) + unit
            outcool = os.path.join(outfolder, '{0}.{1}.cool'.format(out_pre, reslabel))
            command = ['cooler', 'cload', 'pairix', '--assembly', out_assembly, '--nproc', str(nproc_out),
                       '--max-split 12', bin_label, out_path, outcool]
            log.info('Generate contact matrix in cool at {0} ...'.format(reslabel))
            subprocess.check_call(' '.join(command), shell=True)
        else:
            data_folder = os.path.join(os.path.split(pairLiftOver.__file__)[0], 'data')
            juicer_folder = os.path.join(data_folder, 'juicer_tools_1.11.09_jcuda.0.8.jar')
            outhic = os.path.join(outfolder, '{0}.hic'.format(out_pre))
            command = ['java', '-jar', juicer_folder, 'pre', out_path, outhic, out_chroms]
            log.info('Generate contact matrices using juicer ...')
            subprocess.check_call(' '.join(command), shell=True)

        os.remove(out_path)
        os.remove(out_path+'.px2')

    log.info('Done')