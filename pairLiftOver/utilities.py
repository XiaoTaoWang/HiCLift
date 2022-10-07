import subprocess, sys, os, io, logging, cooler, pairLiftOver
from pairLiftOver.liftover import LiftOver
from pairLiftOver.io import open_pairs, _pixel_to_reads, _pairs_write

log = logging.getLogger(__name__)

def get_chrom_order(chroms_file):
    """
    Produce an "enumeration" of chromosomes based on the list
    of chromosomes
    """
    chrom_enum = dict()
    i = 1
    with open(chroms_file, 'r') as f:
        for line in f:
            chrom = line.strip().split()[0]
            chrom_enum[chrom] = i
            i += 1

    return chrom_enum

def make_standard_pairsheader(assembly, chromsizes, columns, shape="upper triangle"):

    header = []
    header.append("## pairs format v1.0.0")
    header.append("#shape: {}".format(shape))

    header.append(
        "#genome_assembly: {}".format(assembly if assembly is not None else "unknown")
    )

    if chromsizes is not None:
        try:
            chromsizes = chromsizes.items()
        except AttributeError:
            pass
        for chrom, length in chromsizes:
            header.append("#chromsize: {} {}".format(chrom, length))
    
    SEP_COLS = " "
    header.append("#columns: " + SEP_COLS.join(columns))

    return header

def get_header(instream, comment_char='#'):
    '''Returns a header from the stream and the remainder of the stream
    with the actual data.

    Parameters
    ----------
    instream : a file object
        An input stream.

    comment_char : str
        The character prepended to header lines (use '@' when parsing sams, 
        '#' when parsing pairsams).

    Returns
    -------
    header : list
        The header lines, stripped of terminal spaces and newline characters.

    remainder_stream : stream/file-like object
        Stream with the remaining lines.
    
    '''
    header = []
    if not comment_char:
        raise ValueError('Please, provide a comment char!')
    comment_byte = comment_char.encode()
    # get peekable buffer for the instream
    inbuffer = instream.buffer
    current_peek = inbuffer.peek()
    while current_peek.startswith(comment_byte):
        # consuming a line from buffer guarantees
        # that the remainder of the buffer starts 
        # with the beginning of the line.
        line = inbuffer.readline()
        # append line to header, since it does start with header
        header.append(line.decode().strip())
        # peek into the remainder of the instream
        current_peek = inbuffer.peek()
    # apparently, next line does not start with the comment
    # return header and the instream, advanced to the beginning of the data
    return header, instream

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
    chain_file, resolution=500, nproc_in=8, nproc_out=8, tmpdir='/tmp', memory='4G', high_res=False):
    
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
    header = make_standard_pairsheader(
        assembly=out_assembly, chromsizes=chromsizes,
        columns=['readID', 'chrom1', 'pos1', 'chrom2', 'pos2', 'strand1', 'strand2'],
        shape='upper triangle'
    )
    if in_assembly != out_assembly:
        header.append('#pairLiftOver: coordinates transformed from {0}'.format(in_assembly))
    else:
        header.append('#pairLiftOver: pure data format conversion')

    outstream.writelines((l+'\n' for l in header))
    outstream.flush()
    
    chrom_index = get_chrom_order(out_chroms)
    if in_format in ['cooler', 'juicer']:
        body_stream = instream
    else:
        _, body_stream = get_header(instream)
    
    # sort command
    command = r'''/bin/bash -c 'export LC_COLLATE=C; export LANG=C; sort -k 2,2 -k 4,4 -k 3,3n -k 5,5n --stable {0} {1} -S {2} {3}'''.format(
        '--parallel={0}'.format(nproc_out), '--temporary-directory={0}'.format(tmpdir), '8G', '--compress-program=lz4c'
    )
    command += "'"

    if in_assembly != out_assembly:
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
        
        log.info('Converting, sorting, and compressing ...')
    else:
        lo = None
        mapping_table = None
        log.info('Dumping contact pairs from {0} ...'.format(in_path))

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
    
    if in_assembly != out_assembly:
        log.info('{0:,} / {1:,} pairs were uniquely mapped to the target genome'.format(mapped_count, total_count))

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
        if out_format == 'cool':
            if high_res:
                outcool = os.path.join(tmpdir, '{0}.{1}.cool'.format(out_pre, '1kb'))
                bin_label = ':'.join([out_chroms, str(1000)])
                log.info('Generate contact matrix using cooler at 2500000,1000000,500000,250000,100000,50000,25000,10000,5000,2000,1000 ...')
                command = ['cooler', 'cload', 'pairix', '--assembly', out_assembly, '--nproc', str(nproc_out),
                           '--max-split 12', bin_label, out_path, outcool]
                subprocess.check_call(' '.join(command), shell=True)
                outmcool = os.path.join(outfolder, '{0}.mcool'.format(out_pre))
                command = ['cooler', 'zoomify', '-p', str(nproc_out), '-r 1000,2000,5000,10000,25000,50000,100000,250000,500000,1000000,2500000',
                           '--balance', '-o', outmcool, outcool]
                subprocess.check_call(' '.join(command), shell=True)
            else:
                outcool = os.path.join(tmpdir, '{0}.{1}.cool'.format(out_pre, '5kb'))
                bin_label = ':'.join([out_chroms, str(5000)])
                log.info('Generate contact matrix using cooler at 2500000,1000000,500000,250000,100000,50000,25000,10000,5000 ...')
                command = ['cooler', 'cload', 'pairix', '--assembly', out_assembly, '--nproc', str(nproc_out),
                           '--max-split 12', bin_label, out_path, outcool]
                subprocess.check_call(' '.join(command), shell=True)
                outmcool = os.path.join(outfolder, '{0}.mcool'.format(out_pre))
                command = ['cooler', 'zoomify', '-p', str(nproc_out), '-r 5000,10000,25000,50000,100000,250000,500000,1000000,2500000',
                           '--balance', '-o', outmcool, outcool]
                subprocess.check_call(' '.join(command), shell=True)
            
            os.remove(outcool)
        else:
            data_folder = os.path.join(os.path.split(pairLiftOver.__file__)[0], 'data')
            juicer_folder = os.path.join(data_folder, 'juicer_tools_1.11.09_jcuda.0.8.jar')
            outhic = os.path.join(outfolder, '{0}.hic'.format(out_pre))
            if high_res:
                command = ['java', '-Xmx'+memory.lower(), '-jar', juicer_folder, 'pre',
                            '-r 2500000,1000000,500000,250000,100000,50000,25000,10000,5000,2000,1000',
                            out_path, outhic, out_chroms]
                log.info('Generate contact matrices using juicer at 2500000,1000000,500000,250000,100000,50000,25000,10000,5000,2000,1000 ...')
            else:
                command = ['java', '-Xmx'+memory.lower(), '-jar', juicer_folder, 'pre',
                            '-r 2500000,1000000,500000,250000,100000,50000,25000,10000,5000',
                            out_path, outhic, out_chroms]
                log.info('Generate contact matrices using juicer at 2500000,1000000,500000,250000,100000,50000,25000,10000,5000 ...')
            subprocess.check_call(' '.join(command), shell=True)

        os.remove(out_path)
        os.remove(out_path+'.px2')

    log.info('Done')