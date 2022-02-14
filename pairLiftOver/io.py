import pipes, subprocess, struct, straw, sys, os, random

def generate_hic_blocks(chromsizes, step=10000000):

    chroms = sort_chromlabels(list(chromsizes))
    for i in range(len(chroms)):
        for j in range(len(chroms)):
            if i > j:
                continue
            c1 = chroms[i]
            c2 = chroms[j]
            for bs_1 in range(0, chromsizes[c1], step):
                for bs_2 in range(0, chromsizes[c2], step):
                    be_1 = min(chromsizes[c1], bs_1 + step - 1)
                    be_2 = min(chromsizes[c2], bs_2 + step - 1)
                    yield c1, bs_1, be_1, c2, bs_2, be_2

def find_digit_parts(chrname):

    collect = []
    for s in chrname[::-1]:
        if s.isdigit():
            collect.append(s)
        else:
            break
    
    if len(collect):
        digit_parts = int(''.join(collect[::-1]))
        return digit_parts
    else:
        return

def sort_chromlabels(chrnames):

    num_table = []
    nonnum_names = []
    for n in chrnames:
        tmp = find_digit_parts(n)
        if tmp is None:
            nonnum_names.append(n)
        else:
            num_table.append((tmp, n))

    num_table.sort()
    sorted_names = [s[1] for s in num_table]

    for s in ['M', 'Y', 'X']:
        for idx, n in enumerate(nonnum_names):
            if n.endswith(s):
                nonnum_names.pop(idx)
                nonnum_names.insert(0, n)
                break
    sorted_names = sorted_names + nonnum_names

    return sorted_names

def readcstr(f):
    buf = ""
    while True:
        b = f.read(1)
        b = b.decode('utf-8', 'backslashreplace')
        if b is None or b == '\0':
            return str(buf)
        else:
            buf = buf + b

def read_hic_header(hicfile):

    if not os.path.exists(hicfile):
        return None  # probably a cool URI

    req = open(hicfile, 'rb')
    magic_string = struct.unpack('<3s', req.read(3))[0]
    req.read(1)
    if (magic_string != b"HIC"):
        return None  # this is not a valid .hic file

    info = {}
    version = struct.unpack('<i', req.read(4))[0]
    info['version'] = str(version)

    masterindex = struct.unpack('<q', req.read(8))[0]
    info['Master index'] = str(masterindex)

    genome = ""
    c = req.read(1).decode("utf-8")
    while (c != '\0'):
        genome += c
        c = req.read(1).decode("utf-8")
    info['Genome ID'] = str(genome)

    nattributes = struct.unpack('<i', req.read(4))[0]
    attrs = {}
    for i in range(nattributes):
        key = readcstr(req)
        value = readcstr(req)
        attrs[key] = value
    info['Attributes'] = attrs

    nChrs = struct.unpack('<i', req.read(4))[0]
    chromsizes = {}
    for i in range(nChrs):
        name = readcstr(req)
        length = struct.unpack('<i', req.read(4))[0]
        if (name.upper() != 'ALL'):
            chromsizes[name] = length

    info['chromsizes'] = chromsizes

    info['Base pair-delimited resolutions'] = []
    nBpRes = struct.unpack('<i', req.read(4))[0]
    for i in range(nBpRes):
        res = struct.unpack('<i', req.read(4))[0]
        info['Base pair-delimited resolutions'].append(res)

    info['Fragment-delimited resolutions'] = []
    nFrag = struct.unpack('<i', req.read(4))[0]
    for i in range(nFrag):
        res = struct.unpack('<i', req.read(4))[0]
        info['Fragment-delimited resolutions'].append(res)

    return info

def read_hic_file(hicfil):

    info = read_hic_header(hicfil)
    chromsizes = info['chromsizes']
    binsize = min(info['Base pair-delimited resolutions'])
    blocks = generate_hic_blocks(chromsizes)
    for c1, bs_1, be_1, c2, bs_2, be_2 in blocks:
        result = straw.straw('NONE', hicfil, '{0}:{1}:{2}'.format(c1, bs_1, be_1), '{0}:{1}:{2}'.format(c2, bs_2, be_2), 'BP', binsize)
        _c1 = 'chr'+c1.lstrip('chr') # assume every chromosome has the prefix "chr"
        _c2 = 'chr'+c2.lstrip('chr')
        for k in range(len(result[0])):
            s1 = result[0][k]
            e1 = min(s1 + binsize, chromsizes[c1])
            s2 = result[1][k]
            e2 = min(s2 + binsize, chromsizes[c2])
            yield _c1, s1, e1, _c2, s2, e2, int(result[2][k])

def open_pairs(path, mode, data_format='pairs', nproc=1):
    
    if data_format == 'cooler':
        pread = subprocess.Popen(['cooler', 'dump', '--join', path],
                              stdout = subprocess.PIPE, bufsize = -1)
        f = pread.stdout
        return f
    elif data_format == 'juicer':
        f = read_hic_file(path)
        return f
    else:
        if path.endswith('.gz'):
            if mode =='w':
                t = pipes.Template()
                t.append('bgzip -c -@ {}'.format(nproc), '--')
                f = t.open(path, 'w')
            elif mode == 'r':
                t = pipes.Template()
                t.append('bgzip -dc -@ {}'.format(nproc), '--')
                f = t.open(path, 'r')

            return f

        elif path.endswith('.lz4'):
            if mode == 'w':
                t = pipes.Template()
                t.append('lz4c -cz', '--')
                f = t.open(path, 'w')
            elif mode == 'r':
                t = pipes.Template()
                t.append('lz4c -cd', '--')
                f = t.open(path, 'r')
            return f

        else:
            return open(path, mode)


def has_correct_order(loci1, loci2, chrom_index):
    
    check = (chrom_index[loci1[0]], loci1[1]) <= (chrom_index[loci2[0]], loci2[1])

    return check

def _core(loci, mapping_table, lo, resolution):

    if mapping_table is None:
        hit = lo.convert_coordinate(loci[0], loci[1])
        if (hit is None) or (len(hit) != 1):
            return
        else:
            return (hit[0][0], hit[0][1])
    else:
        c, s = loci
        s = s // resolution * resolution
        value = mapping_table.get((c, s), None)
        return value

def _pixel_to_reads(outstream, line, chrom_index, mapping_table, lo, resolution, source, total_count, mapped_count):

    if source == 'cooler':
        parse = line.decode().rstrip().split()
        if not len(parse):
            return total_count, mapped_count
        c1_, s1_, e1_, c2_, s2_, e2_, v = parse
        c1_ = 'chr' + c1_.lstrip('chr')
        c2_ = 'chr' + c2_.lstrip('chr')
        s1_, s2_, e1_, e2_, v = int(s1_), int(s2_), int(e1_), int(e2_), int(v)
    else:
        c1_, s1_, e1_, c2_, s2_, e2_, v = line

    total_count += v
    '''
    # strategy 1, only consider the midpoint of each bin
    p1_ = (s1_ + e1_) // 2
    p2_ = (s2_ + e2_) // 2
    hit1 = _core((c1_, p1_), mapping_table, lo, resolution)
    hit2 = _core((c2_, p2_), mapping_table, lo, resolution)
    if (hit1 is None) or (hit2 is None):
        return total_count, mapped_count
    
    if (not hit1[0] in chrom_index) or (not hit2[0] in chrom_index):
        return total_count, mapped_count
    
    if not has_correct_order(hit1, hit2, chrom_index):
        hit1, hit2 = hit2, hit1
    
    mapped_count += v
    for i in range(v):
        cols = ['.', hit1[0], str(hit1[1]), hit2[0], str(hit2[1]), '.', '.']
        outstream.write('\t'.join(cols) + '\n')
    '''
    # strategy 2, liftover start and end coordinates of each bin
    bin_size = e1_ - s1_
    l1 = _core((c1_, s1_), mapping_table, lo, resolution)
    r1 = _core((c1_, e1_), mapping_table, lo, resolution)
    l2 = _core((c2_, s2_), mapping_table, lo, resolution)
    r2 = _core((c2_, e2_), mapping_table, lo, resolution)
    if (l1 is None) or (r1 is None) or (l2 is None) or (r2 is None):
        return total_count, mapped_count
    
    if (not l1[0] in chrom_index) or (not r1[0] in chrom_index) \
        or (not l2[0] in chrom_index) or (not r2[0] in chrom_index):
        return total_count, mapped_count
    
    if (r1[1]-l1[1]!=bin_size) or (r2[1]-l2[1]!=bin_size) or (l1[0]!=r1[0]) or (l2[0]!=r2[0]):
        return total_count, mapped_count
    
    mapped_count += v
    for i in range(v):
        p1 = random.randint(l1[1], r1[1])
        p2 = random.randint(l2[1], r2[1])
        hit1 = (l1[0], p1)
        hit2 = (l2[0], p2)
        if not has_correct_order(hit1, hit2, chrom_index):
            hit1, hit2 = hit2, hit1
        cols = ['.', hit1[0], str(hit1[1]), hit2[0], str(hit2[1]), '.', '.']
        outstream.write('\t'.join(cols) + '\n')

    return total_count, mapped_count

def _pairs_write(outstream, line, chrom_index, mapping_table, lo, resolution, source, total_count, mapped_count):

    parse = line.rstrip().split()
    if not len(parse):
        return total_count, mapped_count
    
    if source=='hic-pro':
        readID, c1_, p1_, strand1, c2_, p2_, strand2 = parse[:7]
    else:
        readID, c1_, p1_, c2_, p2_, strand1, strand2 = parse[:7]
    
    total_count += 1
    p1_, p2_ = int(p1_), int(p2_)
    hit1 = _core((c1_, p1_), mapping_table, lo, resolution)
    hit2 = _core((c2_, p2_), mapping_table, lo, resolution)
    if (hit1 is None) or (hit2 is None):
        return total_count, mapped_count
    
    if (not hit1[0] in chrom_index) or (not hit2[0] in chrom_index):
        return total_count, mapped_count
    
    if not has_correct_order(hit1, hit2, chrom_index):
        hit1, hit2 = hit2, hit1
        strand1, strand2 = strand2, strand1
    
    cols = [readID, hit1[0], str(hit1[1]), hit2[0], str(hit2[1]), strand1, strand2]
    outstream.write('\t'.join(cols) + '\n')
    mapped_count += 1

    return total_count, mapped_count