# Copyright (c) 2014-2025 Sentieon Inc. All rights reserved
import collections
import fnmatch
import io
import itertools
import os
import re
import sys
import tempfile

from . import bgzf
from . import sharder
from . import tabix
from . import tribble
from .compat import *

__all__ = ['VCF', 'Variant']

class Variant(object):
    __slots__ = ('chrom', 'pos', 'id', 'ref', 'alt', 'qual', 'filter',
        'info', 'samples', 'end', 'line')
    def __init__(self, *args):
        for k,v in izip(self.__slots__, args):
            setattr(self, k, v)
    def __str__(self):
        return self.line

class VCF(sharder.Shardable):
    decoders = { 'Integer': int, 'Float': float, 'Flag': bool }
    encoders = { 'Integer': str, 'Float': str }
    kvpat = re.compile(r'(.*?)=(".*?"|.*?)(?:,|$)') 

    def __init__(self, path, mode='r'):
        self.path = path
        self.mode = 'b' not in mode and mode + 'b' or mode
        self.isGVCF = path.endswith('.g.vcf') or path.endswith('.g.vcf.gz')
        self.open(self.path, self.mode)
        if self.mode[0:1] == 'r':
            self.load_header()
            self.parse_header()
        else:
            self.headers = []

    def __getstate__(self):
        odict = self.__dict__.copy()
        if self.mode[0:1] == 'r':
            odict['fp'] = self.fp.tell()
            odict.pop('index')
        else:
            odict['fp'] = None
            odict['index'] = None
        return odict

    def __setstate__(self, ndict):
        path = ndict['path']
        mode = ndict['mode']
        if mode[0:1] == 'r':
            self.open(path, mode)
            self.fp.seek(ndict.pop('fp'))
        self.__dict__.update(ndict)

    def open(self, path, mode):
        self.fp, self.index = None, None
        if path.endswith('.gz'):
            self.fp = bgzf.open(path, mode)
            self.index = tabix.Tabix(path, mode)
        elif path == '-':
            if mode[0:1] == 'r':
                raise RuntimeError('Input vcf cannot be stdin')
            self.fp = os.fdopen(1, mode)
        else:
            self.fp = io.open(path, mode)
            self.index = tribble.TribbleIndex(path, mode)

    def load_header(self):
        self.headers = []
        offset = 0
        line = self.fp.readline().decode()
        while line.startswith('#'):
            self.headers.append(line.rstrip())
            offset = self.fp.tell()
            line = self.fp.readline().decode()
        self.fp.seek(offset)
        self.initoff = offset

    def copy_header(self, src, update=None, remove=None):
        hdrs = collections.OrderedDict()
        pat = re.compile(r'^##([^=]+)=(<ID=([^,]+).*>)?')
        for line in src.headers:
            m = pat.match(line)
            if m is None:
                fld, id = line, None
            else:
                fld, id = m.group(1), m.group(3)
            hdrs.setdefault(fld, collections.OrderedDict())[id] = line
        if update:
            for line in update:
                m = pat.match(line)
                if m is None:
                    continue
                fld, id = m.group(1), m.group(3)
                hdrs.setdefault(fld, collections.OrderedDict())[id] = line
        if remove:
            for line in remove:
                m = pat.match(line)
                if m is None:
                    continue
                fld, id = m.group(1), m.group(3)
                if fld not in hdrs:
                    continue
                if id is None:
                    if fnmatch.fnmatch(hdrs[fld].get(id,''), line):
                        hdrs[fld].pop(id)
                else:
                    for id in fnmatch.filter(hdrs[fld].keys(), id):
                        hdrs[fld].pop(id)
        self.headers = [l for _,v in hdrs.items() for _,l in v.items()]
        self.parse_header()

    def emit_header(self):
        for line in self.headers:
            if line.startswith('#CHROM'):
                continue
            self.fp.write(line.encode() + b'\n')
        cols = '#CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO'.split(',')
        if len(self.samples) > 0:
            cols.append('FORMAT')
            cols += self.samples
        line = '\t'.join(cols)
        self.fp.write(line.encode() + b'\n')
        if self.index is not None:
            maxlen = max(int(t.get('length',0))
                for c,t in iteritems(self.contigs))
            self.index.add(None, maxlen, 0, self.fp.tell())

    def parse_header(self):
        self.contigs = collections.OrderedDict()
        self.alts = {}
        self.filters = {}
        self.infos = {}
        self.formats = {}
        self.samples = []
        for line in self.headers:
            if line.startswith('##contig='):
                d = self.parse_line(line[2:])
                self.contigs[d['ID']] = d
            elif line.startswith('##ALT='):
                d = self.parse_line(line[2:])
                self.alts[d['ID']] = d
            elif line.startswith('##FILTER='):
                d = self.parse_line(line[2:])
                self.filters[d['ID']] = d
            elif line.startswith('##INFO='):
                d = self.parse_line(line[2:])
                self.infos[d['ID']] = d
            elif line.startswith('##FORMAT='):
                d = self.parse_line(line[2:])
                self.formats[d['ID']] = d
            elif line.startswith('#CHROM'):
                self.samples = line[1:].split('\t')[9:]

    @staticmethod
    def parse_kv(kv):
        kv = kv.split('=',1)
        return len(kv) == 2 and kv or [kv[0], True]

    @staticmethod
    def parse_line(line):
        s = line.index('<')
        e = line.index('>')
        return dict(VCF.kvpat.findall(line[s+1:e]))

    @staticmethod
    def parse_field(desc, kv):
        k,v = kv
        d = desc.get(k)
        if d is None:
            return (k, v)
        cvt = VCF.decoders.get(d['Type'], str)
        if v == '.':
            v = None
        elif d['Number'] != '0' and d['Number'] != '1':
            v = v.split(',')
            if all(x == '.' for x in v):
                v = None
            else:
                v = list(map(cvt, v))
        elif cvt:
            v = cvt(v)
        return (k,v)

    def parse_info(self, kv):
        return self.parse_field(self.infos, self.parse_kv(kv))

    def parse_sample(self, kv):
        return self.parse_field(self.formats, kv)

    def parse(self, line):
        vals = line.split('\t')
        if len(vals) < 8 or vals[7] == '.':
            info = {}
        else:
            info = dict(map(self.parse_info, vals[7].split(';')))
        if len(vals) < 9 or vals[8] == '.':
            fmts = []
        else:
            fmts = vals[8].split(':')
        samples = []
        for val in vals[9:]:
            s = dict(map(self.parse_sample, izip(fmts,val.split(':'))))
            samples.append(s)
        vals[1] = int(vals[1])-1
        vals[4] = vals[4].split(',') if vals[4] != '.' else []
        vals[5] = float(vals[5]) if vals[5] != '.' else None
        vals[6] = vals[6].split(';') if vals[6] != '.' else []
        end = info.get('END', vals[1] + len(vals[3]))
        vals = vals[:7] + [info, samples, end, line]
        return Variant(*vals)

    @staticmethod
    def genotypes(a, p):
        if len(a) == 1:
            return (a * p, )
        a, b = a[:-1], a[-1:]
        return [g + b*k for k in xrange(p+1) for g in VCF.genotypes(a, p-k)]

    @staticmethod
    def sort_field(desc, alt, kv):
        k,v = kv
        d = desc.get(k)
        if d is None:
            return kv
        if d['Number'] == '.' and len(v) == len(alt)+1 or d['Number'] == 'R':
            a = ['R'] + alt
            return (k, [e for _,e in sorted(izip(a,v))])
        if d['Number'] == 'A':
            a = alt
            return (k, [e for _,e in sorted(izip(a,v))])
        if d['Number'] == 'G':
            a = ['R'] + alt
            ploidy, glsize = 2, len(a)*(len(a)+1)/2
            if glsize == len(v):
                # diploid
                g = [sorted((y,x)) for i,x in enumerate(a) for y in a[:i+1]]
            elif glsize < len(v):
                # figure out ploidy len(v) = C(ploidy+len(a)-1,len(a)-1)
                while glsize < len(v):
                    glsize = glsize * (ploidy+len(a)) / (ploidy+1)
                    ploidy += 1
                assert glsize == len(v)
                g = VCF.genotypes(a, ploidy)
            else:
                # haploid
                ploidy -= 1
                glsize = len(a)
                assert glsize == len(v)
                g = a
            return (k, [e for _,e in sorted(izip(g,v))])
        return kv

    def sort_info(self, alt, kv):
        return self.sort_field(self.infos, alt, kv)

    def sort_sample(self, alt, kv):
        return self.sort_field(self.formats, alt, kv)

    @staticmethod
    def format_field(desc, kv):
        k,v = kv
        d = desc.get(k)
        if d is None:
            return (k,str(v))
        cvt = VCF.encoders.get(d['Type'], str)
        if v is None:
            return (k,'.')
        elif d['Number'] == '0':
            return (k,None)
        elif d['Number'] == '1':
            return (k,cvt(v))
        return (k,','.join(map(cvt, v)))

    def format_info(self, kv):
        return self.format_field(self.infos, kv)

    def format_sample(self, kv):
        return self.format_field(self.formats, kv)

    def format(self, v):
        flds = [v.chrom, str(v.pos+1), v.id, v.ref, ','.join(v.alt) or '.',
            v.qual is not None and ('%4.2f' % v.qual) or '.',
            ';'.join(v.filter) or '.']
        t = [self.format_info(kv) for kv in sorted(v.info.items())]
        f = lambda kv: kv[1] is None and kv[0] or '='.join(kv)
        flds.append(';'.join(map(f, t)) or '.')
        if len(self.samples) > 0:
            keys = set(itertools.chain(*[iterkeys(s) for s in v.samples]))
            keys.discard('GT'); keys = ['GT'] + sorted(keys)
            flds.append(':'.join(keys))
            for s in v.samples:
                t = [self.format_sample((k,s.get(k)))[1] for k in keys]
                flds.append(':'.join(t))
        t, v.line = v.line, '\t'.join(flds)
        return t

    def __iter__(self):
        self.fp.seek(self.initoff)
        return self

    def __next__(self):
        line = self.fp.readline().decode()
        while line.startswith('#'):
            line = self.fp.readline().decode()
        if len(line) == 0:
            raise StopIteration
        line = line.rstrip()
        try:
            v = self.parse(line)
        except ValueError as e:
            e.args += (line,)
            raise
        return v

    next = __next__

    def range(self, chrom, start=0, end=0x7fffffff):
        return VCFReader(self, chrom, start, end)

    def emit(self, v):
        if v.line is None:
            self.format(v)
        self.fp.write(v.line.encode() + b'\n')
        if self.index is not None:
            self.index.add(v.chrom, v.pos, v.end, self.fp.tell())

    def cmp_variants(self, v1, v2):
        v1_contig_idx = list(self.contigs.keys()).index(v1.chrom)
        v2_contig_idx = list(self.contigs.keys()).index(v2.chrom)
        if not cmp(v1_contig_idx, v2_contig_idx):
            return cmp(v1_contig_idx, v2_contig_idx)
        else:
            return cmp(v1.pos, v2.pos)

    def close(self):
        if self.fp is not None:
            self.fp.close()
        if self.index is not None and self.mode[0:1] == 'w':
            self.index.save()

    def __del__(self):
        self.close()

    def __shard__(self, cse):
        chrom, start, end = cse
        if self.mode[0:1] == 'r':
            return VCFReader(self, chrom, start, end)
        else:
            return VCFWriter(self, chrom, start, end)

    def __accum__(self, tmpf):
        if tmpf is None:
            return
        tfp = io.open(tmpf, 'rb')
        for line in tfp:
            flds = line.decode().rstrip().split('\t')
            chrom = flds[0]
            pos = int(flds[1])-1
            end = pos + len(flds[3])
            if len(flds) > 7:
                info = flds[7]
                i = info.find('END=')
                while i > 0 and info[i-1] != ';':
                    i = info.find('END=', i+4)
                if i >= 0:
                    end = int(info[i+4:].split(';',1)[0])
            self.fp.write(line)
            if self.index is not None:
                self.index.add(chrom, pos, end, self.fp.tell())
        tfp.close()
        os.unlink(tmpf)

class VCFReader(object):
    def __init__(self, vcf, chrom, start, end):
        self.vcf = vcf
        self.chrom = chrom
        self.start = start
        self.end = end

    def __iter__(self):
        c, s, e = self.chrom, self.start, self.end
        self.ranges = list(self.vcf.index.query(c, s, e))
        return self

    def __getattr__(self, key):
        return getattr(self.vcf, key)

    def first(self):
        for i,r in enumerate(self.ranges):
            self.vcf.fp.seek(r[0])
            while True:
                try:
                    v = self.vcf.next()
                except ValueError:
                    o = self.vcf.fp.tell()
                    if o >= r[1]:
                        i += 1
                    else:
                        self.ranges[i] = (o, r[1])
                    self.ranges = self.ranges[i:]
                    raise
                if v.chrom != self.chrom or v.pos >= self.end:
                    raise StopIteration
                if v.end <= self.start:
                    if self.vcf.fp.tell() >= r[1]:
                        break
                    continue
                self.ranges = None
                return v
        raise StopIteration

    def __next__(self):
        if self.ranges is not None:
            return self.first()
        v = self.vcf.next()
        if v.chrom != self.chrom or v.pos >= self.end:
            raise StopIteration
        return v

    next = __next__

    def range(self, chrom, start=0, end=0x7fffffff):
        self.chrom = chrom
        self.start = start
        self.end = end
        return self

class VCFWriter(sharder.ShardResult):
    def __init__(self, vcf, chrom, start, end):
        self.vcf = vcf
        self.fp = tempfile.NamedTemporaryFile(mode='wb', suffix='.vcf',
            delete=False, dir=os.getenv('SENTIEON_TMPDIR'))
        self.path = self.fp.name
        self.chrom = chrom
        self.start = start
        self.end = end

    def __getattr__(self, key):
        return getattr(self.vcf, key)

    def emit_header(self):
        pass

    def emit(self, v):
        if v.chrom != self.chrom or v.pos >= self.end or v.end <= self.start:
            return
        if v.line is None:
            self.vcf.format(v)
        self.fp.write(v.line.encode() + b'\n')

    def close(self):
        self.fp.close()

    def __getdata__(self):
        self.close()
        return self.path

    def range(self, chrom, start=0, end=0x7fffffff):
        self.chrom = chrom
        self.start = start
        self.end = end
        return self

# vim: ts=4 sw=4 expandtab
