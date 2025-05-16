# Copyright (c) 2014-2025 Sentieon Inc. All rights reserved
import collections
import os
import struct
import sys

from . import bgzf
from .compat import *

__all__ = ['Tabix']

class Header(object):
    __slots__ = ('format', 'col_seq', 'col_beg', 'col_end', 'meta', 'skip')
    def __init__(self, *args):
        for k,v in zip(self.__slots__, args):
            setattr(self, k, v)
    def __iter__(self):
        for k in self.__slots__:
            yield getattr(self, k)

class Tabix(object):
    TBI_MAGIC = 0x01494254
    CSI_MAGIC = 0x01495343
    FMT_GENERIC, FMT_SAM, FMT_VCF, FMT_ZERO_BASED = 0, 1, 2, 0x10000

    def __init__(self, path, mode='r'):
        self.path = path
        self.mode = mode
        if mode[0:1] == 'r':
            self.load()
        elif mode[0:1] == 'w':
            self.init()
        else:
            raise IOError('Mode ' + mode + ' not supported')

    def load(self):
        self.indices = collections.OrderedDict()

        if os.path.exists(self.path + '.csi'):
            idxf = self.path + '.csi'
            magic = self.CSI_MAGIC
        else:
            idxf = self.path + '.tbi'
            magic = self.TBI_MAGIC

        with bgzf.open(idxf, 'rb') as fp:
            s4 = struct.Struct('<L')
            s8 = struct.Struct('<Q')
            sh = struct.Struct('<6L')
            data = fp.read(); off = 0
            self.magic, = s4.unpack_from(data, off); off += s4.size
            if self.magic != magic:
                raise RuntimeError('Not a tabix file')
            if self.magic == self.TBI_MAGIC:
                self.min_shift, self.depth = 14, 5
                n_ref, = s4.unpack_from(data, off); off += s4.size
                aux = sh.unpack_from(data, off); off += sh.size
                l_nm, = s4.unpack_from(data, off); off += s4.size
                names = data[off:off+l_nm].split(b'\0'); off += l_nm
            else:
                self.min_shift, = s4.unpack_from(data, off); off += s4.size
                self.depth, = s4.unpack_from(data, off); off += s4.size
                l_aux, = s4.unpack_from(data, off); off += s4.size
                if l_aux < sh.size + s4.size:
                    raise RuntimeError('Invalid header')
                aux = sh.unpack_from(data, off); off += sh.size
                l_nm, = s4.unpack_from(data, off); off += s4.size
                names = data[off:off+l_nm].split(b'\0'); off += l_nm
                off += l_aux - (sh.size + s4.size + l_nm)
                n_ref, = s4.unpack_from(data, off); off += s4.size
            if len(names) != n_ref+1 or len(names[-1]) != 0:
                raise RuntimeError('Header sequence name length mismatch')
            self.header = Header(*aux)
            for i in xrange(n_ref):
                bins = {}
                n_bin, = s4.unpack_from(data, off); off += s4.size
                for _ in xrange(n_bin):
                    bin, = s4.unpack_from(data, off); off += s4.size
                    if self.magic == self.TBI_MAGIC:
                        loffset = 0
                    else:
                        loffset, = s8.unpack_from(data, off); off += s8.size
                    chunks = []
                    n_chunk, = s4.unpack_from(data, off); off += s4.size
                    for _ in xrange(n_chunk):
                        s, = s8.unpack_from(data, off); off += s8.size
                        e, = s8.unpack_from(data, off); off += s8.size
                        chunks.append((s, e))
                    bins[bin] = (loffset, chunks)
                intvs = []
                if self.magic == self.TBI_MAGIC:
                    n_intv, = s4.unpack_from(data, off); off += s4.size
                    for _ in xrange(n_intv):
                        o, = s8.unpack_from(data, off); off += s8.size
                        intvs.append(o)
                    if n_intv == 0:
                        intvs.append(0)
                self.indices[names[i].decode()] = (bins, intvs)
            self.max_shift = self.min_shift + self.depth * 3

    def save(self):
        if self.header is None:
            return
        self.add(None, 0, 0, self.end)

        for ext in ('.tbi', '.csi'):
            f = self.path + ext
            if os.path.exists(f): os.remove(f)

        if self.magic == self.TBI_MAGIC:
            idxf = self.path + '.tbi'
        else:
            idxf = self.path + '.csi'

        with bgzf.open(idxf, 'wb') as fp:
            s4 = struct.Struct('<L')
            s8 = struct.Struct('<Q')
            sh = struct.Struct('<6L')
            nms = b''.join(c.encode()+b'\0' for c,_ in iteritems(self.indices))
            fp.write(s4.pack(self.magic))
            if self.magic == self.TBI_MAGIC:
                fp.write(s4.pack(len(self.indices)))
                fp.write(sh.pack(*self.header))
                fp.write(s4.pack(len(nms)))
                fp.write(nms)
            else:
                fp.write(s4.pack(self.min_shift))
                fp.write(s4.pack(self.depth))
                fp.write(s4.pack(sh.size + s4.size + len(nms)))
                fp.write(sh.pack(*self.header))
                fp.write(s4.pack(len(nms)))
                fp.write(nms)
                fp.write(s4.pack(len(self.indices)))
            for c, (bins, intvs) in iteritems(self.indices):
                fp.write(s4.pack(len(bins)))
                for bin in sorted(bins.keys()):
                    loffset, chunks = bins[bin]
                    fp.write(s4.pack(bin))
                    if self.magic != self.TBI_MAGIC:
                        fp.write(s8.pack(loffset))
                    fp.write(s4.pack(len(chunks)))
                    for s,e in chunks:
                        fp.write(s8.pack(s))
                        fp.write(s8.pack(e))
                if self.magic == self.TBI_MAGIC:
                    fp.write(s4.pack(len(intvs)))
                    for o in intvs:
                        fp.write(s8.pack(o))
        self.header = None

    def query(self, c, s, e):
        ranges = []
        ci = self.indices.get(c)
        if ci is None:
            return ranges
        s = max(s, 0)
        i = s >> self.min_shift
        minoff = ci[1][min(i,len(ci[1])-1)] if ci[1] else 0
        for shift in range(self.max_shift, self.min_shift-3, -3):
            bo = ((1 << self.max_shift - shift) - 1) // 7
            bs = bo + (s >> shift)
            be = bo + (e-1 >> shift)
            if not ci[1]:
                for bi in xrange(bs, bo-1, -1):
                    b = ci[0].get(bi)
                    if b is not None:
                        minoff = max(minoff, b[0])
                        break
            for bi in xrange(bs, be+1):
                b = ci[0].get(bi)
                if b is not None:
                    ranges.extend(b[1])
        if minoff > 0:
            ranges = [(max(s,minoff), e) for s,e in ranges if e > minoff]
        return self.merge(ranges, 16)

    @staticmethod
    def merge(ranges, shift):
        p = None
        for r in sorted(ranges):
            if p is None:
                p = r
            elif r[0] >> shift > p[1] >> shift:
                yield p
                p = r
            else:
                p = (p[0], max(p[1],r[1]))
        if p is not None:
            yield p

    def init(self):
        self.magic = self.TBI_MAGIC
        self.min_shift = 14
        self.depth = 5
        type = list(map(int, os.getenv('VCF_INDEX_TYPE', '1').split(':')))
        if len(type) > 0 and type[0] == 2:
            self.magic = self.CSI_MAGIC
            if len(type) > 1:
                self.min_shift = type[1]
            if len(type) > 2:
                self.depth = type[2]
        self.max_shift = self.min_shift + self.depth * 3
        self.header = Header(self.FMT_VCF, 1, 2, 2, ord('#'), 0)
        self.indices = collections.OrderedDict()
        self.ci = None
        self.pos = 0
        self.end = 0

    def add(self, c, s, e, off):
        if c is None and s > 0:
            # s is the max contig length
            shift = self.min_shift
            limit = 1 << shift
            while s > limit:
                limit <<= 1
                shift += 1
            if shift >= 32:
                raise RuntimeError('Some contigs are too long')
            if shift > self.min_shift + self.depth * 3:
                self.magic = self.CSI_MAGIC;
                self.depth = (shift - self.min_shift + 2) // 3;
                self.max_shift = self.min_shift + self.depth * 3
        if self.ci and self.ci[0] != c:
            self.optimize(self.ci)
            self.ci = None
        if self.ci is None and c is not None:
            self.ci = (c, {}, [])
            self.indices[c] = self.ci[1:]
            self.pos = 0
        if self.ci:
            chrom, bins, intvs = self.ci
            assert chrom == c and s >= self.pos
            be = e-1 >> self.min_shift
            if be >= len(intvs):
                intvs += [self.end] * (be+1 - len(intvs))
            bin = 0
            for shift in range(self.min_shift, self.max_shift+3, 3):
                bs, be = s >> shift, e-1 >> shift
                if bs == be:
                    bo = ((1 << self.max_shift - shift) - 1) // 7
                    bin = bo + bs
                    break
            b = bins.setdefault(bin,[])
            if not b: b.extend((0, []))
            chunks = b[1]
            if chunks and chunks[-1][1] == self.end:
                chunks[-1] = (chunks[-1][0], off)
            else:
                chunks.append((self.end, off))
            self.pos = s
        self.end = off

    def optimize(self, ci):
        bins = ci[1]
        for shift in range(self.min_shift, self.max_shift+3, 3):
            bo = ((1 << self.max_shift - shift) - 1) // 7
            for bin in sorted(bins.keys()):
                if bin < bo:
                    continue
                if bin > bo << 3:
                    break
                b = bins.get(bin)
                if b is None:
                    continue
                chunks = b[1]
                if len(chunks) == 0:
                    del bins[bin]
                    continue
                bs = chunks[0][0] >> 16
                be = chunks[-1][1] >> 16
                if be - bs < 65536 and bo > 0:
                    del bins[bin]
                    bin = bin-1 >> 3
                    b = bins.setdefault(bin,[])
                    if not b: b.extend((0, []))
                    b[1] = list(self.merge(chunks + b[1], 16))
                elif ci[2]:
                    intv = (bin - bo) << (shift - self.min_shift)
                    intv = min(intv, len(ci[2])-1)
                    b[0] = ci[2][intv]

# vim: ts=4 sw=4 expandtab
