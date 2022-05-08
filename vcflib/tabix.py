# Copyright (c) 2014-2021 Sentieon Inc. All rights reserved
import collections
import struct
import sys

from . import bgzf
from .compat import *

__all__ = ['Tabix']

class Header(object):
    __slots__ = ('magic', 'n_ref', 'format',
        'col_seq', 'col_beg', 'col_end', 'meta', 'skip', 'l_nm')
    def __init__(self, *args):
        for k,v in zip(self.__slots__, args):
            setattr(self, k, v)
    def __iter__(self):
        for k in self.__slots__:
            yield getattr(self, k)

class Tabix(object):
    SHIFTS = (14, 17, 20, 23, 26, 29)
    MAXBIN = ((1 << SHIFTS[-1]-SHIFTS[0]+3) - 1) // 7 + 1
    MAGIC = 0x01494254
    FMT_GENERIC, FMT_SAM, FMT_VCF, FMT_ZERO_BASED = 0, 1, 2, 0x10000

    def __init__(self, idxf, mode='r'):
        self.path = idxf
        self.mode = mode
        if mode[0:1] == 'r':
            self.load()
        elif mode[0:1] == 'w':
            self.init()
        else:
            raise IOError('Mode ' + mode + ' not supported')

    def load(self):
        self.indices = collections.OrderedDict()

        with bgzf.open(self.path, 'rb') as fp:
            s4 = struct.Struct('<L')
            s8 = struct.Struct('<Q')
            sh = struct.Struct('<9L')
            data = fp.read()
            off = 0
            h = Header(*sh.unpack_from(data, off)); off += sh.size
            if h.magic != self.MAGIC:
                raise RuntimeError('Not a tabix file')
            self.header = h
            names, l_nm = [], 0
            for i in xrange(h.n_ref):
                eos = data.find(b'\0', off)
                if eos < 0: break
                names.append(data[off:eos].decode())
                l_nm += eos + 1 - off
                off = eos + 1
            if h.l_nm != l_nm:
                raise RuntimeError('Header sequence name length mismatch')
            for i in xrange(h.n_ref):
                bins = {}
                n_bin, = s4.unpack_from(data, off); off += s4.size
                for _ in xrange(n_bin):
                    bin, = s4.unpack_from(data, off); off += s4.size
                    chunks = []
                    n_chunk, = s4.unpack_from(data, off); off += s4.size
                    for _ in xrange(n_chunk):
                        s, = s8.unpack_from(data, off); off += s8.size
                        e, = s8.unpack_from(data, off); off += s8.size
                        chunks.append((s, e))
                    bins[bin] = chunks
                intvs = []
                n_intv, = s4.unpack_from(data, off); off += s4.size
                for _ in xrange(n_intv):
                    o, = s8.unpack_from(data, off); off += s8.size
                    intvs.append(o)
                if n_intv == 0:
                    intvs.append(0)
                self.indices[names[i]] = (bins, intvs)

    def save(self):
        if self.header is None:
            return
        self.add(None, 0, 0, self.end)
        h = self.header
        h.n_ref = len(self.indices)
        nms = b''.join(c.encode()+b'\0' for c,_ in iteritems(self.indices))
        h.l_nm = len(nms)
        with bgzf.open(self.path, 'wb') as fp:
            s4 = struct.Struct('<L')
            s8 = struct.Struct('<Q')
            sh = struct.Struct('<9L')
            fp.write(sh.pack(*h))
            fp.write(nms)
            for c, (bins, intvs) in iteritems(self.indices):
                fp.write(s4.pack(len(bins)))
                for bin in sorted(bins.keys()):
                    chunks = bins[bin]
                    fp.write(s4.pack(bin))
                    fp.write(s4.pack(len(chunks)))
                    for s,e in chunks:
                        fp.write(s8.pack(s))
                        fp.write(s8.pack(e))
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
        i = s >> self.SHIFTS[0]
        minoff = ci[1][i] if i < len(ci[1]) else ci[1][-1]
        for shift in reversed(self.SHIFTS):
            bo = ((1 << 29-shift) - 1) // 7
            bs = bo + (s >> shift)
            be = bo + (e-1 >> shift)
            be = min(be, self.MAXBIN-1)
            for bi in xrange(bs, be+1):
                if bi not in ci[0]:
                    continue
                for chunk in ci[0][bi]:
                    if chunk[1] > minoff:
                        ranges.append(chunk)
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
        h = Header(self.MAGIC, 0, self.FMT_VCF, 1, 2, 2, ord('#'), 0, 0)
        self.header = h
        self.indices = collections.OrderedDict()
        self.ci = None
        self.pos = 0
        self.end = 0

    def add(self, c, s, e, off):
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
            be = e-1 >> self.SHIFTS[0]
            if be >= len(intvs):
                intvs += [self.end] * (be+1 - len(intvs))
            bin = 0
            for shift in self.SHIFTS:
                bs, be = s >> shift, e-1 >> shift
                if bs == be:
                    bo = ((1 << 29-shift) - 1) // 7
                    bin = bo + bs
                    break
            chunks = bins.setdefault(bin,[])
            if chunks and chunks[-1][1] == self.end:
                chunks[-1] = (chunks[-1][0], off)
            else:
                chunks.append((self.end, off))
            self.pos = s
        self.end = off

    def optimize(self, ci):
        bins = ci[1]
        for shift in self.SHIFTS[:-1]:
            bo = ((1 << 29-shift) - 1) // 7
            for bin in sorted(bins.keys()):
                if bin < bo:
                    continue
                if bin > bo << 3:
                    break
                chunks = bins.get(bin)
                if chunks is None:
                    continue
                if len(chunks) == 0:
                    del bins[bin]
                    continue
                bs = chunks[0][0] >> 16
                be = chunks[-1][1] >> 16
                if be - bs < 65536:
                    del bins[bin]
                    bin = bin-1 >> 3
                    chunks += bins.get(bin,[])
                    bins[bin] = list(self.merge(chunks, 16))

# vim: ts=4 sw=4 expandtab
