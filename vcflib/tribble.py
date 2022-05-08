# Copyright (c) 2014-2021 Sentieon Inc. All rights reserved
import bisect
import collections
import heapq
import io
import os
import struct
import time

from .compat import *

__all__ = ['TribbleIndex']

class Header(object):
    __slots__ = ('magic', 'type', 'version', 'filename', 'filesize',
        'timestamp', 'md5', 'flags', 'properties')
    def __init__(self, *args):
        for k,v in zip(self.__slots__, args):
            setattr(self, k, v)

class LinearIndex(object):
    def __init__(self, chrom, off, width, density):
        self.chrom = chrom
        self.end = off
        self.width = width
        self.longest = 0
        self.count = 0
        self.blocks = []
        self.density = density

    def decode(self, data, off):
        chrom,_ = data[off:].split(b'\0',1)
        off += len(chrom) + 1
        self.chrom = chrom.decode()
        s = struct.Struct('<iiiii')
        self.width, bins, self.longest, _, self.count = s.unpack_from(data, off)
        off += s.size;
        s = struct.Struct('<'+'Q'*(bins+1))
        self.blocks = s.unpack_from(data, off)
        off += s.size
        return off

    def encode(self):
        data = bytearray()
        data.extend(self.chrom.encode())
        data.extend(b'\0')
        data.extend(struct.pack('<iiiii', self.width,
            len(self.blocks)-1, self.longest, 0, self.count))
        data.extend(struct.pack('<'+'Q'*len(self.blocks), *self.blocks))
        return data

    def query(self, s, e):
        s = max(s, self.longest) - self.longest
        i = s // self.width
        if i >= len(self.blocks):
            return []
        return [(self.blocks[i], self.blocks[-1])]

    def add(self, s, e, off):
        bin = s // self.width
        if bin >= len(self.blocks):
            self.blocks += [self.end] * (bin+1-len(self.blocks))
        self.longest = max(self.longest, e-s)
        self.count += 1
        self.end = off

    def done(self):
        self.blocks.append(self.end)
        self.optimize()

    def optimize(self):
        maxsize = max(self.blocks[i] - self.blocks[i-1]
            for i in xrange(1, len(self.blocks)))
        fullsize = self.blocks[-1] - self.blocks[0]
        scale = (self.density * fullsize) // (self.count * maxsize)
        if scale > 1:
            bins = (len(self.blocks)-1 + scale-1) // scale
            self.blocks = [self.blocks[i*scale] for i in xrange(bins)]
            self.width *= scale

class IntervalTree(object):
    def __init__(self):
        self.intvls = [(2**32-1, 2**32-1, None)]
        self.splits = None
        self.values = None

    def insert(self, s, e, d):
        self.intvls.append((s, e, d))

    def update(self):
        self.intvls.sort()
        self.splits = []
        self.values = []
        cur, h = 0, []
        for i,v in enumerate(self.intvls):
            while h and h[0][0] <= v[0]:
                self.splits.extend((cur, h[0][0]))
                self.values.append([j for _,j in h])
                cur = heapq.heappop(h)[0]
            if h and cur < v[0]:
                self.splits.extend((cur, v[0]))
                self.values.append([j for _,j in h])
            cur = v[0]
            heapq.heappush(h, (v[1], i))

    def query(self, s, e):
        r = set()
        i = bisect.bisect(self.splits, s) // 2
        while i < len(self.values):
            if e <= self.splits[i*2]:
                break
            r |= set(self.values[i])
            i += 1
        return (self.intvls[i][2] for i in r)

class IntervalTreeIndex(object):
    def __init__(self, chrom, off, density):
        self.chrom = chrom
        self.density = density
        self.tree = IntervalTree()
        self.curr = [0, 0, off, off, 0] # [sloc, eloc, soff, eoff, count]

    def decode(self, data, off):
        chrom,_ = data[off:].split(b'\0',1)
        off += len(chrom) + 1
        self.chrom = chrom.decode()
        s = struct.Struct('<i')
        nitvs, = s.unpack_from(data, off)
        off += s.size
        s = struct.Struct('<iiQi')
        for _ in xrange(nitvs):
            sloc, eloc, boff, size = s.unpack_from(data, off)
            off += s.size
            self.tree.insert(sloc-1, eloc, (boff, boff+size))
        self.tree.update()
        return off

    def encode(self):
        data = bytearray()
        data.extend(self.chrom.encode())
        data.extend(b'\0')
        s = struct.Struct('<i')
        data.extend(s.pack(len(self.tree.intvls)-1))
        s = struct.Struct('<iiQi')
        for sloc, eloc, b in self.tree.intvls:
            if b is None:
                continue
            data.extend(s.pack(sloc+1, eloc, b[0], b[1]-b[0]))
        return data

    def query(self, s, e):
        return TribbleIndex.merge(self.tree.query(s, e), 0)

    def add(self, s, e, off):
        if self.curr[4] == self.density:
            self.tree.insert(self.curr[0], self.curr[1], self.curr[2:4])
            self.curr[0] = s
            self.curr[2] = self.curr[3]
            self.curr[4] = 0
        self.curr[1] = e
        self.curr[3] = off
        self.curr[4] += 1

    def done(self):
        if self.curr[4] > 0:
            self.tree.insert(self.curr[0], self.curr[1], self.curr[2:4])
        self.tree.update()

class TribbleIndex(object):
    MAGIC = 0x58444954
    INDEX_TYPE_LINEAR = 1
    INDEX_TYPE_INTERVAL_TREE = 2
    VERSION = 3
    SEQUENCE_DICTIONARY_FLAG = 0x8000

    DEFAULT_INDEX_BIN_WIDTH = 8000
    GVCF_INDEX_BIN_WIDTH = 128000
    MAX_FEATURES_PER_BIN = 100
    MAX_FEATURES_PER_INTERVAL = 600

    def __init__(self, idxf, mode='r'):
        if not idxf.endswith('.idx'):
            raise ValueError('File name suffix is not .idx')
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

        with io.open(self.path, 'rb') as fp:
            data = fp.read()
            off = 0

            s = struct.Struct('<iii')
            magic, type, version = s.unpack_from(data, off)
            off += s.size
            if (magic != self.MAGIC or version != self.VERSION or
                type != self.INDEX_TYPE_LINEAR and
                type != self.INDEX_TYPE_INTERVAL_TREE):
                raise RuntimeError('Bad magic/type/version')

            file,_ = data[off:].split(b'\0',1)
            off += len(file) + 1
            file = file.decode()

            s = struct.Struct('<QQ')
            size,time = s.unpack_from(data, off)
            off += s.size

            md5,_ = data[off:].split(b'\0',1)
            off += len(md5) + 1

            s = struct.Struct('<ii')
            flags,nprop = s.unpack_from(data, off)
            off += s.size
            if nprop > 0:
                t = data[off:].split(b'\0', 2*nprop)
                if len(t) != 2*nprop+1:
                    raise RuntimeError('Incorrect property count')
                data = t.pop(); off = 0
                t = [s.decode() for s in t]
                properties = zip(*[iter(t)]*2)
            else:
                properties = []

            self.header = Header(magic, type, version, file,
                size, time, md5, flags, properties)

            s = struct.Struct('<i')
            nchrs, = s.unpack_from(data, off)
            off += s.size
            for _ in xrange(nchrs):
                if type ==  self.INDEX_TYPE_LINEAR:
                    ci = LinearIndex('', 0, 0, 0)
                elif type ==  self.INDEX_TYPE_INTERVAL_TREE:
                    ci = IntervalTreeIndex('', 0, 0)
                off = ci.decode(data, off)
                self.indices[ci.chrom] = ci

    def save(self):
        if self.header is None:
            return
        self.add(None, 0, 0, self.end)
        h = self.header
        h.filesize = self.end
        h.timestamp = int(time.time())
        with io.open(self.path, 'wb') as fp:
            fp.write(struct.pack('<iii', h.magic, h.type, h.version))
            fp.write(h.filename.encode()); fp.write(b'\0')
            fp.write(struct.pack('<QQ', h.filesize, h.timestamp))
            fp.write(h.md5); fp.write(b'\0')
            fp.write(struct.pack('<ii', h.flags, len(h.properties)))
            for k,v in h.properties:
                fp.write(k.encode()); fp.write(b'\0')
                fp.write(v.encode()); fp.write(b'\0')
            fp.write(struct.pack('<i', len(self.indices)))
            for k,ci in iteritems(self.indices):
                fp.write(ci.encode())
        self.header = None

    def query(self, c, s, e):
        ci = self.indices.get(c)
        if ci is None:
            return []
        return ci.query(s, e)

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
        type = int(os.getenv('VCF_INDEX_TYPE', '1'))
        if (type != self.INDEX_TYPE_LINEAR and
            type != self.INDEX_TYPE_INTERVAL_TREE):
            type = self.INDEX_TYPE_LINEAR
        self.header = Header(self.MAGIC, type, self.VERSION,
            self.path[:-4], 0, 0, b'', 0, [])
        self.indices = collections.OrderedDict()
        self.ci = None
        self.pos = 0
        self.end = 0

    def add(self, c, s, e, off):
        if self.ci and self.ci.chrom != c:
            self.ci.done()
            self.ci = None
        if self.ci is None and c is not None:
            if self.header.type == self.INDEX_TYPE_LINEAR:
                width = (self.path.endswith('.g.vcf.idx') and
                    self.GVCF_INDEX_BIN_WIDTH or self.DEFAULT_INDEX_BIN_WIDTH)
                density = self.MAX_FEATURES_PER_BIN
                self.ci = LinearIndex(c, self.end, width, density)
            elif self.header.type == self.INDEX_TYPE_INTERVAL_TREE:
                density = self.MAX_FEATURES_PER_INTERVAL
                self.ci = IntervalTreeIndex(c, self.end, density)
            self.indices[c] = self.ci
            self.pos = 0
        if self.ci:
            assert self.ci.chrom == c and s >= self.pos
            self.ci.add(s, e, off)
            self.pos = s
        self.end = off

# vim: ts=4 sw=4 expandtab
