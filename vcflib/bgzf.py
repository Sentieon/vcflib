# Copyright (c) 2014-2021 Sentieon Inc. All rights reserved
import io
import struct
import zlib

from .compat import *

__all__ = ['BGZFile', 'open']

def open(file, mode='rb'):
    return BGZFile(file, mode)

class BGZFile(io.IOBase):
    block_size = 65280
    max_block_size = 65536

    def __init__(self, file, mode=None):
        level = zlib.Z_DEFAULT_COMPRESSION
        if mode is not None:
            for c in mode:
                if '0' <= c <= '9':
                    level = ord(c) - ord('0')
            mode = ''.join(c for c in mode if not '0' <= c <= '9')
        self.level = level

        self.name = None
        self.file = None

        if isinstance(file, basestring):
            self.name = file
            mode = mode and 'b' not in mode and mode + 'b' or mode or 'rb'
            self.file = io.open(file, mode)
        elif file is not None:
            self.file = file
            self.name = getattr(file, 'name', None)
            mode = mode or getattr(file, 'mode', 'rb')
        else:
            raise ValueError('File object cannot be None')

        if mode[0:1] == 'r':
            self.mode = 0
            self.block = None
        elif mode[0:1] == 'w':
            self.mode = 1
            self.zobj = None
            self.block = b''
        else:
            raise IOError('Mode ' + mode + ' not supported')

        self.pos = 0
        self.blkoff = 0
        self.blkend = 0

    def seekable(self):
        return True

    def readable(self):
        return self.mode == 0

    def writable(self):
        return self.mode == 1

    def seek(self, offset, whence=0):
        if whence == 1 and offset == 0:
            return self.pos
        if self.mode != 0 or whence != 0:
            raise ValueError('Invalid seek operation')
        if self.pos >> 16 != offset >> 16:
            self.block = None
        self.pos = offset

    def read(self, size=-1):
        self._checkReadable()
        data = bytearray()
        blk = self.pos >> 16
        off = self.pos & 65535
        while size != 0:
            if self.blkoff != blk or self.block is None:
                if not self._read_block(blk):
                    break
            end = len(self.block)
            if size > 0:
                end = min(end, off + size)
            data += self.block[off:end]
            size -= end - off
            if end < len(self.block):
                off = end
                break
            blk = self.blkend
            off = 0
        self.pos = blk << 16 | off
        return bytes(data)

    def read_until(self, delim, size=-1):
        self._checkReadable()
        data = bytearray()
        blk = self.pos >> 16
        off = self.pos & 65535
        while size != 0:
            if self.blkoff != blk or self.block is None:
                if not self._read_block(blk):
                    break
            end = len(self.block)
            if size > 0:
                end = min(end, off + size)
            eos = self.block.find(delim, off, end)
            if eos >= 0: end = eos+1
            data += self.block[off:end]
            size -= end - off
            if end < len(self.block):
                off = end
                break
            blk = self.blkend
            off = 0
            if eos >= 0:
                break
        self.pos = blk << 16 | off
        return bytes(data)

    def readline(self, size=-1):
        return self.read_until(b'\n', size)

    def write(self, data):
        self._checkWritable()
        ptr, size = 0, len(data)
        while ptr < size:
            n = min(self.block_size - len(self.block), size - ptr)
            self.block += data[ptr:ptr+n]
            if len(self.block) >= self.block_size:
                self._write_block()
            ptr += n
        self.pos = self.blkoff << 16 | len(self.block)

    def flush(self):
        if self.file is None:
            raise ValueError('flush of closed file')
        if self.mode == 1:
            self._write_block()
            self.pos = self.blkoff << 16
        self.file.flush()

    def close(self):
        if self.file is None:
            return
        if self.mode == 1:
            self.flush()
            self._write_eof()
        self.file.close()
        self.file = None

    @property
    def closed(self):
        return self.file is None

    def _read_block(self, offset):
        self.file.seek(offset)
        header = self.file.read(18)
        if len(header) == 0:
            return False
        if len(header) != 18:
            raise IOError('Incorrect header size')
        length = struct.unpack('<H', header[16:18])[0] + 1
        body = self.file.read(length-18-8)
        crc, size = struct.unpack('<LL', self.file.read(8))
        zobj = zlib.decompressobj(-zlib.MAX_WBITS)
        self.block = zobj.decompress(body, self.max_block_size) + zobj.flush()
        if len(self.block) != size:
            raise IOError('Incorrect block size')
        self.blkoff = offset
        self.blkend = offset + length
        zobj = None
        return True

    def _write_block(self):
        zobj = zlib.compressobj(self.level, zlib.DEFLATED, -zlib.MAX_WBITS)
        body = zobj.compress(self.block) + zobj.flush(zlib.Z_FINISH)
        length = len(body)+18+8
        crc = zlib.crc32(self.block) & 0xffffffff
        self.file.write(b'\037\213\010\4\0\0\0\0\0\377\6\0\102\103\2\0')
        self.file.write(struct.pack('<H', length-1))
        self.file.write(body)
        self.file.write(struct.pack('<LL', crc, len(self.block)))
        self.block = b''
        self.blkoff += length
        self.blkend = self.blkoff
        zobj = None

    def _write_eof(self):
        self.file.write(b'\037\213\010\4\0\0\0\0\0\377\6\0\102\103\2\0')
        self.file.write(b'\033\0\3\0\0\0\0\0\0\0\0\0')

# vim: ts=4 sw=4 expandtab
