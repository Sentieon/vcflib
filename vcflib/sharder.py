# Copyright (c) 2014-2021 Sentieon Inc. All rights reserved
from abc import ABCMeta, abstractmethod
import copy
import heapq
import multiprocessing
import operator
import signal

from .compat import *

__all__ = ['Sharder', 'Shardable', 'ShardResult']

class Shardable(object):
    __metaclass__ = ABCMeta

    @abstractmethod
    def __shard__(self, cse):
        return None

    @abstractmethod
    def __accum__(self, cse, data):
        return None

class ShardResult(object):
    __metaclass__ = ABCMeta

    @abstractmethod
    def __getdata__(self):
        return None

def shard(obj, cse):
    if isinstance(obj, Shardable):
        return obj.__shard__(cse)
    return obj

def getdata(obj):
    if isinstance(obj, ShardResult):
        return obj.__getdata__()
    return None

def accum(obj, data):
    if isinstance(obj, Shardable):
        obj.__accum__(data)

def apply(arg):
    idx, shd, func, args, kwargs = arg
    if not isinstance(shd, list) or isinstance(shd[0], basestring):
        shd = [shd]
    ret = []
    for cse in shd:
        av = [shard(o, cse) for o in args]
        kw = dict((k, shard(o, cse)) for k,o in iteritems(kwargs))
        try:
            rv = func(*av, **kw)
        except:
            import traceback
            traceback.print_exc()
            raise
        rv = copy.deepcopy(rv)
        av = [copy.deepcopy(getdata(o)) for o in av]
        kw = dict((k, copy.deepcopy(getdata(o))) for k,o in iteritems(kw))
        ret.append((rv, av, kw))
    return (idx, shd, ret)

class Sharder(object):
    def __init__(self, nproc=None):
        self.nproc = nproc

    @staticmethod
    def prep():
        signal.signal(signal.SIGINT, signal.SIG_IGN)

    def run(self, shards, map_fun, reduce_fun, *args, **kwargs):
        if reduce_fun is None:
            results = None
        elif isinstance(reduce_fun, list):
            reduce_fun, results = lambda s,x: s.append(x) or s, []
        else:
            results = reduce_fun(None)
        tasks = ((i, shd, map_fun, args, kwargs) for i,shd in enumerate(shards))
        q, idx = [], 0
        pool = multiprocessing.Pool(self.nproc, self.prep)
        it = pool.imap_unordered(apply, tasks, chunksize=1)
        while True:
            try:
                heapq.heappush(q, it.next(1))
                while q and q[0][0] == idx:
                    i, shd, ret = heapq.heappop(q)
                    for rv, av, kw in ret:
                        if reduce_fun is not None:
                            results = reduce_fun(results, rv)
                        for o,r in zip(args, av):
                            accum(o, r)
                        for k,o in iteritems(kwargs):
                            accum(o, kw.get(k))
                    idx += 1
            except StopIteration:
                assert not q
                break
            except multiprocessing.TimeoutError:
                continue
            except KeyboardInterrupt:
                pool.terminate()
                pool.join()
                raise
        pool.close()
        pool.join()
        return results

    @staticmethod
    def cut(intvs, step):
        shds = []
        size = 0
        for c,s,e in intvs:
            while s < e:
               n = min(e-s, step-size)
               shds.append((c, s, s+n))
               s += n
               size += n
               if size == step:
                   yield shds
                   shds = []
                   size = 0
        if shds:
           yield shds

if sys.version_info[0] == 2:
    # allow pickling of instance methods
    import copy_reg
    import types

    def _pickle_method(method):
        func_name = method.im_func.__name__
        obj = method.im_self
        cls = method.im_class
        return _unpickle_method, (func_name, obj, cls)

    def _unpickle_method(func_name, obj, cls):
        for cls in cls.mro():
            try:
                func = cls.__dict__[func_name]
            except KeyError:
                pass
            else:
                break
        if obj is None:
            return func
        return func.__get__(obj, cls)

    copy_reg.pickle(types.MethodType, _pickle_method, _unpickle_method)

# vim: ts=4 sw=4 expandtab
