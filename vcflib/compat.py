# Copyright (c) 2014-2021 Sentieon Inc. All rights reserved
import sys

if sys.version_info[0] == 2:
    from itertools import izip

    def iterkeys(d, **kw):
        return d.iterkeys(**kw)
    def iteritems(d, **kw):
        return d.iteritems(**kw)

if sys.version_info[0] == 3:
    cmp = lambda x, y: (x > y) - (x < y)
    basestring = str
    izip = zip
    from functools import reduce
    xrange = range

    def iterkeys(d, **kw):
        return iter(d.keys(**kw))
    def iteritems(d, **kw):
        return iter(d.items(**kw))
