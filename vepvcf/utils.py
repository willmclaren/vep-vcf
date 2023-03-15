#! /usr/bin/env python3

import sys
import gzip
import contextlib
import io


@contextlib.contextmanager
def smart_open(filename=None):
    if filename and filename != "-":
        if filename[-3:] == ".gz":
            fh = gzip.open(filename, "wb")
        else:
            fh = open(filename, "w")
    else:
        fh = sys.stdout

    try:
        yield fh
    finally:
        if fh is not sys.stdout:
            fh.close()


def smart_write(filehandle, string):
    if isinstance(filehandle, (io.RawIOBase, io.BufferedIOBase)):
        filehandle.write(string.encode())
    else:
        filehandle.write(string)
