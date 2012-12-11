#!/usr/bin/env python

from sys import argv
from random import choice
from inspect import getfile, currentframe

def main():
    """outputs a random CGL board in PGM format and a text representation
    of the board for easy inspection.
    """
    if len(argv) < 4:
        file_name = getfile(currentframe())
        print "USAGE: python", file_name, "nrows ncols out_file"
        return 1

    nrows = int(argv[1])
    ncols = int(argv[2])
    out_f = open(argv[3], 'w')
    txt_f = open(argv[3]+'.txt', 'w')

    # write PGM file header
    out_f.write("P5\n%d %d\n255\n" % (nrows, ncols))

    for i in range(nrows):
        for j in range(ncols):
            num = choice([1,2,3])
            txt_f.write("%d " % num)
            out_f.write("%c" % num)
        txt_f.write('\n')

if __name__ == "__main__":
    main()
