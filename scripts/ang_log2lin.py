#!/usr/bin/env python
# encoding: UTF8

import sys
import numpy as np

def main():
    if len(sys.argv) != 4:
        print "Usage: ang_log2lin.py min_logval log_step n_bins"
        exit(1)
    
    (min_logval, log_step, n_bins) = map(float, sys.argv[1:])
    
    e = min_logval
    log_hw = log_step/2.0
    n = 0
    ang_means = []
    ang_widths = []
    log_means = []
    while( n <= n_bins ):
        log_means.append(10.0**e)
        v_min = 10.0**(e-log_hw)
        v_max = 10.0**(e+log_hw)
        v_width = v_max-v_min
        v = v_min + v_width/2.0
        ang_means.append(v)
        ang_widths.append(v_width)
        e += log_step
        n += 1
    
    # print "log_means:"
    # print log_means
    
    print "ang_means:"
    print ang_means
    print "ang_widths:"
    print ang_widths
    


if __name__ == '__main__':
    main()
