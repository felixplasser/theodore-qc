#!/usr/bin/env python3
"""
Script for creating graphs from multiple directories, e.g. potential curves.
"""

from __future__ import print_function, division
from theodore import theo_header, lib_plot
        
def run_plot():
    infilen = 'graph.in'
    
    popt = lib_plot.write_plot_options(infilen)
    ropt = lib_plot.read_plot_options(infilen, False)

    if ropt.init == 0:
        copy = popt.ret_yn('Found %s. Use this file directly rather than performing an interactive input?'%infilen, True)
    else:
        copy = False    
    
    if copy:
        popt.copy(ropt)
    else:        
        popt.plot_input()
        popt.flush()
    
    popt.read_data()
    
    if popt['doplots']: popt.plot()
    if popt['dotxt']:   popt.txt_files()
    if popt['dognu']:   popt.gnu_inp()
    
if __name__ == '__main__':
    theo_header.print_header('Graph plotting')
    
    run_plot()
