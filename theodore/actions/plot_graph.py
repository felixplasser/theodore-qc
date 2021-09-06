#!/usr/bin/env python3
"""
Script for creating graphs from multiple directories, e.g. potential curves.
"""

from __future__ import print_function, division

from .actions import Action
from colt.lazyimport import LazyImportCreator, LazyImporter


with LazyImportCreator() as importer:
    theo_header = importer.lazy_import_as('..theo_header', 'theo_header')
    lib_plot = importer.lazy_import_as('..lib_plot', 'lib_plot')


class PlotGraph(Action):

    name = 'plot_graph'

    _colt_description = 'Graph plotting for potential curves etc.'

    _lazy_imports = LazyImporter({
            '..theo_header': 'theo_header',
            '..lib_plot': 'lib_plot',
    })

    def run():
        theo_header.print_header(title=__class__._colt_description)
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
