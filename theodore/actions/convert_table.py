#!/usr/bin/env python3
"""
Convert the TheoDORE output data to a table in latex or html format.
"""

from __future__ import print_function, division
from .actions import Action
from colt.lazyimport import LazyImportCreator, LazyImporter


with LazyImportCreator() as importer:
    theo_header = importer.lazy_import_as('..theo_header', 'theo_header')
    input_options = importer.lazy_import_as('..input_options', 'input_options')
    lib_file = importer.lazy_import_as('..lib_file', 'lib_file')
    error_handler = importer.lazy_import_as('..error_handler', 'error_handler')



class write_table_options(input_options.write_options):
    def table_input(self):
        print("Convert the output from a TheoDORE run into a latex or html table.\n")
        
        self.read_str("Name of the file to analyze", "ana_file", "tden_summ.txt")
        
        rstr = self.ret_str("Properties of interest (separated by spaces).\n  Leave empty to print all")
        if rstr == '':
            self['prop_list'] = []
        else:
            self.write_list('prop_list', rstr.split(), lformat="'%s'")
        
        self.choose_list('Output format', 'output_format',
                         [  ('latex', 'LaTeX source file'),
                            ('html', 'HTML format')
                         ], 'latex')
        
        if self['output_format'] == 'latex':
            self.write_option('fname', "table.tex")
            self.read_yn('Write state labels as LaTeX formulas', 'lformula', False)
        elif self['output_format'] == 'html':
            self.write_option('fname', "table.html")
            self['lformula'] = False
        
        digs = self.ret_int('Number of decimal digits', 2)
        fformat = '%.' + str(digs) + 'f'
        self.write_option('fformat', fformat)
    
    def write_table(self):
        sfile = lib_file.summ_file(self['ana_file'])
        header = sfile.ret_header()
        ddict  = sfile.ret_ddict()
        state_labels = sfile.ret_state_labels()
        
        if self['prop_list'] == []:
            self['prop_list'] = header[1:]
        
        if self['output_format'] == 'html':
            wfile = lib_file.htmlfile
            wtable = lib_file.htmltable
        elif self['output_format'] == 'latex':
            wfile = lib_file.latexfile
            wtable = lib_file.latextable
        else:
            raise error_handler.ElseError(self['output_format'], 'output_format')
            
        wf = wfile(self['fname'])        
        wf.pre(title='TheoDORE data')
        
        wt = wtable(ncol = len(self['prop_list']) + 1)
        wt.add_row(['State'] + self['prop_list'])
        
        for state in state_labels:
            if not self['lformula']:
                wt.add_el(state)
            else:
                wt.add_el('$%s$'%(state.replace('(', '^').replace(')', '')))
                
            for prop in self['prop_list']:
                try:
                    wt.add_el(self['fformat']%ddict[state][prop])
                except KeyError:
                    wt.add_el('-')
        
        wf.write(wt.ret_table())
        wf.post(lvprt=1)

class read_table_options(input_options.read_options):
    def set_defaults(self):
        self['ana_file'] = 'tden_summ.txt'
        self['output_format'] = 'latex'
        self['lformula'] = False
        self['prop_list'] = []
        self['fname'] = None
        self['fformat'] = '%.2f'


class ConvertTable(Action):

    name = 'convert_table'

    _colt_description = 'Convert the output to latex/html table'

    _lazy_imports = LazyImporter({
            '..theo_header': 'theo_header',
            '..input_options': 'input_options',
            '..lib_file': 'lib_file',
            '..error_handler': 'error_handler',
    })

    def run():
        theo_header.print_header(title=__class__._colt_description)
        infilen = 'table.in'
        
        topt = write_table_options(infilen)    
        ropt = read_table_options(infilen, False)
        
        if ropt.init == 0:
            copy = topt.ret_yn('Found %s. Use this file directly rather than performing an interactive input?'%infilen, True)
        else:
            copy = False    
        
        if copy:
            topt.copy(ropt)
        else:
            topt.table_input()    
        
        topt.write_table()
        
        if not copy:
            topt.flush()
