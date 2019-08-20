from __future__ import print_function, division

from . import error_handler

"""
General file manipulation classes.
"""
class wfile:
    """
    Basic routines for writing a file.
    """
    def __init__(self, fname):
        self.f = open(fname, 'w')
        self.name = fname

    def __str__(self):
        return self.name

    def pre(self, title):
        """
        Initialization.
        """
        pass

    def write(self, wstr):
        """
        Write any string wstr to the file.
        """
        self.f.write(wstr)

    def post(self, lvprt=0):
        """
        Close and write the file.
        """
        self.post_extra()

        self.f.close()

        if lvprt >= 1:
            print("  File %s written."%self.name)

    def post_extra(self):
        """
        Any specific routines for closing the file.
        """
        pass

class wtable:
    """
    Virtual class with routines for creating a general table.
    """
    def __init__(self, ncol=2):
        self.ncol = ncol
        self.icol = 0

        self.str = self.init_extra()

    def init_extra(self):
        return ''

    def add_el(self, el):
        """
        Add an element
        """
        if self.icol == self.ncol:
            self.str += self.new_row()
            self.icol = 0

        self.str += self.new_el(el)

        self.icol += 1

        # Return info if the last column is reached
        return self.icol==self.ncol

    def add_row(self, row_list):
        """
        Add a row at once.
        """
        assert (len(row_list)==self.ncol)

        for el in row_list:
            self.add_el(el)

    def new_row(self):
        raise error_handler.PureVirtualError()

    def new_el(self, el):
        raise error_handler.PureVirtualError()

    def ret_table(self):
        self.str += self.close_table()

        return self.str

    def close_table(self):
        raise error_handler.PureVirtualError()

class asciitable(wtable):
    """
    Routines for creating a simple ASCII table.
    """
    def new_row(self):
        return '\n'

    def new_el(self, el):
        return '% .6f'%el

    def close_table(self):
        return '\n'

class htmlfile(wfile):
    """
    Basic routines for writing an html file.
    """
    def pre(self, title):
        """
        Inizialize the html file.
        """
        self.f.write("<html>\n<head>\n<title>")
        self.f.write(title)
        self.f.write("</title>\n</head>\n<body>\n")

    def post_extra(self):
        """
        Finish up and close the html file.
        """
        self.f.write("</body>\n</html>\n")

class htmltable(wtable):
    """
    Routines for creating an html table.
    """
    def init_extra(self):
        return '<table><tr>'

    def new_row(self):
        return '</tr><tr>\n'

    def new_el(self, el):
        return '<td>%s</td>\n'%el

    def close_table(self):
        return '</tr></table>\n'

class latexfile(wfile):
    """
    Write a file that can be interpreted by LaTeX.
    """
    def pre(self, title=None, graphicx=False, docclass="[a4paper]{article}"):
        self.f.write("\\documentclass%s\n"%docclass)
        if graphicx:
            self.f.write("\\usepackage{graphicx}\n")
            self.f.write("\\newcommand{\\incMO}{\\includegraphics[trim = 1.00cm 1.00cm 1.00cm 1.00cm, clip=true,width=6.00 cm]}\n\n")
            self.f.write("\\newcommand{\\incplot}{\\includegraphics[width=6.00 cm]}\n\n")
        self.f.write("\\begin{document}\n")
        if not title==None: self.f.write("%s\n"%title)

    def post_extra(self):
        self.f.write("\\end{document}\n")

class latextable(wtable):
    """
    Creating a LaTeX table.
    """
    def init_extra(self):
        ret_str  = "\\begin{table}\n"
        ret_str += "\\caption{(Caption)}\n"
        ret_str += "\\begin{tabular}{l%s}\n"%((self.ncol-1)*'r')

        return ret_str

    def new_row(self):
        return "\\\\\n"

    def new_el(self, el):
        ret_str = str(el)
        if not self.icol == self.ncol-1:
            ret_str += ' & '
        return ret_str

    def close_table(self):
        return "\n\\end{tabular}\n\\end{table}\n\n"

class latextabular(wtable):
    """
    Creating a LaTeX table.
    """
    def init_extra(self):
        ret_str = "\\begin{tabular}{l%s}\n"%((self.ncol-1)*'r')

        return ret_str

    def new_row(self):
        return "\\\\\n"

    def new_el(self, el):
        ret_str = str(el)
        if not self.icol == self.ncol-1:
            ret_str += ' & '
        return ret_str

    def close_table(self):
        return "\n\\end{tabular}\n"

class summ_file:
    """
    Class for analyzing the summary files.
    """
    def __init__(self, fname):
        self.ddict = {}
        self.state_labels = []

#       try:
        f = open(fname, 'r')
#        except IOError:
#            print ' WARNING: file %s not found!'%fname
#            self.header = ''
#            return

        self.header = next(f).replace('|','').split()
        next(f)

        while True:
            try:
                line = next(f).replace('|','')
            except StopIteration:
                break

            words = line.split()
            state_label = words[0]
            if state_label in self.ddict:
                errmsg  = "State %s already present.\n"%state_label
                errmsg += "  Please, do not combine tden_summ.txt files here."
                raise error_handler.MsgError(errmsg)
            self.ddict[state_label] = {}
            pdict = self.ddict[state_label]
            self.state_labels.append(state_label)

            pdict['state'] = state_label
            for i, prop in enumerate(self.header[1:]):
                try:
                    pdict[prop] = float(words[i+1])
                except ValueError:
                    pass

        f.close()

    def ret_header(self):
        return self.header

    def ret_ddict(self):
        return self.ddict

    def ret_state_labels(self):
        return self.state_labels
