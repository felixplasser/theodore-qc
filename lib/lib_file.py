"""
General file manipulation classes.
"""

class htmlfile:
    """
    Basic routines for writing an html file.
    """
    def __init__(self, hfname):
        self.hfile = open(hfname, 'w')
        
    def pre(self, title):
        """
        Inizialize the html file.
        """
        self.hfile.write("<html>\n<head>\n<title>")
        self.hfile.write(title)
        self.hfile.write("</title>\n</head>\n<body>\n")
        
    def write(self, wstr):
        """
        Write any string wstr to the file.
        """
        self.hfile.write(wstr)
        
    def post(self):
        """
        Finish up and close the html file.
        """
        self.hfile.write("</body>\n</html>\n")
        self.hfile.close()        
        
class htmltable:
    """
    Routines for creating an html table.
    """
    def __init__(self, ncol=2):
        self.ncol = ncol
        self.icol = 0
        
        self.str = '<table><tr>'
        
    def add_el(self, el):
        """
        Add an element
        """
        if self.icol == self.ncol:
            self.str += '</tr><tr>\n'
            self.icol = 0
        
        self.str += '<td>%s</td>\n'%el
        
        self.icol += 1
        
    def ret_table(self):
        self.str += '</tr></table>'
        
        return self.str
        
    