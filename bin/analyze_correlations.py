"""
Script for analyzing correlations in OmFrag files
Author: Sebastian Mai
"""

from __future__ import print_function, division
from .. import theo_header, input_options, lib_file, error_handler, lib_struc
import numpy
import os
import math


# =======================================================================
try:
    import scipy
    import scipy.cluster.hierarchy as sch
    from scipy.spatial.distance import squareform
except:
    print("scipy/scipy.cluster/scipy.spatial not installed - cannot continue")
    raise

# =======================================================================

try:
    from matplotlib import pyplot as plt
except:
    print("pylab/matplotlib not installed - plotting not possible")
    raise

# =======================================================================

try:
    import openbabel
    OPENBABEL=True
except ImportError:
    print(" *** Warning: python-openbabel not found! ***")
    print(" Using emulation program with limited capabilities ...")
    import OB_repl as openbabel
    OPENBABEL=False
    print(" Cannot draw bonds in LaTeX plots without python-openbabel.")

# =======================================================================
# =======================================================================
# =======================================================================

class AnCorr_options(input_options.write_options):
    """
    Set and store the options for plotting.
    """
    def __init__(self, *args, **kwargs):
        self.OmFrag = []
        self.Omega = {}
        self.maxOm = 0.
        self.numF = 0

        input_options.write_options.__init__(self, *args, **kwargs)

    # =======================================================================

    ## \param fname file with the data produced in a previous analyze_tden.py run
    def read_OmFrag(self, fname='OmFrag.txt'):
        """
        Read the OmFrag.txt file written by analyze_tden.py
        """
        Ofile = open(fname, 'r')
        data = Ofile.readlines()
        Ofile.close()

        self.numF = int(data[0].split()[0])

        matrix = []
        for line in data:

            words = line.split()
            if len(words)<=1:
                continue
            m = [ float(x) for x in words[2:] ]
            matrix.append(m)
            self.Omega[words[0]]=float(words[1])
            self.maxOm = max(self.maxOm, max(m))

        D=numpy.array(matrix)
        a,b=D.shape
        D=D.reshape(a,self.numF,self.numF)
        self.OmFrag=D
        print(('Number of fragments: %i' % self.numF))
        print(('Descriptor data (dimensions: %i obs x %i frags x %i frags):' % D.shape))
        numpy.set_printoptions(linewidth=130,precision=3,suppress=True,threshold=10000)
        #print self.OmFrag
        #print '\n\n'

        # NUMPY indexing:
        # D[0,:,:] is the first Om matrix
        # D[:,0,:] are the first lines from all Om matrices
        # D[:,:,0] are the first rows from all Om matrices
        # D[0,0,:] is the first line of the first Om matrix
        # D[0,:,0] is the first row from the first Om matrix
        # D[:,0,0] are the first matrix elements of all Om matrices

    # =======================================================================

    def AnCorr_input(self):

        # clustering options
        #dist_opts=['sqrt( 0.5*(1-COR) )']
        #ichoice = self.ret_choose_list('Measure used to convert correlation COR into distance?', dist_opts, 1)
        #self.write_option('dist_measure', ichoice)
        self.choose_list(
            'Measure used to convert correlation COR into distance?',
            'dist_measure',
        [
            ('theodore',   'sqrt( 0.5*(1-COR) )')
        ], 'theodore'
        )

        #link_opts=['single','complete','average','weighted']
        #ichoice = self.ret_choose_list('Hierarchical clustering linkage scheme?', link_opts, 3)
        #self.write_option('link_scheme', ichoice)
        self.choose_list(
            'Hierarchical clustering linkage scheme?',
            'link_scheme',
        [
            ('single',   'Simple linkage (nearest distance)'),
            ('complete', 'Complete linkage (farthest distance)'),
            ('average',  'Average linkage (average distance, UPGMA)'),
            ('weighted', 'Weighted-average linkage (WPGMA)')
        ], 'average'
        )

        #div_opts=['distance delta','distance','inconsistent']
        #ichoice = self.ret_choose_list('Method to find best clustering?', div_opts, 1)
        #self.write_option('division_scheme', ichoice)
        self.choose_list(
            'Method to find best clustering?',
            'division_scheme',
        [
            ('distance delta',   'Through largest distance step'),
            ('distance',         'At 0.7*maximum distance')
        ], 'distance delta'
        )

        # interactive plotting options
        self.read_yn('Show interactive plots?', 'doplot', True)


        # Latex diagram options
        self.read_yn('Create LaTeX diagram?', 'dolatex', True)
        if self.opt_dict['dolatex']:

            # get geometry
            self.read_str('File containing geometry?', 'coor_file', default='geom.xyz', autocomp=True)
            try:
                tstruc = lib_struc.structure()
                self['coor_format'] = tstruc.guess_file_type(self['coor_file'])
            except error_handler.MsgError:
                pass
            self.read_str('Format of coordinate file', 'coor_format', self['coor_format'])
            struc = lib_struc.structure()
            struc.read_file(file_path=self['coor_file'], file_type=self['coor_format'])
            self.coor=struc.ret_3xN_matrix()
            self.natom=self.coor.shape[0]
            self.elements=[]
            for i in range(self.natom):
                self.elements.append(struc.ret_symbol(i+1))

            # get the connectivity
            self.connectivity=[]
            if OPENBABEL:
                for curr_at in range(1,1+self.natom):
                    atom = struc.mol.GetAtom(curr_at)
                    for bonded in openbabel.OBAtomAtomIter(atom):
                        bind = bonded.GetIdx()
                        self.connectivity.append( (curr_at,bind) )

            # get at_lists
            self.read_str('File containing fragment definitions?', 'ana_file', default='dens_ana.in', autocomp=True)
            ioptions = input_options.tden_ana_options(self.opt_dict['ana_file'])
            self.at_lists=ioptions.opt_dict['at_lists']

            # latex options
            self.latex_opts={'alpha':90.,               # camera rotation angle 1 for molecule
                             'beta':45.,                # camera rotation angle 2 for molecule
                             'length':10.,              # edge length of matrix
                             'shift':0.3,               # gap between triangle matrices (gap=sqrt(2)*shift)
                             'yaxis':3.,                # height of dendrograms over matrix
                             'yaxis_extra':0.5,         # extra length for arrow tip
                             'cb_length':3,             # length of color bar
                             'cb_width':0.3,            # width of color bar
                             'shift_label':0.5,         # distance of "1.0" tick from matrix edge
                             'tick_length':0.15,        # length of ticks
                             'mat_dend_sep':1.,         # distance of combined analysis from matrix
                             'mol_scale':0.7,           # scale of molecule
                             'font_tick':r'\footnotesize',      # font commands for tick labels
                             'font_title':r'\large',            # font commands for titles
                             'colors':{                         # mapping from matplotlib colors to Latex colors
                                        'b': 'GR',
                                        'g': 'BD',
                                        'r': 'G',
                                        'c': 'RL',
                                        'm': 'T',
                                        'y': 'GL',
                                        'k': 'O',
                                        'w': 'W'
                                        }
                             }

        # debug plot
        #print self.__dict__


    # =======================================================================
    # =======================================================================
    # =======================================================================

    def compute(self):

        # ------------------- hole
        self.output_h={}

        # compute hole correlation matrix
        self.output_h['COR'] = self.calc_corr_matrix(self.OmFrag,userows=False,usecols=True)

        # compute distance matrix
        DIS = self.calc_dist_matrix(self.output_h['COR'])
        self.output_h['DIS'] = numpy.copy(DIS)

        # perform clustering
        Z=squareform( DIS, checks=False )
        Y=sch.linkage( Z, method=self.opt_dict['link_scheme'] )
        cc,C=sch.cophenet(Y,Z)
        I=sch.inconsistent(Y)
        self.output_h['dendro'] = Y
        self.output_h['cophco'] = cc
        self.output_h['incons'] = I

        # analyze clusters
        self.output_h['cluster'], self.output_h['maxd']=self.find_clustering(Y)



        # ------------------- electron
        self.output_e={}

        # compute electron correlation matrix
        self.output_e['COR'] = self.calc_corr_matrix(self.OmFrag,userows=True,usecols=False)

        # compute distance matrix
        DIS = self.calc_dist_matrix(self.output_e['COR'])
        self.output_e['DIS'] = numpy.copy(DIS)

        # perform clustering
        Z=squareform( DIS, checks=False )
        Y=sch.linkage( Z, method=self.opt_dict['link_scheme'] )
        cc,C=sch.cophenet(Y,Z)
        I=sch.inconsistent(Y)
        self.output_e['dendro'] = Y
        self.output_e['cophco'] = cc
        self.output_e['incons'] = I

        # analyze clusters
        self.output_e['cluster'], self.output_e['maxd']=self.find_clustering(Y)



        # ------------------- mixed
        self.output_m={}

        # compute mixed correlation matrix (maximum of hole and electron)
        self.output_m['COR'] = numpy.where(self.output_e['COR']>=self.output_h['COR'],
                                           self.output_e['COR'],
                                           self.output_h['COR'])

        # compute distance matrix
        DIS = self.calc_dist_matrix(self.output_m['COR'])
        self.output_m['DIS'] = numpy.copy(DIS)

        # perform clustering
        Z=squareform( DIS, checks=False )
        Y=sch.linkage( Z, method=self.opt_dict['link_scheme'] )
        cc,C=sch.cophenet(Y,Z)
        I=sch.inconsistent(Y)
        self.output_m['dendro'] = Y
        self.output_m['cophco'] = cc
        self.output_m['incons'] = I

        # analyze clusters
        self.output_m['cluster'], self.output_m['maxd']=self.find_clustering(Y)

        #print self.output


    # =======================================================================
    def calc_corr_matrix(self,D,userows=True,usecols=False):
        a,b,c=D.shape
        if not b==c:
            print(('Data is not shaped correctly: %i,%i,%i' % (a,b,c)))
            exit(1)
        V1=numpy.zeros(b)
        V2=numpy.zeros(b)
        # compute row covariance
        COV=numpy.zeros( (b,b) )
        if userows:
            for M in range(b):
                for N in range(b):
                    for I in range(a):
                        COV[M,N]+=numpy.dot( D[I,:,M],D[I,:,N] )
                    COV[M,N]/=a
                    V1=0.
                    V2=0.
                    for I in range(a):
                        V1+=D[I,:,M]
                        V2+=D[I,:,N]
                    V1/=a
                    V2/=a
                    COV[M,N]-=numpy.dot(V1,V2)
        # compute column covariance
        COV2=numpy.zeros( (b,b) )
        if usecols:
            for M in range(b):
                for N in range(b):
                    for I in range(a):
                        COV2[M,N]+=numpy.dot( D[I,M,:],D[I,N,:] )
                    COV2[M,N]/=a
                    V1=0.
                    V2=0.
                    for I in range(a):
                        V1+=D[I,M,:]
                        V2+=D[I,N,:]
                    V1/=a
                    V2/=a
                    COV2[M,N]-=numpy.dot(V1,V2)
        # average both
        COV+=COV2
        # compute correlation matrix
        COR=numpy.zeros( (b,b) )
        for M in range(b):
            for N in range(b):
                COR[M,N]=COV[M,N]/numpy.sqrt(COV[M,M])/numpy.sqrt(COV[N,N])
        # sanitize data
        COR=numpy.where(COR < -1.,-1.,COR)
        COR=numpy.where(COR > +1.,+1.,COR)
        COR=numpy.nan_to_num(COR)
        numpy.fill_diagonal(COR,1.0)
        # return
        return COR

    # =======================================================================
    def calc_dist_matrix(self,COR):
        if self.opt_dict['dist_measure']=='theodore':
            DIS=numpy.sqrt(0.5*(1.-COR))
        # TODO: here some possible distance measures, but I do not know which are sensible
        #DIS=(1.-COR)**2
        #DIS=1./(1.+COR)-0.5
        #DIS=np.where( COR>0., (0.5*(np.cos(np.pi*COR)+1.))**1, 1.)
        #DIS=np.sqrt(1.-COR**2)
        #DIS=(DIS+DIS.T)/2.
        return DIS

    # =======================================================================
    def find_clustering(self,Y):
        # aliases
        maxclust=self.numF
        if self.opt_dict['division_scheme']=='distance delta':
            # calculate rankings
            ranking={}
            for i in range(2,maxclust+1):
                if i==maxclust:
                    ranking[i]=0.
                else:
                    ranking[i]=(Y[-i+1,2]-Y[-i,2]) * (maxclust-i)/maxclust
            # choose best
            best=max(ranking, key=lambda i: ranking[i])
            if best==maxclust:
                max_d=min(Y[:,2])/2.
            else:
                max_d=(Y[1-best][2]+Y[-best][2])/2.
            #print max_d
            cluster=sch.fcluster(Y,best,criterion='maxclust')
        elif self.opt_dict['division_scheme']=='distance':
            max_d=0.7*max(Y[:,2])
            cluster=sch.fcluster(Y,max_d,criterion='distance')
        cluster=self.get_clusters(cluster)
        # return
        return cluster, max_d

    # =======================================================================
    def get_clusters(self,cluster_list):
        clusters={}
        for i,clu in enumerate(cluster_list):
            if not clu in clusters:
                clusters[clu]=set()
            clusters[clu].add(i+1)
        # flatten
        final=[]
        for clu in clusters:
            final.append(list(clusters[clu]))
        return final

    # =======================================================================
    def print_results(self,output,title):
        print('\n')
        print(('='*80))
        print(('Doing analysis for: ',title))
        print(('='*80))

        ## print correlation matrix
        #COR=output['COR']
        #print '\n&&& Correlation matrix for %s (dimensions: %i x %i ):' % (title,COR.shape[0],COR.shape[1])
        #print numpy.array_str(COR).replace('[',' ').replace(']',' ')

        ## print distance matrix
        #DIS=output['DIS']
        #print '\n&&& Distance matrix for %s (dimensions: %i x %i):' % (title,DIS.shape[0],DIS.shape[1])
        #print numpy.array_str(DIS).replace('[',' ').replace(']',' ')

        #print distances
        self.print_dendro(output['dendro'])

        # print cophenetic coefficient
        print(('\nThe cophenetic coefficient is %6.3f.' % output['cophco']))
        print('If this value is close to 1.0 then the dendrogram is \na good representation of the elementwise distances.')

        # print main output
        print('\nResult of the clustering procedure:')
        print('----------------------------------')
        print(('The %i original fragments from input' % (self.numF)))
        print('can be merged optimally in the following way:')
        for i,c in enumerate(output['cluster']):
            s='Cluster #% 3i  =  Fragments' % (i+1)
            for x in c:
                s+=' % 3i' % x
            print(s)
        print(('This means that instead of %i fragments, only %i fragments might suffice.' % (self.numF,len(output['cluster']))))

        # Plot
        if self.opt_dict['doplot']:
            #sch.dendrogram(output['dendro'])
            #dendro=output['dendro'].astype(numpy.float64)
            plt.figure(figsize=(20, 8))
            plt.title('Hierarchical Clustering Dendrogram -- Clustering of atoms for %s' % title)
            plt.xlabel('Fragment index (counting starts at one)')
            plt.ylabel('Distance')
            self.fancy_dendrogram(output['dendro'],
                leaf_rotation=90.,  # rotates the x axis labels
                leaf_font_size=8.,  # font size for the x axis labels
                max_d=output['maxd'],  # plot a horizontal cut-off line
                annotate_above=0.1,  # useful in small plots so annotations don't overlap
                labels=[ i+1 for i in range(self.numF)]
            )
            plt.show()

        ## Write data for latex figure generation
        #print '\nData for plot generation:'
        #print '-------------------------'
        #print '&&& Best clustering for %s:' % title
        #print output['cluster']
        #print '&&& Threshold for %s:' % title
        #print output['maxd']
        #print '&&& Dendrogram for %s:' % title
        #print sch.dendrogram(output['dendro'],
                             #no_plot=True,
                             #color_threshold=output['maxd'],
                             #labels=[ '%i'%(i+1) for i in range(self.numF)])


    # =======================================================================
    def print_dendro(self,dendro):
        n=len(dendro)+1
        print('\n iclust   merge[1]  merge[2]  distance  nfrags  nclust')
        print('-----------------------------------------------------')
        #for i in range(n):
        print(('% 3i-% 3i  elementary fragment  % 8.3f  % 6i  % 6i' % (1,n,0.,1,n)))
        for i in range(n-1):
            print(('% 7i   % 8i  % 8i  % 8.3f  % 6i  % 6i' % (n+i+1,
                                                             dendro[i][0]+1,
                                                             dendro[i][1]+1,
                                                             dendro[i][2],
                                                             dendro[i][3],
                                                             n-i-1)))

    # =======================================================================

    def fancy_dendrogram(self, *args, **kwargs):
        max_d = kwargs.pop('max_d', None)
        if max_d and 'color_threshold' not in kwargs:
            kwargs['color_threshold'] = max_d
        annotate_above = kwargs.pop('annotate_above', 0)

        ddata = sch.dendrogram(*args, **kwargs)

        if not kwargs.get('no_plot', False):
            for i, d, c in zip(ddata['icoord'], ddata['dcoord'], ddata['color_list']):
                x = 0.5 * sum(i[1:3])
                y = d[1]
                if y > annotate_above:
                    plt.plot(x, y, 'o', c=c)
                    plt.annotate("%.3g" % y, (x, y), xytext=(0, -5),
                                textcoords='offset points',
                                va='top', ha='center')
            if max_d:
                plt.axhline(y=max_d, c='k')
        return ddata

    # =======================================================================
    # =======================================================================
    # =======================================================================

    def create_latex(self,output_h, output_e, output_m):

        # LaTeX header
        header='''
\\documentclass[10pt]{article}

\\usepackage[utf8x]{inputenc}
\\usepackage[T1]{fontenc}
\\usepackage{tikz,pgfplots,sfmath}

\\usepackage[active,tightpage]{preview}
\\PreviewEnvironment{tikzpicture}
\\setlength\PreviewBorder{0.1pt}
\\setlength{\\textwidth}{8.25cm}

\\pgfplotsset{every tick label/.append style={/pgf/number format/1000 sep={}}}
\\pgfplotsset{compat=1.3}
\\usetikzlibrary{plotmarks,fit,positioning,calc}

\\usepackage{tikz-3dplot}

\\usepackage{xcolor}
\definecolor{W}{HTML}{FFFFFF}   %% white
\definecolor{K}{HTML}{000000}   %% black
\definecolor{GR}{HTML}{909090}  %% grey

\definecolor{BL}{HTML}{0099e6}  %% blue light
\definecolor{B}{HTML}{006699}   %% blue medium
\definecolor{BD}{HTML}{004466}  %% blue dark
\definecolor{BR}{HTML}{0080FF}  %% blue alternative

\definecolor{RW}{HTML}{f5d0df}  %% red very light
\definecolor{RL}{HTML}{d65454}  %% red light
\definecolor{R}{HTML}{893636}   %% red medium
\definecolor{RD}{HTML}{ad1737}  %% red dark

\definecolor{GW}{HTML}{d0f5df}  %% green mint
\definecolor{GL}{HTML}{10d05d}  %% green light
\definecolor{G}{HTML}{009933}   %% green medium
\definecolor{GD}{HTML}{004d1a}  %% green dark

\definecolor{Y}{HTML}{E0E040}   %% yellow
\definecolor{O}{HTML}{FFA500}   %% orange
\definecolor{TW}{HTML}{d9f1ff}  %% cyan light
\definecolor{T}{HTML}{4DC0FF}   %% cyan
\definecolor{OL}{HTML}{FFAA80}  %% orange light
\definecolor{P}{HTML}{DD88BB}   %% pink

\\begin{document}
  \\tdplotsetmaincoords{%f}{%f}
  \\begin{tikzpicture}[
    %%tdplot_main_coords,
    >=stealth, 
    every node/.style={circle, inner sep=0pt, outer sep=0pt, minimum size=0pt},
    thick,
    font=\sffamily,
    scale=1.0
    ]

%%  \\draw[ultra thick, -stealth, black!25] (0,0,0) -- (4,0,0) node[anchor=west] {X};
%%  \\draw[ultra thick, -stealth, black!25] (0,0,0) -- (0,4,0) node[anchor=west] {Y};
%%  \\draw[ultra thick, -stealth, black!25] (0,0,0) -- (0,0,4) node[anchor=west] {Z};

''' % (self.latex_opts['alpha'],self.latex_opts['beta'])

        footer=r'''
  \end{tikzpicture}
\end{document}
'''

        # fetch data
        geom=[]
        for i,el in enumerate(self.elements):
            geom.append( [el]+self.coor[i].tolist() )
        natom=len(geom)
        connectivity=self.connectivity
        nbonds=len(connectivity)
        frags=self.at_lists
        nfrags=len(frags)

        DIS_hole=output_h['COR']
        cluster_hole=output_h['cluster']
        nclust_hole=len(cluster_hole)
        maxd_hole=output_h['maxd']
        dendro_hole=sch.dendrogram(output_h['dendro'],
                                  no_plot=True,
                                  color_threshold=output_h['maxd'],
                                  labels=[ '%i'%(i+1) for i in range(self.numF)])

        DIS_electron=output_e['COR']
        cluster_electron=output_e['cluster']
        nclust_electron=len(cluster_electron)
        maxd_electron=output_e['maxd']
        dendro_electron=sch.dendrogram(output_e['dendro'],
                                  no_plot=True,
                                  color_threshold=output_e['maxd'],
                                  labels=[ '%i'%(i+1) for i in range(self.numF)])

        DIS=output_m['DIS']
        cluster=output_m['cluster']
        nclust=len(cluster)
        maxd_molecule=output_m['maxd']
        dendro_all=sch.dendrogram(output_m['dendro'],
                                  no_plot=True,
                                  color_threshold=output_m['maxd'],
                                  labels=[ '%i'%(i+1) for i in range(self.numF)])


        # fetch options
        coltrans        =self.latex_opts['colors']                              # colors for dendrogram branches
        length          =self.latex_opts['length']
        shift           =self.latex_opts['shift']
        yaxis           =self.latex_opts['yaxis']
        cb_length       =self.latex_opts['cb_length']
        cb_width        =self.latex_opts['cb_width']
        yaxis_extra     =self.latex_opts['yaxis_extra']
        shift_label     =self.latex_opts['shift_label']
        tick_length     =self.latex_opts['tick_length']
        mat_dend_sep    =self.latex_opts['mat_dend_sep']
        mol_scale       =self.latex_opts['mol_scale']
        font_tick       =self.latex_opts['font_tick']
        font_title      =self.latex_opts['font_title']
        alpha           =self.latex_opts['alpha']
        beta            =self.latex_opts['beta']

        colors_all=['K','R','RD','RL','RW','W','TW','T','BL','BD','K']          # color map for matrix




        # create string
        latex_string=header+'\n\n\n'


        # build stuff for hole
        # triangle
        l=length*(nfrags+1)/nfrags
        unit=length/nfrags
        latex_string+='  \draw[K,ultra thick] (0,%f) -- (%f,%f) -- (0,%f) -- (0,%f);\n' % (-shift, 
                                                                                            l,
                                                                                            -l-shift,
                                                                                            -l-shift,
                                                                                            -shift)
        # axis
        latex_string+='  \draw[K,ultra thick, -stealth] (0,%f) -- (%f,%f) node[anchor=south] {%s $r$};\n' % (-shift-unit,
                                                                                                            -yaxis-yaxis_extra,
                                                                                                            -shift-unit,
                                                                                                            font_tick)
        labels=[1.0,0.8,0.5,0.0]
        for i in labels:
            xco=-(yaxis*math.sqrt( 0.5*(1.-i) )+shift_label)
            latex_string+='  \draw[K,thick,rectangle] (%f,%f) -- (%f,%f) node[anchor=south] {%s %.1f};\n'     % (xco,
                                                                                                                -shift-unit-tick_length,
                                                                                                                xco,
                                                                                                                -shift-unit+tick_length,
                                                                                                                font_tick,
                                                                                                                i)
        latex_string+='  \\node[rotate=90,anchor=center,rectangle] at (%f,%f) {%s Hole analysis};\n' % (-yaxis-yaxis_extra,-shift-unit-l/2.,font_title)

        # colorbar
        k=len(colors_all)-1
        colors=colors_all
        for i in range(k):
            latex_string+='  \shade[left color=%s, right color=%s] (%f,0.25) rectangle (%f,%f);\n' % ( colors[i],
                                                                                                colors[i+1],
                                                                                                -0.5*shift_label-cb_length*(1.-float(i)/k),
                                                                                                -0.5*shift_label-cb_length*(1.-float(i+1)/k),
                                                                                                +cb_width+0.25)
        latex_string+='  \draw[K,ultra thick] (%f,0.25) rectangle (%f,%f);\n' % (-0.5*shift_label-cb_length,
                                                                            -0.5*shift_label,
                                                                            +cb_width+0.25)
        latex_string+='  \\node[anchor=south, rectangle] at (%f,%f) {Correlation $r$};' % (-0.5*shift_label-0.5*cb_length,0.25+cb_width+tick_length+0.4)
        labels=[1.0,0.5,0.0,-0.5,-1.0]
        for i in labels:
            xco=-0.5*shift_label-0.5*cb_length +0.5*cb_length*i
            latex_string+='  \draw[K,thick,rectangle] (%f,%f) -- (%f,%f) node[anchor=south,rectangle] {%s %.1f};\n' % (xco,+cb_width+0.25,xco,tick_length+cb_width+0.25,font_tick,i)

        # matrix
        for i in range(nfrags):
            x=dendro_hole['leaves'][i]
            for j in range(nfrags):
                if i>j:
                    continue
                y=dendro_hole['leaves'][j]
                value=DIS_hole[x][y]
                if value==1.:
                    value=0.9999
                col=self.colormap(value, colors_all)
                latex_string+='  \\fill[%s,draw=K,ultra thin] (%f,%f) rectangle (%f,%f);\n' % (col,
                                                                                    i*unit,
                                                                                    -shift-unit-j*unit,
                                                                                    (i+1)*unit,
                                                                                    -shift-unit-(j+1)*unit )

        # labels
        for i in range(nfrags):
            label=dendro_hole['ivl'][i]
            latex_string+='  \\node[anchor=center,inner sep=1pt] at (%f,%f) {\\scriptsize %s};\n' % (-0.5*shift_label,
                                                                                    -shift-1.5*unit-i*unit,
                                                                                    label)
            latex_string+='  \\node[anchor=center,inner sep=1pt] at (%f,%f) {\\scriptsize \\scalebox{0.8}[1]{%s}};\n' % ((0.5+i)*unit,
                                                                                    -shift-l-0.5*shift_label,
                                                                                    label)

        # dendrogram
        for i in range(len(dendro_hole['icoord'])):
            xco=[ -yaxis*x-shift_label  for x in dendro_hole['dcoord'][i] ]
            yco=[ -x/10.*unit-shift-unit for x in dendro_hole['icoord'][i] ]
            color=coltrans[ dendro_hole['color_list'][i]]
            latex_string+='  \\draw[%s,line width=3pt] (%f,%f) -- (%f,%f) -- (%f,%f) -- (%f,%f);\n' % (color,xco[0],yco[0],xco[1],yco[1],xco[2],yco[2],xco[3],yco[3])
        latex_string+='  \draw[K,ultra thick] (%f,%f) -- (%f,%f);\n' % (-yaxis*maxd_hole-shift_label,
                                                                        -shift-unit-tick_length,
                                                                        -yaxis*maxd_hole-shift_label,
                                                                        -shift-l+tick_length)








        # build stuff for electron
        # triangle and axis
        latex_string+='  \draw[K,ultra thick] (%f,0) -- (%f,%f) -- (%f,0) -- (%f,0);\n' % (shift, 
                                                                                        l+shift,
                                                                                        -l,
                                                                                        l+shift,
                                                                                        shift)
        latex_string+='  \draw[K,ultra thick, -stealth] (%f,0) -- (%f,%f) node[anchor=east,] {%s $r$};\n' % (+shift+unit,
                                                                                                            +shift+unit,
                                                                                                            +yaxis+yaxis_extra,
                                                                                                            font_tick)
        labels=[1.0,0.8,0.5,0.0]
        for i in labels:
            xco=+(yaxis*math.sqrt( 0.5*(1.-i) )+shift_label)
            latex_string+='  \draw[K,thick,rectangle] (%f,%f) -- (%f,%f) node[anchor=east] {%s %.1f};\n'     % (shift+unit+tick_length,
                                                                                                            xco,
                                                                                                            shift+unit-tick_length,
                                                                                                            xco,
                                                                                                            font_tick,
                                                                                                            i)
        latex_string+='  \\node[anchor=center,rectangle] at (%f,%f) {%s Electron analysis};\n' % (shift+unit+l/2.,yaxis+yaxis_extra,font_title)

        # matrix
        for i in range(nfrags):
            x=dendro_electron['leaves'][i]
            for j in range(nfrags):
                if i<j:
                    continue
                y=dendro_electron['leaves'][j]
                value=DIS_electron[x][y]
                if value==1.:
                    value=0.9999
                col=self.colormap(value,colors_all)
                latex_string+='  \\fill[%s,draw=K,ultra thin] (%f,%f) rectangle (%f,%f);\n' % (col,
                                                                                    shift+unit+i*unit,
                                                                                    -j*unit,
                                                                                    shift+unit+(i+1)*unit,
                                                                                    -(j+1)*unit, )

        # labels
        for i in range(nfrags):
            label=dendro_electron['ivl'][i]
            latex_string+='  \\node[anchor=center,inner sep=1pt] at (%f,%f) {\\scriptsize \\scalebox{0.8}[1]{%s}};\n' % (shift+1.5*unit+i*unit,
                                0.5*shift_label,
                                label)
            latex_string+='  \\node[anchor=center,inner sep=1pt] at (%f,%f) {\\scriptsize %s};\n' % (shift+l+0.5*shift_label,
                                -(0.5+i)*unit,
                                label)

        # dendrogram
        for i in range(len(dendro_electron['icoord'])):
            xco=[ x/10.*unit+shift+unit for x in dendro_electron['icoord'][i] ]
            yco=[ +yaxis*x+shift_label  for x in dendro_electron['dcoord'][i] ]
            color=coltrans[ dendro_electron['color_list'][i]]
            latex_string+='  \\draw[%s,line width=3pt] (%f,%f) -- (%f,%f) -- (%f,%f) -- (%f,%f);\n' % (color,xco[0],yco[0],xco[1],yco[1],xco[2],yco[2],xco[3],yco[3])
        latex_string+='  \draw[K,ultra thick] (%f,%f) -- (%f,%f);\n' % (+shift+unit+tick_length,
                                                                        +yaxis*maxd_electron+shift_label,
                                                                        +shift+l-tick_length,
                                                                        +yaxis*maxd_electron+shift_label)




        # build stuff for combined analysis
        # triangle and axis
        scaling=l+shift+yaxis+yaxis_extra-shift_label
        scaling=l/2.+shift
        xshift=shift+l+shift_label+mat_dend_sep
        yshift=-shift-l+shift_label
        yshift=yaxis+yaxis_extra-scaling
        latex_string+='  \draw[K,ultra thick] (%f,%f) -- (%f,%f);\n' % (xshift,
                                                                        yshift-shift_label,
                                                                        xshift+length,
                                                                        yshift-shift_label)
        latex_string+='  \draw[K,ultra thick, -stealth] (%f,%f) -- (%f,%f) node[anchor=east] {%s $R^2$};\n' % (xshift,
                                                                                                            yshift-shift_label,
                                                                                                            xshift,
                                                                                                            yshift+scaling,
                                                                                                            font_tick)
        if l>=6.:
            labels=[0.2*i for i in range(-2,6)]
        else:
            labels=[0.25*i for i in range(-1,5)]
        for i in labels:
            xco=yshift+scaling*math.sqrt( 0.5*(1.-i) )
            latex_string+='  \draw[K,thick,rectangle] (%f,%f) -- (%f,%f) node[anchor=east] {%s %.1f};\n' % (xshift+tick_length,
                                                                                                            xco,
                                                                                                            xshift-tick_length,
                                                                                                            xco,
                                                                                                            font_tick,
                                                                                                            i)
        if l<6.:
            title='H+E analysis'
        else:
            title='Combined (H+E) analysis'
        latex_string+='  \\node[anchor=center,rectangle] at (%f,%f) {%s %s};\n' % (xshift+length/2.,
                                                                                    yshift+scaling,
                                                                                    font_title,
                                                                                    title)
        latex_string+='  \\node[anchor=center,rectangle] at (%f,%f) {%s Automatic fragmentation};\n' % (xshift+length/2.,
                                                                                                        yshift-shift_label*2,
                                                                                                        font_title)


        # dendrogram
        colors_of_frags={}
        for i in range(len(dendro_all['icoord'])):
            yco=[ yshift+scaling*x  for x in dendro_all['dcoord'][i] ]
            xco=[ xshift+unit*x/10. for x in dendro_all['icoord'][i] ]
            color=coltrans[ dendro_all['color_list'][i]]
            latex_string+='  \\draw[%s,line width=3pt] (%f,%f) -- (%f,%f) -- (%f,%f) -- (%f,%f);\n' % (color,xco[0],yco[0],xco[1],yco[1],xco[2],yco[2],xco[3],yco[3])

            if dendro_all['dcoord'][i][0]==0.:
                x=int( (dendro_all['icoord'][i][0]-5.)/10. )
                x=dendro_all['leaves'][x]
                colors_of_frags[x]=color
            if dendro_all['dcoord'][i][3]==0.:
                x=int( (dendro_all['icoord'][i][3]-5.)/10. )
                x=dendro_all['leaves'][x]
                colors_of_frags[x]=color

        for i in range(len(dendro_all['icoord'])):
            yco=[ yshift+scaling*x  for x in dendro_all['dcoord'][i] ]
            xco=[ xshift+unit*x/10. for x in dendro_all['icoord'][i] ]
            color=coltrans[ dendro_all['color_list'][i]]
            if color=='GR':
                for j,jx in enumerate(dendro_all['leaves'][:-1]):
                    x1=xshift+0.5*unit+j*unit
                    x2=xshift+0.5*unit+(j+1)*unit
                    #sys.stderr.write('%f %f %f %f\n' % (x1,xco[0],xco[3],x2))
                    if x1<=xco[0]<x2 and yco[0]<yshift+scaling*maxd_molecule:
                        color=colors_of_frags[jx]
                        #sys.stderr.write('%f %f %f %s\n' % (x1,xco[0],x2,color))
                        latex_string+='  \\draw[%s,line width=3pt] (%f,%f) -- (%f,%f);\n' % (color,xco[0],yco[0],xco[0],yshift+scaling*maxd_molecule)
                    if x1<=xco[3]<x2 and yco[3]<yshift+scaling*maxd_molecule:
                        color=colors_of_frags[jx]
                        #sys.stderr.write('%f %f %f %s\n' % (x1,xco[3],x2,color))
                        latex_string+='  \\draw[%s,line width=3pt] (%f,%f) -- (%f,%f);\n' % (color,xco[3],yco[3],xco[3],yshift+scaling*maxd_molecule)

        # labels
        for i in range(nfrags):
            label=dendro_all['ivl'][i]
            latex_string+='  \\node[anchor=center,inner sep=1pt] at (%f,%f) {\\scriptsize \\scalebox{0.8}[1]{%s} };\n' % (xshift+0.5*unit+i*unit,
                                                                                    yshift-0.5*shift_label,
                                                                                    label)

        latex_string+='  \draw[K,ultra thick] (%f,%f) -- (%f,%f);\n' % (xshift+tick_length,
                                                                        yshift+scaling*maxd_molecule,
                                                                        xshift+length-tick_length,
                                                                        yshift+scaling*maxd_molecule)




        # center molecule to origin
        com=[0.,0.,0.]
        for atom in geom:
            for xyz in range(3):
                com[xyz]+=atom[1+xyz]
        for xyz in range(3):
            com[xyz]/=natom
        for iatom in range(natom):
            for xyz in range(3):
                geom[iatom][xyz+1]-=com[xyz]


        # calculate depth information
        alpha=math.radians(alpha)
        beta=math.radians(beta)
        normvec=[ math.sin(alpha)*math.sin(beta),-math.sin(alpha)*math.cos(beta),math.cos(alpha) ]

        zatom={}
        for iatom in range(natom):
            zatom[iatom+1]=sum( [ normvec[i]*geom[iatom][i+1] for i in range(3) ] )
        zmin=min( [ zatom[i] for i in zatom ] )
        zmax=max( [ zatom[i] for i in zatom ] )

        zbond={}
        for con in connectivity:
            zbond[con]=sum( [ normvec[i]*  (geom[con[0]-1][i+1] + geom[con[1]-1][i+1])/2. for i in range(3) ] )

        zboth=zatom.copy()
        zboth.update(zbond)
        zboth=sorted(zboth, key=lambda q: zboth[q])


        # make latex
        latex_string+='''



        \\begin{scope}[tdplot_main_coords,xshift=%fcm,yshift=%fcm,scale=%f]

        ''' % (xshift+length/2.,yshift-1.5*shift_label-(shift+l+yshift+shift_label)/2.,mol_scale)

        for iatom in range(natom):
            latex_string+='  \\node (x%i) at (%8.4f,%8.4f,%8.4f) {};\n' % (iatom,geom[iatom][1],geom[iatom][2],geom[iatom][3])
        latex_string+='\n\n\\tikzset{every node/.append style={very thick,circle, inner sep=0pt, draw, fill=white, font=\sffamily\\footnotesize}}\n\n'

        # make correlation lines
        thres=0.316
        for i in range(nfrags):
            iatom=frags[i][0]-1
            for j in range(nfrags):
                jatom=frags[j][0]-1
                if i!=j and DIS[i][j]<=thres:
                    cor=1.-2.*DIS[i][j]**2
                    latex_string+='  \\draw[%s, opacity=%5.3f,line width=%5.3fpt] (x%i) -- (x%i);\n' % ('O', min(1.,4.*(cor-0.75)),max(0.,10.*(cor-0.7)),iatom,jatom)


        boosts={'h': 0.6,
                're': 1.9,
                'mn': 1.7,
                'tc': 1.8,
                'ru': 1.8}


        # make bonds and atoms
        #colors=['RL','B','G',  'O','T','GL',  'P','BD','green',  'W','GR','black'   'Y']
        for obj in zboth:
            if isinstance(obj,int):
                iatom=obj
                z=100.*(zatom[iatom]-zmin)/(zmax-zmin)/2.+100./2.
                fill=z
                width=2./math.pi*math.atan2(1.,(100.-z)/50.)   *25.
                ## find color
                for ifrag,frag in enumerate(frags):
                    if iatom in frag:
                        break
                #for iclust,clust in enumerate(cluster):
                    #if ifrag+1 in clust:
                        #break
                #color=colors[iclust]
                ## find color from dendro_all
                element=geom[iatom-1][0].lower()
                if element in boosts:
                    boost=boosts[element]
                else:
                    boost=1.0
                color=colors_of_frags[ifrag]
                latex_string+='  \\node[fill=%s!%i, minimum width=%fpt] (x%i) at (x%i) {%s};\n' % (color,fill,width*boost,iatom-1,iatom-1,ifrag+1)
            else:
                con=obj
                latex_string+='  \\draw[very thick] (x%i) -- (x%i);\n' % (con[0]-1,con[1]-1)



        latex_string+='''



        \\end{scope}
        '''


        latex_string+=footer

        filename='diagram.tex'
        f=open(filename,'w')
        f.write(latex_string)
        f.close()



    # =============================================
    def colormap(self,value, colors_all):
        # value can be between -1 and +1 (red to white to blue)
        value=(value+1.)*10./2.
        big=int(value)
        small=value-big
        string='%s!%i!%s' % (colors_all[big+1],int(small*100),colors_all[big])
        return string

# =======================================================================
# =======================================================================
# =======================================================================

def run_plot():
    anOm = AnCorr_options('ancorr.in')
    anOm.read_OmFrag()
    anOm.AnCorr_input()
    anOm.compute()

    anOm.print_results(anOm.output_h,'Hole')
    anOm.print_results(anOm.output_e,'Electron')
    anOm.print_results(anOm.output_m,'Mixed')

    if anOm.opt_dict['dolatex']:
        anOm.create_latex(anOm.output_h, anOm.output_e, anOm.output_m)


def analyze_correlations():
    theo_header.print_header('Correlation analysis and clustering of Omega matrices ')
    run_plot()
