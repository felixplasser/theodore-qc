from .actions import ActionFactory
from .analyze_tden import AnalyzeTden, AnalyzeTdenUnr
from .analyze_sden import AnalyzeSden
from .theoinp import TheodoreInput
from .plot_vist import PlotVist
from .jmol_MOs import JMolMOs
from .babel import Babel
from .cc_opt import CCOpt
from .cc_check import CCCheck
from .dgrid_prep import DGridPrep
from .draw_moments import DrawMoments
from .extract_molden import ExtractMolden
from .fcd import FCD
from .jmol_vibs import JMolVibs
from .parse_libwfa import ParseLibwfa
from .plot_OmFrag import PlotOmFrag
from .plot_Om_bars import PlotOmBars
from .plot_frag_decomp import PlotFragDecomp
from .plot_graph import PlotGraph
from .plot_graph_nx import PlotGraphNx
from .spectrum import Spectrum
from .tden_OV import TDenOv
from .vmd_plots import VMDPlots
from .analyze_tden_soc import AnalyzeTdenSoc
from .convert_table import ConvertTable


from .. import theo_header


settings = {
            'description': theo_header.ret_header(),
            'logo': None,
            'error_order': ['logo', 'description', 'args', 'usage', 'space', 'comment', 'space', 'error'],
            'arg_format': {
                'name': 20,
                'comment': 50,
                'seperator': ' | ',
            },
            'subparser_args': {
                'title': 'Actions: %s',
            },
            'subparser_format': {
                'name': 20,
                'comment': 50,
                },
}

run = ActionFactory.from_commandline(description=settings, as_parser=True)
