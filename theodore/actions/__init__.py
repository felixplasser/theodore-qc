from .actions import ActionFactory
from .analyze_tden import AnalyzeTden, AnalyzeTdenUnr
from .analyze_sden import AnalyzeSden
from .theoinp import TheodoreInput
from .plot_vist import PlotVist
from .. import theo_header


def run():
    settings = {
            'description': theo_header.ret_header(),
            'logo': None,
            'main_order': ['logo', 'description', 'pos_args', 'opt_args', 'usage', 'space'],
            'error_order': ['logo', 'usage', 'space', 'error'],
            'arg_format': {
                'name': 12,
                'comment': 40,
                'seperator': 2,
            },
    }
    ActionFactory.from_commandline(description=settings)
