from .actions import ActionFactory
from .analyze_tden import AnalyzeTden
from .analyze_sden import AnalyzeSden
from .theoinp import TheodoreInput
from .plot_vist import PlotVist
from .. import theo_header


def run():
    ActionFactory.from_commandline(description=theo_header.ret_header())
