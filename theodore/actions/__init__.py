from .actions import ActionFactory
from .analyze_tden import AnalyzeTden
from .analyze_sden import AnalyzeSden


def run():
    ActionFactory.from_commandline()
