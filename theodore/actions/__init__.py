from .actions import ActionFactory
from .analyze_tden import AnalyzeTden, AnalyzeTdenUnr
from .analyze_sden import AnalyzeSden
from .theoinp import TheodoreInput
from .plot_vist import PlotVist
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
