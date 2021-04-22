"""
One script to rule them all
"""
import argparse

from collections import namedtuple


TheoAction = namedtuple("TheoAction", ("name", "run", "add_arguments", "help"))


class TheodoreParser:

    actions_name = '__theo_actions__'

    def __init__(self, actions, **kwargs):
        self._parser, self._subparser = self._setup(**kwargs)
        self._actions = actions
        self._add_subparser(self._subparser, actions)

    def _add_subparser(self, subparser, actions):
        for action in actions:
            parser = subparser.add_parser(action.name, help=action.help)
            action.add_arguments(parser)

    def _setup(self, **kwargs):
        parser = argparse.ArgumentParser(**kwargs)
        subparsers = parser.add_subparsers(title="actions")
        subparsers.required = True
        subparsers.dest = self.actions_name
        return parser, subparsers
  
    def run(self): 
        args = self._parser.parse_args()
        action = getattr(args, self.actions_name)
        return self._actions[action].run(args)


if __name__ == '__main__':
    ...
