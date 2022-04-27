from colt import Plugin


class ActionFactory(Plugin):
    _is_plugin_factory = True
    _plugins_storage = '_actions'

    _user_input = """
    options =  :: str
    """

    @classmethod
    def _extend_user_input(cls, questions):
        questions.generate_cases("options", {action.name: action.colt_user_input
                                             for action in cls.plugins.values()})
    @classmethod
    def from_config(cls, config):
        config = config['options']
        for action in cls.plugins.values():
            if action.name == config.value:
                return action.from_config(config)
        raise Exception(f"Action '{config[config.value]}' unknown")


class Action(ActionFactory):

    _register_plugin = False

    @classmethod
    def from_config(cls, config):
        if not isinstance(cls.run, staticmethod):
            cls.run = staticmethod(cls.run)
        return cls.run(**config)
         
    def run(**options):
        raise NotImplementedError
