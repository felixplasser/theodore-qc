Usage
-----

.. colt_commandline:: theodore run

    main_order = usage, args, space
    alias = theodore

    [arg_format]
    name = 20
    comment = 60

    [subparser_format]
    name = 25
    comment = 60

.. colt_commandline:: theodore run
   :subparsers: options(*)
   :header: True

    main_order = comment, usage, pos_args, opt_args, subparser_args, space
    alias = theodore

    [arg_format]
    name = 20
    comment = 60

    [subparser_format]
    name = 25
    comment = 60
