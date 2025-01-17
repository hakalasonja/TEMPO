#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Utility file meant to encapsulate custom exception classes. These custom exceptions
are used in testing to ensure that the user gets feedback on common errors to 
differentiate between a library issue and a user input issue.

If you find one of these exceptions in your project, please
refer to the documentation for troubleshooting.
When using iterables, please ensure that you strictly use the supported types, as checks
are only done for whether or not you provide iterables. 

I.e. while Iterable support is mandatory in most cases, 
this does not mean that simply using an iterable object will guarantee that the code works. 
Using lists is the safest approach to guaranteeting the success of the code. 

If you do NOT receive one of these exceptions, but something is going wrong, please
open an issue on Github --> Unless the issue occurs within your own custom provided functions!


@author: georgew79
"""

class TEMPO_ImproperInputException(Exception):
    def __init__(self, message):
        super().__init__(message)

class TEMPO_NullValueException(Exception):
    def __init__(self, message):
        super().__init__(message)

def print_warning(message: str):
    """
    Overridable function meant to dictate how / where warnings / other messages are printed.

    Parameters
    ----------
    message : str
        The message to print
    """

    print(message)

def castable(obj, type, allow_strs=False):
    # strings are NOT allowed as castable usually.
    try:
        x = type(obj)
        if isinstance(obj, str) and not allow_strs:
            return False
        return True
    except ValueError:
        return False
