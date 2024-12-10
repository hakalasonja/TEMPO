#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Tue Nov 19 

Utility file meant to encapsulate custom exception classes. These custom exceptions
are used in testing, and ensuring that the user gets feedback on common errors to 
differentiate between a library issue and a user input issue.

If you find one of these exceptions have been thrown in your project, please
refer to the docs for the class you are trying to use to make sure everything is right.
ALSO ensure that when using iterables you strictly use supported types. I.e. checks
are only done for whether or not you provide iterables, however many type hints
include Iterables of - say - type Pulse. These iterables cannot be strings, therefore,
they must be lists of Pulse objects. Please ensure that you correctly use lists then. 

I.e. while Iterable support is the only thing that is mandatory in most cases, this does not 
mean that the use of an iterable object will guarantee the code works. Using lists is the wisest
approach to guaranteeting the success of the provided code. 

If you do NOT receive one of these exceptions, but something is going wrong, please
open an issue on github --> Unless the issue occurs within your own custom provided functions!

Essentially, these exception definitions will catch if there is something 
obvious wrong with the user's input. 

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