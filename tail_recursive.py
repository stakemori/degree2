# -*- coding: utf-8; mode: sage -*-
# Borrowed from http://code.activestate.com/recipes/496691/
class tail_recursive(object):

    def __init__(self, func):
        self.func = func
        self.firstcall = True
        self.CONTINUE = object()

    def __call__(self, *args, **kwd):
        if self.firstcall:
            func = self.func
            CONTINUE = self.CONTINUE
            self.firstcall = False
            try:
                while True:
                    result = func(*args, **kwd)
                    if result is CONTINUE: # update arguments
                        args, kwd = self.argskwd
                    else: # last call
                        return result
            finally:
                self.firstcall = True
        else: # return the arguments of the tail call
            self.argskwd = args, kwd
            return self.CONTINUE
