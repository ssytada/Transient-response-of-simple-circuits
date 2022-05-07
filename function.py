# -*- coding: utf-8 -*-

import symengine

def fundamental_function(t):
    f = 1 - symengine.exp(-(t))
    return f