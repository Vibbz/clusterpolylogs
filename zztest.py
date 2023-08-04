from sympy import Symbol
from sympy.abc import a,b,c
import re
import pickle
import numpy as np
from wolframclient.evaluation import WolframLanguageSession
from wolframclient.language import wl, wlexpr
'''
print(session.evaluate('Range[5]'))

expr=2*a+6*b-c

print([expr.coeff(var) for var in ['a','b','c']])

session.terminate()
'''
