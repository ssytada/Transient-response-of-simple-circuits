# -*- coding: utf-8 -*-

import math
from typing import Any
import symengine
import sympy
import configparser

import function

config_ini = configparser.ConfigParser()
config_ini.read('config.ini', encoding='utf-8')

sec = int(config_ini['PARAMETERS']['sec'])
div = int(config_ini['PARAMETERS']['div'])
ps_turn = int(config_ini['PARAMETERS']['ps_turn'])
a = int(config_ini['PARAMETERS']['a'])
p = int(config_ini['PARAMETERS']['p'])

t,s = sympy.symbols('t,s')

def laplace_tranform(t,s):
    f = sympy.laplace_transform(function.fundamental_function(t),t,s)[0]
    return f
lf = laplace_tranform(t,s)

def image_function(s,l,w) -> complex:
    f= symengine.Subs(l,s,w)
    return f

def theory_solution(a, sec, div):
    print("===== Theory =====")
    for i in range(1,sec * div + 1 ,1):
        time = float(i) / div
        at = a / time
        for j in range(1,sec * div + 1 ,1):
            s = function.fundamental_function(time)
        print(time,s)

def numerical_solution(a, p, ps_turn, sec, div):
    print("===== Numerial =====")
    eu = [ 0 for i in range( p+2 ) ]
    teu = [ 0 for i in range( p+2 ) ]
    teu[0] = 1.0
    eu[p+1] = 0.0

    for i in range( 1, p+1, 1 ):
        teu[i] = teu[i-1] * ( p + 1 - i ) / i
    for i in range( p, -1, -1 ):
        eu[i] = eu[i+1] + teu[i]

    for i in range (1,sec * div + 1, 1):
        time = float(i)/ div
        at = a / time
        x = 1.0
        fs = 0.0
        fse = 0.0

        for j in range(1, ps_turn, 1):
            cs = complex(at, (x - 0.5) * math.pi / time)
            cs = image_function(s,lf,cs)
            fs = fs + (-1) ** x * cs.imag
            x += 1.0
        l = 1
        fse = fs * eu[0]
        
        while True :
            if p < l :
                break
            
            cs = complex( at, (x - 0.5) * math.pi / time )
            cs = image_function(s,lf,cs)
            fse = fse + (-1) ** x * cs.imag * eu[l]
            x += 1
            l += 1
        
        fse = math.exp(a) * ( fse / eu[0] ) / time
        print(time,fse)

print("=============================")
print("function  : %s" % (symengine.expand(function.fundamental_function(t))))
print("parameter : a = %d" % (a))
print("            p = %d" % (p))
print("            UpperSec = %d" % (sec))
print("            TimeWidth = %d" % (div))
print("            PartialSum = %d" % (ps_turn))
print("=============================")

theory_solution(a, sec, div)
numerical_solution(a, p, ps_turn, sec, div)
