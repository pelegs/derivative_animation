#!/usr/bin/env python3
# -*- coding: iso-8859-15 -*-

import numpy as np
from math import sqrt
from sympy import *
from sympy.solvers import solve
import PyGnuplot as gp
from subprocess import call


def loop(frame, extremum=False):
    # Calculate values
    m = f_1.subs(x, x0)
    y0 = f.subs(x, x0)
    s = m*x + f.subs(x, x0) - m*x0
    Dx = R/sqrt(1+m**2)
    dx1, dx2 = x0-Dx/2.0, x0+Dx/2.0

    # Gnuplot stuff
    gp.c('set output "frames/{:03d}.png"'.format(frame))
    gp.c('unset arrow')
    gp.c('unset label')
    gp.c('set object 2 rect from graph 0,1 to graph 0.5,0.81 front fc rgb "white"')
    gp.c('set label "f\'(x)={:0.2f}" at graph 0.05, 0.9 tc rgb "blue" front font "Hershey/Symbol_Math, 40"'.format(float(m)))
    if extremum:
        y_2 = f_2.subs(x, x0)
        if y_2 < 0:
            type = 'maximum'
        elif y_2 == 0:
            type = 'inflection'
        else:
            type = 'minimum'
        gp.c('set label "{}" at graph 0.05, 0.83 tc rgb "red" front font "Hershey/Symbol_Math, 40"'.format(type))
    gp.c('set arrow from {},{} to {},{} lw 1 dt 1 lc rgb "gray" nohead'.format(x0, ymin, x0, ymax))
    gp.c('set arrow from {},{} to {},{} lw 1 dt 1 lc rgb "gray" nohead'.format(xmin, y0, xmax, y0))
    gp.c('plot {} notitle lw 3 lc rgb "#FF0000",\\'.format(f))
    gp.c('[{}:{}] {} notitle lw 3 lc rgb "#0000FF",\\'.format(dx1, dx2, s))
    gp.c('"-" u 1:2:3 with circles fill solid lc rgb "#0000FF" notitle')
    gp.c('{} {} 0.02'.format(x0, y0))
    gp.c('e')
    gp.c('unset label')


# Cleaning stuff
call('rm -f deriv.gif frames/*.png', shell=True)

# Function stuff
x = symbols('x')
f = (x**3 - x) * exp(-x**2/2) * x**2/2
f_1 = diff(f, x)
f_1_zeros = sorted([zero.evalf() for zero in solve(f_1, x)])
f_2 = diff(f_1, x)

# Parameters
xmin, xmax = -3, 3
ymin, ymax = -3, 3
N = 300
R = 0.5
pause_length = 50

# Gnuplot stuff
gp.c('set term png size 900, 900')
gp.c('set sample 1500')
gp.c('set xrange[{}:{}]'.format(xmin, xmax))
gp.c('set yrange[{}:{}]'.format(ymin, ymax))
gp.c('set xlabel "x" font "Hershey/Symbol_Math, 20"')
gp.c('set ylabel "y" font "Hershey/Symbol_Math, 20" rotate by 0')
gp.c('set size square')

# Generate frames
zero = f_1_zeros.pop(0)

xs = np.linspace(xmin, xmax, N)
dxs = xs[1]-xs[0]
frame = 0
for x0 in xs:
    # Check if point is, or closest to, a zero of the derivative
    if x0 <= zero <= x0+dxs:
        x0 = zero
        if f_1_zeros:
            zero = f_1_zeros.pop(0)
        for j in range(pause_length):
            frame += 1
            loop(frame, extremum=True)
    else:
        loop(frame)
    frame += 1
print('')

# Convert to animation
call(['convert', '-loop', '0', 'frames/*.png', 'deriv.gif'])
call(['eog', 'deriv.gif'])
