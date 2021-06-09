"""
	visualize a realization of a 2D random field
	created: April 2020
	author©: Alois Pichler
"""

using Gnuplot; Gnuplot.options.gpviewer= true;	# external viewer
include("KernelsRKHS.jl")

function Kernel(x, y)
#	MaternKernel(x,y; ℓ= .1, ν= 1.5)
	SigmoidKernel(x,y; ℓ= .2)
end


#	│	test plots
#	╰────────────────────────────────────────────────────
println("───────────── ", time())
x= range(0, stop= 1, length= 50)
y= range(0, stop= 1, length= 51)
f2DRKHS=RKHS([[x,y] for y in y for x in x], kernel= Kernel)	#	new realization

@gp :- "reset"
@gp "set terminal wxt 0" :-
@gp :- "set title 'Gaussian random field'" :-
field= rand(MvNormal(0.0001I+ Gram([[x,y] for y in y for x in x], Kernel)))
@gsp :- "set size square" "set auto fix" :-
@gsp :- "unset xtics" "unset ytics" "unset ztics" "unset border" :-
@gsp :- x y reshape(field, length(x), :) "w pm3d title 'Random field'"

@gsp "set terminal windows 0 title 'function'" "reset" :-		# or wxt
@gsp :- "set title 'RKHS function'" "set xlabel 'x'" "set ylabel 'y'" :-
@gsp :- "unset xtics" "unset ytics" "unset ztics" "unset border" :-
@gsp :- x y reshape(f2DRKHS.([[x,y] for y in y for x in x]), length(x), :) "w pm3d title ''"

time()	