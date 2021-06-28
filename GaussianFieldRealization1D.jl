"""
	display samples of a Gaussian random field
	created: 2021, May
	author©: Alois Pichler
"""

using Gnuplot; Gnuplot.options.gpviewer= true;	# external viewer
include("KernelsRKHS.jl")

function Kernel(x, y)
#	MaternKernel(x, y; ℓ= .2, ν= 1.5)
#	LaplaceKernel(x,y; ℓ= 2.)
	GaussKernel(x,y; ℓ= .2)
#	BrownianKernel(x,y; HurstH=0.2)
#	BrownianBridge(x,y; ℓ= .2)
		SigmoidKernel(x,y; ℓ= .2)
#	multiquadraticKernel(x,y)
#	StudentTKernel(x,y)
end

println("────────────────────── ", time())

copies= 5
x= range(0, 1, length=1000)
@gp "set terminal wxt 0; reset; set title 'Gaussian random field'" :-
@gp :- "set title 'RKHS realizations'" "set xlabel 'x'" "set border 1" :-
for i= 1: copies
	fRKHS= RKHS(collect(x), kernel= Kernel)	#	new realization
	@gp :- x fRKHS.(x) "w line title 'realization $i' lw ($i<2?3:1)"
end
@gp

time()	