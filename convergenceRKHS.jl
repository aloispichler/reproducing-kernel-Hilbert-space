"""
	convergence of RKHS in a single point; display samples of a Gaussian random field
	created: 2021, May
	author©: Alois Pichler
"""

using Revise, Distributions
using Gnuplot; Gnuplot.options.gpviewer= true	# use external viewer
include("KernelsRKHS.jl")
include("../statistics/plotHistogram.jl")

function Kernel(x, y)
	GaussKernel(x,y; ℓ= .4)
end

println("───────────── ", time())
#f0= RKHSbyWeights([1.0], [1.], Kernel)	# the function f0(x)= E[f|X=x]
#f0(x)= 0.9* Kernel(x,0.3)- 0.3* Kernel(x,0.7)

#	│	choose an initial functions
#	╰────────────────────────────────────────────────────
f0(x)= max(0.0, x*(1.0-x))	#f0(x)= max(0.0, x*(1.0-x)), 1.0/(1+x^2)
#f0(x)= 1 < x < 2 ? 1 : 0	#Treppenfunktion
xSupport= range(-1.1, 2.1, length= 101)				# plot and integrate display area
@gp "set terminal wxt 0 raise; reset" :-
@gp :- xlabel="x domain" ylabel="y" title = "RKHS regression" "set key top left"
@gp :- xSupport f0.(xSupport) "w line lw 2 lc 'red' title 'f0'"		# f0 function

#	│	set parameters and go!
#	╰────────────────────────────────────────────────────
λ= .1			#	regression parameter
n= 100; 		#	sample values on the x-axis, tends to ∞
copies= 1000	#	samples repetition
samplesFrom= Normal(0.5, 1.0)		#Uniform(-1, 1), Normal(0.5, 1.), Exponential(.3)
Xi= zeros(n); fi= zeros(n)
fHat= RKHS
fλ= RKHS(f0, 100, samplesFrom, Kernel; λ= λ)
@gp :- xSupport fλ.(xSupport) "w line lw 3 lc 'purple' title 'fλ'"

f̅= zeros(length(xSupport)); normfλ= zeros(copies)
for i= 1:copies
	Xi= rand(samplesFrom, n); sort!(Xi)				# Xi-Werte, sort inplace
	fi= f0.(Xi)+ rand(Normal(0.0, .2), n)
	fHat= RKHS(Xi, fi; λ= λ, kernel= Kernel)
#	fHat= RKHSbyWeights(Xi, (fi- fλ(Xi))/λ, Kernel)	# Dommel's guess
	f̅ += fHat.(xSupport)/ copies					 # compute f̅ f᷉ f̂ ḟ x̅ X̅  x̃ x̄
	normfλ[i]= RKHSnorm(fHat - fλ)
end

@gp :- Xi fi "lc '#A69F00' pt 2 lw .1 tit 'samples'"
@gp :- xSupport f̅ "w line lw 2 lc 'green' title 'f̅'"
@gp :- xSupport fHat.(xSupport) "w line lc 'blue' title 'f̂|sample'"	# sample f f᷉

@gp "set terminal wxt 1 raise; reset" :-
@gp :- Xi fHat.w "w line lc 'blue' title 'f̂.weights'" :-
@gp :- "set y2tics; set ytics nomirror" :-
@gp :- fλ.x fλ.w "title 'fλ.weights' with line lw 2 lc'purple' axis x1y1" :-
#@gp :- (fλ.x)' (fλ.w-fHat.w) "title 'fλ.weights' with line lw 3 lc'red' axis x1y1" :-
@gp :- Xi (λ* fHat.w- fi+ fλ.(Xi)) "w line lw 1  lc 'green' title 'f̂.weights (right)' axis x1y2"


# histogram and density
plotHistogram(n* normfλ.^2; title= "n|f̂ - f0|^2, λ=$λ", windowNr= 2)
time()