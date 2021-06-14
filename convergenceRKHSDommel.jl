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
#	│	choose an initial functions
#	╰────────────────────────────────────────────────────
#f0= RKHSbyWeights([1.0], [1.], Kernel)	# the function f0(x)= E[f|X=x]

fL2(x)= max(0.0, x*(1.0-x))
xSupport= range(-1.1, 1.2, length= 101)				# plot and integrate display area
#	│	set parameters and go!
#	╰────────────────────────────────────────────────────
λ= .1			#	regression parameter
n= 1000; 		#	sample values on the x-axis, tends to ∞
copies= 100		#	samples repetition
samplesFrom= Uniform(-1, 1)		#Uniform(-1, 1), Normal(0.5, 1.), Exponential(.3)
f0= RKHSbyWeights([0.0], [1.0], Kernel)

Xi= zeros(n); fi= zeros(n)
fHat= RKHS; fHat2= RKHS; fHat0= RKHS
fλ= RKHS(f0, 100, samplesFrom, Kernel; λ= λ/n)		#	Erwartungswert/ Grenzfunktion
@gp "set terminal wxt 0 raise; reset" :-
@gp :- xlabel="x domain" ylabel="y" title = "RKHS regression" "set key top left"
@gp :- xSupport f0.(xSupport) "w line lw 2 lc 'red' title 'f0'"		# f0 function
@gp :- xSupport fλ.(xSupport) "w line lw 3 lc 'purple' title 'fλ'"

f̅= zeros(length(xSupport)); normf0= zeros(copies)
for i= 1:copies
	Xi= rand(samplesFrom, n); sort!(Xi)				# Xi-Werte, sort inplace
	fi= f0.(Xi)+ rand(Normal(0.0, .2), n)
	fHat0= RKHS(Xi, fi; λ= λ/n, kernel= Kernel)
	fHat= RKHSbyWeights(Xi, (fi- fλ.(Xi))/λ, Kernel)	# Dommel's guess
	gλ= RKHS(f0, Xi, samplesFrom, Kernel; λ= λ/n)	#	Erwartungswert/ Grenzfunktion
	fHat2= RKHSclean(fHat+ gλ)	# Dommel's guess
	f̅ += fHat2.(xSupport)/ copies					 # compute f̅ f᷉ f̂ ḟ x̅ X̅  x̃ x̄
	normf0[i]= RKHSnorm(fHat0 - f0)
end

@gp :- Xi fi "lc 'blue' pt 2 lw .1 tit 'samples'" 		# #A69F00
@gp :- xSupport f̅ "w line lw 2 lc 'green' title 'f̅'"
@gp :- xSupport fHat0.(xSupport) "w line lc 'blue' title 'f̂0|sample'"	# sample f f᷉
@gp :- xSupport fHat2.(xSupport) "w line lc 'turquoise' title 'f̂2|sample'"	# sample f f᷉

@gp "set terminal wxt 1 raise; reset" :-
@gp :- Xi fHat0.w "w line lc 'blue' title 'f̂0.weights'" :-
@gp :- Xi fHat2.w "w line lc 'turquoise' title 'f̂2.weights'" :-
@gp :- "set y2tics; set ytics nomirror" :-
@gp :- fλ.x fλ.w "title 'fλ.weights' with line lw 2 lc'purple' axis x1y1" :-
#@gp :- (fλ.x)' (fλ.w-fHat.w) "title 'fλ.weights' with line lw 3 lc'red' axis x1y1" :-
@gp :- Xi (λ* fHat.w- fi+ fλ.(Xi)) "w line lw 1  lc 'green' title 'f̂.weights (right)' axis x1y2"


# histogram and density
plotHistogram(n* normf0.^2; title= "n|f̂ - f0|^2, λ=$λ", windowNr= 2)
time()