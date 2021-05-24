"""
	generate and display RKHS functions
	created: 2020, July
	author©: Alois Pichler
"""

using Gnuplot; Gnuplot.options.gpviewer= true;	# external viewer
using Test, Revise
include("KernelsRKHS.jl")

@generated function ≂(x, y)		# compare structs element by element. (type with \\eqsim)
    if !isempty(fieldnames(x)) && x == y
        mapreduce(n -> :(x.$n == y.$n), (a,b)->:($a && $b), fieldnames(x))
    else
        :(x == y)
end	end

function Kernel(x, y)
	GaussKernel(x,y; ℓ= .1)
end

println("───────────── ", time())
f0= RKHSbyWeights([0.1, .7, 1.1], [.6, 0.6, 3.0], Kernel; fillGram= false)
f1= RKHSbyWeights([0.1, 1.1, 0.7, 0.3, 1.1], [1., 2, 1, 0., 3], Kernel)
println("Test: without Gram")
@show f1c= RKHSclean(f1; fillGram= false); 
@test f1c ≂ f0

println("Test: with Gram")
f3= RKHSclean(f1; fillGram= true);
f4= RKHSbyWeights(f1.x, f1.w, f1.kernel; fillGram= true); f4c= RKHSclean(f4)
@test f3 ≂ f4c

println("Test: with Gram")
f0w= RKHSbyWeights(f0.x, f0.w, f0.kernel; fillGram= true)
f1c= RKHSclean(f4)
@test f0w ≂ f1c

f5= RKHSbyWeights([0.1, .3, 0.1, 0.3], [1., 2, 0., 3], Kernel)
xSupport= range(-1.1, 2.1, length= 101)				# plot and integrate display area
@gp "set terminal wxt 0 raise" "reset" :-
@gp :- xlabel="x domain" ylabel="y" title = "Test"
@gp :- xSupport f1.(xSupport) "w line lw 3 lc 'red' title 'f1'"			# f1 function
@gp :- xSupport f5.(xSupport) "w line lw 3 lc 'blue' title 'f5'"		# f5 function
f6= f5-f1
@gp :- xSupport f6.(xSupport) "w line lw 3 lc 'green' title 'f5-f1'"	# f5-f1 function

@show RKHSnorm(f0-f5)
@show RKHSnorm(f0,f5)

fNull= RKHSbyWeights([0.], [0.], Kernel)
fEins= RKHSbyWeights([0.], [1.], Kernel)
@show RKHSnorm(fEins-fNull)
@show RKHSnorm(fNull,fEins)

f3= RKHSbyWeights([0.1, .7, 1.1], [.6, 0.6, 3.0], Kernel; fillGram= true)
f3.GramMatrix[2,3]= NaN
f3c= RKHSclean(f3, fillGram= true)
@test f3c ≂ f3
time()