"""
	reproducing kernel Hilbert space RKHS – struct definitions
	created: 2021, May
	author©: Alois Pichler
"""

using SpecialFunctions, LinearAlgebra, Distributions

#	│	RKHS: define the struct: Reproducing Kernel Hilbert space
#	╰────────────────────────────────────────────────────
struct RKHS{T} 				# RKHS vector supported by T, describing the RKHS function
	x::Vector{T}						# array of testpoints
	w::Vector{Float64}					# corresponding weights
	GramMatrix::AbstractArray{Float64}	# Gram matrix, symmetric
	kernel								# the kernel function
end

#	│	functor to evaluate the RKHS function f at x
#	╰────────────────────────────────────────────────────
function (f::RKHS{T})(x::T) where T	# callable functor
	length(f.x) == length(f.w) || throw(DimensionMismatch("Support must match length of weights."))
	length(f.w) > 0 || return 0.0		# manage zero element
	fx= 0.0
	for j= 1:length(f.w)				# evaluate the kernel at x
		fx+= f.w[j]* f.kernel(x, f.x[j])
	end
	return fx/ length(f.w)
end

#	│	outer constructors
#	╰────────────────────────────────────────────────────
function RKHSbyWeights(x::Vector{T}, weights, kernel= GaussKernel; fillGram= false)	where T	#	function
	if fillGram
		RKHS(x, weights, Gram(x, kernel), kernel)
	else
		RKHS(x, weights, Array{Float64}(undef, 0, 0), kernel)
end	end

function RKHS(f0::RKHS; λ= 0.0, x= f0.x, kernel= f0.kernel)		#	λ-projection
	GramM= Gram(x, kernel); n= length(f0.w)
	RKHS(x, (λ*I*n+ GramM)\ (n* f0(x)), GramM, kernel)
end

function RKHS(x::Vector{T}, fx::Vector{Float64}; λ= 0.0, kernel= GaussKernel)::RKHS{T} where T	#	regression
	GramM= Gram(x, kernel); n= size(GramM,1)
	RKHS(x, (λ*I*n+ GramM)\ (n* fx), GramM, kernel)
end

function RKHS(x::Vector{T}; kernel= GaussKernel)::RKHS{T} where T	#	Gaussian random field
	GramM= Gram(x, kernel)
	RKHS(x, rand(Distributions.MvNormalCanon( 0.0000001I+ GramM)), GramM, kernel)
end

function RKHS(f0, nPoints::Int64, distD::UnivariateDistribution, kernel= GaussKernel; λ= 0.0)	# employ a distribution
	RKHS(f0, unique!(quantile(distD, range(1/2.0/nPoints, 1-1/2.0/nPoints, length= nPoints))), distD, kernel; λ=λ)
	# 	x= unique!(quantile(distD, range(1/2.0/nPoints, 1-1/2.0/nPoints, length= nPoints)))		# clever choice of sample points :-)
	# 	GramM= Gram(x, kernel); nPoints= length(x)
	# 	RKHS(x[:,:], (λ*I+ GramM/ n)\ f0.(x), GramM, kernel)
end

#	fλ= RKHS(f0, collect(xSupport), samplesFrom, Kernel; λ= λ)
function RKHS(f0, x::Vector{Float64}, D::UnivariateDistribution, kernel= GaussKernel; λ= 0.0)	# employ a distribution
	unique!(sort!(x))	# remove duplicates
	GramM= Gram(x, kernel); n= size(GramM, 1)
	fx= Array{Float64}(undef, length(x)); p= Array{Float64}(undef, length(x)); tmp= 0.0
	for (i,y) in enumerate(x)
		fx[i]= f0(y)
		if i< length(x)		# adjust the distribution
			tmp2= cdf(D, (x[i]+x[i+1])/ 2); p[i]= tmp2- tmp; tmp= tmp2
		else
			p[i]= 1.0- tmp
	end end
	RKHS(x, (λ*I+ Diagonal(p)*GramM)\ (Diagonal(p)* fx)* n, GramM, kernel)
end


#	│	clean the RKHS object
#	╰────────────────────────────────────────────────────
function RKHSclean(f::RKHS{T}; fillGram= true)::RKHS{T} where T
	length(f.x) == length(f.w) || throw(DimensionMismatch("Support must match length of weights."))
	ind= sortperm(f.x);					# permutation of the sorted x
	fw= copy(f.w);						# don't alter the weights vector
	if iszero(fw[ind[1]]); ind[1]= 0; end
	for i = 2:length(fw)				# scan for duplicates and removeables
		if ind[i-1] > 0 && f.x[ind[i]] == f.x[ind[i-1]]
			fw[ind[i]]+= fw[ind[i-1]]	# roll forward the weight of duplicate
			ind[i-1]= 0					# remove duplicate predecessor
		end
		if iszero(fw[ind[i]])			# weight 0 encountered
			ind[i]= 0					# remove actual
	end	end
	filter!(!iszero, ind)				# sweep out now all redundant elements
	if fillGram
		if size(f.GramMatrix) == (length(f.x),length(f.x))
			for i in ind				# scan all entries
				for j in ind
					if isnan(f.GramMatrix[i,j])			# provide missing entry
						if isnan(f.GramMatrix[j,i])		# employ symmetry
							f.GramMatrix[i,j]= f.GramMatrix[j,i]= f.kernel(f.x[i], f.x[j])
						else
							f.GramMatrix[i,j]= f.GramMatrix[j,i]
			end	end	end	end
			return RKHS(f.x[ind], fw[ind]* length(ind)/ length(f.w), Symmetric(f.GramMatrix[ind,ind]), f.kernel)
		else							# compute from scratch
			GramMatrix= Gram(f.x[ind], f.kernel)
			return RKHS(f.x[ind], fw[ind]* length(ind)/ length(f.w), Symmetric(GramMatrix), f.kernel)
		end
	else	# don't fill the matrix
		return RKHS(f.x[ind], fw[ind]* length(ind)/ length(f.w), Array{Float64}(undef, 0, 0), f.kernel)
end	end

#	│	- (minus) operation
#	╰────────────────────────────────────────────────────
function Base.:-(f1::RKHS{T}, f2::RKHS{T})::RKHS{T} where T		# minus operation
	f1.kernel == f2.kernel || throw(ArgumentError("RKHS - (minus): the kernel functions differ"))
	n1= length(f1.w); n2= length(f2.w)
	if length(f1.GramMatrix) > 0 || length(f2.GramMatrix) > 0
		GramM= fill(NaN, n1+n2, n1+n2)			# fill with NaN
		if length(f1.GramMatrix)== n1; GramM[1:n1,1:n1]= f1.GramMatrix; end				# rescue old
		if length(f2.GramMatrix)== n2; GramM[n1+1:n1+n2,n1+1:n1+n2]= f2.GramMatrix; end	# rescue old
		return RKHS([f1.x;f2.x], [(n1+n2)*f1.w/n1; -(n1+n2)*f2.w/n2], GramM, f1.kernel)
	else
		return RKHS([f1.x;f2.x], [(n1+n2)*f1.w/n1; -(n1+n2)*f2.w/n2], Array{Float64}(undef,0,0), f1.kernel)
end	end

#	│	inner product of RKHS functions
#	╰────────────────────────────────────────────────────
function RKHSinner(f1::RKHS, f2::RKHS)
	if f1.x== f2.x && f1.kernel== f2.kernel && size(f1.GramMatrix) == (length(f1.w),length(f1.w))
		return f1.x * f1.GramMatrix * f2.w/ length(f1.w)/ length(f2.w)
	elseif f1.kernel==f2.kernel
		sum= 0.0;
		for i= 1:length(f1.x) 
			for j= 1:length(f2.x)
				sum+= f1.w[i]* f1.kernel(f1.x[i], f2.x[j])
		end; end
		return sum/ length(f1.w)/ length(f2.w)
	else		#	the kernels differ
		@error "RKHSinner: the kernel functions differ." maxlog= 1
		return NaN
end	end

#	│	RKHS norm
#	╰────────────────────────────────────────────────────
function RKHSnorm(f::RKHS)
	ftmp= RKHSclean(f; fillGram= true)	# square root of the innter product
	sqrt(ftmp.w'* ftmp.GramMatrix* ftmp.w/ length(ftmp.w)/ length(ftmp.w))
end

#	Gram matrix
function Gram(ξ::Vector{T}, kernel) where T
	n= length(ξ); GramMatrix= Array{Float64}(undef, n, n)
	for i= 1:n
		for j= 1:i	# filled by symmetry
			GramMatrix[j,i]= GramMatrix[i,j]= kernel(ξ[i], ξ[j])
	end; end
	return GramMatrix
end

#	│	finally some Kernel functions
#	╰────────────────────────────────────────────────────
function LaplaceKernel(x, y; ℓ=1.)	#	Laplace covariance
	return exp(-norm(y-x)/ ℓ)
end

function MaternKernel(x, y; ℓ=1., ν=1.5)	# Matérn covariance
	z= sqrt(2*ν)* norm(y-x)/ ℓ				#	fractinoal smoothness
	z< 0.01 ? gamma(ν)* 2^(ν-1) : z^ν* besselk(ν, z)
end

function SigmoidKernel(x, y; ℓ=1.)	#	Sigmoid kernel
	return 1/(exp(norm(y-x)/ ℓ) + exp(-norm(y-x)/ ℓ))
end

function GaussKernel(x, y; ℓ=1.)	#	Gauss kernel
	return exp(-(norm(y-x)/ ℓ)^2)
end

function BrownianKernel(x, y; ℓ=1.)	#	Wiener process
	HurstH= 0.5						#	fractional Brownian motion
	return (x.^(2*HurstH)+ y.^(2*HurstH)- abs.(x.-y).^(2*HurstH))/2
#	return min.(x,y)[1]
end

function BrownianBridge(x, y; ℓ=1.)	#	Brownian bridge
	return min.(x,y)[1]- x[1]*y[1]
end

function multiquadraticKernel(x, y; ℓ=1., β= -.5)
	(1+ (norm(y-x)/ ℓ)^2)^ β
end

function StudentTKernel(x, y; ℓ=1., β= 1.0)
	1.0/ (1.0+ (norm(y-x)/ ℓ)^β)	#	StudentT/ Cauchy
end