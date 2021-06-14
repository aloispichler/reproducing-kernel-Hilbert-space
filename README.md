# Reproducing-kernel-Hilbert-space

Provides a framework in _Julia_ for working with functions from kernel Hilbert spaces.

The RKHS object consists of supporting points xi, i=1,..., n, of any type, a user-sepcified kernel function and weights wi. The file _KernelsRKHS_ collects the main definitions.
The function evaluation x↦Σi  wi k(x, xi) is implemented as a functor.

The package includes multiple constructors, the norm of an RKHS function, addition, multiplication and so on. Even constructors involving general probability distributions are provided and various, exemplary kernel functions.
