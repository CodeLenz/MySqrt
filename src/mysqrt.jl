using LinearAlgebra

# Julia version for the algorith in 
# https://mathsfromnothing.au/matrix-square-root/?i=1
function Sqrt(A::AbstractMatrix)
  
	# Dimension
	n,m = size(A) 

	# Check if it is square
	m==n || throw("Sqrt:: Matrix should be square")

    # For some reason, very small numbers cause a catastrofic
	# error in sqrt and Sqrt. It should be something related
	# to schur(A)
	A[abs.(A).<sqrt(eps(1.0))].=zero(eltype(A))

	# Check if A is diagonal...we cannot miss this oportunity :o)
	if isdiag(A)
		return sqrt.(A)
	end

	# Schur decomposition of A
	S,U = Schur{Complex}(schur(A))

 	# Diagonal of S
	D = diag(S)

    # Check if matrix has sqrt
	# by looking for zeros in D
	any(isapprox.(D,0)) && throw("Sqrt:: matrix has no sqrt ")
	
	# sqrt of main diagonal of S
	X = diagm(sqrt.(D))

    if !isdiag(S)
  		@inbounds for j=2:n
    		@inbounds for i=j-1:-1:1
				k = i+1:j-1
				X1 = @view X[i,k]
				X2 = @view X[k,j]
				somat = transpose(X1)*X2
				denom = X[i,i]+X[j,j]
				isapprox(abs(denom),0) && throw("Sqrt::division by zero")	
    		  	X[i,j] = (S[i,j]-somat)/(denom)
   			end
  		end
	end

    return U*X*adjoint(U)

end

