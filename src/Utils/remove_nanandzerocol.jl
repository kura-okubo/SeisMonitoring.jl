"""
	remove_nanandzerocol(A::AbstractArray)

	Remove any nan and all zero column from CorrData.
"""
function remove_nanandzerocol(A::AbstractArray)

	N = size(A, 2)
	nancol = ones(Int64, N)
	for i = 1:N
		if any(isnan.(A[:, i])) || all(iszero, A[:,i])
			# this has NaN in its column
			nancol[i] = 0
		end
	end

	nancol=convert.(Bool, nancol)

	#NOTE: you cannot unbind entire array, so remove_nanandzerocol! is not used here.
	return A[:, nancol], nancol

end

remove_nanandzerocol_t(C, nancol) = (return C.t[nancol])

"""
	remove_nanandzerocol!(C::CorrData)

Remove any nan and all zero column from CorrData.
Modify C.t as well.
"""

function remove_nanandzerocol!(C::CorrData)
	C.corr, nancol = remove_nanandzerocol(C.corr)
	C.t = remove_nanandzerocol_t(C, nancol)
end
