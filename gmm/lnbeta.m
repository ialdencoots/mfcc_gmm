% Returns the natural logarithm of the beta distribution at alpha
function b = lnbeta(alpha)

	b = sum(lgamma(alpha)) - lgamma(sum(alpha));

end
