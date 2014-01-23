% Returns the natural logarithm of the value of the
% p-variate gamma distribution at a
function g = lnmultivargamma(p,a)

	g = p*(p-1)/4*log(pi) + sum(lgamma(a .+ (1.-(1:p))./2));

end
