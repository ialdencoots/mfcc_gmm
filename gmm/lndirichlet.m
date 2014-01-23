% Returns the natural logarithm of the value of the dirichlet distribution with
% dirichlet parameter vector alpha at x
function d = lndirichlet(x, alpha)

	d = sum((alpha.-1).*log(x)) - lnbeta(alpha);

end
