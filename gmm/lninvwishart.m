% w = lninvwishart(S, nu, Sigma)
%
% Returns the natural logarithm of the value of the Inverse-Wishart distribution
% with parameters S,nu at Sigma
function w = lninvwishart(S, nu, Sigma)

	p = size(S,1);
	w = nu/2 * log(det(S)) - (nu*p/2) * log(2) - lnmultivargamma(p,nu/2) - ((nu+p+1)/2) * log(det(Sigma)) - .5*trace(S*inv(Sigma));
	
end
