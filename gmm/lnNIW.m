%w = lnNIW(m_0, kappa, S_0, nu, mu, Sigma)
%
% Returns a row vector of natural log value of the Normal-inverse-Wishart distribution
% with parameters m_0, kappa, S_0, nu at each Mu, SIGMA
function w = lnNIW(m_0, kappa, S_0, nu, Mu, SIGMA)

	M = size(Mu,2);
	w = zeros(1,M);

	for k = 1:M
		w(k) = log(multivargaussian(Mu(:,k), m_0, (1/kappa).*SIGMA(:,:,k))) + lninvwishart(S_0, nu, SIGMA(:,:,k));
	end
end	
