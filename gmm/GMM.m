% p = GMM(X, pi, Mu, SIGMA)
%
% Returns p(x|pi,Mu,SIGMA) for each column(data point) of X for a
% Gaussian Mixture Model with weights pi, mean vectors Mu,
% and covariance matrices SIGMA.
function p = GMM(X, pi, Mu, SIGMA)

	if nargin ~= 4
		error('Usage: GMM(pi, Mu, SIGMA).')
		return
	elseif size(SIGMA,3) ~= size(Mu,2) || size(SIGMA,3) ~= length(pi) 
		error('Each Gaussian in the mixture must have a full set of parameters.')
		return
	elseif size(SIGMA,1) ~= size(Mu,1)
		error('Covariance matrix and mean vector dimensions must match.')
		return
	elseif size(SIGMA,1) ~= size(SIGMA,2)
		error('Covariance matrices in SIGMA must be square.')
		return
	elseif size(Mu,1) ~= size(X,1)
		error('Data must be of the same dimensionality as the model.')
		return
	end

	N = size(X,2); % # of data samples
	p = zeros(1,N);

	% Compute pdf value
	for k = 1:length(pi)
		p += pi(k) * multivargaussian(X, Mu(:,k), SIGMA(:,:,k));	
	end

end
