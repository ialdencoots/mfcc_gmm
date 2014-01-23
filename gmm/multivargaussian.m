% p = multivargaussian(X, mu, Sigma)
%
% Returns p(x|mu,Sigma) for each column (data point) of X for
% a multivariate gaussian distribution with mean vector mu and 
% covariance matrix Sigma.
function p = multivargaussian(X, mu, Sigma)

	if nargin ~= 3
		error('Usage: multivargaussian(X, mu, Sigma).')
		return
	elseif size(Sigma,1) ~= size(Sigma,2)
		error('Covariance matrix Sigma must be square.')
		return
	elseif size(mu,2) ~= 1
		error('Mean vector mu must be a column vector.')
		return
	elseif size(Sigma,1) ~= length(mu)
		error('Covariance matrix and mean vector dimensions must match.')
		return
	elseif length(mu) ~= size(X,1)
		error('Data must be of the same dimensionality as the model.')
		return
	end

	d = length(mu); % dimensionality
	N = size(X,2); % # of data samples

	% Compute pdf value
	normconst = 1 / ((2*pi)^(d/2) * sqrt(det(Sigma)));
	Diff = X .- mu;
	p = zeros(1,N);

	for i = 1:N
		gauss = exp(-0.5 * Diff(:,i)' * inv(Sigma) * Diff(:,i));
		p(i) = normconst * gauss;
	end

end
