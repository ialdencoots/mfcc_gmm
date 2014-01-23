% [pi Mu SIGMA] = EMMLE(X,M)
% Returns M sets of D-variate gaussian parameters for an M component GMM
% obtained through maximum likelihood estimation expectation maximization
% X is a DxN matrix of N D-dimensional input data points (feature vectors)
% pi(mixture weights) is 1xM, Mu(means) is DxM, and SIGMA(variances) is DxDxM
function [pi Mu SIGMA] = EMMLE(X, M)

	D = size(X,1);
	N = size(X,2);

	% Initialize parameters

	% Responsibility matrix
	R = ones(M,N)./M;
	%
	% Need to init Mu, SIGMA
	
	%%%%
	% M step
	%%%%

	r = sum(R,2);

	for k = 1:M
		Mu(:,k) = sum(R(k,:) .* X,2) ./ r(k);
	end

	pi = r./N;

	%%%%
	% E Step
	%%%%

	for k = 1:M
		R(k,:) = pi(k) .* multivargaussian(X, Mu(:,k), SIGMA(:,:,k))	
	end
	R = R ./ sum(R);

end
