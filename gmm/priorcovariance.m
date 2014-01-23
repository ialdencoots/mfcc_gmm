% S = priorcovariance(X, M)
%
% Returns a covariance matrix to be used as a prior for the 
% DxN dataset matrix X of N D-dimensional data points.
% M is the number of gaussian components used in the GMM.
function S = priorcovariance(X, M)

	D = size(X,1);
	N = size(X,2);

	xj = mean(X,2);

	s = (1/N) .* sum((X .- xj).^2, 2);

	S = diag((1/(M^(1/D)) .* s.^2));

end
