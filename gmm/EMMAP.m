% [pi Mu SIGMA] = EMMAP(X, alpha, I, m_0, kappa, S_0, nu)
%
% Returns M sets of D-variate gaussian parameters for an M component GMM
% obtained through maximum a posteriori estimation expectation maximization
% run for I iterations.
% pi(mixture weights) is 1xM, Mu(means) is DxM, and SIGMA(variances) is DxDxM
%
% X is a DxN matrix of N D-dimensional input data points (feature vectors)
% m_0, kappa, S_0, and nu are Normal-inverse-Wishart distribution parameters
% m_0 must be Dx1 and S_0 must be DxD 
% alpha is an M element row vector of Dirichlet parameters
%
% The algorithm used is an implmentation of MAP estimation for GMMs as given
% in Kevin P. Murphy's "Machine Learning, A Probabilistic Perspective."
function [pi Mu SIGMA] = EMMAP(X, alpha, I, m_0, kappa, S_0, nu)

	if nargin ~= 3 && nargin ~= 7
		error('Usage: EMMAP(X, alpha, I, [m_0, kappa, S_0, nu]).')
		return
	end

	% Matrix dimension variables
	D = size(X,1); % feature dimensionality
	N = size(X,2); % # of samples
	M = length(alpha); % # of mixture components

	if nargin == 2
		nu = D + 2;
		S_0 = priorcovariance(X,M);
		m_0 = mean(X,2);
		kappa = .01;
	end

	% Initialize outputs
	Mu = zeros(D,M);
	SIGMA = zeros(D,D,M);

	x = zeros(D,M);
	S = zeros(D,D,M);

	% Keep track of log likelihood
	loglikelihood = 0;
	thresh = .0000001; % Stopping criterion
	iter = 1;

	% Responsibility matrix, uniformly initialized
	% R(k,i) is the responsibility component k has for data point i
	R = ones(M,N)./M;

	plotnum = 1;
	graphpts = [];
	f = figure;
	while (true)

		%%%%
		% M step
		%%%%
	
		% r(k) is the responsibility of component k
		r = sum(R, 2);
	
		for k = 1:M
	
			% Calculate x (mean)
			x(:,k) = sum(R(k,:) .* X,2) ./ r(k);
	
			% Calculate S (covariance)
			diff = X .- x(:,k);
			Sk = zeros(D,D);
			for i = 1:N
				Sk += R(k,i) .* diff(:,i) * diff(:,i)';
			end
			S(:,:,k) = Sk;
	
			% Calculate posterior mean estimate using priors
			if (kappa)
				Mu(:,k) = (r(k) .* x(:,k) + kappa .* m_0) / (r(k) + kappa);
			else
				Mu(:,k) = x(:,k);
			end
	
			% Calculate posterior covariance estimates using priors
			if (kappa)
				prefac = (kappa * r(k)) / (kappa + r(k));
				meandiff = x(:,k) - m_0;
				denom = kappa + r(k) + D + 2;
				SIGMA(:,:,k) = (S_0 + S(:,:,k) + prefac .* meandiff * meandiff') / denom;
			else
				SIGMA(:,:,k) = (S_0 + S(:,:,k)) / (nu + r(k) + D + 2);
			end
		end
	
		% Calculate mixture weights
		alphasum = sum(alpha);
		pi = (r' + alpha - 1) / (N + alphasum - D);
	
		%%%%
		% E Step
		%%%%

		% Update responsibility matrix and compute log likelihood

		likes = zeros(M,N);
		for k = 1:M
			likes(k,:) = multivargaussian(X, Mu(:,k), SIGMA(:,:,k));
		end
		R = pi' .* likes;
		R = R ./ sum(R);

		loglikelihood(iter) = (
	   		sum(sum(R.*log(pi')))
		+	sum(sum(R.*log(likes)))
		+	lndirichlet(pi,alpha)
		+	sum(lnNIW(m_0, kappa, S_0, nu, Mu, SIGMA))
			);

		% Check if stopping criterion has been met
		if (iter == I)
		%	if (loglikelihood(iter) - loglikelihood(iter-1) < thresh)
		%	if (iter > 10000)
				break
		%	end
		end

%		if iter == 1 || ~mod(iter,10)
%			graphpts(plotnum) = iter;
%			subplot(4,4,plotnum)
%			scatter(X(1,:),X(2,:),'y')
%			xlabel('C_0')
%			ylabel('C_1')
%			for i = 1:length(pi)
%				plot_gaussian_ellipsoid(Mu(:,i),SIGMA(:,:,i))
%			end
%			subplot(4,4,[13:16])
%    		plot(loglikelihood)
%			hold on
%			scatter(graphpts, loglikelihood(graphpts),30, 'r','d')
%			hold off
%			xlabel('# of iterations')
%			ylabel('Log likelihood + log prior')
%			drawnow
%			if (plotnum == 12)
%		%		saveas(f,'hallevikc0c1','png');
%			return
%			end
%			plotnum++;
%		end

		iter++;
	end


end
