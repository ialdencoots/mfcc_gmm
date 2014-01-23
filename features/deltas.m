% D = deltas(C, n)
% Returns an L x K-2n matrix of delta coefficients
% from the L x K cepstrum coefficient matrix C
% where L is the number of MFCCs in each of the K signal frames
% n is the number of frames before and after the target frame used to calculate deltas
% If n is excluded, it is assumed to be 2
function D = deltas(C, n)
	
	if nargin < 1 || nargin > 2
		error('Usage: deltas(C, [n]).')	
		return
	elseif nargin == 1
		n = 2;
	end

	L = size(C,1);
	K = size(C,2);
	i = 1:n;
	D = zeros(L,K);
	
	% Calculate deltas
	for k = 1+n : K-n
		D(:,k) = sum(i.*(C(:,k+i) - C(:,k-i)),2);
	end

	D = D(:,1+n : K-n);
	denom = 2 * sum((1:n).^2);
	D = D ./ denom;

end
