% F = featurevectors(C, n)
%
% Returns an 3L x K-4n matrix of feature coefficient vectors from 
% MFCC matrix C, where L is the number of MFCCs in each of the 
% K signal frames, n is the number of frames before and after the
% target frame used to calculate deltas.
% The first and last 2n frames are excluded from the result
% Rows 1-L are MFCCs, rows L+1-2L are delta, and 2L+1-3L are 
% delta-delta coefficients. If n is excluded, it is assumed to be 2
function F = featurevectors(C, n)

	if nargin < 1 || nargin > 2
		error('Usage: featurevectors(C, [n]).')	
		return
	elseif nargin == 1
		n = 2;
	end

	K = size(C,2);

	% Calculate deltas
	D = deltas(C,n);
	% Calculate delta-deltas
	DD = deltas(D,n);

	% Append all coefficients to output vectors
	F = [C(:,1+2*n : K-2*n) ; D(:,1+n : K-3*n) ; DD];

end
