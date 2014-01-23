% C = cepstrumcoefficients(S, L)
% Returns a LxK matrix of mel frequency cepstrum coefficients
% from the K framed signal matrix of log energy values, S
% L is the number of cepstrum coefficients to be returned for each frame
% If L is not given, it is assumed to be 13
function C = cepstrumcoefficients(S, L)

	if nargin < 1 || nargin > 2
		error('Usage: cepstrumcoefficients(S, [L]).')
		return
	elseif nargin == 1
		L = 13;
	end

	C = dct(S)(1:L,:);
	%C = dct(S,L);

end
