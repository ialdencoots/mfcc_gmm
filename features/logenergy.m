% E = logenergy(P, H)
%
% Returns a MxK matrix of log energy filter coefficients
% P is an NxK matrix of K power spectrum frames of length N
% H is an NxM filterbank matrix of M filters for frames of length N
function E = logenergy(P, H)

	if nargin ~= 2
		error('Usage: logenergy(P, H).')
		return
	end

	warning("off","Octave:broadcast");

	M = size(H,2); % num filters
	K = size(P,2); % num frames
	E = zeros(M,K);

	% Sum the power values in each filter for each frame
	% Integration of power over frequency gives energy
	for k = 1:K
		E(:,k) = sum(H .* P(:,k))';
	end

	E = log(E);

end
