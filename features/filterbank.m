% H = filterbank(minfreq, maxfreq, M, freqs)
%
% Returns an N x M matrix H of filterbank multipliers
% extending between minfreq and maxfreq
% where N is the number of samples in the signal, i.e. length(freqs)
% and M is the number of triangular filters in the filterbank
% freqs is a vector of the frequency values of the input signal
function H = filterbank(minfreq, maxfreq, M, freqs)

	if nargin ~= 4
		error('Usage: filterbank(minfreq, maxfreq, M, freqs).')
		return
	end

	% evenly spaced mel-value filter boundaries
	melb = linspace(freq2mel(minfreq), freq2mel(maxfreq), M+2);
	% convert boundaries to frequencies
	fb = mel2freq(melb);

	N = length(freqs);
	H = zeros(N,M);

	% Calculate filterbank multipliers
	f = 1;
	for m = 1:M
		f += 1; % index into fb for center of filter m
		for k = 1:N
			if freqs(k) > fb(f+1) % filter is 0-valued after filter peak
				break;
			elseif freqs(k) < fb(f-1) % 0-valued before peak
				continue
			elseif freqs(k) < fb(f)
				H(k,m) = (freqs(k)-fb(f-1)) / (fb(f) - fb(f-1));
			else
				H(k,m) = (freqs(k)-fb(f+1)) / (fb(f) - fb(f+1));
			end
		end
	end

end
