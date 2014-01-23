% F = extractfeatures(file)
%
% Extract featurevector matrix from the given .wav file with default 
% function parameters
% Feature vectors are 39 coefficients long: 13 MFCC, 13 delta,
% and 13 delta-delta.
function F = extractfeatures(file)

	if nargin ~= 1
		error('Usage: extractfeatures("file").')
		return
	end

	[signal, fs, bps] = wavread(file);
	signal = signal(:,1);

	S = framesignal(signal, fs);
	[P freqs] = powerspectrum(S, fs);
	H = filterbank(300, 8000, 20, freqs);
	E = logenergy(P, H);
	C = cepstrumcoefficients(E);
	F = featurevectors(C);

end
