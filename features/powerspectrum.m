% [P freqs] = powerspectrum(S, fs, K)
%
% Returns the power spectrum P of framed signal S with sampling rate fs
% And its corresponding frequency values (independent variable)
% K is the number of points in the FFT (performs best when K is a power of 2)
% If K is not given, it is assumed to be 512
% P = (1/N)|FFT(S)|^2  This is the periodogram estimate of the power spectrum
% where N is the number of samples in a signal frame
% P is then truncated to K/2
function [P freqs] = powerspectrum(S, fs, K)

	if nargin < 2 || nargin > 3
		error('Usage powerspectrum(S, fs, [K]).')
		return
	elseif nargin == 2
		K = 512;
	end

	N = size(S,1);
	% Compute FFT and take half the values
	Sft = fft(S,K)(1:floor(K/2),:);
	% Compute power spectrum
	P = (1/N) .* abs(Sft).^2;
	%P = abs(Sft).^2;

	% Generate the frequency values corresponding to power spectrum points
	mult = fs/K;
	freqs = linspace(0,floor(K/2), floor(K/2)) .* mult;

end
