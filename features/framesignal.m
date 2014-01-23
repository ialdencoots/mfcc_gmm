% S = framesignal(signal, fs, framelength, steplength)
%
% Splits the signal with sampling rate fs into frames framelength ms long
% spaced steplength ms apart and then applies a hamming window to each frame.
% Each column in the resultant matrix S is a frame
% If not given, framelength is assumed to be 25ms and steplength 10ms
function S = framesignal(signal, fs, framelength, steplength)

	if nargin < 1 || nargin == 3
		error ('Usage: framesignal(signal, fs, [framelength, steplength]).')
		return
	elseif nargin == 2
		framelength = 25;
		steplength = 10;
	end
	
	sigsize = length(signal);

	% the number of samples per frame
	framesize = floor((framelength/1000) * fs);
	% the number of samples between the start of adjacent frames 
	stepsize = floor((steplength/1000) * fs);

	% Count the number of frames
	numframes = 0;
	indx = 1;
	while (true)	
		numframes += 1;
		if (indx + framesize >= sigsize)
			break;
		end
		indx += stepsize;
	end

	% Initialize S
	S = zeros(framesize, numframes);

	% Populate S with values from the signal
	indx = 1 - stepsize;
	for (i = 1:numframes)
		indx += stepsize;
		if (i == numframes)
			len = length(signal(indx:end));
			S(1:len,i) = signal(indx : end);
		else
			S(:,i) = signal(indx : indx + framesize - 1);
		end
	end

	% Initialize hamming window filter
	w = hamming(framesize);

	% Filter results with hamming window
	S = S .* w;

end
