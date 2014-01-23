% f = mel2freq(m)
% Converts mels to frequency (Hz)
function f = mel2freq(m)
	
	% Formula from Douglas O'Shaughnessy's 1987 book "Speech Communication: human and machine"
	% converted from log_10 to log_e and inverted
	f = 700 .* (exp(m./1127) .- 1);
	
end
