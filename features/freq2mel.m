% m = freq2mel(f)
% Converts frequency in Hz to mels
function m = freq2mel(f)

	% Formula from Douglas O'Shaughnessy's 1987 book "Speech Communication: human and machine"
	% converted from log_10 to log_e
	m = 1127 .* log(1 + f./700);

end
