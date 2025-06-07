moufflon = getgenbank('AB060288'); %Retrieve the nucleotide sequence of the prion protein (PrP) for the moufflon from the GenBank database.
takin = getgenbank('AB060290'); %Retrieve the nucleotide sequence of the prion protein (PrP) for the golden takin from GenBank.
Matches = seqdotplot(moufflon,takin,1,1);
Matches1 = seqdotplot(moufflon,takin,11,7); %The window size is set to 11 and the step size is set to 7 for the comparison.
Matches2 = seqdotplot(moufflon,takin,15,11); %The window size is set to 15 and the step size is set to 11 for the comparison.

example = getgenbank('M10051');