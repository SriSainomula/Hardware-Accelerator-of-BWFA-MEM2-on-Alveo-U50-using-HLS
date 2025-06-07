% Dot Matrix Sequence Alignment
% Define the sequences
%seq1 = 'ACGTGACCTG';
%seq2 = 'ACGTCGACCT';

seq1 = 'AGCTTAGCTA';
seq2 = 'CGTTAGCTAG';

%seq1 = 'abcdaefghijbklcmnopd';
%seq2 = 'abcdefghijklmnopqrst';

% Convert sequences to uppercase (to avoid case sensitivity)
seq1 = upper(seq1);
seq2 = upper(seq2);

% Get lengths of the sequences
len1 = length(seq1);
len2 = length(seq2);

% Initialize a dot matrix
dotMatrix = zeros(len2, len1);

% Fill the dot matrix with matches
for i = 1:len1
    for j = 1:len2
        if seq1(i) == seq2(j)
            dotMatrix(j, i) = 1; % Mark a match
        end
    end
end

% Plot the dot matrix
figure;
imshow(dotMatrix, 'InitialMagnification', 'fit');
colormap(flipud(gray));
title('Dot Matrix for Sequence Alignment');
xlabel('Sequence 1');
ylabel('Sequence 2');
%set(gca, 'XTick', 1:len1, 'XTickLabel', seq1, 'YTick', 1:len2, 'YTickLabel', seq2);
axis on;
grid on;

