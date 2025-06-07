% Define the DNA sequences
%seq1 = 'abcdaefghijbklcmnopd';
%seq2 = 'abcdefghijklmnopqrst';

seq2 = 'AGCTTAGCTA';
seq1 = 'CGTTAGCTAG';

% Create a dot matrix
dotMatrix = zeros(length(seq1), length(seq2));
% Fill the dot matrix
for i = 1:length(seq1)
    for j = 1:length(seq2)
        if seq1(i) == seq2(j)
            dotMatrix(i, j) = 1;
        end
    end
end
% Plot the dot matrix
figure;
imagesc(dotMatrix);
colormap(flipud(gray));
axis equal;
xlabel('Sequence 1');
ylabel('Sequence 2');
title('Dot Matrix Comparison');
% Initialize variables to store information about the longest contiguous line
maxLineLength = 0;
maxLineStart = [0, 0];
maxLineEnd = [0, 0];
% Check horizontal lines
for i = 1:size(dotMatrix, 1)
    lineLength = 0;
    for j = 1:size(dotMatrix, 2)
        if dotMatrix(i, j) == 1
            lineLength = lineLength + 1;
            if lineLength > maxLineLength
                maxLineLength = lineLength;
                maxLineStart = [i, j - lineLength + 1];
                maxLineEnd = [i, j];
            end
        else
            lineLength = 0;
        end
    end
end
% Check vertical lines
for j = 1:size(dotMatrix, 2)
    lineLength = 0;
    for i = 1:size(dotMatrix, 1)
        if dotMatrix(i, j) == 1
            lineLength = lineLength + 1;
            if lineLength > maxLineLength
                maxLineLength = lineLength;
                maxLineStart = [i - lineLength + 1, j];
                maxLineEnd = [i, j];
            end
        else
            lineLength = 0;
        end
    end
end
% Check diagonal lines (from bottom-left to top-right)
for k = -(size(dotMatrix, 1) - 1):(size(dotMatrix, 2) - 1)
    lineLength = 0;
    for i = max(1, k+1):min(size(dotMatrix, 1), size(dotMatrix, 2) + k)
        j = i - k;
        if dotMatrix(i, j) == 1
            lineLength = lineLength + 1;
            if lineLength > maxLineLength
                maxLineLength = lineLength;
                maxLineStart = [i - lineLength + 1, j - lineLength + 1];
                maxLineEnd = [i, j];
            end
        else
            lineLength = 0;
        end
    end
end
% Check diagonal lines (from top-left to bottom-right)
for k = 1:(size(dotMatrix, 1) + size(dotMatrix, 2) - 1)
    lineLength = 0;
    for i = max(1, k-size(dotMatrix, 2)+1):min(k, size(dotMatrix, 1))
        j = k - i + 1;
        if dotMatrix(i, j) == 1
            lineLength = lineLength + 1;
            if lineLength > maxLineLength
                maxLineLength = lineLength;
                maxLineStart = [i - lineLength + 1, j + lineLength - 1];
                maxLineEnd = [i, j];
            end
        else
            lineLength = 0;
        end
    end
end
% Highlight the longest contiguous line with maximum matches
hold on;
plot([maxLineStart(2) - 0.5, maxLineEnd(2) + 0.5], [maxLineStart(1) - 0.5, maxLineEnd(1) - 0.5], 'r', 'LineWidth', 2);
hold off;
% Display the longest contiguous line information
fprintf('Longest contiguous line with maximum matches:\nStart: (%d, %d)\nEnd: (%d, %d)\nLength: %d\n', ...
    maxLineStart(1), maxLineStart(2), maxLineEnd(1), maxLineEnd(2), maxLineLength);