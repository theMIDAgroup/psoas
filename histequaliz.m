function out = histequaliz(im) 
%% HISTEQUALIZ
% IN 
% - im: original image
% OUT
% - out: image with normalized histogram
%
%
% 
% by Daniela Schenone. Creative Commons (CC) License, 
% Attribution - NonCommercial 4.0  International (CC BY-NC 4.0)
% https://creativecommons.org/licenses/by-nc/4.0/legalcode

%% CODE
maxim = max(max(im));
minim = min(min(im));

histim = zeros(maxim- minim+1,1);
for jj = minim:maxim
    histim(jj - minim + 1) = numel(find(im==jj)); 
end

nhist = histim/sum(sum(histim)); %Perchè sommare (due volte)? In teoria la somma è sempre il n. elementi
out = zeros(size(im));
%nhistcum = cumsum(nhist);

for jj = minim:maxim
    temp = sum(nhist(1:jj - minim + 1));
    out(im == jj) = floor(maxim*temp);
end