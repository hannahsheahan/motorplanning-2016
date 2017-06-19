function [ colours ] = ColourSelect( filename, kcolors )
%% function:  [colours] = ColourSelect(filename, kcolors)
% This function takes a .jpg image filename e.g. 'coloursample.jpg', and
% the number of colours that you wish to extract from it (kcolors), and uses the
% kmeans++ algorithm to cluster the image and find what colours to use for
% plotting. This is stored in the rgb output, 'colors' [n x 3 matrix].
% Issues: N/A
% Notes:  If wishing to select very pale colours (~white), adjust
%         variable 'whitethresh'
% Author: Hannah Sheahan, sheahan.hannah@gmail.com
% Year:   2017
%------------------------------------------------------------

rng(1);                         % reproduce the same colour order each time

rgb = imread(filename);
[nRow,nCol,nChannel] = size(rgb);
A = zeros(nRow*nCol,nChannel+2);

for i = 1:nRow
     for j = 1:nCol
         z = nCol*(i-1)+j;
         A(z,1) = i;
         A(z,2) = j;
         A(z,3:end) = rgb(i,j,:);
     end
end

[~,I,~]=unique(A(:,3:end),'rows');
[idx, Centres] = kmeans(A(:,3:5), kcolors+1);

% Centres will hold the major pixel values at the centre of each cluster,
% but note that one of these pixel values will be white space so get rid of
% this if all the elements in the colour are greater than whitethresh.
whitethresh = 256 - 5;
for i=1:(kcolors+1)
    tmp = (Centres(i,:) >= whitethresh);
    iswhite(i) = tmp(1)*tmp(2)*tmp(3);
end

Centres = Centres(~iswhite, :);  % get rid of the whitespace value
colours = Centres./256;  % set to MATLABs normalised pixel values

% An optional check that it's selecting the right colours:
%{
subplot(2,1,1);
imagesc(rgb)
subplot(2,1,2);
C=ones(nRow,nCol,nChannel);
%A(I,1:2) contains the unique rgb pixel locations
%make such pixels black color
C(A(I,1),A(I,2),:)=0;
imagesc(C)
%}

end

