%Script takes 4D multi-tifs with channel order: red, green, blue, phase.
%outputs folder with spores and mothercells highlighted
%and background subtracted fluorescence values for mothercell/forespore
%pairs as well as the mc to fs fluorescence ratio
clc; 
close all;
clear all;
fprintf('Running mc_fs_fluorescence.m\n'); 
workspace;
imtool close all;

[baseFileName,PathName,FilterIndex] = uigetfile({'*.tif';'*.tiff'},'Select the .tif image file');

fullFileName = fullfile(PathName, baseFileName);
[pathstr,name,ext] = fileparts(fullFileName);
numframes = numel(imfinfo(fullFileName));
Image = imread(fullFileName,1);

nImage=im2double(Image);
nImage=nImage/max(max(nImage));
minImage=min(min(nImage));
exImage = zeros(size(nImage,1),size(nImage,2));
for i = 1:size(nImage,1)
    for j = 1:size(nImage,2)
        exImage(i,j)=nImage(i,j)-minImage*((-1/(1-minImage))*nImage(i,j)+(1/(1-minImage)));
    end
end

binaryImageSpore=im2bw(exImage,graythresh(exImage));
binaryImageMC=im2bw(exImage,graythresh(exImage)/5);
binaryImageMCInvert=imcomplement(binaryImageMC);
binaryImageSporeInvert=imcomplement(binaryImageSpore);

binaryImageSporetemp = imfill(binaryImageSpore, 'holes');
binaryImageMCtemp = imfill(binaryImageMC, 'holes');

binaryImageMCInvert(binaryImageMCtemp==0)=0;
binaryImageSporeInvert(binaryImageSporetemp==0)=0;
binaryImageSpore=binaryImageSporeInvert;

labeledImageSpore = bwlabel(binaryImageSpore, 8); 
coloredLabelsSpore = label2rgb (labeledImageSpore, 'hsv', 'k', 'shuffle');
labeledImageMC = bwlabel(binaryImageMCInvert, 8); 
coloredLabelsMC = label2rgb (labeledImageMC, 'hsv', 'k', 'shuffle');
M_spore = regionprops(labeledImageSpore, exImage, 'all');
M_mother = regionprops(labeledImageMC, exImage, 'all');

load('sporeFit.mat')
X=cat(2,[M_spore.Area]',[M_spore.MajorAxisLength]',[M_spore.MinorAxisLength]',...
    [M_spore.Eccentricity]',[M_spore.ConvexArea]',[M_spore.Solidity]',[M_spore.Extent]');
temp=predict(sporeFitV4,[M_spore.Area]');
filterIndexSpore = find(temp);


load('motherFit.mat')

Y=cat(2,[M_mother.Area]',[M_mother.MajorAxisLength]',[M_mother.MinorAxisLength]',...
    [M_mother.Eccentricity]',[M_mother.ConvexArea]',[M_mother.Solidity]',[M_mother.Extent]');
tempMC=predict(motherFitV6,Y);
filterIndexMC = find(tempMC);

ctrSpore = zeros(length(filterIndexSpore),2);
for j = 1 : length(filterIndexSpore)
        ctrSpore(j,1) = M_spore(filterIndexSpore(j)).Centroid(1);
        ctrSpore(j,2) = M_spore(filterIndexSpore(j)).Centroid(2);
end
for i = 1 : length(filterIndexMC)
    ctr=M_mother(filterIndexMC(i)).Centroid;
    theta=M_mother(filterIndexMC(i)).Orientation;
    lengthMajor=M_mother(filterIndexMC(i)).MajorAxisLength;
    xMajor = ctr(1)  +  [ -1 +1 ] * lengthMajor*cosd(-theta)/2;
    yMajor = ctr(2)  +  [ -1 +1 ] * lengthMajor*sind(-theta)/2;
    for j = 1 : length(filterIndexSpore)
        evalCen(i,j)=pdist([ctr(1),ctr(2);ctrSpore(j,1),ctrSpore(j,2)],'Euclidean');
        evalLine(i,j)=point_to_line([ctrSpore(j,1),ctrSpore(j,2),0],[xMajor(1),yMajor(1),0],[xMajor(2),yMajor(2),0]);
    end
end
for j = 1 : length(filterIndexSpore)
    for i = 1 : length(filterIndexMC)
        eval(i,j)=evalCen(i,j)+evalLine(i,j);
    end
    pair(j)=find(eval(:,j)==min(eval(:,j)));
end

for j = 1 : length(filterIndexSpore)
    pairCen(j)=evalCen(pair(j),j);
    pairLine(j)=evalLine(pair(j),j);
end
pairMatrix=cat(1,pairCen,pairLine)';
load('filterSnakes1.mat');
load('filterSnakes2.mat');
tempEval1=predict(filterSnakes,pairMatrix);
tempEval2=predict(filterSnakesV2,pairCen');
tempEval1(tempEval2==0)=0;
sporeNum=1:1:length(pair);
pairFilter=pair(tempEval1'==1);
sporeFilter=sporeNum(tempEval1'==1);
mkdir(name)
for s = 1 : length(pairFilter)
    filterImageMC = ismember(labeledImageMC, filterIndexMC(pairFilter(s)));
    filterImageSpore = ismember(labeledImageSpore, filterIndexSpore(sporeFilter(s)));
    subImage(s).MC=filterImageMC;
    subImage(s).Spore=filterImageSpore;
end

Image2 = imread(fullFileName,2);
oImage=exImage;
nImage=im2double(Image2);
nImage=nImage/max(max(nImage));
minImage=min(min(nImage));
exImage = zeros(size(nImage,1),size(nImage,2));
for i = 1:size(nImage,1)
    for j = 1:size(nImage,2)
        exImage(i,j)=nImage(i,j)-minImage*((-1/(1-minImage))*nImage(i,j)+(1/(1-minImage)));
    end
end
BWB=im2bw(oImage,graythresh(oImage)/4);
BWclose=imclose(BWB,strel('disk',5));
BWB=imfill(BWclose,'holes');
BWB2=imcomplement(BWB);
background=mean(Image2(BWB2));
Image2=Image2-background;
Image2(Image2<0)=0;


M_spore = regionprops(labeledImageSpore, exImage, 'all');
for s = 1 : length(pairFilter)
    sporePixels = M_spore(filterIndexSpore(sporeFilter(s))).PixelIdxList;
    MCPixels = M_mother(filterIndexMC(pairFilter(s))).PixelIdxList;
    meanSpore(s) = mean(Image2(sporePixels));
    meanMC(s) = mean(Image2(MCPixels));
    ratio(s)=meanSpore(s)/meanMC(s);
end

averageMC=mean(meanMC);
stdMC=std(meanMC);
secondFilter=0;
for s = 1 : length(pairFilter)
    if meanMC(s) < averageMC-2*stdMC
        ratio(s)=NaN;
        meanSpore(s)=NaN;
        meanMC(s)=NaN;
        secondFilter(s)=0;
    else
        secondFilter(s)=1;
    end
end

[uniqueValues,uniqueFilter,uniqueFilterBack]=unique(pairFilter(secondFilter==1));
pairFilterApplied=pairFilter(secondFilter==1);
sporeFilterApplied=sporeFilter(secondFilter==1);

convertIdx=0;
j=1;
for i = 1 : length(secondFilter)
    if secondFilter(i)==1
        convertIdx(j)=i;
        j=j+1;
    end
end
        
filterImageMC = ismember(labeledImageMC, filterIndexMC(pairFilterApplied(uniqueFilter)));
filterImageSpore = ismember(labeledImageSpore, filterIndexSpore(sporeFilterApplied(uniqueFilter)));

BWO=im2bw(exImage,graythresh(exImage));
labeledBWO=bwlabel(BWO,8);
BWmeasure=regionprops(labeledBWO,exImage,'all');

filterImageMC = ismember(labeledImageMC, filterIndexMC(pairFilterApplied(uniqueFilter)));
filterImageSpore = ismember(labeledImageSpore, filterIndexSpore(sporeFilterApplied(uniqueFilter)));


textFontSize = 12;	
labelShiftX = -5;
figure(length(pairFilterApplied(uniqueFilter))+1)
imshow(exImage);
hold on;
boundaries1 = bwboundaries(filterImageMC);
boundaries2 = bwboundaries(filterImageSpore);
boundaries=cat(1,boundaries1,boundaries2);
numberOfBoundaries = size(boundaries, 1);
for k = 1 : numberOfBoundaries
	thisBoundary = boundaries{k};
	plot(thisBoundary(:,2), thisBoundary(:,1), 'g', 'LineWidth', 2);
end
for k = 1 : length(sporeFilterApplied)
    sporeCentroid = M_spore(filterIndexSpore(sporeFilterApplied(k))).Centroid;	
    sporeCentroidX(k,1)=sporeCentroid(1);
    sporeCentroidY(k,1)=sporeCentroid(2);
    text(sporeCentroidX(k,1) + labelShiftX, sporeCentroidY(k,1), num2str(k), 'FontSize', textFontSize, 'FontWeight', 'Bold', 'color', 'r');

end
print(['-f' num2str(length(pairFilterApplied(uniqueFilter))+1)],[name '/' name '_all'],'-dpng','-r0')
hold off;

sporeAreas=regionprops(filterImageSpore,'Area');

ratio=ratio(uniqueFilter);
MC=meanMC(uniqueFilter);
Spore=meanSpore(uniqueFilter);
output=cat(1,ratio,MC,Spore,[sporeAreas.Area]);
mean(ratio,'omitnan')
std(ratio,'omitnan')
sum(~isnan(ratio))
outFileName=[name 'results.txt'];
writetable(cell2table(num2cell(output)),outFileName,'delimiter','\t')

figure(9000)
imshow(BWB2)
print(['-f' num2str(9000)],[name '/' name '_background'],'-dpng','-r0')


function d = point_to_line(pt, v1, v2)
      a = v1 - v2;
      b = pt - v2;
      d = norm(cross(a,b)) / norm(a);
end

