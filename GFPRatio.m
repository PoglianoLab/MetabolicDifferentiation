clc;
clearvars;
close all;
clear all;
fprintf('Running GFPRatio.m...\n');

format long g;
format compact;
captionFontSize = 7;
outFileName='test.tif';

% Check that user has the Image Processing Toolbox installed.
hasIPT = license('test', 'image_toolbox');
if ~hasIPT
	message = sprintf('Sorry, but you do not seem to have the Image Processing Toolbox.\nDo you want to try to continue anyway?');
	reply = questdlg(message, 'Toolbox missing', 'Yes', 'No', 'Yes');
	if strcmpi(reply, 'No')
		return;
	end
end

baseFileName = 'ImageFileHere.tif';
folder = fileparts(which(baseFileName)); % Determine where file folder is (works with all versions).
fullFileName = fullfile(folder, baseFileName);
[pathstr,name,ext] = fileparts(fullFileName);
if ~exist(fullFileName, 'file')
	if ~exist(baseFileName, 'file')
		warningMessage = sprintf('Error: the input image file\n%s\nwas not found.\nClick OK to exit the program.', fullFileName);
		uiwait(warndlg(warningMessage));
		fprintf(1, 'Finished running getTheseMfSnakesoffMyMfPlane.m.\n');
		return;
	end
	fullFileName = baseFileName;
end

% If we get here, we should have found the image file.
numframes = numel(imfinfo(baseFileName));
originalImage = imread(fullFileName,1);

%normalize image over [0,1] and convert to double
doubleOri=im2double(originalImage);
doubleOri=doubleOri/max(max(doubleOri));
minOri=min(min(doubleOri));
exOri = zeros(size(doubleOri,1),size(doubleOri,2));
for i = 1:size(doubleOri,1)
    for j = 1:size(doubleOri,2)
        exOri(i,j)=doubleOri(i,j)-minOri*((-1/(1-minOri))*doubleOri(i,j)+(1/(1-minOri)));
    end
end

%use two thresholds on same image to get high (spores) and low (mothercells) intensity
binaryImageSpore=im2bw(exOri,graythresh(exOri));
binaryImageMC=im2bw(exOri,graythresh(exOri)/5);
binaryImageMCInvert=imcomplement(binaryImageMC);
binaryImageSporeInvert=imcomplement(binaryImageSpore);

binaryImageSporetemp = imfill(binaryImageSpore, 'holes');
binaryImageMCtemp = imfill(binaryImageMC, 'holes');

binaryImageMCInvert(binaryImageMCtemp==0)=0;
binaryImageSporeInvert(binaryImageSporetemp==0)=0;
binaryImageSpore=binaryImageSporeInvert;


% Label each object so we can make measurements of it (color for visual clarity)
labeledImageSpore = bwlabel(binaryImageSpore, 8); 
coloredLabelsSpore = label2rgb (labeledImageSpore, 'hsv', 'k', 'shuffle');
labeledImageMC = bwlabel(binaryImageMCInvert, 8); 
coloredLabelsMC = label2rgb (labeledImageMC, 'hsv', 'k', 'shuffle');
blobMeasurements = regionprops(labeledImageSpore, exOri, 'all');
numberOfBlobs = size(blobMeasurements, 1);
blobMeasurementsMC = regionprops(labeledImageMC, exOri, 'all');
numberOfBlobsMC = size(blobMeasurementsMC, 1);

%load spore fit tree classifier from data on hundreds of spores previously measured
load('sporeFitV4.mat')
allBlobIntensities = [blobMeasurements.MeanIntensity];
allBlobAreas = [blobMeasurements.Area];
allBlobCentroids = [blobMeasurements.Centroid];
A=[blobMeasurements.Area];
B=[blobMeasurements.MajorAxisLength];
C=[blobMeasurements.MinorAxisLength];
D=[blobMeasurements.Eccentricity];
E=[blobMeasurements.ConvexArea];
F=[blobMeasurements.Solidity];
G=[blobMeasurements.Extent];
X=cat(2,A',B',C',D',E',F',G');
temp=predict(sporeFitV4,A');
keeperIndexesSpore = find(temp);
%keep only objects which pass spore classifier
keeperBlobsImageSpore = ismember(labeledImageSpore, keeperIndexesSpore);


%load mothercell fit tree classifier from data on hundreds of mothercells previously measured
load('motherFitV6.mat')
allBlobIntensitiesMC = [blobMeasurementsMC.MeanIntensity];
allBlobAreasMC = [blobMeasurementsMC.Area];
allBlobCentroidsMC = [blobMeasurementsMC.Centroid];
A1=[blobMeasurementsMC.Area];
B1=[blobMeasurementsMC.MajorAxisLength];
C1=[blobMeasurementsMC.MinorAxisLength];
D1=[blobMeasurementsMC.Eccentricity];
E1=[blobMeasurementsMC.ConvexArea];
F1=[blobMeasurementsMC.Solidity];
G1=[blobMeasurementsMC.Extent];

X1=cat(2,A1',B1',C1',D1',E1',F1',G1');
%keep only objects which pass mothercell classifier
tempMC=predict(motherFitV6,X1);
keeperIndexesMC = find(tempMC);

%for debugging
% textFontSize = 9;	
% labelShiftX = -5;
% figure(length(pairFilter)+1)
% imshow(exOri);
% hold on;
% boundaries = bwboundaries(labeledImageMC);
% numberOfBoundaries = size(boundaries, 1);
% for k = 1 : numberOfBoundaries
% 	thisBoundary = boundaries{k};
% 	plot(thisBoundary(:,2), thisBoundary(:,1), 'g', 'LineWidth', 2);
%     blobCentroid = blobMeasurementsMC(k).Centroid;	
%     blobCentroidX(k,1)=blobCentroid(1);
%     blobCentroidY(k,1)=blobCentroid(2);
%     text(blobCentroidX(k,1) + labelShiftX, blobCentroidY(k,1), num2str(k), 'FontSize', textFontSize, 'FontWeight', 'Bold', 'color', 'r');
% 
% end
% hold off;


keeperBlobsImageMC = ismember(labeledImageMC, keeperIndexesMC);

%pair each spore with its corresponding mothercell (strict to reduce false positives at the expense of missing some true positives)
for j = 1 : length(keeperIndexesSpore)
        ctrSpore(j,1) = blobMeasurements(keeperIndexesSpore(j)).Centroid(1);
        ctrSpore(j,2) = blobMeasurements(keeperIndexesSpore(j)).Centroid(2);
end
for i = 1 : length(keeperIndexesMC)
    ctr=blobMeasurementsMC(keeperIndexesMC(i)).Centroid;
    theta=blobMeasurementsMC(keeperIndexesMC(i)).Orientation;
    lengthMajor=blobMeasurementsMC(keeperIndexesMC(i)).MajorAxisLength;
    xMajor = ctr(1)  +  [ -1 +1 ] * lengthMajor*cosd(-theta)/2;
    yMajor = ctr(2)  +  [ -1 +1 ] * lengthMajor*sind(-theta)/2;
    for j = 1 : length(keeperIndexesSpore)
        evalCen(i,j)=pdist([ctr(1),ctr(2);ctrSpore(j,1),ctrSpore(j,2)],'Euclidean');
        evalLine(i,j)=point_to_line([ctrSpore(j,1),ctrSpore(j,2),0],[xMajor(1),yMajor(1),0],[xMajor(2),yMajor(2),0]);
    end
end
for j = 1 : length(keeperIndexesSpore)
    for i = 1 : length(keeperIndexesMC)
        eval(i,j)=evalCen(i,j)+evalLine(i,j);
    end
    pair(j)=find(eval(:,j)==min(eval(:,j)));
end

for j = 1 : length(keeperIndexesSpore)
    pairCen(j)=evalCen(pair(j),j);
    pairLine(j)=evalLine(pair(j),j);
end
X2=cat(1,pairCen,pairLine)';

%load classifiers to eliminate false positive spore/mothercell pairs
load('filterSnakesV1.mat');
load('filterSnakesV2.mat');
tempEval=predict(filterSnakes,X2);
tempEvalV2=predict(filterSnakesV2,pairCen');
tempEval(tempEvalV2==0)=0;
sporeNum=1:1:length(pair);
pairFilter=pair(tempEval'==1);
sporeFilter=sporeNum(tempEval'==1);
mkdir(name)
for s = 1 : length(pairFilter)
    keeperBlobsImageMC = ismember(labeledImageMC, keeperIndexesMC(pairFilter(s)));
    keeperBlobsImageSpore = ismember(labeledImageSpore, keeperIndexesSpore(sporeFilter(s)));
    subImage(s).MC=keeperBlobsImageMC;
    subImage(s).Spore=keeperBlobsImageSpore;
    figure(s)
    imshow(keeperBlobsImageMC+keeperBlobsImageSpore)
    print(['-f' num2str(s)],[name '/' name '_object' num2str(s)],'-dpng','-r0')
end

originalImage1 = imread(fullFileName,2);
oriImg=exOri;


%normalize 2nd channel image over [0,1]
doubleOri=im2double(originalImage1);
doubleOri=doubleOri/max(max(doubleOri));
minOri=min(min(doubleOri));
exOri = zeros(size(doubleOri,1),size(doubleOri,2));
for i = 1:size(doubleOri,1)
    for j = 1:size(doubleOri,2)
        exOri(i,j)=doubleOri(i,j)-minOri*((-1/(1-minOri))*doubleOri(i,j)+(1/(1-minOri)));
    end
end
BWB=im2bw(oriImg,graythresh(oriImg)/4);
BWclose=imclose(BWB,strel('disk',5));
BWB=imfill(BWclose,'holes');
BWB2=imcomplement(BWB);
background=mean(originalImage1(BWB2));
originalImage1=originalImage1-background;
originalImage1(originalImage1<0)=0;


blobMeasurements = regionprops(labeledImageSpore, exOri, 'all');
for s = 1 : length(pairFilter)
    sporePixels = blobMeasurements(keeperIndexesSpore(sporeFilter(s))).PixelIdxList;
    MCPixels = blobMeasurementsMC(keeperIndexesMC(pairFilter(s))).PixelIdxList;
    meanSpore(s) = mean(originalImage1(sporePixels));
    meanMC(s) = mean(originalImage1(MCPixels));
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
        
keeperBlobsImageMC = ismember(labeledImageMC, keeperIndexesMC(pairFilterApplied(uniqueFilter)));
keeperBlobsImageSpore = ismember(labeledImageSpore, keeperIndexesSpore(sporeFilterApplied(uniqueFilter)));


%remove tricky ones
BWOri=im2bw(exOri,graythresh(exOri));
labeledBWOri=bwlabel(BWOri,8);
BWmeasure=regionprops(labeledBWOri,exOri,'all');

%need to quantify how much overlap in pixels (if not enough, then dont
%count it
% sporeObj=0;
% MCObj=0;
% for j = 1 : size(pairFilter,2)
%     sporeIdxPixels = blobMeasurements(keeperIndexesSpore(sporeFilter(j))).PixelIdxList;
%     MCIdxPixels = blobMeasurementsMC(keeperIndexesMC(pairFilter(j))).PixelIdxList;
%     for k = 1 : size(BWmeasure,1)
%         combineSpore=cat(1,BWmeasure(k).PixelIdxList,sporeIdxPixels);
%         combineMC=cat(1,BWmeasure(k).PixelIdxList,MCIdxPixels);
%         if length(find(ismember(BWmeasure(k).PixelIdxList,sporeIdxPixels)))==length(sporeIdxPixels)
%             sporeObj=k;
%         end
%         if length(find(ismember(BWmeasure(k).PixelIdxList,MCIdxPixels)))>0.8*length(MCIdxPixels)
%             MCObj=k;
%         end
%     end
%     if sporeObj~=MCObj && sporeObj~=0 && MCObj~=0
%         secondFilter(j)=0;
%         ratio(j)=NaN;
%     end
%     sporeObj=0;
%     MCObj=0;
% end

keeperBlobsImageMC = ismember(labeledImageMC, keeperIndexesMC(pairFilterApplied(uniqueFilter)));
keeperBlobsImageSpore = ismember(labeledImageSpore, keeperIndexesSpore(sporeFilterApplied(uniqueFilter)));

%print boundaries for quality control
textFontSize = 12;	
labelShiftX = -5;
figure(length(pairFilterApplied(uniqueFilter))+1)
imshow(exOri);
hold on;
boundaries1 = bwboundaries(keeperBlobsImageMC);
boundaries2 = bwboundaries(keeperBlobsImageSpore);
boundaries=cat(1,boundaries1,boundaries2);
numberOfBoundaries = size(boundaries, 1);
for k = 1 : numberOfBoundaries
	thisBoundary = boundaries{k};
	plot(thisBoundary(:,2), thisBoundary(:,1), 'g', 'LineWidth', 2);
end
for k = 1 : length(sporeFilterApplied)
    blobCentroid = blobMeasurements(keeperIndexesSpore(sporeFilterApplied(k))).Centroid;	
    blobCentroidX(k,1)=blobCentroid(1);
    blobCentroidY(k,1)=blobCentroid(2);
    text(blobCentroidX(k,1) + labelShiftX, blobCentroidY(k,1), num2str(k), 'FontSize', textFontSize, 'FontWeight', 'Bold', 'color', 'r');

end
print(['-f' num2str(length(pairFilterApplied(uniqueFilter))+1)],[name '/' name '_all'],'-dpng','-r0')
hold off;

%output GFP ratio and save in results text file
ratio=ratio(uniqueFilter);
MC=meanMC(uniqueFilter);
Spore=meanSpore(uniqueFilter);
output=cat(1,ratio,MC,Spore);
mean(ratio,'omitnan')
std(ratio,'omitnan')
sum(~isnan(ratio))
outFileName=[name 'resultsV5.txt'];
writetable(cell2table(num2cell(output)),outFileName,'delimiter','\t')

figure(9000)
imshow(BWB2)
print(['-f' num2str(9000)],[name '/' name '_background'],'-dpng','-r0')


function d = point_to_line(pt, v1, v2)
      a = v1 - v2;
      b = pt - v2;
      d = norm(cross(a,b)) / norm(a);
end

fprintf('finished sporeMothercellGFPRatio.m...\n');