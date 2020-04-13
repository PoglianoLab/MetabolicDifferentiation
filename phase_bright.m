%Script takes 4D multi-tifs with channel order: red, green, blue, phase.
%outputMatrix=[numCompleted,numWorking,numTotal,percentCompleted,percentWorking,totalSporulating];

clc; 
clearvars;
close all;
clear all;
fprintf('Running phase_bright.m...\n');
workspace;
imtool close all;

[baseFileName,PathName,FilterIndex] = uigetfile({'*.tif';'*.tiff'},'Select the .tif image file');

fullFileName = fullfile(PathName, baseFileName);
[pathstr,name,ext] = fileparts(fullFileName);
numframes = numel(imfinfo(fullFileName));
Image = imread(fullFileName,1);

FM = imread(fullFileName,1);
G = imread(fullFileName,2);
DAPI = imread(fullFileName,3);
Phase = imread(fullFileName,4);

FMdouble=im2double(FM);
FMdouble=FMdouble/max(max(FMdouble));
Gdouble=im2double(G);
Gdouble=Gdouble/max(max(Gdouble));
DAPIdouble=im2double(DAPI);
DAPIdouble=DAPIdouble/max(max(DAPIdouble));
Phasedouble=im2double(Phase);
Phasedouble=Phasedouble/max(max(Phasedouble));
minI(1)=min(min(FMdouble));
minI(2)=min(min(Gdouble));
minI(3)=min(min(DAPIdouble));
minI(4)=min(min(Phasedouble));
exFM = zeros(size(FMdouble,1),size(FMdouble,2));
exG = zeros(size(Gdouble,1),size(Gdouble,2));
exDAPI = zeros(size(DAPIdouble,1),size(DAPIdouble,2));
exPhase = zeros(size(Phasedouble,1),size(Phasedouble,2));
for i = 1:size(FMdouble,1)
    for j = 1:size(FMdouble,2)
        exFM(i,j)=FMdouble(i,j)-minI(1)*((-1/(1-minI(1)))*FMdouble(i,j)+(1/(1-minI(1))));
        exG(i,j)=Gdouble(i,j)-minI(2)*((-1/(1-minI(2)))*Gdouble(i,j)+(1/(1-minI(2))));
        exDAPI(i,j)=DAPIdouble(i,j)-minI(3)*((-1/(1-minI(3)))*DAPIdouble(i,j)+(1/(1-minI(3))));
        exPhase(i,j)=Phasedouble(i,j)-minI(4)*((-1/(1-minI(4)))*Phasedouble(i,j)+(1/(1-minI(4))));
    end
end


exGadjust=imadjust(exG);
BWG=im2bw(exGadjust,graythresh(exGadjust));
BWDAPI=im2bw(exDAPI,graythresh(exDAPI)/2);
BWGfill=imfill(BWG,'holes');
BWGerode=imerode(BWGfill,strel('disk',3));
BWspore=BWGerode-BWDAPI;
BWspore(BWspore<0)=0;

BWFM=im2bw(exFM,graythresh(exFM)/3);
BWFMdilate=imdilate(BWFM,strel('disk',1));
BWFMfill=imfill(BWFMdilate,'holes');
BWbackground=imcomplement(BWFMfill);
averageBackground=mean(exPhase(BWbackground));
exPhaseAdjust=exPhase-averageBackground;
exPhaseAdjust(exPhaseAdjust<0)=0;

totalCellsImage=BWFMfill-(BWFM+BWbackground);
totalErode=imerode(totalCellsImage,strel('disk',3));

labeledCells = bwlabel(totalErode, 8); 
cellMeasurements = regionprops(labeledCells, exPhase, 'all');
Ac=[cellMeasurements.Area];
Bc=[cellMeasurements.MinorAxisLength];
Cc=[cellMeasurements.MajorAxisLength];
tempCells= Ac > 10 & Ac < 300;
filteredCells=find(tempCells);
filteredCells=ismember(labeledCells,filteredCells);


labeledSpore = bwlabel(BWspore, 8); 
sporeMeasurements = regionprops(labeledSpore, exPhaseAdjust, 'all');
allSporeIntensities = [sporeMeasurements.MeanIntensity];
A=[sporeMeasurements.Area];
B=[sporeMeasurements.MajorAxisLength];
C=[sporeMeasurements.MinorAxisLength];
D=[sporeMeasurements.Perimeter];
E=[sporeMeasurements.EulerNumber];
F=[sporeMeasurements.ConvexArea];
Z=F-A;
temp = A > 60 & A < 160 & B > 10 & B < 30 & C > 7 & C < 14 & D > 22 & D < 65 & E ==1 & Z < 50;

keeperIndexesSpore = find(temp);
keeperSpores = ismember(labeledSpore, keeperIndexesSpore);
sporeMask=imclearborder(keeperSpores);

phaseIntensities=allSporeIntensities(temp);

exDAPIadjust=imadjust(exDAPI);
BWDAPI=im2bw(exDAPIadjust,graythresh(exDAPIadjust));
BWDAPIerode=imerode(BWDAPI,strel('disk',2));

labeledEngulf=bwlabel(BWDAPIerode, 8);
engulfMeasurements = regionprops(labeledEngulf,exPhase,'all');
load('engulfmdl.mat')
A1=[engulfMeasurements.Area];
B1=[engulfMeasurements.MajorAxisLength];
C1=[engulfMeasurements.MinorAxisLength];
D1=[engulfMeasurements.Eccentricity];
E1=[engulfMeasurements.ConvexArea];
F1=[engulfMeasurements.Solidity];
G1=[engulfMeasurements.Extent];
H1=[engulfMeasurements.MeanIntensity];
X1=cat(2,A1',B1',C1',D1',E1',F1',G1',H1');
temp1=predict(membranes,X1);

keeperEngulf=ismember(labeledEngulf,find(temp1));

numCompleted=sum(temp);
numWorking=sum(temp1);
numTotal=sum(tempCells);
percentCompleted=numCompleted/(numCompleted+numWorking);
percentWorking=numWorking/(numCompleted+numWorking);
totalSporulating=(numCompleted+numWorking)/numTotal;

outputMatrix=[numCompleted,numWorking,numTotal,percentCompleted,percentWorking,totalSporulating];

outFileName=[name '_intensities.txt'];
writetable(cell2table(num2cell(phaseIntensities)),outFileName,'delimiter','\t')

% figure(1)
% imshow(keeperSpores+exFM)
% 
% figure(2)
% imshow(keeperEngulf+exFM)
% 
% figure(3)
% imshow(filteredCells+exFM)
% 
% figure(4)
% histogram(phaseIntensities,10)
% 
% figure(5)
% imshow(BWFM)