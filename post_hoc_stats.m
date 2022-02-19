% script to do statistics post hoc on data
%%
% load in data
dataFile = fullfile(getenv('OUTPUT_DIR_SERVER'),'absDer_jan_19_2022','absDer.mat');
load(dataFile);

%%
numChannelsVec = [];
meanWidthVec = [];
numChannelsVecSub = [];
meanWidthVecSub = [];
%
for jj = 1:length(permResultsCell)
    width = permResultsCell{jj}.meanWidth;
    widthSub = permResultsCell{jj}.meanWidthSub;
    num = permResultsCell{jj}.numSigChannels;
    numSub = permResultsCell{jj}.numSigChannelsSub;
    
    meanWidthVec = [meanWidthVec width];
    meanWidthVecSub = [meanWidthVecSub widthSub];
    numChannelsVec = [numChannelsVec num];
    numChannelsVecSub = [numChannelsVecSub numSub];
    
end
%% plot
figure
histogram(meanWidthVec,BinWidth=10)
hold on
histogram(meanWidthVecSub,BinWidth=10)
legend({'8mm','4mm'})
title('LL40 Mean Width')


figure
histogram(numChannelsVec,BinWidth=1)
hold on
histogram(numChannelsVecSub,BinWidth=1)
legend({'8mm','4mm'})
title('LL40 Number of Channels')

%% stats
[pNum,hNum,statsNum] = signrank(numChannelsVec,numChannelsVecSub);

[pWidth,hWidth,statsWidth] = signrank(meanWidthVec,meanWidthVecSub);



