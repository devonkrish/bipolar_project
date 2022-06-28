% script to do statistics post hoc on data
%%
savePlots = 1;

folderDataBase = '/Users/davidcaldwell/Box/KLEENLAB/David/Results/June results/';
folderFiguresCell = {fullfile(folderDataBase,'LL20'),fullfile(folderDataBase,'LL40'),fullfile(folderDataBase,'LL100'),fullfile(folderDataBase,'absDer')};
saveName = {'LL20','LL40','LL100','absDer'};
statsStruct = struct;

subjCell = [];
for jj = 1:length(permResultsCell)
    subjCell{end+1} = permResultsCell{jj}.subj;
end

numEls = length(permResultsCell);
c= cmocean('thermal',numEls);

% load in data
for jj = 1:length(saveName)
    processedInt = saveName{jj};
    folderFigures = folderFiguresCell{jj};
    dataFile  = fullfile(folderFigures,[processedInt '.mat']);
    load(dataFile);
    
    %%
    numChannelsVec = [];
    meanWidthVec = [];
    meanWidthSID = {};
    numChannelsVecSub = [];
    meanWidthVecSub = [];
    numChannelsSID = {};
    %
    for jjj = 1:length(permResultsCell)
        
        % if sum(strcmp(patientsVetted{jj},patientsInt))
        width = permResultsCell{jjj}.meanWidth;
        widthSub = permResultsCell{jjj}.meanWidthSub;
        num = permResultsCell{jjj}.numSigChannels;
        numSub = permResultsCell{jjj}.numSigChannelsSub;
        subj = permResultsCell{jjj}.subj;
        
        if ~isnan(width) & ~isnan(widthSub)
            meanWidthVec = [meanWidthVec width];
            meanWidthVecSub = [meanWidthVecSub widthSub];
            meanWidthSID{end+1} = subj;
        end
        
        if ~isnan(num) & ~isnan(numSub)
            numChannelsVec = [numChannelsVec num];
            numChannelsVecSub = [numChannelsVecSub numSub];
            numChannelsSID{end+1} = subj;
        end
        %  end
        
    end
    %% plot
    
    colorHist = cmocean('thermal',2);
    colorParallel = cmocean('matter',numEls);
    
    histFig = figure;
    tiledlayout(2,2,'TileSpacing','Compact','Padding','Compact');
    histFig.Position = [839 109 1408 1229];
    nexttile
    histogram(meanWidthVec,BinWidth=10)
    hold on
    histogram(meanWidthVecSub,BinWidth=10)
    legend({'High density','Sub-sampled'})
    title([processedInt ' Mean Width Histogram'])
    xlim([0 max([meanWidthVec(:);meanWidthVecSub(:)])+10])
    set(gca,'FontSize',14)
    
    nexttile
    histogram(numChannelsVec,BinWidth=1)
    hold on
    histogram(numChannelsVecSub,BinWidth=1)
    legend({'High density','Sub-sampled'})
    title([processedInt ' Number of Channels Histogram'])
    set(gca,'FontSize',14)
    xlim([0 max([numChannelsVec(:);numChannelsVecSub(:)])+3])
    
    ax= nexttile;
    ax.ColorOrder  = colorParallel;
    grid(ax,'on')
    xlim([0.5 2.5])
    ylim([0 max([meanWidthVec(:);meanWidthVecSub(:)])+10])
    hold on
    p1=parallelcoords([meanWidthVec;meanWidthVecSub]','group',[1:length(meanWidthVec)],'Labels',{'High density','Sub-sampled'},'LineWidth',2,'Marker','o');
    xlabel('High density vs. Sub-sampled')
    ylabel('Mean width')
    title([processedInt ' mean spike subject wise comparison'])
    set(gca,'FontSize',14)
    ax.XTickLabel ={'','High density','','Sub-sampled',''};
    
    ax = nexttile;
    ax.ColorOrder  = colorParallel;
    xlim([0.5 2.5])
    grid(ax,'on')
    ylim([0 max([numChannelsVec(:);numChannelsVecSub(:)])+3])
    hold on
    p2=parallelcoords([numChannelsVec;numChannelsVecSub]','group',[1:length(numChannelsVec)],'Labels',{'High density','Sub-sampled'},'LineWidth',2,'Marker','o');
    ylabel('Number of channels with significant spikes')
    xlabel('High density vs. Sub-sampled')
    title([processedInt ' number spike channels subject wise comparison'])
    set(gca,'FontSize',14)
    ax.XTickLabel ={'','High density','','Sub-sampled',''};
    
    if savePlots
        exportgraphics(histFig,fullfile(folderFigures,[processedInt '_2x2.png']),'Resolution',600)
        exportgraphics(histFig,fullfile(folderFigures,[processedInt '_2x2.eps']))
    end
    
    
    %% stats
    [pNum,hNum,statsNum] = signrank(numChannelsVec,numChannelsVecSub);
    
    [pWidth,hWidth,statsWidth] = signrank(meanWidthVec,meanWidthVecSub);
    
    statsStruct.name{jj} = processedInt;
    statsStruct.meanWidthSID{jj} = meanWidthSID;
    statsStruct.numChannelsSID{jj} = numChannelsSID;
    statsStruct.pNum{jj} = pNum;
    statsStruct.pWidth{jj} = pWidth;
    statsStruct.meanWidthVec{jj} = meanWidthVec;
    statsStruct.meanWidthVecSub{jj} = meanWidthVecSub;
    statsStruct.numChannelsVec{jj} = numChannelsVec;
    statsStruct.numChannelsVecSubj{jj} = numChannelsVecSub;
    statsStruct.meanWidthDiff{jj} = meanWidthVec - meanWidthVecSub;
    statsStruct.numChannelsDiff{jj} = numChannelsVec - numChannelsVecSub;
    
end
%%
figure
tiledlayout(2,1,'TileSpacing','Compact','Padding','Compact');
nexttile
hold on
for jj = 1:length(saveName)
    subjSpecific = statsStruct.meanWidthSID{jj};
    colormapSpecific = [];
    for jjj = 1:length(subjSpecific)
        ind = find(strcmp(subjCell,subjSpecific(jjj)));
        colormapSpecific = [colormapSpecific; c(ind,:)];
    end
    swarmchart(jj,statsStruct.meanWidthDiff{jj},[],colormapSpecific,'filled');
end
title('Differences in Mean Spike Duration by Condition')
xticks([0 1 2 3 4 5])
xticklabels({'','LL20','LL40','LL100','Absolute Derivative',''})

nexttile
hold on
for jj = 1:length(saveName)
    subjSpecific = statsStruct.numChannelsSID{jj};
    colormapSpecific = [];
    for jjj = 1:length(subjSpecific)
        ind = find(strcmp(subjCell,subjSpecific(jjj)));
        colormapSpecific = [colormapSpecific; c(ind,:)];
    end
    swarmchart(jj,statsStruct.numChannelsDiff{jj},[],colormapSpecific,'filled');
end
title('Differences in Number of Channels with Spikes by Condition')
xticks([0 1 2 3 4 5])
xticklabels({'','LL20','LL40','LL100','Absolute Derivative',''})
xlabel('condition')

tempFig = gcf;
if savePlots
    exportgraphics(tempFig,fullfile(folderDataBase,['across_conditions.png']),'Resolution',600)
    exportgraphics(tempFig,fullfile(folderDataBase,['across_conditions.eps']))
end


%%
saveNameSpecific = 'spikeStats.mat';
save(fullfile(folderDataBase,saveNameSpecific),'statsStruct');




