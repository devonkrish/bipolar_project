

% BIPOLAR PAIR ANALYSIS: EACH VS. ALL

pt='EC183'; % EC175 and EC183 both have intact 16x16 square grids (channel #s 1:256)
nchtocheck=256; %***ideally, should eventually convert this into an index of channels to check (STG, IFG, etc) which could be defined before the spectra loop to save computational time... 
windowstocheck=1:100;
binsz=3; % bin size in mm

onlygrids=true;
onlydepths=false;
onlystrips=false;

cm=cool(6); cm(1,:)=[0 0 0]; 
datadir='/Volumes/KLEEN_DRIVE/David/Bipolar project/baseline-high-density-data/bandpassfiltered/';
cd([datadir])
load('/Volumes/KLEEN_DRIVE/David/Bipolar project/taggedspikes.mat')
sfx=512;
frxrange=[2 200]; %frequency range to examine
  ft=[2 5 10 20 50 100 200]; ftl=cellstr(num2str(ft')); %frequency labels for plots
 
u=dir; uptbl={}; for i=1:length(u); uname=u(i).name; uptbl{i,1}=uname(1:end-23); end; uptbl(1:2)=[]; clear i u uname

p=find(strcmpi(pts,pt)); %patient number ID
pblocks=strfind(uptbl,pts{p}); 
for i=1:length(pblocks); 
    isbl(i,1)=~isempty(pblocks{i}); 
end
ptbl=find(isbl); if ~isempty(ptbl); disp(['Loading ' pts{p} ' blocks...']); end

% load all blocks for this patient and stack their baseline windows together
d=[]; nwind=0;
for b=1:length(ptbl); disp(uptbl{ptbl(b)})
    % load using using "_jk" versions of baseline windows, updated 2/2022
    load([datadir uptbl{ptbl(b)} '_baselineWindows_jk.mat'])
    % get rid of baseline windows containing spikes or artifact
    spksarti=hasspk | hasarti;
    nonspks_windows(:,spksarti)=[];
    hasstim(spksarti)=[]; %update indices for which windows overlap with stimuli/speech
    hasspeech(spksarti)=[]; 
    clear hasspkvec hasspk hasartivec hasarti spksarti % now clear spike- and artifact-related variables from workspace

    % convert to 3D matrix, combine all windows from consecutive blocks for each patient
    for i=1:size(nonspks_windows,2); 
        d(:,:,i+nwind)=nonspks_windows{2,i}'; 
    end
    nwind=size(d,3);

    clear nonspks_windows info
end; clear b

nch=size(d,2); 

% load electrode component infor (grid/strip/depth and how many linear contacts they have in a row
[bpN,bpT]=xlsread(['/Users/davidcaldwell/code/high_density_ecog/AN_ElectrodeInfoTDT.xlsx'],pts{p});
[em,eleclabels,anatomy]=getelecs(pts{p},2);

cm=cool(6); cm(1,:)=[0 0 0]; 
datadir='/Volumes/KLEEN_DRIVE/David/Bipolar project/baseline-high-density-data/bandpassfiltered/';
cd([datadir])
load('/Volumes/KLEEN_DRIVE/David/Bipolar project/taggedspikes.mat')
sfx=512;
frxrange=[2 200]; %frequency range to examine
  ft=[2 5 10 20 50 100 200]; ftl=cellstr(num2str(ft')); %frequency labels for plots
  

% if wanting to only look at grids, strips, or depths, then nan the others
if onlygrids||onlystrips||onlydepths;
  for r=1:size(bpT,1)
    if any(strcmpi(bpT(r,2),{'grid','minigrid'})) && ~onlygrids  || ...
       strcmpi(bpT(r,2),'strip')              && ~onlystrips || ...
       strcmpi(bpT(r,2),'depth')              && ~onlydepths;
     d(:,bpN(r,1):bpN(r,2),:)=nan;
    end
  end; clear r
end


% bad channels
badchI=isnan(mean(mean(d,1),3))'; %any channel with nans in any window will be a "bad channel"
badchI(badchidx{p})=1; %follow up with all marked as bad channels in original preprocessing
x=[]; xbch=true(size(badchI)); for i=1:size(bpT,1); x=[x bpN(i,1):bpN(i,2)]; end; xbch(x)=false; 
badchI(xbch)=1; %find and nan empty channels unaccounted for by component rows
d(:,badchI,:)=nan;
okc=~badchI; clear x xbch


%*at this point, isolate speech or stim or non-speech/stim

%*at this point, isolate STG or IFG channels
[STG]=getelecs_region(pt,'stg',2);
[IFG]=getelecs_region(pt,{'po','pt'},2);


%% ALL PAIRS (each vs. all others) analysis and example plot
 d=d(:,:,windowstocheck); clear Straces_allch; %free up RAM by getting rid of whatever won't be used (only using first ___ number of windows)
                 % ***hint hint: opportunity here to select speech or stim windows!
 [M,Mbpdist,frx]=bpspectra_EachVsAll(d,sfx,frxrange,em,nchtocheck);
 
 
 %% mean across windows
 M=sq(mean(M,4)); 
 
 
 %% unstack and line up into 2D matrix to create channels^2 X frequencies for easier indexing-->binning
 Mflat=[];
 Mflat_bpdist=[];
 for c2=1:nchtocheck; 
     Mflat=[Mflat; sq(M(:,c2,:))];               
     Mflat_bpdist=[Mflat_bpdist; Mbpdist(:,c2)]; % corresponding distance index
 end

 binz=0:binsz:85; % can change binsz PRN at this point to test different bin resolutions on the plots below
 clear mx my mz
 nbinz=length(binz)-1;
 mz=[]; 
 %bin by distance and take the mean, creating: frequency X distance (binned)
 for i=1:nbinz 
     mz(:,i)=nanmean(Mflat(Mflat_bpdist>=binz(i) & Mflat_bpdist<binz(i+1),:),1);
     binindex_min_ltnmax(i,:)=[binz(i) binz(i+1)];
 end
 
% Notes:
% mz is frequency X distance (binned), after having averaged across windows
% binindex_min_ltmax tells you for each column of mz what was the 
%         1] minimum (>0) distance, and
%         2] the "less than max" (<) distance
%         that were used to index bipolar pairs for that bin
%         Note: first bin includes zero, which corresponds to
%         the bin containing the referential channels
 
 % zscore the log transformed power according to frequency
 nfrx=length(frx);
 for i=1:nfrx; nns=~isnan(mz(i,:)); 
     mz_z(i,nns)=zscore(log(mz(i,nns))); 
 end; 

 %%
 figure('color','w','position',[277 223 1095 963]); colormap(parula); 

% Plot confusion matrix of bipolar distances between all pairs
 subplot(8,2,1:2:5); pcolor(Mbpdist); shf; hold on; plot([1 1; 1 256]',[1 256;256 256]','k-','linewidth',1); plot([1 1 1 128 256],[1 128 256 256 256],'k.'); %plot([1 128 256 256 256]/2,[1 1 1 128 256]/2,'k.'); 
    axis equal off; set(gca,'ydir','reverse','xdir','reverse'); colormap(gca,flipud(cmocean('deep'))); colorbar('fontsize',12); title('Bipolar pair distance','fontsize',14,'fontweight','normal')

% Histogram of number of pairs per bin
 subplot(8,2,7); histogram(make1d(Mbpdist),0:binsz:85,'facecolor',.5*[1 1 1]); set(gca,'fontsize',12); xlabel('Binned bipolar distance (mm)','fontsize',14); 
 ylabel('Counts/bin','fontweight','normal'); axis tight; grid on; cb=colorbar; set(cb,'visible','off');

 % line plot of same as above, zscored
 subplot(2,2,2); hold on; cmo=fliplr(cmocean('thermal',nfrx+round(nfrx)*.1)); colormap(gca,cmo); cb=colorbar; set(cb,'ticks',0:.25:1,'ticklabels',0:50:200,'fontsize',12)
 for i=nfrx:-1:1; nns=~isnan(mz(i,:));
     plot(binz(2:size(mz_z,2)+1),mz_z(i,:),'-','color',cmo(i,:),'linewidth',1); 
 end; grid on; axis tight; 
 set(gca,'fontsize',14); ylabel('ln(power)'); xlabel('Bipolar distance (mm)','fontsize',14); 
 text(max(xlim)+diff(xlim)/4,mean(ylim),'Frequency (Hz)','fontsize',12,'rotation',90,'horizontalalignment','center')

 subplot(2,2,3); 
 pcolor(binz(2:size(mz_z,2)+1),frx,mz_z); shading flat; set(gca,'ydir','normal'); ylabel('Frequency (Hz)'); xlabel('Distance (mm)'); set(gca,'fontsize',14); colorbar; 
 text(max(xlim)+diff(xlim)/4,mean(ylim),'ln(power)','fontsize',12,'rotation',90,'horizontalalignment','center')
 title('(z-scored by frequency)','fontweight','normal')
 set(gca,'yscale','log','ytick',ft,'yticklabel',ftl);

 sp(2,2,4)
 [mx,my]=meshgrid(binz(2:size(mz_z,2)+1)+binsz/2,frx);
 surf(mx,my,(mz_z)); xlabel('Bipolar distance (mm)','fontsize',14); ylabel('Frequency (Hz)','fontsize',14); zlabel('ln (Power)','fontsize',14)
 view(140,25); set(gca,'xdir','reverse','ydir','reverse','ytick',ft,'yticklabel',ftl,'fontsize',14); xlim([binsz/2 85]);  axis tight; 
 set(gca,'yscale','log'); 
  title('(z-scored by frequency)','fontweight','normal')%set(gca,'yscale','linear'); set(gca,'zscale','log'); set(gca,'zscale','linear')
 %ylim([0 30]); %zlim([.035 3]); caxis([.035 3]); %low freqs only
 %ylim([30 200]); %zlim([0 0.035]); caxis([0 0.035]); %high freqs only
 ylim([0 200]); 
