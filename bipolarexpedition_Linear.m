% BIPOLAR PAIR ANALYSIS: LINEAR

onlygrids=true;
onlydepths=false;
onlystrips=false;

cm=cool(6); cm(1,:)=[0 0 0]; 
datadir='/Volumes/KLEEN_DRIVE/David/Bipolar project/baseline-high-density-data/'; %bandpassfiltered/';
cd([datadir])
load('/Volumes/KLEEN_DRIVE/David/Bipolar project/taggedspikes.mat')
sfx=512;
frxrange=[2 200]; %frequency range to examine
  ft=[2 5 10 20 50 100 200]; ftl=cellstr(num2str(ft')); %frequency labels for plots

%which patients are ok to do
okpt=false(1,length(pts)); 
    okpt([4 12:16 19:23])=1;
% how many bipolar steps to loop through
maxbpd=5;

% for i=1:length(Spt); ptblall{i,1}=[Spt{i} '_' Sblock{i}]; end %seems like Devon forgot to include EC181_B7-9 so will use the below method instead
% uptbl=unique(ptblall); 

%% LINEAR PAIRS analysis and plots
figure(1); set(gcf,'color','w','position',[372 1 1297 1337]); 
u=dir; uptbl={}; for i=1:length(u); uname=u(i).name; uptbl{i,1}=uname(1:end-28); end; uptbl(1:2)=[]; clear i u uname
hasmat=false(maxbpd+1,length(pts));
for bpd=0:maxbpd; %bipolar distance (# of electrodes to subsample)
    
    for p=find(okpt) %[4 12:23]
        pblocks=strfind(uptbl,pts{p}); 
        for i=1:length(pblocks); 
            isbl(i,1)=~isempty(pblocks{i}); 
        end
        ptbl=find(isbl); if ~isempty(ptbl); disp(['Loading ' pts{p} ' blocks...']); end

        d=[]; nwind=0;
        for b=1:length(ptbl); disp(uptbl{ptbl(b)})
        % load using using "_jk" versions of baseline windows, updated 2/2022
        load([datadir uptbl{ptbl(b)} '_baselineWindows_fromraw.mat'])
        % get rid of baseline windows containing spikes or artifact
        spksarti=hasspk | hasarti;
        nonspks_windows(:,spksarti)=[];
        hasstim(spksarti)=[];
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

        [bpN,bpT]=xlsread(['/Users/davidcaldwell/code/high_density_ecog/AN_ElectrodeInfoTDT.xlsx'],pts{p});
        [em,eleclabels,anatomy]=getelecs(pts{p},2);

        % d is a matrix of samples by channels by trials consisting of referential intracranial EEG data
        bpdist=nan(nch,1); %will fill in euclidean distance in 3D space for each bipolar pair
        if bpd>0
         % linear grid rows, strips, depths
         for jj = 1:size(d,3) %windows
            for r=1:size(bpT,1) %each row of the sheet is a component (grid, strip, or depth)
                if any(strcmpi(bpT(r,2),{'grid','minigrid'})) %grids (2-D)
                    N=bpN(r,3)-bpd; % number of new consecutive bipolar contacts of distance "bpd"
                    if N>0
                      for i=bpN(r,1):N+bpd:bpN(r,2); %every grid row
                        c1=[i:i+N-1];
                        c2=[i:i+N-1]+bpd;
                        d(:,c1,jj)=d(:,c1,jj)-d(:,c2,jj);
                        d(:,i+N:i+N+bpd-1,jj)=NaN; %last channel in the line will be NaNs
                        for i=1:length(c1); bpdist(c1(i))=distance3D(em(c1(i),:),em(c2(i),:)); end
                      end
                    else %if bipolar spacing is longer than #electrodes in component, nan the component electrodes
                        d(:,bpN(r,1):bpN(r,2),:)=nan; bpdist(:,bpN(r,1):bpN(r,2))=nan;
                    end
                else; i=bpN(r,1); %strips, depths (1-D)
                    N=diff(bpN(r,1:2))+1-bpd; % number of new consecutive bipolar contacts of distance "bpd"
                    if N>0
                        c1=[i:i+N-1];
                        c2=[i:i+N-1]+bpd;
                        d(:,c1,jj)=d(:,c1,jj)-d(:,c2,jj);
                        d(:,i+N:i+N+bpd-1,jj)=NaN; %last channel in the line will be NaNs
                        for i=1:length(c1); bpdist(c1(i))=distance3D(em(c1(i),:),em(c2(i),:)); end
                    else %if bipolar spacing is longer than #electrodes in component, nan the component electrodes
                        d(:,bpN(r,1):bpN(r,2),:)=nan; bpdist(:,bpN(r,1):bpN(r,2))=nan;
                    end
                end; clear N i
            end
         end; clear r jj
        end

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


      if any(okc)
        [trm,frx,s]=bpspectra_Linear(d,sfx,frxrange,okc);

        % Matrices (bipolarDistance X patient X frequency) aggregated for each patient 
        TRM(bpd+1,p,:)=mean(trm(okc,:),1); % Mean
        TRSD(bpd+1,p,:)=std(trm(okc,:),[],1); %Standard deviation
        TRSE(bpd+1,p,:)=TRSD(bpd+1,p,:)/sqrt(sum(okc)); %SEM
        TRbpdist{bpd+1,p}=bpdist; %actual eucidean distances for each bp pair 

        SS{bpd+1,p}=s; 
        hasmat(bpd+1,p)=1;
        Sokc{bpd+1,p}=okc; 

        figure(1); sp(5,5,p); hold on; %each patient in their own plot
         ribbons(frx,trm,cm(bpd+1,:),.5,'sem',0,0); 
         grid on; title(pts{p}); drawnow; axis tight; %set(gca,'xscale','log'); 
         if p==4; ylabel('ln(power)'); xlabel('Frequency (Hz)'); 
             if     onlygrids;  legend({'referential','', '4mm','', '8','','12','','16','','20'},'location','ne'); 
             elseif onlydepths; legend({'referential','', '5mm','','10','','15','','20','','25'},'location','ne'); 
             elseif onlystrips; legend({'referential','','10mm','','20','','30','','40','','50'},'location','ne'); 
             end
         end
      end

      if p==4 % ECoG trace plots for increasing bipolar spacing (example patient)
        figure(5); set(gcf,'color','w','position',[1638 1 668 1344])
        sp(6,2,(bpd+1)*2-1); hold on; 
        chtoplot=49:64; %example channels [EC143-49:64, EC175-]
        windowtoplot=12; %example window
        ts=1/sfx:1/sfx:1; %timestamps for 1-sec window
        eegplotbytime2021(d(:,chtoplot,windowtoplot),512,300,[],0,[.3 .3 .3],1);
    %             for c=1:length(chtoplot); plot(ts,-c*1000+d(:,chtoplot(c),windowtoplot),'color',[0 .6 .6],'linewidth',1); end
        if ~exist('yl','var'); yl=ylim; end; ylim(yl);
        axis off
        sp(2,2,2); hold on; 
            for c=chtoplot
                [specttoplot(c-chtoplot(1)+1,:),frx]=spectrogramjk_chronuxmtfft(sq(d(:,c,windowtoplot)),sfx,frxrange,[.5,1],0);
            end; 
            ribbons(frx,log(specttoplot),cm(bpd+1,:),.3,'sem',0,0); grid on; set(gca,'xlim',frxrange) %,'xscale','log','xtick',ft,'xticklabel',ftl
            clear c specttoplot
      end; %set(gcf,'position',[1198 785 498 481]) %to resize spectra for figure



    end % patient loop
    clear d isbl ptbl pblocks s trm badchI okc nch nwind
end % bipolar spacing loop



% additional plots for Linear bipolar analysis: across all patients
% fill nans for empty (non-analyzed) patients
TRM(:,~okpt,:)=nan; 
TRSD(:,~okpt,:)=nan; 
TRSE(:,~okpt,:)=nan; 
% TRbpdist(:,~okpt)=nan;  %********

%% plots aggregated across patients
figure; set(gcf,'color','w'); %,'position',[1055 216 1217 826]
hold on
ps=nan(maxbpd,length(pts),length(frx));
for bpd=[1:maxbpd 0]
    h=find(hasmat(bpd+1,:));
    for p=h
        ps_=SS{bpd+1,p};
        ps(bpd+1,p,:)=mean(mean(log(ps_(:,Sokc{bpd+1,p},:)),3),2);
        %ps(bpd+1,p,:)=nanmean(nanmean(log(ps_),3),2);

%         ps_=nanmean(nanmean(log(ps_),3),2); ps_=(ps_-mean(ps_))/std(ps_); %z-score
%             ps(bpd+1,p,:)=ps_;
    end


    %ribbons(frx,sq(ps(bpd+1,h,:)),cm(bpd+1,:),.3,'sem',0,0); grid on; set(gca,'xlim',frxrange)
    plot(frx,mean(sq(ps(bpd+1,h,:)),1),'color',cm(bpd+1,:),'linewidth',2); 
    hold on
end; grid on; set(gca,'xlim',frxrange); ylabel('ln(power)'); xlabel('Frequency (Hz)'); set(gca,'xscale','linear')

%%
figure('color','w','position',[1000 517 354 821]); colormap(cmocean('thermal')); x=4*(0:maxbpd);
sp(4,1,4); f=frx>0&frx<=20; imagesc(x,frx(f),sq(nanmean(TRM(:,:,f),2))'); set(gca,'ydir','normal','xtick',x,'xticklabel',cellstr(num2str(x'))'); xlabel('Bipolar distance (mm)'); ylabel('Frequency (Hz)')
sp(4,1,3); f=frx>20&frx<=50; imagesc(x,frx(f),sq(nanmean(TRM(:,:,f),2))'); set(gca,'ydir','normal','xtick',[]);
sp(4,1,2); f=frx>50&frx<=100; imagesc(x,frx(f),sq(nanmean(TRM(:,:,f),2))'); set(gca,'ydir','normal','xtick',[]);
sp(4,1,1); f=frx>100&frx<=170; imagesc(x,frx(f),sq(nanmean(TRM(:,:,f),2))'); set(gca,'ydir','normal','xtick',[]); title('Average across patients')



