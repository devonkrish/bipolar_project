

% paths
upath=regexp(userpath,'/','split'); upath=['/' upath{2} '/' upath{3} '/'];
drivecheck=regexp(upath,'/','split'); if strcmp(drivecheck{end-1},'jonkleen'); jk='Drive1/'; else; jk=''; end %get volume address if using laptop
cdn=cd; cd([upath 'Desktop/'])
datapath=['/Volumes/' jk 'KLEEN_DRIVE/AN/TRANSFORMED DATA/']; %
rawdatapath=['/Users/jkleen/Documents/ChangLabServer/data_store1/human/prcsd_data/']; %
paramspec='SPKS h referential/';
savepath='/Volumes/KLEEN_DRIVE/David/Bipolar project/';

[tinfo,tdata,tinfo_vars,tdata_vars,T,TRi]=get_ANTask_exceldata;


pts=unique(tinfo(:,1)); 
    npt=length(pts);
    types={'AN','REC','REP'};
  
    Spt={};
    Sblock={};
    Stype={};
    spikes=[];
    spikesCHIDX=[];
    Stimeofpeakinblock=[];
    Stimeofpeakinblock_at_sfx=[];
    trSelem=[];
    trSelemIDX=[];
    trSnum=[];
    tinfoS=[];
    tdataS=[];
    
    Str=[];
    Selem=[];
    Ssamplesfromelem=[];
    Stinfo=[];
    Strelements=[];
    
    Straces_allch={};
    
    sfx=512;
    [blp,alp]=butter(2,[255]/(sfx/2),'low'); %lowpass below nyquist for downsampling target 512 Hz new sampling frequency
    
for p=1:length(pts); disp(['processing ' pts{p}])
    ptdir=[datapath paramspec pts{p}];
    
    tridx  =  strcmp(tinfo(:,1),pts{p});
    tinfo_pttyp=tinfo(tridx,:);
    tdata_pttyp=tdata(tridx,:);
    
    spikesCHIDX=[];
    StracesAllCh=[];
    badCHIDX=[];
    
    for t=1:3
        typ=types{t};
      if exist([ptdir '/' typ '/'])==7; 
        load([ptdir '/' typ '/' pts{p} '_' typ '_allch.mat'])
        
        nblocks=length(q.Sdata);
        for b=1:nblocks; ptbl=[pts{p} '_' q.blocks{b}]; disp(ptbl)
                    
                %%%%%%% pull raw data
                    rawfolder=[rawdatapath pts{p} '/' pts{p} '_' q.blocks{b} '/RawHTK'];
                    if rawfolder==0; d=[]; dsfx=[]; foldr=0;return; end
                        if strcmpi(ptbl,'EC131_B13'); rawfolder='/Users/jkleen/Documents/ChangLabServer/data_store1/human/prcsd_data/EC131/EC131_13/RawHTK'; end %this one is mis-named on server, missing the "B" in B13
                        if strcmpi(ptbl,'EC135_B5'); rawfolder='/Users/jkleen/Documents/ChangLabServer/data_store1/human/prcsd_data/EC135/EC135_b5/RawHTK'; end %this one is mis-named on server, small "b"
                    cd(rawfolder)
                    files=dir(rawfolder);
                    files={files.name};
                    w=0;
                    for i=1:length(files)
                        if ~isempty(strfind(files{i},'.htk'))
                            w= w+1;
                            wavs(w)=files(i);
                        end
                    end
                    nchan=w;
                    pathsplit=regexp(rawfolder,'/','split');
                    dirname=[]; audfold=[];
                    for i=2:length(pathsplit)-2
                        dirname=[dirname '/' pathsplit{i}];
                    end
                    dirname=[dirname '/'];
                    audfold=[dirname,pathsplit{i+1},'/Analog/'];
                    recordingblocksize=64; if length(wavs)<64; recordingblocksize=16; end
                    nrecordingblocks=nchan/recordingblocksize;
                    i=0; clear raw
                        for n=1:nrecordingblocks %4 for humans
                            for k=1:recordingblocksize %64 for humans
                                i=i+1;i
                                if i==129;
                                    1;
                                end
                                rawpath=fullfile(rawfolder,['Wav' num2str(n) num2str(k) '.htk']); pns{i}=['Wav' num2str(n) num2str(k) '.htk'];
                                if ~logical(exist(rawpath)); qw=strfind(rawpath,'Wav'); if ~isempty(qw); rawpath(qw:qw+2)='Raw'; end; end
                                RawECOG=readhtk(rawpath)*1e6;
                                if i==1; 
                                    [~,dsfx]=readhtk(rawpath); 
                                    gh=1:length(RawECOG);
                                    gh=gh(1:(dsfx/512):end); %downsample to 512 Hz
                                    raw=nan(nchan,length(gh)); clear gh
                                end
                                RawECOG = filtfilt(blp,alp,RawECOG); %lowpass below nyquist
                                RawECOG=RawECOG(1:(dsfx/512):end); %downsample to 512 Hz
                                raw(i,:)=RawECOG-mean(RawECOG); %zero center
                                %raw(i,:)=notchjk(raw(i,:),512,60,0); %notch filter
                                clear RawECOG
                            end
                        end
                    
                        
                [em,~,~]=getelecs(pts{p},2);
                nch=size(em,1);
                if size(raw,1)>nch; raw(1+size(em,1):end,:)=[]; end; %clip empty channels
                
                raw=raw'; %transpose
                
                raw(:,q.blocks_badch{b})=nan; %make bad channels NaNs
                    
                    if strcmp(ptbl(1:8),'EC129_B3'); raw(:,[193:256 65:128])=raw(:,[65:128 193:256]); end %Wav2 and Wav4 were switched during recording for B34 and B35 (looks for 30-something)
                    
                    % CAR median
                    if ~strcmp(pts{p},'EC130'); CAR=nanmedian(raw,2); raw=raw-repmat(CAR,1,size(raw,2)); end %CAR median unless EC130 because not enough channels (risk for more artifact not less)
            
                    clear rawfolder files w wavs nchan pathsplit dirname audfold recordingblocksize nrecordingblocks si rawpath n k i RawECOG
                %%%%%%%
                
                
            pt={}; bl={}; ty={}; stracesallch={};
            nspks=size(q.Sdata(b).spikes,1);
            for i=1:nspks; 
                pt{i,1}=pts{p}; 
                bl{i,1}=q.blocks{b}; 
                ty{i,1}=q.typ; 
                %stracesallch{i,1}=squeeze(q.Sdata(b).StracesAllCh(i,:,:));
                
                
                if q.Sdata(b).Stimeofpeak_at_sfx(i)-256>=1 & q.Sdata(b).Stimeofpeak_at_sfx(i)+256<=size(raw,1) 
                   stracesallch{i,1}=raw(q.Sdata(b).Stimeofpeak_at_sfx(i)-256:q.Sdata(b).Stimeofpeak_at_sfx(i)+256,:)';
                else
                    stracesallch{i,1}=nan(sfx+1,nch)'; % % if spike and surrounding window are beyond edge of data, just make this spike into NaNs
                end
            end
            
            
            Spt=[Spt;pt];
            Sblock=[Sblock;bl]; 
            Stype=[Stype;ty];
            
            spikes=             [spikes              ;   q.Sdata(b).spikes];
            Stimeofpeakinblock=        [Stimeofpeakinblock         ;   q.Sdata(b).Stimeofpeak];
            Stimeofpeakinblock_at_sfx= [Stimeofpeakinblock_at_sfx  ;   q.Sdata(b).Stimeofpeak_at_sfx];
            trSelem=            [trSelem             ;   q.Sdata(b).trSelem];
            trSelemIDX=         [trSelemIDX          ;   q.Sdata(b).trSelemIDX];
            trSnum=             [trSnum              ;   q.Sdata(b).trSnum];
            tinfoS=             [tinfoS              ;   q.Sdata(b).tinfoS];
            tdataS=             [tdataS              ;   q.Sdata(b).tdataS];
            
            spikesCHIDX=        [spikesCHIDX         ;   q.Sdata(b).spikesCHIDX];
            StracesAllCh=       [StracesAllCh        ;   stracesallch];
            badCHIDX=           [badCHIDX                q.blocks_badch{b}];
            
            Str=                [Str              ;   q.Sdata(b).Str];
            Selem=              [Selem            ;   q.Sdata(b).Selem];
            Ssamplesfromelem=   [Ssamplesfromelem ;   q.Sdata(b).Ssamplesfromelem];
            Stinfo=             [Stinfo           ;   q.Sdata(b).Stinfo];
            Strelements=        [Strelements      ;   q.Sdata(b).Strelements];
            SelemCOLUMNINFO=    q.Sdata(b).SelemCOLUMNINFO;
            
            
            %% load and pack baseline windows directly from original files
            %load(['/Volumes/' jk 'KLEEN_DRIVE/AN/DATA/' pts{p} '/' pts{p} '_' q.blocks{b} '/MAT files/' pts{p} '_' q.blocks{b} '_notched.mat']);
            
            
            
            
            tdata_block  =  tdata(strcmp(tinfo(:,1),pts{p}) & strcmp(tinfo(:,2),q.blocks{b}),:);
                tdata_block(:,7:20)=round(tdata_block(:,7:20)*sfx); 
                ntr=size(tdata_block,1);
            % load artifact, if occurred in block
            fn_arti=['/Volumes/' jk 'KLEEN_DRIVE/AN/DATA/' pts{p} '/' pts{p} '_' q.blocks{b} '/' pts{p} '_' q.blocks{b} '_ARTI.mat'];
            if exist(fn_arti)==2; load(fn_arti); arti_block=round(Events(:,1:2)*(sfx/dsfx)); else; arti_block=[];
            end
            ntp=size(raw,1); %number of timepoints
            hasspkvec=   false(1,ntp); hasspk=[];
            hasartivec=  false(1,ntp); hasarti=[];
            hasstimvec=  false(1,ntp); hasstim=[];
            hasspeechvec=false(1,ntp); hasspeech=[];
            
            for i=1:nspks; hasspkvec(q.Sdata(b).spikes(i,1):q.Sdata(b).spikes(i,2))=true; % spikes
            end
            for i=1:size(arti_block,1); 
                hasartivec(arti_block(i,1):arti_block(i,2))=true; % artifact
            end
            for i=1:ntr; if t~=2; hasstimvec(tdata_block(i,7):tdata_block(i,18))=true; end % during audio stimulus
                                  hasspeechvec(tdata_block(i,19):tdata_block(i,20))=true;  % during speech
            end
            
            nonspks_windows={};
            windowsize=sfx;
            idx=1; 
            for i=1:windowsize:ntp-windowsize %mod(ntp,windowsize)+1
                nonspks_windows{1,idx}=i; %timepoint where the window starts
                nonspks_windows{2,idx}=raw(i:i+windowsize-1,:)';
                hasspk   (idx)=any(hasspkvec   (i:i+windowsize-1));
                hasarti  (idx)=any(hasartivec  (i:i+windowsize-1));
                hasstim  (idx)=any(hasstimvec  (i:i+windowsize-1));
                hasspeech(idx)=any(hasspeechvec(i:i+windowsize-1));
                idx=idx+1;
            end
            info=[pts{p} '_' q.blocks{b} '_baselineWindows.mat'];
            save(['/Volumes/KLEEN_DRIVE/David/Bipolar project/baseline-high-density-data/' info(1:end-4) '_fromraw.mat'],'nonspks_windows','info','hasspkvec','hasspk','hasartivec','hasarti','hasstimvec','hasstim','hasspeechvec','hasspeech','-v7.3');
            clear d dsfx fn tdata_block ntr ntp windowsize idx i nonspks_windows info hasspkvec hasspk hasstimvec hasstim hasspeechvec hasspeech hasartivec hasarti
            
            %%
            clear nspks pt bl ty stracesallch
            
        end % by block
      end %
    end % by type (of task)
    
    
    spikes_chidx{p}=spikesCHIDX; clear spikesCHIDX %has different numbers of channels for each patient
    Straces_allch=[Straces_allch; StracesAllCh]; clear StracesAllCh 
    badchidx{p}=unique(badCHIDX); clear badCHIDX 
    
end % by patient

save([savepath 'taggedspikes'],'Spt','Sblock','Stype','spikes','Stimeofpeakinblock','Stimeofpeakinblock_at_sfx','Str','Selem','Ssamplesfromelem','Stinfo','Strelements','SelemCOLUMNINFO','trSelem','trSelemIDX','trSnum','tinfoS','tdataS','spikes_chidx','Straces_allch','badchidx','pts','-v7.3')
