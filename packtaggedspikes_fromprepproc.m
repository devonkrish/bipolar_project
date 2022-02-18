

% paths
upath=regexp(userpath,'/','split'); upath=['/' upath{2} '/' upath{3} '/'];
drivecheck=regexp(upath,'/','split'); if strcmp(drivecheck{end-1},'jonkleen'); jk='Drive1/'; else; jk=''; end %get volume address if using laptop
cdn=cd; cd([upath 'Desktop/'])
datapath=['/Volumes/' jk 'KLEEN_DRIVE/AN/TRANSFORMED DATA/']; %
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
        for b=1:nblocks; disp(q.blocks{b})
                    pt={}; bl={}; ty={}; stracesallch={};
            nspks=size(q.Sdata(b).spikes,1);
            for i=1:nspks; 
                pt{i,1}=pts{p}; 
                bl{i,1}=q.blocks{b}; 
                ty{i,1}=q.typ; 
                stracesallch{i,1}=squeeze(q.Sdata(b).StracesAllCh(i,:,:));
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
            load(['/Volumes/' jk 'KLEEN_DRIVE/AN/DATA/' pts{p} '/' pts{p} '_' q.blocks{b} '/MAT files/' pts{p} '_' q.blocks{b} '_notched.mat']);
            tdata_block  =  tdata(strcmp(tinfo(:,1),pts{p}) & strcmp(tinfo(:,2),q.blocks{b}),:);
                tdata_block(:,7:20)=round(tdata_block(:,7:20)*sfx); 
                ntr=size(tdata_block,1);
            % load artifact, if occurred in block
            fn_arti=['/Volumes/' jk 'KLEEN_DRIVE/AN/DATA/' pts{p} '/' pts{p} '_' q.blocks{b} '/' pts{p} '_' q.blocks{b} '_ARTI.mat'];
            if exist(fn_arti)==2; load(fn_arti); arti_block=round(Events(:,1:2)*(512/dsfx)); else; arti_block=[];
            end
            
            hasspkvec=   false(1,size(d,1)); hasspk=[];
            hasartivec=  false(1,size(d,1)); hasarti=[];
            hasstimvec=  false(1,size(d,1)); hasstim=[];
            hasspeechvec=false(1,size(d,1)); hasspeech=[];
            
            for i=1:nspks; hasspkvec(q.Sdata(b).spikes(i,1):q.Sdata(b).spikes(i,2))=true; % spikes
            end
            for i=1:size(arti_block,1); 
                hasartivec(arti_block(i,1):arti_block(i,2))=true; % artifact
            end
            for i=1:ntr; if t~=2; hasstimvec(tdata_block(i,7):tdata_block(i,18))=true; end % during audio stimulus
                                  hasspeechvec(tdata_block(i,19):tdata_block(i,20))=true;  % during speech
            end
            
            nonspks_windows={};
            windowsize=512;
            idx=1; 
            for i=1:512:size(d,1)-windowsize %mod(size(d,1),windowsize)+1
                nonspks_windows{1,idx}=i; %timepoint where the window starts
                nonspks_windows{2,idx}=d(i:i+windowsize-1,:)';
                hasspk   (idx)=any(hasspkvec   (i:i+windowsize-1));
                hasarti  (idx)=any(hasartivec  (i:i+windowsize-1));
                hasstim  (idx)=any(hasstimvec  (i:i+windowsize-1));
                hasspeech(idx)=any(hasspeechvec(i:i+windowsize-1));
                idx=idx+1;
            end
            info=[pts{p} '_' q.blocks{b} '_baselineWindows.mat'];
            %save(['/Volumes/KLEEN_DRIVE/David/Bipolar project/baseline-high-density-data/' info(1:end-4) '_jk.mat'],'nonspks_windows','info','hasspkvec','hasspk','hasartivec','hasarti','hasstimvec','hasstim','hasspeechvec','hasspeech','-v7.3');
            clear d dsfx fn sfx tdata_block ntr windowsize idx i nonspks_windows info hasspkvec hasspk hasstimvec hasstim hasspeechvec hasspeech hasartivec hasarti
            
            %%
            clear nspks pt bl ty stracesallch
            
        end
      end
    end 
    
    
    spikes_chidx{p}=spikesCHIDX; clear spikesCHIDX %has different numbers of channels for each patient
    Straces_allch=[Straces_allch; StracesAllCh]; clear StracesAllCh 
    badchidx{p}=unique(badCHIDX); clear badCHIDX 
    
end

save([savepath 'taggedspikes'],'Spt','Sblock','Stype','spikes','Stimeofpeakinblock','Stimeofpeakinblock_at_sfx','Str','Selem','Ssamplesfromelem','Stinfo','Strelements','SelemCOLUMNINFO','trSelem','trSelemIDX','trSnum','tinfoS','tdataS','spikes_chidx','Straces_allch','badchidx','pts','-v7.3')
