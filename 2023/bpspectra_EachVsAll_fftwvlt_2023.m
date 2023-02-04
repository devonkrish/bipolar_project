function [M,Mrefave,Mbpdist,frx]=bpspectra_EachVsAll_fftwvlt(d,sfx,frxrange,em,nchtocheck,fft0wvlt1)
% for each window in d (samples x channels x windows), take channel c1
% minus channel c2 and computes spectrum on that bipolar signal, placing
% into a confusion matrix. Also does referential spectra (on the diagonal).
% OUTPUTS:
%   M is spectra for all bipolar pairs: channels X channels X frequencies X windows 
%   Mrefave is average of the two individual spectra for each pair
%   Mbpdist is a confusion matrix of the euclideandistances for every bipolar pair (channel x channel)
%   frx is the frequency index for maz

if ~exist('fft0wvlt1','var'); fft0wvlt1=0; end % does chronux multitaper spectra if 0, and wavelet spectra if 1

%just getting frequency index
if     fft0wvlt1==0
 [~,frx]=spectrogramjk_chronuxmtfft(zeros(1,size(d,1)),sfx,frxrange,[.5,1],0); 
elseif fft0wvlt1==1
 [~,~,frx]=zWPHZ(zeros(1,size(d,1)),sfx,frxrange,0); 
end
 
 nwindtocheck=size(d,3);
 
 % Create matrix M:  channels X channels X frequencies X windows 
 M=nan(nchtocheck,nchtocheck,length(frx),nwindtocheck); 
 Mrefave=M; %copy
 MMbpdist=nan(nchtocheck,nchtocheck); % will log the euclidean distance for each of these pairs
 tic
 for w = 1:nwindtocheck; %disp([num2str(w/nwindtocheck) '%'])
     parfor c1=1:nchtocheck %parfor here
         for c2=1:nchtocheck 
             trc=sq(d(:,c1,w))-sq(d(:,c2,w));
             if ~any(isnan(trc))
                 if fft0wvlt1==0
                     M(c1,c2,:,w)=                    spectrogramjk_chronuxmtfft(trc,sfx,frxrange,[.5,1],0);
                     Mrefave(c1,c2,:,w)=exp(mean([log(spectrogramjk_chronuxmtfft(sq(d(:,c1,w)),sfx,frxrange,[.5,1],0))...
                                                  log(spectrogramjk_chronuxmtfft(sq(d(:,c2,w)),sfx,frxrange,[.5,1],0));],   2));
                 elseif fft0wvlt1==1
                     spct=zWPHZ(trc,sfx,frxrange,0); spct=squeeze(mean(spct(sfx/4:sfx-sfx/4,:)));
                       M(c1,c2,:,w)=spct; 
                     spct1=log(zWPHZ(sq(d(:,c1,w)),sfx,frxrange,0)); spct1=squeeze(mean(spct1(sfx/4:sfx-sfx/4,:)));
                     spct2=log(zWPHZ(sq(d(:,c2,w)),sfx,frxrange,0)); spct2=squeeze(mean(spct2(sfx/4:sfx-sfx/4,:)));
                       Mrefave(c1,c2,:,w)=exp(mean([spct1...
                                                    spct2;],   2));
                 end
             end  
         end
     end; toc
 end; clear c1 c2 trc spct spct1 spct2
 % include referential signal power on the diagonal for storage/plotting
 for w = 1:nwindtocheck 
     for c1=1:nchtocheck
             trc=sq(d(:,c1,w));
             if ~any(isnan(trc))
                 if     fft0wvlt1==0
                    M(c1,c1,:,w)=spectrogramjk_chronuxmtfft(trc,sfx,frxrange,[.5,1],0);
                 elseif     fft0wvlt1==1
                     spctref=zWPHZ(trc,sfx,frxrange,0); spctref=squeeze(mean(spctref(sfx/4:sfx-sfx/4,:)));
                     M(c1,c1,:,w)=spctref
                 end
             end  
     end;
 end; clear spctref

 % euclidean distance for each pair
 for c1=1:nchtocheck; for c2=1:nchtocheck; Mbpdist(c1,c2)=distance3D(em(c1,:),em(c2,:)); end; end
    Mbpdist(Mbpdist==0)=nan; %remove spurious zero values, then
    for c1=1:nchtocheck; Mbpdist(c1,c1)=0; end %referential signal denoted as "0" distance for coding/indexing purposes below
 % remove mirror image in confusion matrix
 for c1=1:nchtocheck; for c2=c1+1:nchtocheck; M(c1,c2,:)=nan; Mbpdist(c1,c2)=nan; end; end

      
