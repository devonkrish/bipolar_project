function [trm,frx,s]=bpspectra_Linear_2023(d,sfx,frxrange,okc)
% trm is mean of natural log transform of power across windows/trials
% frx is index of frequencies
% is all spectra data (frx X channel X window) that was used to log-transform then mean to create trm

[~,nch,nwind]=size(d);
%just getting frequency index
[~,frx]=spectrogramjk_chronuxmtfft(sq(d(:,1,1)),sfx,frxrange,[.5,1],0); 

% FFT
D=d(:,okc,:);
nfrx=length(frx); %number of anticipated frequencies
spect=nan(nfrx,size(D,2),size(D,3));
L=size(D,2);
tic
parfor w=1:nwind
    W=sq(D(:,:,w))'; 
    for c=1:L %c=find(~badchAll(b,:))
        spect(:,c,w)=spectrogramjk_chronuxmtfft(W(c,:),sfx,frxrange,[.5,1],0);
    end
end; clear D W c w L
toc
s=nan(nfrx,nch,nwind);
s(:,okc,:)=spect; clear spect 
ln_s=s;
ln_s(:,okc,:)=log(s(:,okc,:)); %natural log transform of power
%             ln_s=s; %or just skip transform

trm=sq(mean(ln_s,3))'; %trial mean




