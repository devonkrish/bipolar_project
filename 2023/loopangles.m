
clear F; f=1;
F(f)=getframe(gcf); f=f+1; 
anglemin=-180;
anglemax=anglemin+45;

m=M; Md=Mbp_distance; Ma=Mbp_angle; %make a copy (only do this once) in case you want to repeat later with different angle range (start with next line)
    
angstep=5; %angle step in degrees

for i=1:360/angstep

    M=m; Mbp_distance=Md; Mbp_angle=Ma; 
    ww=0;
    for c1=1:size(Mbp_distance,1);
    for c2=1:size(Mbp_distance,2);
      if c1==c2; 
          Mbp_distance(c1,c2)=nan; 
          Mbp_angle   (c1,c2)=nan;
      end %remove diagonal (data for referential recording) 
      if ~([rad2ang(Mbp_angle(c1,c2))>(anglemin    ) && rad2ang(Mbp_angle(c1,c2))<(anglemax    )] ||...
           [rad2ang(Mbp_angle(c1,c2))>(anglemin-360) && rad2ang(Mbp_angle(c1,c2))<(anglemax-360)] ||...
           [rad2ang(Mbp_angle(c1,c2))>(anglemin+360) && rad2ang(Mbp_angle(c1,c2))<(anglemax+360)]      );
        Mbp_distance(c1,c2)=nan;
        Mbp_angle(c1,c2)=nan;
        disp([num2str(c1) ' ' num2str(c2)])
        ww=ww+1;
      end
    end
    end
    ww
    
%% Unstack and line up into 2D matrix to create channels^2 X frequencies for easier indexing-->binning
mb=[]; 
mARb=[]; 
binz=0:binsz:85; % can change binsz PRN at this point to test different bin resolutions on the plots below
clear mx my mz  
nbinz=length(binz)-1;
parfor w=windowstocheck; disp(num2str(w)); % parfor here
    Mflat=[];
      Mflatbp_distance=[];
      Mflatbp_angle=[];
    Mrefaveflat=[];
    for c2=1:nchtocheck; 
        Mflat=[Mflat; squeeze(M(:,c2,:,w))];               
        Mflatbp_distance=[Mflatbp_distance; Mbp_distance(:,c2)]; % corresponding distance index
        Mflatbp_angle   =[Mflatbp_angle;    Mbp_angle(:,c2)]; % corresponding angle index
        Mrefaveflat=[Mrefaveflat; squeeze(Mrefave(:,c2,:,w))];              
    end
    
    %bin by distance and take the mean, creating: frequency X binned distance
    for i=1:nbinz 
        mb(:,i,w)=nanmean(Mflat(Mflatbp_distance>=binz(i) & Mflatbp_distance<binz(i+1),:),1);
   %      binindex_min_ltnmax(i,:)=[binz(i) binz(i+1)];
        mARb(:,i,w)=nanmean(Mrefaveflat(Mflatbp_distance>=binz(i) & Mflatbp_distance<binz(i+1),:),1);
    end
    disp([num2str(round(windowstocheck(w)/windowstocheck(end)*100,1)) '% of windows'])
end; disp('Done')



% Notes:
% mb is frequency X binned distance X windows, and you can average across windows
% binindex_min_ltmax tells you for each column of mb what was the 
%         1] minimum (>0) distance, and
%         2] the "less than max" (<) distance
%         that were used to index bipolar pairs for that bin
%         Note: first bin includes zero, which corresponds to
%         the bin containing the referential channels




%% zscore the log transformed power according to frequency
mbm=squeeze(mean(log(mb),3)); 
mARbm=squeeze(mean(log(mARb),3)); 
 nfrx=length(frx);
 mbm_z=mbm; mARbm_z=mARbm;
 for i=1:nfrx; nns=~isnan(mbm(i,:)); 
     mbm_z(i,nns)=zscore((mbm(i,nns))); 
     mARbm_z(i,nns)=zscore((mARbm(i,nns))); 
 end; 

 %
 figure('color','w','position',[[54 223 1907 1102]]); colormap(parula); 
 toplot=mean((mbm_z),3);  % mb or mb_z
 cm_distance=flipud(cmocean('deep',1+ceil(max(max(Md)))));

% Plot confusion matrix of bipolar distances between all pairs
 subplot(8,3,1:3:7); pcolorjk(Mbp_distance); shf; hold on; %plot([1 1; 1 256]',[1 256;256 256]','k-','linewidth',1); plot([1 1 1 128 256],[1 128 256 256 256],'k.'); %plot([1 128 256 256 256]/2,[1 1 1 128 256]/2,'k.'); 
    axis equal off; set(gca,'ydir','normal','xdir','reverse'); colormap(gca,cm_distance); colorbar('fontsize',12); title('Bipolar pair distance (mm)','fontsize',14,'fontweight','normal')
    caxis([0 80]); 
% Plot confusion matrix of bipolar distances between all pairs
 subplot(8,3,2:3:8); pcolorjk(rad2deg(Mbp_angle)); shf; hold on; %plot([1 1; 1 256]',[1 256;256 256]','k-','linewidth',1); plot([1 1 1 128 256],[1 128 256 256 256],'k.'); %plot([1 128 256 256 256]/2,[1 1 1 128 256]/2,'k.'); 
    axis equal off; set(gca,'ydir','normal','xdir','reverse'); colormap(gca,hsv(360)); caxis([-180 180]); colorbar('fontsize',12); title('Bipolar pair angle (deg)','fontsize',14,'fontweight','normal')

% Histogram of number of pairs per bin
 subplot(8,2,7); histogram(make1d(Mbp_distance),0:binsz:85,'facecolor',.5*[1 1 1]); set(gca,'fontsize',12); xlabel('Binned bipolar distance (mm)','fontsize',14); 
 ylabel('Counts/bin','fontweight','normal'); axis tight; grid on; cb=colorbar; set(cb,'visible','off'); xlim(xldist)

% Histogram of angles
 subplot(8,6,22); angs=make1d(Mbp_angle); nns=~isnan(angs);
 polarhistogram(angs(nns),-pi:2*pi*(5/360):pi,'facecolor',.5*[1 1 1]); title('Bipolar angle distribution (5^o steps)','fontsize',14,'fontweight','normal'); 
 set(gca,'rtick',[min(get(gca,'rtick')) maxrt/2 maxrt],'ThetaTick',[0 90 180 270],'fontSize',9); 

%  % line plot of same as below, zscored
%  subplot(2,3,3); hold on; cmo=fliplr(cmocean('thermal',nfrx+round(nfrx)*.1)); colormap(gca,cmo); cb=colorbar; set(cb,'ticks',0:.25:1,'ticklabels',0:50:200,'fontsize',12)
%  for i=nfrx:-1:1; nns=~isnan(toplot(i,:));
%      plot(binz(2:size(toplot,2)+1),toplot(i,:),'-','color',cmo(i,:),'linewidth',1); 
%  end; grid on; axis tight; xlim(xldist)
%  set(gca,'fontsize',14); ylabel('ln(power)'); xlabel('Bipolar distance (mm)','fontsize',14); xlim(xldist)
%  text(max(xlim)+diff(xlim)/4,mean(ylim),'Frequency (Hz)','fontsize',12,'rotation',90,'horizontalalignment','center')

 % plot illustration of location of bipolar pairs, colored by same distance scale
 subplot(2,3,3); hold on; 
 elecsbrain(pt,0,[1:nchtocheck],[0 0 0],'l',0,5,2); alpha 0.05; litebrain('r',0); zoom(1.5)
 for c1=1:size(Mbp_angle,1)
   for c2=1:size(Mbp_angle,2)
       if ~isnan(Mbp_angle(c1,c2))
            plot3([em(c1,1) em(c2,1)],[em(c1,2) em(c2,2)],[em(c1,3) em(c2,3)],'-','color',cm_distance(1+round(Mbp_distance(c1,c2)),:),'LineWidth',.5)%,'LineWidth',.75)
       end
   end
 end
 
 

 subplot(2,2,3); 
 pcolorjk(binz(2:size(toplot,2)+1),frx,toplot); shading flat; set(gca,'ydir','normal'); ylabel('Frequency (Hz)'); xlabel('Distance (mm)'); set(gca,'fontsize',14); colorbar; 
 %text(max(xlim)+diff(xlim)/4,mean(ylim),'ln(power)','fontsize',12,'rotation',90,'horizontalalignment','center')
 title({'ln(power), z-scored by frequency',''},'fontweight','normal')
 set(gca,'yscale','log','ytick',ft,'yticklabel',ftl); xlim(xldist); caxis([-1 1]*(max(abs(caxis)))); colormap(gca,cmocean('balance')); %cm=cmocean('balance',100); colormap(gca,cm(:,[2 1 3]))
 caxis([-4.5 4.5]);
%  subplot(2,2,3); 
%  pcolorjk(frx,binz(2:size(toplot,2)+1),toplot'); shading flat; set(gca,'xdir','normal'); xlabel('Frequency (Hz)'); ylabel('Distance (mm)'); set(gca,'fontsize',14); colorbar; 
%  text(max(ylim)+diff(ylim)/4,mean(xlim),'ln(power)','fontsize',12,'rotation',90,'horizontalalignment','center')
%  title('(z-scored by frequency)','fontweight','normal')
%  set(gca,'xscale','log','xtick',ft,'xticklabel',ftl); ylim(xldist)

 subplot(2,2,4)
 [mx,my]=meshgrid(binz(2:size(toplot,2)+1)+binsz/2,frx);
 surf(mx,my,(toplot)); xlabel('Bipolar distance (mm)','fontsize',14); ylabel('Frequency (Hz)','fontsize',14); zlabel({'ln(power),','z-scored by frequency'},'fontsize',14)
 view(140,25); set(gca,'xdir','reverse','ydir','reverse','ytick',ft,'yticklabel',ftl,'fontsize',14); 
 colormap(gca,cmocean('balance')); caxis([-1 1]*(max(abs(caxis)))); 
 xlim([binsz/2 85]);  axis tight; 
 set(gca,'yscale','log'); 
 xlim(xldist)
 %set(gca,'yscale','linear'); set(gca,'zscale','log'); set(gca,'zscale','linear')
 %ylim([0 30]); %zlim([.035 3]); caxis([.035 3]); %low freqs only
 %ylim([30 200]); %zlim([0 0.035]); caxis([0 0.035]); %high freqs only
 ylim([0 200]); 
 zlim([-4.5 4.5]);
%  if doanglerange; 
     title([num2str(anglemin) ' to ' num2str(anglemax)]); 
%  end

F(f)=getframe(gcf); f=f+1; 
pause(3)
close

anglemin=anglemin+angstep; 
anglemax=anglemax+angstep;


end

vidfilename=['/Users/jkleen/Desktop/loopangles_video_' pt];
v=VideoWriter(vidfilename,'MPEG-4'); 
v.FrameRate = 15; 
open(v); 
writeVideo(v,F); 
close(v); 
