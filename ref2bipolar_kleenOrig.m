
% d is a matrix of samples by channels consisting of referential intracranial EEG data

pt='EC219';
[bipolarN,bipolarT]=xlsread([upath 'Box/KLEENLAB/AN Code/_AN_2021/AN_ElectrodeInfoTDT.xlsx'],pts{p});

for r=1:size(bipolarT,1) %each row of the sheet is a component (grid, strip, or depth)
    if any(strcmpi(bipolarT(r,2),{'grid','minigrid'})) %grids (2-D)
        N=bipolarN(r,3)-1; 
      for i=bipolarN(r,1):N+1:bipolarN(r,2)
        d(:,i:i+N-1)=d(:,i:i+N-1)-d(:,i+1:i+N); 
        d(:,i+N)=NaN; %last channel in the line will be NaNs
      end
    else; i=bipolarN(r,1); N=bipolarN(r,2)-i; %strips, depths (1-D)
        d(:,i:i+N-1)=d(:,i:i+N-1)-d(:,i+1:i+N); 
        d(:,i+N)=NaN; %last channel in the line will be NaNs
    end; clear N i
end
