function r=plotmap()
load 'm_coasts.mat';
%ind=find(ncst(:,1)<0);ncst(ind,1)=ncst(ind,1)+360;ind=find(ncst(:,1)>359);ncst(ind,1)=NaN;
hold on;plot(ncst(:,1),ncst(:,2),'k','linewidth',1);