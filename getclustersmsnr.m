function [clusters,clusterprbs,pmo,pkmap,chd] = getclustersmsnr(pksin,thp,pkgp,samprate)

load pkprbmats prbpkt hvals

numchan=length(pksin(:,1));
clustersc=[];
clusters=[];
clusterprbs=[];
prbu=prbpkt(numchan,:);
ht=hvals(numchan,:);
prbd=[100 10 1 zeros(1,40)];
pmo=pkgp*ones(100,numchan-1);
pksh=pksin(find(pksin));
dum=hist(pksh,[1:5000]);
pksh(find(pksh > 870 & pksh < 920));
smh=cumsum(dum);
dd=smh(pkgp:end)-smh(1:end-pkgp+1);
chd=zeros(1,numchan);
ddd=dd(2:end)-dd(1:end-1);
gds=find(ddd);
dd=dd(gds);
pkmap=zeros(5000,1);
cla=find(dd>thp);
for ii=1:length(cla)
  if dd(min(cla(ii)+1,length(dd)))>dd(cla(ii))...
    &&  gds(min(cla(ii)+1,length(dd))) - gds(cla(ii)) < (pkgp-1)
    bad1=gds(cla(ii:ii+1));
    continue
  end
  if dd(cla(max(ii-1,1)))>thp && dd(cla(max(ii-1,1)))>dd(cla(ii)) ...
    &&  gds(cla(ii)) - gds(cla(max(ii-1,1))) < (pkgp-1)
    bad2=gds(cla(ii));   
    continue
  end
  if cla(ii)<2
    continue
  end
  cll=gds(cla(ii));
  clr=max(gds(cla(ii)-1)+pkgp,cll+pkgp);
  if (cll-gds(cla(max(ii-1,1)))) < pkgp-1
    cll=gds(cla(max(ii-1,1)));
  end
  pksc=pksh(find(pksh > cll-1 & pksh < clr+1));
  if length(pksc) > 1
    
    chans=ismember(pksin,pksc);
    [chans,cs]=find(chans);
    chd=chd+ismember([1:numchan],chans);
    pkscsrt=sort(pksc);
    dp=pkscsrt(2:end)'-pkscsrt(1:end-1)';

    if length(pkscsrt) > numchan
      [dum,minds]=sort(abs(dp));
      dp=dp(minds(1:numchan-1));
      
    end
    dp=sort(dp);
    dp(end+1:numchan-1)=pkgp;
    pval=sum(prbd(dp+1));
    pval=find(ht==pval);
    clusterprbs(end+1)=prbu(pval); 
    pmo(length(clusterprbs),1:length(dp))=dp';

  else
      clusterprbs(end+1)=0;
  end
  clusters(end+1)=round((1000/samprate)*mean(pksc));
  pkmap(round((1000/samprate)*pksc))=length(clusters);

end

allpks=round((1000/samprate)*pksh);
restpks=setdiff(allpks,find(pkmap));  %shouldn't be any

if ~isempty(restpks)
  addmap=length(clusters)+[1:length(restpks)];
  clusters(end+[1:length(restpks)])=restpks;
  pkmap(restpks)=addmap;
  clusterprbs(end+[1:length(restpks)])=0;
end
pmo=pmo(1:length(clusterprbs),:);

endfunction
