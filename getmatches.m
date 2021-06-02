function [mtch,ddo] = getmatches(pksf,tmn,poset,samprate)
  
  if ~exist('samprate','var')
    samprate=1000;
  end
  if ~exist('poset','var')
    poset=[0:50];
  end
    
  penfunc=[ones(1,round(15*samprate/1000)) zeros(1,1000)];
  mind=find(tmn> 0 & tmn < 5000);
  tmnu=tmn(mind);
  mxmtch=0;
  for kk=[1:length(poset)]
  for ii=1:length(pksf(:,1))
    gdsl=find(pksf(ii,:) > 0 & ~isnan(pksf(ii,:)));

    mtchl=0;
    for jj=1:length(gdsl)

      ddr=pksf(ii,gdsl(jj))-int32(tmnu-poset(kk));
      [dd,mind]=min(abs(ddr));
      dd=max(1,dd);
      mtchl=mtchl+penfunc(min(dd,length(penfunc)));
 
    end
    mtch(kk+1,ii)=mtchl;
    ddo(kk+1,ii)=ddr(mind);
  end
    
end
[dum,ind]=max(sum(mtch,2));
mtch=mtch(ind,:);
ddo=ddo(ind,:);
endfunction