function [tmsall,rrs,lhood,nmseq,skips] = getchanseqsnr(pksin,pksgl,rru,srchtol,rthold,rrlim)
tmsall=[];
lhood=[];
rrs=[];
nmseq=[];
skips=[];
pksgl=pksgl(find(~isnan(pksgl)));
  
 if rru > 2000 %shouldn't be here
   rru1=rru/2;
   pksout1 = fillgapsmc(pksgl,pksin,rru1);
   pksout=pksout1; 
 else
   pksout = fillgapsmc(pksgl,pksin,rru);  %find peaks that "fit" within gaps of parent sequence
  end

  if ~isempty(pksout)
    seqadd=permsearch(pksout,rrlim(1),srchtol,0);
    pksglr=repmat(pksgl,[length(seqadd(:,1)),1]);
    tmsall=sort([pksglr seqadd],2);

  else
    return
  end
  
    nmseq=sum(tmsall > 0,2);
    for ii=1:length(tmsall(:,1))

      [lhood(ii),rrs(ii),skips(ii)]=getmaxscore(tmsall(ii,:),rru);
      %getmaxscore computes temporal regularity measure for each sequence

    end

    gds=find(lhood > rthold(nmseq));
    tmsall=tmsall(gds,:);
    lhood=lhood(gds);
    rrs=rrs(gds);
    nmseq=nmseq(gds);
    skips=skips(gds);


endfunction
