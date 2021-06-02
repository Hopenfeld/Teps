function snrb = getpairqual(pksx,pksxa,d2mr,rrlim)

pinds=find(ismember(pksxa,pksx));
d2mrp=abs(d2mr(pinds));


for kk=1:length(pksx)-1
    
   opks = pksx(kk+1:end);
   dp = opks - pksx(kk);
   gd = find(dp > rrlim(1) & dp < rrlim(2));
   if isempty(gd)
    snrb(kk,:)=0;
   end
   for jj = 1:length(gd)
      pksina = find(pksxa > pksx(kk) + 60 & pksxa < opks(gd(jj)) - 60);
      pksin = kk+1:kk+gd(jj)-1;
      betwamps=sort(abs(d2mr(pksina)),'descend');
      betwamps=betwamps(1:min(5,length(betwamps)));
      m1=mean(betwamps);
      prat=d2mrp(kk)/d2mrp(kk+gd(jj));
      if m1 < d2mrp(kk) && m1 < d2mrp(kk+gd(jj)) && prat < 10 && prat > 1/10
        snrl=min([d2mrp(kk),d2mrp(kk+gd(jj))])/m1;
        snrb(kk,kk+gd(jj))=min(snrl,4);
        
      else
        
        snrb(kk,kk+gd(jj))=0;

      end
         
   end
 
    
end

[nr,nc]=size(snrb);
snrb=[snrb zeros(nr,length(pksx)-nc)];
  
