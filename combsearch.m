function [seqout,rro] = combsearch(seqin,rrin,chtol)

stpk=seqin(1);
  frctol = linspace(0.1,0.3,7);
  rru=rrin;
  stseq=2;
  stp=0;
  stppk=6;
  ctr=0;
while ~stp
  seqin(stseq:end);
  lrat=(seqin(stseq:end)-stpk)/rru;  
  frc=round(lrat);
  rfrc=abs(lrat-frc);


  %[seqin(2:end)' rfrc' frctol(max(1,frc))'] 
  frctol(min(max(1,frc),length(frctol)));

  gds=find(rfrc < frctol(min(max(1,frc),length(frctol))) & frc > 0 & frc < stppk);
  
  if max(frc) > stppk
    agflg=1;
  else
    agflg=0;
  end
  
  if isempty(gds)
    seqr=[];
    return
  end
  
  frc=frc(gds);
  nf=length(frc);
  mf=max(frc);
  inds=sub2ind([mf nf],frc,[1:nf]);
  mp=zeros(mf,nf);
  mp(inds)=1;
  sequ=seqin(stseq-1+gds);
  nm=sum(mp == 1,2);
  nm=nm(find(nm));
  ivec=[0 cumsum(nm)'];
  nr=prod(nm);
  pm=1;
  seqr=[];
  for kk=1:length(nm)
    
    cv=sequ(ivec(kk)+[1:nm(kk)]);
    nkk=length(cv);
    madd=repmat(cv,pm,1);
    madd=repmat(madd(:),nr/length(madd(:)),1);
    seqr(:,kk)=madd(:);
    pm=pm*nkk;
    
  end

  

dd=seqr(:,2:end)-seqr(:,1:end-1);
pt1=dd(:,2:end);
pt2=dd(:,1:end-1);
rmat=pt1./pt2;
rmat(find(isnan(rmat)))=1;
frcs=[1/4 1/3 1/2 2/3 3/4 4/3 3/2];
rmatint=round(rmat);
da=abs(rmatint-rmat);
da1=abs(rmat-frcs(1));
da2=abs(rmat-frcs(2));
da3=abs(rmat-frcs(3));
da4=abs(rmat-frcs(4));
da5=abs(rmat-frcs(5));
da6=abs(rmat-frcs(6));
skips=sum(max(0,rmatint-1),2);
rs=prod(da < chtol | da1 < chtol | da2 < chtol | da3 < chtol | da4 < chtol | da5 < chtol | da6 < chtol,2);
gdrs=find(rs);
seqr=seqr(gdrs,:);
scr=10*skips(gdrs)+sum(da(gdrs,:),2);
[dum,mind]=min(scr);

if length(seqr(:,1)) > 1
  seqr = seqr(:,~all(isnan(seqr)));
end

dd=seqr(:,2:end)-seqr(:,1:end-1); 

 if agflg
   seqa=seqr(mind,:);
   seqa=seqa(find(seqa));
   stpk=seqa(end);
   stseq=gds(end)+1;
   rru=seqa(end)-seqa(end-1);
   rru=rru/round(rru/rrin);
   stp=0;
   ctr=ctr+1;
 else 
   stp=1;
   if exist('seqa','var')
     seqout=[seqa seqr(1,:)];
   else
     seqout=seqr(1,:);
   end
   rro=seqa(end)-seqa(end-1);
   rro=rro/round(rro/rrin);
 end
 
end  
  
endfunction
return
%greedy algo search; tmslast is high quality prior segment sequence
    rrp=tmslast(end)-tmslast(end-1);

      rrtst=(pksxadjn-pksxadjn(1));
      seqtst=abs(rrtst-rrp);
      seqp=[];

      kpg=1;
      while kpg
        [mt,mtind]=min(seqtst);
        if mt < 40
          seqp(end+1)=pksxadjn(mtind);
          rrp=rrtst(mtind);
          pksxadj(1:mtind)=NaN;
          rrtst=(pksxadjn-pksxadjn(mtind));
          seqtst=abs(rrtst-rrp);
        else
          kpg=0;
        end
      end