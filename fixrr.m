function [rrfixo,stdv,ns,tmsout] = fixrr(rruse,rrin,tmseqin)
stdv=0;
rrfixo=0;
ns=0;
tmsout=[];


    if length(rrin) <2 || rruse == 0
        return
    end
    if rruse >= 1100 && rruse < 1500
        rrfixo=rruse/2;
    elseif rruse >= 1500
        q1 = getqualrr(rruse/2,rrin);
        q2 = getqualrr(rruse/3,rrin);
        if q2 < q1
            rrfixo=rruse/3;
        else
            rrfixo=rruse/2;
        end
    elseif rruse >= 1000 && rruse < 1100
            q1 = getqualrr(rruse,rrin);
            q2 = getqualrr(rruse/2,rrin);
            if q2 < q1
                rrfixo=rruse/2;
            else
                rrfixo=rruse;
            end
     else
            rrfixo=rruse;
    
    end


if exist('tmseqin','var') 

  tmsout=find(tmseqin);
  
numseq=length(rrin)+1;
qrat=rrin/rrfixo;
clints=round(qrat);
ds=abs(clints-qrat);
sds=qrat(1:end-1)+qrat(2:end);
odso=ds > 0.15;
ods=odso(1:end-1) .* odso(2:end);
sds=abs(round(sds)-sds);
sdsind=sds < 0.05;
bds=find(sdsind);
bds1=find(ods);
bds=intersect(bds,bds1)+1;
if odso(1)==1 && sdsind(1) == 0
  bds=[bds 1];
end
if odso(end)==1 && sdsind(end) == 0
  bds=[bds numseq];
end
tmsout=tmseqin(setdiff(tmsout,bds));
end

    
    
function qual = getqualrr(rruse, rrin)

  lrat=rrin/rruse;
  qual=sum(abs(lrat-round(lrat)));
  