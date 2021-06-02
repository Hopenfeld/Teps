function [seqr,rrs,pknos] = permsearch(seq,rrmin,chtol,minpks)
  
[seqr,rrs] = permsearchloc(seq(find(seq <= 5000)),rrmin,chtol,minpks);

if minpks > 0
  sseg=find(seq > 5000);
  if ~isempty(sseg)
    [seqra,rrsa] = permsearchloc(seq(sseg),rrmin,chtol,minpks-1);
    [seqr,rrs]=mergeseqs(seqr,seqra,rrs,rrsa);
  end

end

  ns=sum(seqr > 0,2);
  bds=find(ns < minpks+1 & rrs <1000);
  seqr(bds,:)=[];
  rrs(bds)=[];

seqmap(seq)=[1:length(seq)];
seqmap(max(seq)+1)=0;
seqt=seqr;
seqt(find(isnan(seqt) | seqt ==0))=max(seq)+1;
pknos=seqmap(seqt);

function [seqr,rrs,pknos] = permsearchloc(seq,rrmin,chtol,minpks);

ns=length(seq);
pvm=bsxfun(@(a,b) bitand(a,b) ~= 0,2.^[0:ns-1],[1:2^ns-1]');
nr=length(pvm(:,1));
seqr=pvm .* repmat(seq,nr,1);
%[rr,cc]=find(seqr);
%vals=seq(cc);
%seqr1=sparse(rr,
for ii=1:nr
  ng=find(seqr(ii,:));
  seqr(ii,1:length(ng))=seqr(ii,ng);
  seqr(ii,length(ng)+1:end)=NaN;
endfor

ns=sum(seqr > 0,2);
seqr=seqr(find(ns>minpks),:);
dd=seqr(:,2:end)-seqr(:,1:end-1);

[rr,cc]=find(dd<rrmin);
seqr(unique(rr),:)=[];
dd=seqr(:,2:end)-seqr(:,1:end-1);
pt1=dd(:,2:end);
pt2=dd(:,1:end-1);
rmat=pt1./pt2;
rmat(find(isnan(rmat)))=1;
frcs=[1/4 1/3 1/2 2/3 3/4 4/3 3/2]; 
da=abs(round(rmat)-rmat);
da1=abs(rmat-frcs(1));
da2=abs(rmat-frcs(2));
da3=abs(rmat-frcs(3));
da4=abs(rmat-frcs(4));
da5=abs(rmat-frcs(5));
da6=abs(rmat-frcs(6));

rs=prod(da < chtol | da1 < chtol | da2 < chtol | da3 < chtol | da4 < chtol | da5 < chtol | da6 < chtol,2);
seqr=seqr(find(rs),:);
seqr = seqr(:,~all(isnan(seqr)));
dd=seqr(:,2:end)-seqr(:,1:end-1);
rrs=min(dd,[],2);
lwrt=find(rrs > 1500 & rrs <= 3000);
rrs(lwrt)=rrs(lwrt)/2;

function [seqr,rrs]=mergeseqs(seqr,seqra,rrs,rrsa);

lws=find(rrs >= 1000);
lwsa=find(rrsa >= 1000);
seqra=seqra(lwsa,:);
rrsa=rrsa(lwsa);

addseq=[];
addrr=[];
for ii=1:length(lws)
  ru=seqr(lws(ii),:);
  gv=find(~isnan(ru) & ru >0);
  sequ=seqr(lws(ii),gv);
  rat=(seqra(:,1)-sequ(end))/(sequ(end)-sequ(end-1));
  gds = find(abs(round(rat)-rat) < 0.3);
  if ~isempty(gds)

    seqrep=repmat(sequ,length(gds),1);
    seqrepa=[seqrep seqra(gds,:)];
    seqr(end+[1:length(gds)],1:length(seqrepa(1,:)))=seqrepa;
    rrs(end+[1:length(gds)])=rrs(lws(ii));
  end
  
end
seqr(lws,:)=[];
rrs(lws)=[];



