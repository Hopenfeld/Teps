function stats = domultisnr(sigs,trpks,samprate,moset,coradj)
%
load pkprbmats prbsnr

pksall=[]; nmsall=[]; regscall=[]; compall=[]; snrall=[]; rrsall=[];stats=[];

d=6; %number of samples between second difference points
cluslim=[4 20]; %minimun and maximum peak widths (in samples)
pkgp=6; %temporal alignment limit (in samples)
%rthold=[0 0 -2000 -4000 -6000 -8000 -10000 -12000*ones(1,10)]; %temporal regularity thresholds as function 
                                                               %based on sequence length 
rtv=[1:15];
rthold=-1000*(rtv-2);
cfs=[-0.5437 0.002004 0.001716];   %score coefficients for skips, temporal regularity and SNR/prominence
cfsp=[cfs 1]; %with pk timing probability
rrlim=[400 1300];  %rr interval limits in ms
                                                          
numchan=length(sigs(:,1));
nchs=11+5*(max(0,2-numchan)); %number of high quality of peaks in multichannel search
numcand=max(25,40-5*(4-numchan)); %number of "best" candidate peaks in each channel
datoutsc=[];
snrmat=[];
trpksadj=round(trpks*samprate/1000);  %adjust true peaks for sampling rate
regpen=inline("5*min(0,x+700) + x");
snradjrr=inline("min(0,x-550)/300 + 1");
regpenrr=inline("max(0,550-x)/300 + 1");
regbonrr=inline("1-min(0.7,max(0,x-800)/200)");


for ii=1:numchan  %process each channel

  siggu=sigs(ii,:);
  trpksgcp=trpksadj-round(moset*samprate/1000);
  [pksx,snrb,pksten,snrba,ddd,nslvlu,statspks] = getcandpks(siggu,d,samprate,cluslim,numcand,trpksgcp,rrlim);
  pksxn=round(1000*pksx/256);
  
  statsg=getmatches(pksxn,trpks,moset);
  
  %pksten
  %pksx
  %round(pksx*1000/256)
  %pksx - numcand best peaks in channel
  %snrb - matrix of SNR/prominence between peak pairs in pksx
  %pksten - possibly very high quality peaks within pksx
  %snrba - SNR/prominence matrix associated with pksten
  %ddd - second difference of channel raw signal siggu
  %nslvlu - estimated noise level of channel
  %statspks - various statistics
  nslvls(ii)=nslvlu;
  if 1==1
    [tmsallsc,nregsc,rrssc,skipssc,nmseqsc,snrssc] = singchan(pksten,snrba,rthold);

    %tmsallsc - rows are candidate sequences
    %nregsc - average temporal regularity score for each sequence
    %rrssc - average rr interval corresponding to each sequence
    %skipssc - number of missed peaks in each sequence
    %nmseqsc - number of peaks in each sequence
    %snrssc - SNR/prominence score for each sequence
    if ~isempty(nregsc)
      snradj=snradjrr(rrssc);
      regadjrr=regpenrr(rrssc).*regbonrr(rrssc);
      snrssckp=round(100*snradj(:).*snrssc);
      sfu=round(100*snrssc);
      rlhs=regpen(regadjrr.*nregsc);
      smat=[skipssc(:) rlhs snrssckp(:)]';
      gdh=find(nregsc > -800);
      scr=round(100*cfs*smat); %score the sequences
      [dum,sinds]=sort(scr,'descend');
      kp=sinds(1:min(3,length(sinds))); %keep 3 best sequences
      datadd=[rrssc(:)(kp) scr(:)(kp) rlhs(kp) skipssc(:)(kp) tmsallsc(kp,:)]; %store rr intervals, scores and sequences
      datoutsc(end+[1:length(datadd(:,1))],[1:length(datadd(1,:))])=datadd; %will be added to best multichannel sequences
      if dum(1) > 300
        stats.dat=datoutsc;
        return        
      end
    
    end

  end
  
  pksa(ii,1:length(pksx))=pksx';
  snrmat=[snrmat;snrb]; %assemble global SNR matrix in row/col/value format and 
                        %adjust according to channel correlation

end

if numchan > 1
  [spksa,pkprb,~,pkmap]=getclustersmsnr(pksa,0,pkgp,samprate);
  pkprb=pkprb*coradj;
%get global peaks spksa and associated Bayesian timing coherence probability pkprb by clustering peak times (pksa) 
%across channels; pkmap is the mapping from local channel peaks to the global peaks
else
  spksa=round(pksx*1000/samprate);
  pkprb=0*spksa;
  pkmap(spksa)=[1:length(spksa)];

end
pkprobo=pkprb;
rr=pkmap(snrmat(:,1)); %create global SNR matrices snrmat and snrmatglobal, the latter indexed by global peak time
cc=pkmap(snrmat(:,2));
mapinv([1:length(spksa)])=spksa;
snrmat=full(sparse(rr,cc,snrmat(:,3)));
[rr,cc,vv]=find(snrmat);
rrb=mapinv(rr);
ccb=mapinv(cc);
snrmatglobal=sparse(rrb(:), ccb(:), vv(:),5000, 5000);

mxc=max(snrmat,[],2); %get peak SNR/prominence values mxc and mxr from snrmat
mxr=max(snrmat);
mns=mean(nslvls)-2.5; %heuristic factor assuming estimated noise lvl too low
if isnan(mns)
  mns=0;
end

snrind=round(mns/2.5)+5;
snrprb=squeeze(prbsnr(max(numchan,2),snrind,:)); %probability for given SNR values as a function
                                          %of the number of channels
npad=length(pkprb)-length(mxc);
mxc(end+[1:npad])=0;
npad=length(pkprb)-length(mxr);
mxr(end+[1:npad])=0;
snrprbs=snrprb(max(min(length(snrprb),round(mxc(:)+mxr(:))),1)); %the the SNR probabilities
pkprb=pkprb(:);
pkqms=pkprb+snrprbs-pkprb.*snrprbs; %compute peak probabilities as union of timing probabilities
                                    %and SNR probabilities

[cval,inds]=sort(pkqms,'descend');
pkmapo=pkmap;
pksb=find(pkmap);
indmap(inds)=[1:length(inds)];
pksbr=indmap(pkmap(pksb));
pkmap(pksb)=pksbr;
spksa=spksa(inds);
pkprb=pkprb(inds);
pkqms=pkqms(inds);
prbmap(spksa)=pkprb;
snrprbs=snrprbs(inds);
pkptr=[1:nchs]; 
[pksx,inds]=sort(spksa(pkptr)); %get nchs highest probability global peaks
rpks=sort(spksa(pkptr(end)+1:end)); %rpeaks is the rest of the global peaks
prbgpks=pkprb(inds);
[tmsall,rrs,pknos] = permsearch(pksx,rrlim(1),0.25,3);  %find parent sequences based on best globabl peaks
[rr,~,gv]=find(pknos);
glprbs=accumarray(rr,prbgpks(gv)',[],@sum); %sum of global probablitites for each parent sequence; not used
nseq=sum(tmsall > 0,2); %sequence length of parent sequence

for ii=1:length(tmsall(:,1))
  
    [regsc(ii),rrs(ii),skips(ii)]=getmaxscore(tmsall(ii,:),rrs(ii)); %temp regularity, rr int. and skips for 
                                                               %each mother sequence
end

gds=find(regsc > rthold(nseq));
bds=find((rrs < 625 & nseq <5) | (rrs < 500 & nseq < 6));  %need more sequence entries 
                                                           %for lower RR to keep search reasonable
gds=setdiff(gds,bds);
if ~isempty(gds)
  tmsall=tmsall(gds,:);
  regsc=regsc(gds);
  rrs=rrs(gds);
  skips=skips(gds);
  glprbs=glprbs(gds);
  nseq=nseq(gds);
elseif ~isempty(datadd)
  tmsall=datoutsc(:,5:end);
  gdtms=find(tmsall > 0 & ~isnan(tmsall));
  tmsinds=0*tmsall;
  tmsinds(gdtms)=pkmapo(tmsall(gdtms));
  pkscmap=pkmapo(tmsall(gdtms));
  pkprsc(pkscmap)=pkprobo(pkscmap);
  tmsinds(find(isnan(tmsinds)))=0;
  [rr,~,gv]=find(tmsinds);
  glprbs=accumarray(rr,pkprsc(gv)',[],@sum);
  tmsall(gdtms)=mapinv(pkmapo(tmsall(gdtms)));
  rrs=datoutsc(:,1);
  skips=datoutsc(:,4);
  nseq=sum(tmsall>0,2);
  rpks=sort(spksa);
  datoutsc=[]; 
end
if ~isempty(tmsall)
  datout=[];
  scra=[];
  sra1a=[];
  mxsc=-1000;

  for ii=1:length(tmsall(:,1))  %loop through parent sequences
  
    tmsloc=tmsall(ii,:);
    tmsloc=tmsloc(find(tmsloc >0 & ~isnan(tmsloc)));
 
    if skips(ii) > 0 
      [tmsalla,rrsa,regsca,nmseqa,skipsa] = getchanseqsnr(rpks,tmsloc,rrs(ii),0.1,rthold,rrlim);
      %if the parent sequence has at least one skipped peak, fill in the gaps between 
      %widely spaced peaks and generate new candidate sequences tmsalla based on parent
      %sequence and peaks that fill gaps
    else
      tmsalla=tmsloc; rrsa=rrs(ii); regsca=regsc(ii); nmseqa=nseq(ii); skipsa=skips(ii);
    end
    if ~isempty(tmsalla)
      snrvs = getsnrfrtmsall(tmsalla,snrmatglobal); %for sequences that span segments, 
      snradj=snradjrr(rrsa);                                                  %need to fix gap between snrs between segs
      snrvs=round(100*snradj(:).*snrvs);
      regadjrr=regpenrr(rrsa) .* regbonrr(rrsa);
      lhs=double((regadjrr(:).*regsca(:)))./(nmseqa(:)-2);
      lhs=regpen(lhs);
      smat=[skipsa(:) lhs snrvs(:) glprbs(ii)*ones(length(lhs),1)]';
      scr1=cfsp*smat;    %peak timing coherence glprbs not used here but improve detection in certain cases
      [sscrs,sinds]=sort(scr1,'descend');
      kp=sinds(1:min(length(sinds),10)); %keep 10 best sequences per parent sequence
      scr1=round(100*scr1(kp))';
      datadd=[ rrsa(:)(kp) scr1 lhs(kp) skipsa(:)(kp) tmsalla(kp,:)];
      datout(end+[1:length(datadd(:,1))],[1:length(datadd(1,:))])=datadd;
    end
  end
end
if ~isempty(datoutsc)
 datout(end+[1:length(datoutsc(:,1))],[1:length(datoutsc(1,:))])=datoutsc;
end
if ~isempty(datout)
  datout = cullseqsloc(datout);
  %remove the lower scoring sequences in closely related groups
  stats.dat=datout;
else
  stats.dat=[];
end
endfunction

function snrvs = getsnrfrtmsall(pksin,snrmatloc)

[nr,nc]=size(pksin);
rs=pksin(:,1:end-1)';
cs=pksin(:,2:end)';
bds=find(isnan(rs) | rs ==0 | isnan(cs) | cs ==0);
rs(bds)=1; cs(bds)=1;  %map "bad" peak pairs to snrmatloc(1,1), which is always 0
inds=sub2ind(size(snrmatloc),rs,cs);
snrvs=full(snrmatloc(inds));
snrvs=reshape(snrvs,nc-1,nr);
snrvs=sum(snrvs',2);

endfunction

function [tmsall,nreg,rrs,skips,nmseq,snrs] = singchan(pksten,snrba,rthold)
pksave=[]; tmsall=[]; nreg=[]; skips=[]; nmseq=[]; snrs=[]; 
lhood=[];,rrshsnr=[];skipshsnr=[];
[tmsall,rrs,pknos] = permsearch(pksten,300,0.15,3);
if ~isempty(rrs)
  snrs=getsnrfrtmsall(pknos,snrba);
  nmseq=sum(tmsall>0,2);
    
  for ik=1:length(tmsall(:,1))
    [lhood(ik),rrs(ik),skips(ik)]=getmaxscore(tmsall(ik,:),rrs(ik));
    %getmaxscore computes temporal regularity measure for each sequence
  end
    
  gds=find(lhood > rthold(nmseq));
  nreg=double(lhood(:))./(nmseq(:)-2);
 
end

endfunction

function datout = cullseqsloc(datin)

datout=[];
rrs=unique(datin(:,1));

%tmsallcs=unique(tmsallcs(find(tmsallcs > 0 & ~isnan(tmsallcs))));
for ii = 1:length(rrs)
  clmp=find(datin(:,1)==rrs(ii) & datin(:,2) > -500);
  tmsallcs=datin(clmp,5:end);
  scrs=datin(clmp,2);
  skps=datin(clmp,4);
  hs=find(scrs > 400);
  if ~isempty(hs)
    [dum,mind]=max(scrs(hs));
    tmsmx=tmsallcs(hs(mind),:);
    tmsmx=tmsmx(find(tmsmx > 0 & ~isnan(tmsmx)));
    hsa=hs(setdiff([1:length(hs)],mind));
    tmsc=tmsallcs(hsa,:);
    nsc=sum(tmsc>0,2);
    nmtch=getmatches(tmsc,tmsmx,0,1000);
    bds=find(nsc(:)==nmtch(:));
    scrs(hsa(bds))=-Inf;
  end
    
    [dum,minds]=sort(scrs,'descend');
    mout=minds(1:min(3,length(minds)));
  

  datout(end+[1:length(mout)],:)=datin(clmp(mout),:);
  
end


  tmsallcs=datout(:,5:end);
  scrs=datout(:,2);
  hs=find(scrs > 400);
  if ~isempty(hs)
    [dum,mind]=max(scrs(hs));
    tmsmx=tmsallcs(hs(mind),:);
    tmsmx=tmsmx(find(tmsmx > 0 & ~isnan(tmsmx)));
    hsa=hs(setdiff([1:length(hs)],mind));
    tmsc=tmsallcs(hsa,:);
    nsc=sum(tmsc>0,2);
    nmtch=getmatches(tmsc,tmsmx,0,1000);
    bds=hsa(find(nsc(:)==nmtch(:)));
    datout(bds,:)=[];

    
    
  end
    

  

endfunction
