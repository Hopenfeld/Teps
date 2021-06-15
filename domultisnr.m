function stats = domultisnr(sigs,trpks,samprate,moset)
%
load pkprbmats prbsnr

pksall=[]; nmsall=[]; regscall=[]; compall=[]; snrall=[]; rrsall=[];stats=[];

nchs=11; %number of high quality of peaks in multichannel search
d=6; %number of samples between second difference points
cluslim=[4 20]; %minimun and maximum peak widths (in samples)
pkgp=6; %temporal alignment limit (in samples)
rthold=[0 0 -2000 -4000 -6000 -8000 -10000 -12000*ones(1,10)]; %temporal regularity thresholds as function 
                                                               %based on sequence length 
cfs=[-0.5437 0.002004 0.001716];   %score coefficients for skips, temporal regularity and SNR/prominence
rrlim=[400 1300];  %rr interval limits in ms
                                                          
numchan=length(sigs(:,1));
numcand=max(20,30-5*(4-numchan)); %number of "best" candidate peaks in each channel

datoutsc=[];
snrmat=[];
trpksadj=round(trpks*samprate/1000);  %adjust true peaks for sampling rate

for ii=1:numchan  %process each channel

  siggu=sigs(ii,:);
  [pksx,snrb,pksten,snrba,ddd,nslvlu,statspks] = getcandpks(siggu,d,samprate,cluslim,numcand,trpksadj,rrlim);
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
      snrssckp=round(100*snrssc); 
      smat=[skipssc(:) nregsc snrssckp(:)]';
      gdh=find(nregsc > -800);
      %[tmsallsc(gdh,:) round(nregsc(gdh)) snrssckp(gdh)(:) rrssc(gdh)(:)]
      scr=round(100*cfs*smat); %score the sequences
      [dum,sinds]=sort(scr,'descend');
      kp=sinds(1:min(3,length(sinds))); %keep 3 best sequences
      datadd=[rrssc(:)(kp) scr(:)(kp) tmsallsc(kp,:)]; %store rr intervals, scores and sequences
      datoutsc(end+[1:length(datadd(:,1))],[1:length(datadd(1,:))])=datadd; %will be added to best multichannel sequences
    end

  end
  
  pksa(ii,1:length(pksx))=pksx';
  snrmat=[snrmat;snrb]; %assemble global SNR matrix in row/col/value format and 
                        %adjust according to channel correlation

end

if numchan > 1
  [spksa,pkprb,~,pkmap]=getclustersmsnr(pksa,0,pkgp,samprate);
%get global peaks spksa and associated Bayesian timing coherence probability pkprb by clustering peak times (pksa) 
%across channels; pkmap is the mapping from local channel peaks to the global peaks
else
  spksa=round(pksx*1000/256);
  pkprb=0*spksa;
  pkmap(spksa)=[1:length(spksa)];

end

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
tmsall=tmsall(gds,:);
regsc=regsc(gds);
rrs=rrs(gds);
skips=skips(gds);
glprbs=glprbs(gds);
nseq=nseq(gds);
datout=[];
scra=[];
sra1a=[];
mxsc=-1000;

for ii=1:length(tmsall(:,1))  %loop through parent sequences
  
  tmsloc=tmsall(ii,:);
 
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
                                                  %need to fix gap between snrs between segs
    lhs=double(regsca(:))./(nmseqa(:)-2); 
    snrvs=round(100*snrvs);
    smat=[skipsa(:) lhs snrvs(:)]';
    scr1=cfs*smat;    %peak timing coherence glprbs not used here but improve detection in certain cases
    [sscrs,sinds]=sort(scr1,'descend');
    kp=sinds(1:min(length(sinds),10)); %keep 10 best sequences per parent sequence
    scr1=round(100*scr1(kp))';
    datadd=[ rrsa(:)(kp) scr1 tmsalla(kp,:)];
    datout(end+[1:length(datadd(:,1))],[1:length(datadd(1,:))])=datadd;
  end
end
if ~isempty(datoutsc)
 datout(end+[1:length(datoutsc(:,1))],[1:length(datoutsc(1,:))])=datoutsc;
end

datout = cullseqsloc(datout);
%remove the lower scoring sequences in closely related groups

stats.dat=datout;

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

for ii = 1:length(rrs)
  clmp=find(datin(:,1)==rrs(ii) & datin(:,2) > -500);
  scrs=datin(clmp,2);
  [dum,minds]=sort(scrs,'descend');
  mout=minds(1:min(3,length(minds)));
  datout(end+[1:length(mout)],:)=datin(clmp(mout),:);
  
end
endfunction
