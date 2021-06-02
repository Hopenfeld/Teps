function [pks,snrb,pksbst,snrbt,ddd,nslvlu,stats]  = getcandpks(sigin,d,samprate,cluslim,numcand,trpks,rrlim)

load nnetprominence
load pkprbmats qualmat
[nrq,ncq]=size(qualmat);
omap=[];
pksx=[];
cscr=[];
ds=[];
lr=[];
rsa=[];numpks=[];tmsall=[];bnew=[];scrs=[];rrsa=[];po=[];dd=[];xdd=[];
snrb1=[]; probsmb=[]; stats=[]; pksten=[]; snrb1t=[]; snrb1o=[]; snrb=[];
pkparams=[]; nslvlu=[];

sampadj=1000/samprate;
sigin=sigin;
dorig=d;
direc=1;

stopd=0;
qmax=0;
segquallst=-10;  %for adaptive peak selection; not implemented

while ~stopd  %this loop is set up for adaptive peak selection based on changing the 
              %magnitude of the differencing distance d; this is not implemented so the
              %loop is only run once
  
  dd=sigin(d:end)-sigin(1:end-d+1);
  ddd=dd(d:end)-dd(1:end-d+1);
  [pkval,pospks] = findpeaks(ddd,"DoubleSided");
  pospks=pospks(find(pkval > 0));
  [pkval,negpks] = findpeaks(ddd,"DoubleSided","MaxPeakWidth",5); %constrain negative peak width
  negpks=negpks(find(pkval < 0));
  [dt,pr,d2m,d2p] = processpksloc(ddd,negpks,pospks,sampadj);
  %dt is number of samples between major negative peak (in 2nd difference signal) and previous 
  %positive peak
  %d2m is the magnitude of major negative peak
  d2mr=d2m;
  pksxa=round(negpks * sampadj);
  gdt=find(dt(:,1)>cluslim(1) & dt(:,1) < cluslim(2));  %filter out too wide/narrow peaks
  dt=dt(gdt,:);
  d2m=d2m(gdt);  
  negpks=negpks(gdt);
  pr=pr(gdt,:);
  d2p=d2p(gdt);
  d2d=d2m(2:end)-d2m(1:end-1);
  d2dd=d2d(2:end)-d2d(1:end-1); %peak prominence (2nd difference) of second difference ddd peaks
  d2dd=[2*d2d(1); d2dd; -2*d2d(end)]; %add starting and ending peaks; correct for lack of 2 neighbors
  [srtd2,srt]=sort(d2dd,'ascend');
  srtd2(end:20)=0;
  [nrvalso,numtokeep] = getbestpks(srtd2(1:20),10); %get potentially high quality peaks by looking for big
                                                    %drops in d2dd
  segqual = min(max(1,round(100*sim(net,nrvalso(1:20)))),ncq); %estimate noise level with neural net
  nslvl=[1:nrq]*qualmat(:,segqual)/sum(qualmat(:,segqual));

  if 1==1 %segqual > qmax %condition if adaptive peak selection is invoked; 
    qmax=segqual;
    dmax=d;
    srtu=srt;
    negpksu=negpks;
    srtd2u=srtd2;
    d2mu=d2m;
    d2mru=d2mr;
    valsu=nrvalso;
    numpku=numtokeep;
    dddu=ddd;
    pksxau=pksxa;
    nslvlu=nslvl;
  end
  dq=segqual-segquallst;
  if segqual > 0.9 || d == dorig %segual condition if adaptive peak selection is invoked; "or" as
                                 %implemented here is always satisfied therefore removing adaptive peak selection
    stopd = 1;
  else
    segquallst=segqual;
    if d == dorig + 2
      direc=-1;
    end
    d=d+direc;
  end
  
end

if length(d2dd) > 3
  srt=srtu(1:min(numcand,length(srtu))); %select the numcand best peaks
  npksa=length(srtd2u);
  stats=segqual;
  [pks,srt1]=sort(negpksu(srt));  %put best peaks in sequential order
  srtb=srtu(1:min(numpku,length(srtu)));  %highest quality peaks
  [pksbst,srt1bst]=sort(negpksu(srtb));
end

d2m=d2mu(srt(srt1));
pksx = round(pks * sampadj);
pksbst=round(pksbst * sampadj);
snrbt = getpairqual(pksbst,pksxau,d2mru,rrlim); %SNR prominence matrices
snrb = getpairqual(pksx,pksxa,d2mru,rrlim);
[rr,cc,vals]=find(snrb);
snrb=[pksx(rr)' pksx(cc)' vals];

function [dt,pr,d2m,d2p] = processpksloc(ddd,negpks,pospks,sampadj)

pr=zeros(length(negpks),2);
dt=zeros(length(negpks),2);
d2m=zeros(length(negpks));
d2p=zeros(length(negpks));
for kk=1:length(negpks)

  dp=negpks(kk)-pospks(find(pospks < negpks(kk)));
  [dum,mind]=min(dp);

  if ~isempty(mind)
    pospr=pospks(mind);
    cpv=abs(ddd(negpks(kk)));

    if kk > 1 
      if pospr > negpks(kk-1) 
        dt(kk,1)=negpks(kk)-pospr;
        
        d2m(kk)=cpv; 
        d2p(kk)=ddd(pospr);

        if mind < length(pospks)
          nxtp=pospks(mind+1);
          if ddd(nxtp)  > 0 & nxtp < negpks(min(kk+1, length(negpks)))
            dt(kk,2)=nxtp-negpks(kk);
            d2p(kk)=max(ddd(pospr),ddd(nxtp));
          end
        end
        
        pa = find(negpks(1:kk-1) > negpks(kk)-500/sampadj);
        if ~isempty(pa)
          svals=sort(abs(ddd(negpks(pa))),'descend');
          pr(kk,1)=cpv/mean(svals(1:min(5,length(svals))));
        end
        if kk < length(negpks)
          pa = find(negpks(kk+1:end) < negpks(kk)+500/sampadj);
          svals=sort(abs(ddd(negpks(pa))),'descend');
          if ~isempty(pa)
            pr(kk,2)=cpv/mean(svals(1:min(5,length(svals))));
          end

        end       
        
      end
    else
      dt(kk,1)=negpks(kk)-pospr;
      d2m(kk)=abs(ddd(negpks(kk)));
    end
    
  end
 
end

function [nrvalso,numtokeep] = getbestpks(nrvalso,nmlim)

nrvalso=nrvalso/mean(nrvalso(3:6));
dv6=nrvalso(4:end-6)-nrvalso(10:end);
    [drp,mind]=max(dv6);
    mind=mind+3;
    dr1=nrvalso(1:end-1)-nrvalso(2:end);
    dr1s=dr1(mind+[0:5]);
    drps=sort(dr1s','descend')/drp;
    abv=nrvalso(mind)-0.8*drp;
    fst=find(nrvalso < abv);
    numtokeep=min(fst(1)-1,nmlim);