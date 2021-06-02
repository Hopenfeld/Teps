
function [scoresn,rrso] = meshsegs(tmsnxt,tmscur,rrs,scores,scoresnxt)
  %if the last part of a sequence from the previous segment (tmslast) temporally 
  %meshes with the start of a sequence from the current segment (tmscur)
  %add a fraction of the previous sequence's score to the present sequence's score
  
  scoresn=scores;
  rrso=zeros(1,length(tmscur(:,1)));

    for jj=1:length(tmscur(:,1))
      ns=0;
      tmb=tmscur(jj,find(tmscur(jj,:)));
      tmb=tmb(find(tmb>2500)); 
      if length(tmb) < 2
        continue
      end
      
      rrn=median(tmb(2:end)-tmb(1:end-1));

      for ii=1:length(tmsnxt(:,1)) 
        tma=tmsnxt(ii,find(tmsnxt(ii,:)))+5000;
        tma=tma(find(tma<7500));
        if length(tma) < 2
          continue
        end
        seql=[tmb tma];
        rrd=median(tma(2:end)-tma(1:end-1));
        slst=scoresnxt(ii);
        if abs(rrn - rrd) > 200
          continue
        end
        rrsm=rrd;
        [regsc,rrl,skips,seql] = getmaxscore(seql,rrsm,1);
        regsc=regsc/(length(seql)-2);
        ovap=~isempty(find(seql<=5000)) & ~isempty(find(seql >5000));  
        
        if  (regsc > -700 & ovap) 

          ns=max(ns,int16(max(0,slst)*min(1,(rrsm+double(regsc))/350)));

        end

      end
      scoresn(jj)=scoresn(jj)+ns;
      rrso(jj)=rrn;

  end
  
  
endfunction



