function pksout = fillgapsmc(pksgl,pksin,rru)

if rru < 625
  tol=40;
elseif rru <1000
  tol=60;
else
  tol=150;
end

bdall=[];
pksm=[0 pksgl 5000];
pksout=[];
numdo=length(pksm)-1;
  for ii=1:numdo
    
    pku=pksm(ii);
    if ii == 1
      rrloc=pksm(3)-pksm(2);
    elseif ii == numdo
      rrloc=pku-pksm(ii-1);
    else     
      rrloc=pksm(ii+1)-pku;
    end
    rrloc=rrloc/(round(rrloc/rru));
    pksr=pksin(find(pksin > pku & pksin < pksm(ii+1)));
    if isempty(pksr)
      continue
    end
    
 
    if ii==1
        dd=pksm(2)-pksr;
    elseif ii < length(pksm)-1
        dd=min(pksr-pku,pksm(ii+1)-pksr);
    else
    
        dd=pksr-pku;
    endif
    %if abs(rrloc/rru) < 1.5 && ii < length(pksm)-1
     %  bd=pksr;
    %else
      gd=find(dd > rrloc-tol & dd < rrloc+tol);
      gd1=find(dd > 2*(rrloc-tol) & dd < 2*(rrloc+tol));
      gd2=find(dd > 3*(rrloc-tol) & dd < 3*(rrloc+tol));
      bd=setdiff(pksr,pksr([gd gd1 gd2]));

    %end

    bdall=[bdall bd];
    
  endfor

  pksout=setdiff(pksin,bdall);
  %pksout=sort([pksout pksgl]);
  pksout=sort([pksout ]);
  return
  pksout=unique(pksout(find(pksout)));
  pksout=pksout(find(~isnan(pksout)));
  pkas=pksout-pksgl(1);
  pkas=round(pkas/rru);
  pkas=pkas-min(pkas)+1
  [pksu,pinds]=unique(pkas)
  pksb=unique(sort(pkas(setdiff([1:length(pkas)],pinds))))
  pksb=ismember(pkas,pksb)
  fv=find(pkas .* pksb)
  pksout=pksout(fv)





  