  getcandpks:
  
  if 1==0
  f=figure
  plot(ddd,'b-o')
  hold on
  scatter(negpks,ddd(negpks),'k','filled')
  trpks=round(trpks)-4;
  scatter(trpks,1.2*ddd(trpks),'r','filled')
  scatter(pospks,ddd(pospks),'g','filled')

end

%[srtd2(:) round(1000*negpksu(srtu(:))'/256)]
