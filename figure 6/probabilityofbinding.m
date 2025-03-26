function [F] = probabilityofbinding1(concs,tfcombinations,intvect,Q,kmx,nobms)


con=repelem(concs,nobms);
zoffc=kmx.*con;
zonc=zoffc.*Q;
Zon=0;
Zoff=0;


for i=1:1:length(tfcombinations)
    Zoff=Zoff+prod(zoffc(tfcombinations{i}))*intvect(i);
    Zon=Zon+prod(zonc(tfcombinations{i}))*intvect(i);
end


Zoff = Zoff+1;



F=Zon/(Zon+Zoff);

end