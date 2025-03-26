function [F] = probabilityofbinding(concs,tfcombinations,intvect,Q,kmx,nobms,be)

% calculates binding probability using (1) in He et al.

con=repelem(concs,nobms);
zoffc=kmx.*con;
zonc=zoffc.*Q;
Zon=0;
Zoff=0;


for i=1:1:length(tfcombinations)
    Zoff=Zoff+prod(zoffc(tfcombinations{i}))*intvect(i);
    Zon=Zon+prod(zonc(tfcombinations{i}))*intvect(i);
end

Zon=Zon+be;  % modified to include basal expression (be)
Zoff = Zoff+1;



F=Zon/(Zon+Zoff);

end