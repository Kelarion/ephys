function [fits,prms] = fitTheseModRF(mdl,rf,xPos,yPos,init_std,nStoch)
% [fits,prms] = fitTheseModRF(mdl,rf,xPos,yPos,init_std,nStoch)
%
% run iteratively through all the models in mdl

nMod = length(mdl);
nX = length(xPos); nY = length(yPos);

if any([mdl(:).fitIpsi]) % fix this when you have time!!! it'll break easily!
    ipsiMod = [mdl(:).fitIpsi]; % assumes only one model has the ipsilateral fit
    saveThisOne = (~ipsiMod) & contains({mdl(:).String},mdl(ipsiMod).String);
%     bounds = setBounds(mdl(useMod),rf,xPos,yPos,init_std);
%     [x, fit_mod] = stochfit(mdl(useMod),rf,xPos,yPos,bounds,nStoch);
else
    saveThisOne = false(nMod,1);
end

fits = zeros(nY,nX,nMod);
prms = cell(nMod,1);
for iMod = 1:nMod
    if mdl(iMod).fitIpsi
        ipsmdl = mdl(iMod);
        switch mdl(iMod).ipsitype
            case 'gaussian'
                ipsmdl.String = 'gaussian';
                ipsmdl.func = @D2GaussFunctionRot;
        end
        ips_bounds = setBounds(ipsmdl,rf,-xPos,yPos,1.5); % (other hemifield)
        [x_ipsi, ipsi_mod] = stochfit(ipsmdl,rf-fit_mod,xPos,yPos,ips_bounds,nStoch);
        full_mod = fit_mod+ipsi_mod;
        x = {x,x_ipsi};
    else
        bounds = setBounds(mdl(iMod),rf,xPos,yPos,init_std);
        [x, full_mod] = stochfit(mdl(iMod),rf,xPos,yPos,bounds,nStoch);
        if saveThisOne(iMod)
            fit_mod = full_mod;
        end
    end
    
    fits(:,:,iMod) = full_mod;
    prms{iMod} = x;
end

end