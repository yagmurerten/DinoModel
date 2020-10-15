function extmort_new=flight(start_bin,extmort_coef,extmort)

% the function to reduce the extrinsic mortality with a certain coefficient
% the moment flight evolves

% for now we assume flight must have evolved somewhere between
% 210 Mya (roughly start of the body shrinkage) and 150 Mya (Archaeopteryx)

% extmort_coef = how flight changes the extrinsic mortality
% e.g. if it is 2, flying things live 2 times longer, hence their extrinsic
% mortality is halved
% if 0, extmort remains the same as the biggest animals, i.e. it is
% decoupled from the body size completely
extmort_new=extmort;
if extmort_coef~=0
    extmort_new(start_bin:end)=extmort_new(start_bin:end)./extmort_coef;
end

end