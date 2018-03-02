function Call_FoF_extrapolate(i1,i2,i3)

for model=1:6
    %clearvars -except  i1 i2 i3 model
   % dataset = i1; %1,9
   % randtry = i2; %1,2,3,4,5
   % ttt = i3; %1,2,3,4
   % FoF_extrapolate(dataset,randtry,ttt,model);
    FoF_extrapolate(i1,i2,i3,model);
end