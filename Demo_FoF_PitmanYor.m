%parpool(2)
for model=8
    for randtry=1
        parfor ttt=1:8
            for dataset=1 %[1,11,9]
                [ttt,dataset,model,randtry]
                FoF_extrapolate_PitmanYor(dataset,randtry,ttt,model);
                close all;
            end
        end
    end
end


for model=1
    for randtry=1
        parfor ttt=1:8
            for dataset=1 %[1,11,9]
                [ttt,dataset,model,randtry]
                %FoF_extrapolate_PitmanYor(dataset,randtry,ttt,model);
                FoF_extrapolate(dataset,randtry,ttt,model);
                close all;
            end
        end
    end
end