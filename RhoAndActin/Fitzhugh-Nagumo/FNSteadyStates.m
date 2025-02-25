a=0.4:0.4:5.2;
b=0.4:0.4:5.2;
Iext=0:0.4:20;
tau=3:3:60;
%Behavior = zeros(length(Iext),length(tau));
Behaviors=[];
Params=[];
for iA=1:length(a)
    for iB=1:length(b)
        for iI=1:length(Iext)
            for iT=1:length(tau)
                [rts,Stab,spiral]=findFNRoots(a(iA),b(iB),Iext(iI),tau(iT),0,0,0);
                nrts=length(rts(:,1));
                if (nrts > 1)
                    ThisBehavior=3;
                elseif (nrts==1&&Stab>0)
                    ThisBehavior=-1;
                else
                    ThisBehavior=1;
                end
                if (ThisBehavior>0)
                    Behaviors=[Behaviors ThisBehavior];
                    Params = [Params [a(iA);b(iB);Iext(iI);tau(iT)]];
                end
            end
        end
    end
end
