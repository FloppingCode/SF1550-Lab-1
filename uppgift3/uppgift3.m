
% Kod för uppgift 3...
% a)
load eiffel1.mat
figure(1)
trussplot(xnod,ynod,bars)
b=zeros(2*261,1); b(1*2-1)=1;
x = A\b
xbel = xnod + x(1:2:end); ybel = ynod + x(2:2:end);
hold on
trussplot(xbel,ybel,bars,"c")
hold on
plot(xbel(1),ybel(1),'r*')
hold off
% b)
t = []
N = []
for i = 1:4
    if i == 1
        load eiffel1.mat
    elseif i == 2
        load eiffel2.mat
    elseif i == 3
         load eiffel3.mat
    else
        load eiffel4.mat
    end
    
    n = size(A);
    N = [N n(1)];
    timesum = 0;
    for sample = 1:100
        b = randn(N(end),1);
        tic;
        x = A\b;
        temp = toc;
        timesum = timesum+temp;
    end
    timesum/100;
    t = [t timesum];
end
loglog(N,t,'ro',LineStyle='--')
%% Tidtabell
% Skapa en 4x4-matris T som innehåller beräkningstiderna.
% Raderna ska motsvara de olika modellerna (eiffel1-eiffel4) och
% kolumnerna de olika metoderna, ordnade som "Naiv", "LU",
% "Gles" och "Gles LU".
% Följande kod skapar en snygg tabell med resultaten:
%tab=array2table(T,'VariableNames',{'Naiv' 'LU' 'Gles' 'Gles LU'},'RowNames',
%{'eiffel1' 'eiffel2' 'eiffel3' 'eiffel4'});
%disp(tab);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [jmin,jmax]=kanslighet(A,metod)
% Indata:
%
% A - matrisen
% metod - villken metod som används:
% 1 = Naiv metod
% 2 = LU-faktorisering
%
% Utdata:
%
% jmin - index för minst känsliga nod
% jmax - index för mest känsliga nod
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% kod...
end
