clear all

benchmark = 2017;
ff = 19;
runTime = 1;
D = 30;
popuSize = 100;
FES = 1 * 10^4 * D;

F=0.5;
CR=0.9;
%%
basePath = 'result\DEpin';
imgPath = basePath+"\imgs";
dataPath = basePath+"\data";
netPath = basePath+"\network";
mkdir(imgPath)
mkdir(dataPath)
mkdir(netPath)
mkdir(basePath);
Vertices_max = FES;
Arcs_max = Vertices_max*4;
FF_sum = {};
B_sum = {};

for f=ff

for run_id = 1 : runTime

    bestfitmix=[];

    [popu, lu, rgo, o, ~, ~, a, alpha, b] = Initialization(popuSize, D, benchmark, f);
    fit = Evaluation(popu, benchmark, f);
    fit=fit';
    bestfit=min(fit);
    nFES=popuSize;

    %%
n=1;m=1;
Vertices = zeros(Vertices_max,D);
Arcs = zeros(Arcs_max,3);
arc_num = 0;
k = popuSize;
Vertices(1:k,:) = popu;
SN = [1: popuSize]';
SN_temp = [1: popuSize]';
CV = popuSize;

while nFES < FES
            
            [~,bee]=min(fit);
            bestfitmix=[bestfitmix;bestfit];
            U = popu;
            % F = ones(popuSize, 1) - 2 .* rand(popuSize, 1) .* ones(popuSize, 1);
            % Get indices for mutation
            [r1, r2, r3] = getindex(popuSize);
            % Implement DE/rand/1 mutation
            bestpop = repmat(popu(bee),popuSize,1);
            V = bestpop + F .* (popu(r1, :) - popu(r2, :));
            % Check whether the mutant vector violates the boundaries or not
            [V] = BoundaryDetection(V,lu);
            % Implement binomial crossover
            for i = 1:popuSize
                j_rand = floor(rand * D) + 1;
                t = rand(1, D) < CR;
                t(1, j_rand) = 1;
                t_ = 1 - t;
                U(i, :) = t .* V(i, :) + t_ .* popu(i, :);
            end
            % evaluate population
            fit_U = Evaluation(U, benchmark, f);
            fit_U=fit_U';
            if min(fit_U)<bestfit
                bestfit=min(fit_U);
            end

%%
            for i = 1:popuSize
                if fit_U(i, :) <= fit(i, :)
                    %Record Vertices and Edges/Arcs
                    k = k + 1;
                    Vertices(k,:) = U(i, :);
                    CV = CV +1;
                    
                    arc_num = arc_num + 1;
                    Arcs(arc_num,:) = [SN(i) CV 0];
                    arc_num = arc_num + 1;
                    Arcs(arc_num,:) = [SN(r1(i)) CV 0];

                    SN_temp(i) = CV;
                end
            end

            SN = SN_temp;


            % get bsf_fit_var
            for i = 1 : popuSize
                if fit_U(i, :) <= fit(i, :)
                   popu(i, :) = U(i, :);
                   fit(i, :) = fit_U(i, :);
                end
            end
                
            %%
nFES = nFES + popuSize;
if mod(m,100)==0 || nFES >= FES
PIN(Arcs,arc_num,k,n,run_id,imgPath,dataPath);
n=n+1;
end
m=m+1;
if nFES > FES; break; end

end
xlswrite('poisson-close',bestfitmix);
end
end

function [r1, r2, r3]  = getindex(popuSize)
r1 = zeros(1, popuSize);
r2 = zeros(1, popuSize);
r3 = zeros(1, popuSize);
for ii = 1 : popuSize
    sequence = 1 : popuSize;
    sequence(ii) = [];
    temp = floor(rand * (popuSize - 1)) + 1;
    r1(ii) = sequence(temp);
    sequence(temp) = [];
    temp = floor(rand * (popuSize - 2)) + 1;
    r2(ii) = sequence(temp);
    sequence(temp) = [];
    temp = floor(rand * (popuSize - 3)) + 1;
    r3(ii) = sequence(temp);
end
end

function PIN(Arcs,arc_num,k,n,time,imgPath,dataPath)

FF_sum=[];B_sum=[];
Arcs = Arcs(1:arc_num,:);
A=Arcs;
B=[];
for i=1:k
    b=length(find(A==i));
    B=[B;b];
end
E=[];
for i=1:max(B)
    d=length(find(B==i));
    E=[E;d];
end
FF=E/k;
P= cumsum(FF)/(sum(FF));
FF_sum = [FF_sum, FF];
B_sum = [B_sum, B]; 

figure
plot(FF,'o')
xlabel('Degree','FontName','Times New Roman','FontSize',18);
ylabel('PDF','FontName','Times New Roman','FontSize',18)
title(['(F) F',num2str(n)],'FontName','Times New Roman','FontSize',18)
set(gcf,'visible','off'); % 不显示图片
saveas(gcf,imgPath+'\PDF_DE_'+num2str(n)+ '_time_'+ num2str(time)+ '.png');
saveas(gcf,imgPath+'\PDF_DE_'+num2str(n)+ '_time_'+ num2str(time)+ '.fig');

figure
plot(P,'o-')
xlabel('Degree','FontName','Times New Roman','FontSize',18);
ylabel('CDF','FontName','Times New Roman','FontSize',18)
title(['(F) F',num2str(n)], 'FontName','Times New Roman','FontSize',18)
set(gcf,'visible','off'); % 不显示图片
saveas(gcf,imgPath+'\CDF_DE_'+num2str(n)+ '_time_'+ num2str(time)+ '.png');
saveas(gcf,imgPath+'\CDF_DE_'+num2str(n)+ '_time_'+ num2str(time)+ '.fig');
save(dataPath+'\Data_'+num2str(n)+'_'+num2str(time)+'.mat', 'Arcs', 'E', 'B')

fprintf('DE %5.0f |time%5.0f',n, time)
end
