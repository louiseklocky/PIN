clear all
runTime=1;
ifSaveImg=1;
node_num=8;
numbb=10;

for problem_i=1
FF_sum={};
B_sum={};
nihePath='LSHADEpin';

for time_i=1:runTime
    data_path=nihePath+"/data/Data_"+num2str(problem_i)+"_"+num2str(time_i)+".mat";
    load(data_path);
    k=max(Arcs(:,1));
    FF=E/k;
    FF_sum=[FF_sum,FF];
    B_sum=[B_sum,B];
end
[sse_sum,r2_sum,nihe_result_sum]=new_nihe(runTime,problem_i,FF_sum,B_sum,nihePath,node_num,ifSaveImg,numbb);
end
