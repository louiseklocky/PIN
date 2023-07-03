function [sse_a, r2_a, nihe_result_a]=nihe_analysis(sse_sum, r2_sum, nihe_result_sum)
    % 后两行为mean和std
    method_num = size(sse_sum, 2);
    nihe_result_m = zeros([2, method_num]);

    sse_ms = [mean(sse_sum);std(sse_sum)];
    r2_ms = [mean(r2_sum);std(r2_sum)];
    nihe_result_m(1,:) = mean(nihe_result_sum);
    nihe_result_m(2,:) = nihe_result_m(1,:)>=0.5;

%     sse_sum=sse_sum';
%     r2_sum=r2_sum';

    sse_a = [sse_sum;sse_ms];
    r2_a = [r2_sum;r2_ms];
    nihe_result_a = [nihe_result_sum; nihe_result_m];
end
