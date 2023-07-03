function [sse_sum, r2_sum, nihe_result_sum]=new_nihe(totalTime, func, FF_sum, B_sum, basePath, node_num, ifSaveImg,numbb)
    % 测试数据拟合程度
    % 1:Normal  2:Gamma  3:Possion  4:Exp  5:rayleigh
    % 6:Weibull 7:Chi-square 8:Logistic 9:Powerlaw
    method_num = 3;
    sse_sum = zeros([totalTime method_num]);
    r2_sum = zeros([totalTime method_num]);
    nihe_result_sum = zeros([totalTime method_num]);
    nihe_name = ["Normal", "Gamma", "Possion", "Exp", "Rayleigh" ...
                "Weibull", "Chi-square"];
    % save_path = basePath+"/nihe";
    % save_path_imgs = save_path + "/" + num2str(func);
    % mkdir(save_path_imgs)
    time = 1;
    while time <= totalTime
        FF = FF_sum{time};
        pf = cumsum(FF)/(sum(FF));
        B = B_sum{time};
        i=[1:max(B)]';

        % Normal
%         [mu,sigma]=normfit(B);
%         p1=normcdf(i,mu,sigma);
%         p1_1=normpdf(i,mu,sigma);
%         
%         %sse = sum((YReal - YPred).^2);
%         sse_sum(time, 1) = sum((pf - p1).^2);
%         %r2 = 1 - (sum((YPred - YReal).^2) / sum((YReal - mean(YReal)).^2));
%         r2_sum(time, 1) = 1 - (sum((pf - p1).^2) / sum((p1 - mean(p1)).^2));
%         %r2 = max(0,1 - sum((y(:)-f(:)).^2)/sum((y(:)-mean(y(:))).^2));
%         %r2 = max(0,1 - sum((p1(:)-A1(:)).^2)/sum((p1(:)-mean(p1(:))).^2))
%         %r2 = 1 - (sum((p1 - A1).^2) / sum((A1 - mean(A1)).^2))
%         
%         [H1,s1]=kstest2(FF, p1_1);
%         if H1==0
%             %disp('Data obey Normal Distribution');
%             nihe_result_sum(time, 1) = 1;
%         %else
%         %     disp('Data donot obey Normal Distribution');
%         end
        
        % Gamma
%         phat=gamfit(B);
%         p2=gamcdf(i,phat(1),phat(2));
%         p2_1=gampdf(i,phat(1),phat(2));
% 
%         sse_sum(time, 2) = sum((pf - p2).^2);
%         r2_sum(time, 2) = 1 - (sum((pf - p2).^2) / sum((p2 - mean(p2)).^2));
%         [H2,s2]=kstest2(FF, p2_1);
%         if H2==0
%             nihe_result_sum(time, 2) = 1;
%         end

        % Possion
        lamda=poissfit(B);
        p3=poisscdf(i,lamda);
        p3_1=poisspdf(i,lamda);

        sse_sum(time, 1) = sum((pf - p3).^2);
        r2_sum(time, 1) = 1 - (sum((pf - p3).^2) / sum((p3 - mean(p3)).^2));
        [H3,s3]=kstest2(FF, p3_1);
        if H3==0
            nihe_result_sum(time, 1) = 1;
        end
        
        % Exp
%         mu=expfit(B);
%         p4=expcdf(i,mu);
%         p4_1=exppdf(i,mu);
% 
%         sse_sum(time, 4) = sum((pf - p4).^2);
%         r2_sum(time, 4) = 1 - (sum((pf - p4).^2) / sum((p4 - mean(p4)).^2));
%         [H4,s4]=kstest2(FF, p4_1);
%         if H4==0
%             nihe_result_sum(time, 4) = 1;
%         end
        
%         % rayleigh
%         [phat] = raylfit(B);
%         p5=raylcdf(i,phat);
%         p5_1=raylpdf(i,phat);
% 
%         sse_sum(time, 5) = sum((pf - p5).^2);
%         r2_sum(time, 5) = 1 - (sum((pf - p5).^2) / sum((p5 - mean(p5)).^2));
%         [H5,s5]=kstest2(FF, p5_1);
%         if H5==0
%             nihe_result_sum(time, 5) = 1;
%         end
        
%         Weibull
%         [phat] = wblfit(B);
%         p6=wblcdf(i,phat(1),phat(2));
%         p6_1=wblpdf(i,phat(1),phat(2));
% 
%         sse_sum(time, 6) = sum((pf - p6).^2);
%         r2_sum(time, 6) = 1 - (sum((pf - p6).^2) / sum((p6 - mean(p6)).^2));
%         [H6,s6]=kstest2(FF, p6_1);
%         if H6==0
%             nihe_result_sum(time, 6) = 1;
%         end
        
        % Uniform
    %     [ahat,bhat] = unifit(B);
    %     p6=unifcdf(i,ahat,bhat);
    %     p6_1=unifpdf(i,ahat,bhat);
    %     [H6,s6]=kstest2(FF, p6_1);
    % 
    %     sse_sum(time, 6) = sum((pf - p6).^2);
    %     r2_sum(time, 6) = 1 - (sum((pf - p6).^2) / sum((p6 - mean(p6)).^2));
    %     if H6==0
    %         nihe_result_sum(time, 6) = 1;
    %     end
        
%         % Chi-square
%         p7=chi2cdf(i,7);
%         p7_1=chi2pdf(i,7);
%         [H7,s7]=kstest2(FF, p7_1);
%         sse_sum(time, 7) = sum((pf - p7).^2);
%         r2_sum(time, 7) = 1 - (sum((pf - p7).^2) / sum((p7 - mean(p7)).^2));
%         if H7==0
%             nihe_result_sum(time, 7) = 1;
%         end

        % Logistic
        pd_Logistic = fitdist(B,'Logistic');
        p8 = cdf(pd_Logistic, i);
        p8_1 = pdf(pd_Logistic, i);
        [H8,s8]=kstest2(FF, p8_1);
        sse_sum(time, 2) = sum((pf - p8).^2);
        r2_sum(time, 2) = 1 - (sum((pf - p8).^2) / sum((p8 - mean(p8)).^2));
        if H8==0
            nihe_result_sum(time, 2) = 1;
        end

        % Powerlaw
        x_start_num = node_num+1;
        x=x_start_num:max(B);
        W1=FF;
        W1(1:node_num,:)=[];
        pf = cumsum(W1)/(sum(W1));
        if length(FF)>=x_start_num+numbb && FF(x_start_num+numbb)~=0
            alpha=(log(FF(x_start_num))-log(FF(x_start_num+numbb)))/(log(x_start_num+numbb)-log(x_start_num));
        elseif length(FF)>=x_start_num+numbb+1 && FF(x_start_num+numbb+1)~=0
            alpha=(log(FF(x_start_num))-log(FF(x_start_num+numbb+1)))/(log(x_start_num+numbb)-log(x_start_num));
        elseif length(FF)>=x_start_num+numbb-1 && FF(x_start_num+numbb-1)~=0
            alpha=(log(FF(x_start_num))-log(FF(x_start_num+numbb-1)))/(log(x_start_num+numbb)-log(x_start_num));
        elseif FF(length(FF))~=0
            alpha=(log(FF(x_start_num))-log(FF(length(FF))))/(log(x_start_num+numbb)-log(x_start_num));  
        elseif FF(length(FF)-1)~=0
            alpha=(log(FF(x_start_num))-log(FF(length(FF)-1)))/(log(x_start_num+numbb)-log(x_start_num));  
        end
        C=FF(x_start_num)/(x_start_num^(-alpha));
        x_p = x';
        y1=C*x_p.^(-alpha);
        y_powerlaw = cumsum(y1)/(sum(y1));

        sse_sum(time, 3) = sum((pf - y_powerlaw).^2);
        r2_sum(time, 3) = 1 - (sum((pf - y_powerlaw).^2) / sum((y_powerlaw - mean(y_powerlaw)).^2));
        

        if ifSaveImg
            figure
            plot(pf,'o-')
            hold on
            % plot(p1,'-r*')
            % hold on
            % plot(p2,'-g+')
            % hold on
            % plot(p3,'-bx')
            % hold on
%             plot(p4,'-cs')
%             hold on
%             plot(p5,'-yd')
%             hold on
%             plot(p6,'-m<')
%             hold on
%             plot(p7,'-bx')
%             hold on           
            plot(y_powerlaw,'-r^')
            hold on
            % plot(p8,'-cs')
            % hold on
            
            %clf(fig)
            % set(gcf,'visible','off'); % 不显示图片
            %legend('Original data','Normal','Gama','Possion','Exponential','Rayleigh','Weibull','Chi-square','Logistic','Power law')
            title('CDF to Power law: F9 iteration 28','Interpreter','latex');
            legend('LSHADE','Possion','Power law')
            
            xlabel('Degree')
            ylabel('CDF')
    
            % saveas(gcf, save_path_imgs+ '/time_'+ num2str(time)+ '.png');
            % saveas(gcf, save_path_imgs+ '/time_'+ num2str(time)+ '.fig');
        end

        
        
        time = time + 1;
    end

    % [sse_a , r2_a, nihe_result_a]=nihe_analysis(sse_sum, r2_sum, nihe_result_sum);
    % save_path_datas = save_path+ '/func_'+ num2str(func)+ '.mat';
    % save(save_path_datas, 'sse_a', 'r2_a', 'nihe_result_a','sse_sum','r2_sum');
    %
    %[phat] = betafit(B)
    %p6=betacdf(i,phat(1),phat(2));
    %