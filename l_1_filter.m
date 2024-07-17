function [MObje_max MX Alpha_opt MIG_opt y] = l1_filter(data)

MObje_max=[];MX=[];Alpha_opt=[]; MIG_opt=[];
for i0=1:size(data,2)
    y=data(:,i0);
    n = length(y);
    e = ones(n,1);
    D = spdiags([e -2*e e], 0:2, n-2, n);
    
    lambda_max=l1tf_lambdamax(y);
    
    Obje=[];X=[];FIG=[];
    for alpha=0:0.00005:0.01
        Order=[0:0.00005:0.01];
        
        lambda=alpha*lambda_max;
        
        [x,status] = l1tf(y,lambda);
        
        
        for i1=1:length(x)-1
            k(i1)=x(i1+1)-x(i1);
        end
        
        [~,loca_k]=unique(k);
        loca_sort=sort(loca_k);
        kk=k(loca_sort);
        seg_num=[loca_sort([2:end],:)-loca_sort([1:end-1],:)+1; length(x)-loca_sort(end)+1];
        
        loca_sort=[loca_sort' length(x)];
        
        for i2=1:length(seg_num)%%Calculate trend features
            k_slope(i2)=(x(loca_sort(i2+1))-x(loca_sort(i2)))/(seg_num(i2)-1);
            b(i2)=x(loca_sort(i2+1))-k_slope(i2)*seg_num(i2);
            delta(i2)=sqrt(sum((y(loca_sort(:,i2):loca_sort(:,i2+1))'-(k_slope(i2).*[1:seg_num(i2)]+b(i2))).^2)/seg_num(i2));
        end
        
        LG=[k_slope' b' delta' seg_num];
        
        LG=mat2cell(LG,length(seg_num),4);
        for i3=1:length(seg_num)
            memb(i3)=sum(exp(-((y(loca_sort(:,i3):loca_sort(:,i3+1))'-(k_slope(i3).*[1:seg_num(i3)]+b(i3))).^2)./(2*delta(i3)^2)))/seg_num(i3);
        end
        
        obje=sum(memb)/length(seg_num);%%objective function
        
        Obje=[Obje obje];X=[X x];FIG=[FIG LG];
        clearvars -except D lambda_max y  n data Obje X FIG MObje_max MX Alpha_opt  MIG_opt Order
    end
    
    Obje_max=max(Obje);
    xx=X(:,find(Obje==Obje_max));IG_opt=FIG(find(Obje==Obje_max));
    
    alpha_opt=Order(find(Obje==Obje_max));
    
    MObje_max=[MObje_max Obje_max];MX=[MX xx];Alpha_opt=[Alpha_opt alpha_opt]; MIG_opt=[MIG_opt IG_opt];
end