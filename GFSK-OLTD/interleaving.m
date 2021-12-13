function [alpha,ck] = interleaving(bk, state)
   L_total = length(bk) ;
 for jj = 1:100
    a = (1:1:L_total);
    Mat = repmat(a,L_total,1);
    for ii = 1:L_total
        Mat(ii,:) = abs(ii*ones(1,L_total)-a)>=state;
    end

    init = randperm(L_total,1); %随机选择第一个元素
    alpha = [init];
    k = 1 ;
    n2 = init;

    while k<L_total
        Mat(:,alpha) = 0; %被选择过的元素置0
        row = Mat(n2,:); %取出当前行
        notZero = find(row);
        if length(notZero)>=1
            numnorzero = randperm(length(notZero));
            n2 = notZero(numnorzero(1)); %随机选出下一个元素
        end

        alpha = [alpha,n2];
        k=k+1;
    end
    if alpha(end)-alpha(end-1) >= state
        break
    end
end
ck = bk(alpha(1:L_total));
end
