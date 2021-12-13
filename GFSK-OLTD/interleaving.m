function [alpha,ck] = interleaving(bk, state)
   L_total = length(bk) ;
 for jj = 1:100
    a = (1:1:L_total);
    Mat = repmat(a,L_total,1);
    for ii = 1:L_total
        Mat(ii,:) = abs(ii*ones(1,L_total)-a)>=state;
    end

    init = randperm(L_total,1); %���ѡ���һ��Ԫ��
    alpha = [init];
    k = 1 ;
    n2 = init;

    while k<L_total
        Mat(:,alpha) = 0; %��ѡ�����Ԫ����0
        row = Mat(n2,:); %ȡ����ǰ��
        notZero = find(row);
        if length(notZero)>=1
            numnorzero = randperm(length(notZero));
            n2 = notZero(numnorzero(1)); %���ѡ����һ��Ԫ��
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
