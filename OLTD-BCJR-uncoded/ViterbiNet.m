function akHat = ViterbiNet(Net, y, algotype, M)
if M==4
    yTest = [real(y), imag(y)].';
    v_fYcat = cell(size(yTest, 1),1);
    for ii=1:length(yTest)
        v_fYcat{ii} = yTest(:,ii);
    end

    Likelihood = predict(Net, reshape(yTest,[1,2,1,length(yTest)])); 
    Likelihood1 = Likelihood.*rand(size(y));
    if algotype == "Viterbi"
        akHat = Apply_LANNViterbi(y, Likelihood1, M);
    elseif algotype == "BCJR"
        akHat = Apply_LANNBCJR(y, Likelihood1);
    end
elseif M == 2
    yTest = y.';
    v_fYcat = cell(size(yTest, 1),1);
    for ii=1:length(yTest)
        v_fYcat{ii} = yTest(:,ii);
    end
    Likelihood = predict(Net, reshape(yTest,[1,1,1,length(yTest)])); 
    if algotype == "Viterbi"
        akHat = Apply_LANNViterbi(y, Likelihood, M);
    elseif algotype == "BCJR"
        akHat = Apply_LANNBCJR(y, Likelihood, M);
    end
end



end