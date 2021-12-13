function akHat = ApplyLANN(Net, y, algotype, M)
yTest = [real(y), imag(y)];
Likelihood = predict(Net, yTest);

if algotype == "Viterbi"
    akHat = Apply_LANNViterbi(y, Likelihood, M);
elseif algotype == " BCJR"
    akHat = Apply_LANNBCJR(y, Likelihood);
end




end