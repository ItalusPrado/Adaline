clc
clear all
close all

%Setting initial values
samples = 200;
deviation = 0.4;
a = 3;
b = 5;
bias = 1;
learningRate = 0.5;
seasons = 500;
tolerance = 0.1;

%Generating random points
X = linspace(-2,2,samples)';

noise = deviation*randn(samples,1);
Yd = a*X + b + noise;

%Normalization
X = ( (X-min(X)) / (max(X)-min(X)));
Yd = ( (Yd-min(Yd)) / (max(Yd)-min(Yd)));

%Add Bias to X and creating weights
X_Bias(1:samples, 1) = bias;
X = [X_Bias, X];
N = size(X,2);
weights = zeros(1,N);

%Adaline function
A = [];
Vetor_Erros = [];
for i = 1:seasons
    for j = 1:samples
        x = X(j,1:2);
        result = sum(weights.*x);
        error = Yd(j)-result;
        weights = weights+learningRate*error*x;
    end
    quadraticError = 0;
    for k = 1:samples
        x = X(k,1:2);
        result = sum(weights.*x);
        error = Yd(k)-result;
        quadraticError = quadraticError+error^2;
    end
    Vetor_Erros = [Vetor_Erros quadraticError];
    %Ploting
    clf
    plot(X(:, 2), Yd(:,1),'r.')
    grid on
    hold on
    Reta1 = weights(2)*X + weights(1);
    
    A = [A,Reta1];
    hold on
    for m = 1:2:i*2
        
        if m == 1
            xdata =  eval(sprintf('plot(X,A(:,%.0f:%.0f),''r'')', m,m+1));
        elseif m < i*2-1
            xdata =  eval(sprintf('plot(X,A(:,%.0f:%.0f),''g'')', m,m+1));
        else
            xdata =  eval(sprintf('plot(X,A(:,%.0f:%.0f),''b'',''LineWidth'',2)', m,m+1));
        end
        
    end
    hold off
    if quadraticError <= tolerance
        break;
    end
    % Pausa para melhor vizualização da evolução da reta (aproximação)
%     %pause(0.25)
end

plot(Vetor_Erros)
