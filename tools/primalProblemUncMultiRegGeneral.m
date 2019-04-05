function [f,g] = primalProblemUncMultiRegGeneral(w,X,XX,A,AA,p,b,ref,exv,alpha1,alpha2,ro,numSamples,sessionNum,numEigs)
	f = 0;
    g = [];
    for i = 1:sessionNum
        wi{i} = w(numEigs*(i-1)+1:numEigs*i,:);
        ui{i} = X{i}' * wi{i}; us_i{i}=ui{i}.^2; ex_i{i}=exp(-us_i{i}/2); gaussi{i} =  ui{i}.*ex_i{i};  
        
        dif_i{i} = sum(ex_i{i}-exv);
        
        F{i} = repmat(2 * dif_i{i} / numSamples,numEigs,1) .* (X{i} * gaussi{i} / numSamples) + ...
               2 * alpha1 * (XX{i}*wi{i}-X{i}*ref');
           
        f = -sum((dif_i{i}/numSamples).^2) + alpha1 * sum(sum((wi{i}'*X{i}-ref).^2));
        g = [g;F{i}];
    end 
    
    f = f + ro/2 * sum(sum((w - p + b).^2)) + alpha2 * sum(sum((A*w).^2));
    g = g + ro * (w-(p-b)) + 2 * alpha2 * AA * w;
    
%     w1 = w(1:numEigs,:);
%     w2 = w(numEigs+1:numEigs*2,:);
%     w3 = w(numEigs*2+1:numEigs*3,:);
%     
%     u1 = X1' * w1; u_1=u1.^2; ex1=exp(-u_1/2); gauss1 =  u1.*ex1; 
%     u2 = X2' * w2; u_2=u2.^2; ex2=exp(-u_2/2); gauss2 =  u2.*ex2; 
%     u3 = X3' * w3; u_3=u3.^2; ex3=exp(-u_3/2); gauss3 =  u3.*ex3; 
%     
%     dif1 = sum(ex1-exv); dif2 = sum(ex2-exv);dif3 = sum(ex3-exv);
%     
%     F1 = repmat(2 * dif1 / numSamples,numEigs,1) .* (X1 * gauss1 / numSamples) + 2 * alpha * (XX1*w1-X1*ref1');
%     F2 = repmat(2 * dif2 / numSamples,numEigs,1) .* (X2 * gauss2 / numSamples) + 2 * alpha * (XX2*w2-X2*ref2'); 
%     F3 = repmat(2 * dif3 / numSamples,numEigs,1) .* (X3 * gauss3 / numSamples) + 2 * alpha * (XX3*w3-X3*ref3');
%          
%     f = -sum((dif1/numSamples).^2) - sum((dif2/numSamples).^2) - sum((dif3/numSamples).^2) +...
%         alpha * (sum(sum((w1'*X1-ref1).^2)) + ...
%                  sum(sum((w2'*X2-ref2).^2)) + ...
%                  sum(sum((w3'*X3-ref3).^2))) + ...
%                  ro/2 * sum(sum((w - p + b).^2)) + ...
%                  alpha2 * sum(sum((A*w).^2));
% 
%     g = [F1;F2;F3] + ro * (w-(p-b)) + 2 * alpha2 * AA * w;
                      
end