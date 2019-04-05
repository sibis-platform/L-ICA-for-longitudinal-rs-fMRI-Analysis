function [f,g] = primalProblemUncMultiGeneral(w,X,XX,A,bp,AA,Abp,p,b,ref,exv,alpha,ro1,ro2,numSamples,sessionNum,numEigs)
	f = 0;
    g = [];
    for i = 1:sessionNum
        wi{i} = w(numEigs*(i-1)+1:numEigs*i,:);
        ui{i} = X{i}' * wi{i}; us_i{i}=ui{i}.^2; ex_i{i}=exp(-us_i{i}/2); gaussi{i} =  ui{i}.*ex_i{i};  
        
        dif_i{i} = sum(ex_i{i}-exv);
        
        F{i} = repmat(2 * dif_i{i} / numSamples,numEigs,1) .* (X{i} * gaussi{i} / numSamples) + ...
               2 * alpha * (XX{i}*wi{i}-X{i}*ref');
           
        f = -sum((dif_i{i}/numSamples).^2) + alpha * sum(sum((wi{i}'*X{i}-ref).^2));
        g = [g;F{i}];
    end 
    
    
    %u2 = X2' * w2; u_2=u2.^2; ex2=exp(-u_2/2); gauss2 =  u2.*ex2; 
    %u3 = X3' * w3; u_3=u3.^2; ex3=exp(-u_3/2); gauss3 =  u3.*ex3; 
    %dif2 = sum(ex2-exv);dif3 = sum(ex3-exv);
    %F2 = repmat(2 * dif2 / numSamples,numEigs,1) .* (X2 * gauss2 / numSamples) + 2 * alpha * (XX2*w2-X2*ref2'); 
    %F3 = repmat(2 * dif3 / numSamples,numEigs,1) .* (X3 * gauss3 / numSamples) + 2 * alpha * (XX3*w3-X3*ref3');
         
    %f = -sum((dif1/numSamples).^2) - sum((dif2/numSamples).^2) - sum((dif3/numSamples).^2) +...
    %    alpha * (sum(sum((w1'*X1-ref1).^2)) + sum(sum((w2'*X2-ref2).^2)) + sum(sum((w3'*X3-ref3).^2))) + ...
    %    ro1/2 * sum(sum((A*w+b1).^2)) + ro2/2 * sum(sum((w-(p-b2)).^2));

    %g = [F1;F2;F3] + ro1 * (AA*w+Ab1)+ ro2 * (w-(p-b2));
      
    f = f + ro1/2 * sum(sum((A*w+bp).^2)) + ro2/2 * sum(sum((w-(p-b)).^2));
    g = g + ro1 * (AA*w+Abp)+ ro2 * (w-(p-b));
    
end