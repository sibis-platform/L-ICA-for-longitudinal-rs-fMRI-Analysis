function [f,g] = primalProblemUncMulti(w,X1,X2,X3,XX1,XX2,XX3,A,b1,AA,Ab1,p,b2,ref1,ref2,ref3,exv,alpha,ro1,ro2,numSamples,numEigs)
    w1 = w(1:numEigs,:);
    w2 = w(numEigs+1:numEigs*2,:);
    w3 = w(numEigs*2+1:numEigs*3,:);
    
    u1 = X1' * w1; u_1=u1.^2; ex1=exp(-u_1/2); gauss1 =  u1.*ex1; 
    u2 = X2' * w2; u_2=u2.^2; ex2=exp(-u_2/2); gauss2 =  u2.*ex2; 
    u3 = X3' * w3; u_3=u3.^2; ex3=exp(-u_3/2); gauss3 =  u3.*ex3; 
    
    dif1 = sum(ex1-exv); dif2 = sum(ex2-exv);dif3 = sum(ex3-exv);
    
    F1 = repmat(2 * dif1 / numSamples,numEigs,1) .* (X1 * gauss1 / numSamples) + 2 * alpha * (XX1*w1-X1*ref1');
    F2 = repmat(2 * dif2 / numSamples,numEigs,1) .* (X2 * gauss2 / numSamples) + 2 * alpha * (XX2*w2-X2*ref2'); 
    F3 = repmat(2 * dif3 / numSamples,numEigs,1) .* (X3 * gauss3 / numSamples) + 2 * alpha * (XX3*w3-X3*ref3');
         
    f = -sum((dif1/numSamples).^2) - sum((dif2/numSamples).^2) - sum((dif3/numSamples).^2) +...
        alpha * (sum(sum((w1'*X1-ref1).^2)) + sum(sum((w2'*X2-ref2).^2)) + sum(sum((w3'*X3-ref3).^2))) + ...
        ro1/2 * sum(sum((A*w+b1).^2)) + ro2/2 * sum(sum((w-(p-b2)).^2));

    g = [F1;F2;F3] + ro1 * (AA*w+Ab1)+ ro2 * (w-(p-b2));
                      
end