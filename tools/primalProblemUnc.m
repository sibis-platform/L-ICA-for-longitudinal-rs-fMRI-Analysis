function [f,g] = primalProblemUnc(w,X1,X2,X3,XX1,XX2,XX3,A,b1,AA,Ab1,p,b2,ref1,ref2,ref3,exv,alpha,ro1,ro2,numSamples,numEigs)
    w1 = w(1:numEigs);
    w2 = w(numEigs+1:numEigs*2);
    w3 = w(numEigs*2+1:numEigs*3);
    
    u1 = X1' * w1; u_1=u1.^2; ex1=exp(-u_1/2); gauss1 =  u1.*ex1; 
    u2 = X2' * w2; u_2=u2.^2; ex2=exp(-u_2/2); gauss2 =  u2.*ex2; 
    u3 = X3' * w3; u_3=u3.^2; ex3=exp(-u_3/2); gauss3 =  u3.*ex3; 
            
    F1 = 2 * sum(ex1-exv) / numSamples * X1 * gauss1 / numSamples + 2 * alpha * (XX1*w1-X1*ref1') + ro1 * (AA(1:numEigs,:)*w   +Ab1(1:numEigs))   + ro2 * (w1-(p(1:numEigs)-b2(1:numEigs)));
    F2 = 2 * sum(ex2-exv) / numSamples * X2 * gauss2 / numSamples + 2 * alpha * (XX2*w2-X1*ref2') + ro1 * (AA(numEigs+1:numEigs*2,:)*w +Ab1(numEigs+1:numEigs*2)) + ro2 * (w2-(p(numEigs+1:numEigs*2)-b2(numEigs+1:numEigs*2)));
    F3 = 2 * sum(ex3-exv) / numSamples * X3 * gauss3 / numSamples + 2 * alpha * (XX3*w3-X1*ref3') + ro1 * (AA(numEigs*2+1:numEigs*3,:)*w+Ab1(numEigs*2+1:numEigs*3))+ ro2 * (w3-(p(numEigs*2+1:numEigs*3)-b2(numEigs*2+1:numEigs*3)));
    
%     F1 = 2 * sum(ex1-exv) / numSamples * X1 * gauss1 / numSamples + 2 * alpha * (w1-ref1) + ro1 * (AA(1:numEigs,:)*w   +Ab1(1:numEigs))   + ro2 * (w1-(p(1:numEigs)-b2(1:numEigs)));
%     F2 = 2 * sum(ex2-exv) / numSamples * X2 * gauss2 / numSamples + 2 * alpha * (w2-ref2) + ro1 * (AA(numEigs+1:numEigs*2,:)*w +Ab1(numEigs+1:numEigs*2)) + ro2 * (w2-(p(numEigs+1:numEigs*2)-b2(numEigs+1:numEigs*2)));
%     F3 = 2 * sum(ex3-exv) / numSamples * X3 * gauss3 / numSamples + 2 * alpha * (w3-ref3) + ro1 * (AA(numEigs*2+1:numEigs*3,:)*w+Ab1(numEigs*2+1:numEigs*3))+ ro2 * (w3-(p(numEigs*2+1:numEigs*3)-b2(numEigs*2+1:numEigs*3)));
       
    f = -(sum(ex1-exv)/numSamples)^2 - (sum(ex2-exv)/numSamples)^2 - (sum(ex3-exv)/numSamples)^2 +...
        alpha * (sum((w1'*X1-ref1).^2) + sum((w2'*X2-ref2).^2) + sum((w3'*X3-ref3).^2)) + ...
        ro1/2 * sum((A*w+b1).^2) + ro2/2 * sum((w-(p-b2)).^2);
%     f = -(sum(ex1-exv)/numSamples)^2 - (sum(ex2-exv)/numSamples)^2 - (sum(ex3-exv)/numSamples)^2 +...
%         alpha * (sum((w1-ref1).^2) + sum((w2-ref2).^2) + sum((w3-ref3).^2)) + ...
%         ro1/2 * sum((A*w+b1).^2) + ro2/2 * sum((w-(p-b2)).^2);
    g = [F1;F2;F3];
end