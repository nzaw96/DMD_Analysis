%solving Burger's using DMD
function [Phi, Lambda, b] = DMDsolution(X, Xprime, r, t)
[U, Sig, V] = svd(X,'econ');
%dt = 0.025;
Ur = U(:,1:r);
Sigr = Sig(1:r,1:r);
Vr = V(:,1:r);
Atilde = Ur'*Xprime*Vr*inv(Sigr); %bringing A down to reduced dimension
[W, Lambda] = eig(Atilde);
%Phi = Xprime*Vr*inv(Sigr)*W; %going back up to original dimension. These are DMD modes
Phi = Xprime*Vr*inv(Sigr)*W;
x1 = X(:,1);
b = Phi\x1;

    



