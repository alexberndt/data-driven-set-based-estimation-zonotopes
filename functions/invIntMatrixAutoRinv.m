function [intMatOut] = invIntMatrixAutoRinv(intMatIn)
%INVMATRIX Summary of this function goes here
%   Detailed explanation goes here

Acenter = (intMatIn.Sup + intMatIn.Inf)/2;
Adelta = (intMatIn.Sup - intMatIn.Inf)/2;


M = pinv(eye(size(Adelta*pinv(Acenter))) - Adelta*abs(pinv(Acenter)));
mu = diag(M);
Tmu = diag(mu);
Tv =  pinv(2*Tmu - eye(size(Tmu)));
Blower = -abs(pinv(Acenter))*M +( pinv(Acenter) + abs(pinv(Acenter)))*Tmu ;
Bupper = abs(pinv(Acenter))*M + (pinv(Acenter) - abs(pinv(Acenter)))*Tmu ;

Blowerlower =min(Blower,Blower*Tv);
Bupperupper =max(Bupper,Bupper*Tv);

intMatOut = intervalMatrix((Bupperupper+Blowerlower)/2,(Bupperupper-Blowerlower)/2);

% M = pinv(eye(size(Acenter)) - abs(pinv(Acenter))*Adelta);
% mu = diag(M);
% Tmu = diag(mu);
% Tv =  pinv(2*Tmu - eye(size(Acenter)));
% Blower = -M*abs(pinv(Acenter)) +Tmu*( pinv(Acenter) + abs(pinv(Acenter))) ;
% Bupper = M*abs(pinv(Acenter)) +Tmu* (pinv(Acenter) - abs(pinv(Acenter))) ;
% 
% Blowerlower =min(Blower,Tv*Blower);
% Bupperupper =max(Bupper,Tv*Bupper);
end

