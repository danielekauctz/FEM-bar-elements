% -------------------------------------------------------------------------
% NATURAL FREQUENCIES AND MODE SHAPES
% Daniele Kauctz Monteiro (2022)
% danielekauctz@hotmail.com
% -------------------------------------------------------------------------
% Parameters:
% fn: natural frequencies (Hz)
% wn: natural frequencies (rad/s)
% phi: mode shapes
% -------------------------------------------------------------------------
function [fn,phi] = modal_analysis(Kr,Mr)
    
    [phi,w2] = eig(Kr, Mr);

    [w2,iw] = sort(diag(w2));
    wn = sqrt(real(w2));
    fn = wn/(2*pi);
    phi = phi(:,iw);

    for kk = 1:length(phi)
       phi(:,kk) = (abs(phi(:,kk))./max(abs(phi(:,kk)))).*sign(real(phi(:,kk)));
    end

end