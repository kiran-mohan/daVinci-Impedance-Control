function cb_conf(hs,~,hs1,hs2,hs3,hs4,hs5,hs6)

% Disable the sliders once Confirm has been clicked.
% This is done to avoid setting the sliders to different positions as
% computations are being handled
set(hs1,'enable','off');
set(hs2,'enable','off');
set(hs3,'enable','off');
set(hs4,'enable','off');
set(hs5,'enable','off');
set(hs6,'enable','off');

end