function cb_txtupdate(hs,~,htv)
% Update the text boxes with the values of the sliders that have been set.
set(htv,'string',[num2str(get(hs,'value')),' N']);

end