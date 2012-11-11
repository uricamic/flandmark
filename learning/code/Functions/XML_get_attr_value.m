function val = XML_get_attr_value( el, key )
%
% [val] = XML_get_attr_value( el, key )
%
% Read value 'val' of attribute 'key' from DOM element 'el'.

% Author: Vit Zyka, vit.zyka@seznam.cz
% 2007-07-19 v0.1 new
% 2009-07-30 v0.1.1 used method getNamedItem instead of searching in 'for'

val = [];
attrs = el.getAttributes;
attr  = attrs.getNamedItem(key);
if isempty(attr), return; end
val   = attr.getValue;
% for i=1:attrs.getLength
%   attr = attrs.item(i-1);
%   if( strcmp(key,attr.getName) )
%     val = attr.getValue;
%   end
% end

return

%% end of XML_get_attr_value.m
