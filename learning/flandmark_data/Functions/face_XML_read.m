function [ T ] = face_XML_read( data )
%
% [ T ] = face_XML_read( data )
%
% Reads FACE data structure T from XML file or DOM object
%
% Input:
%   data - XML filename | DOM object
%
% Output (face structure):
%
% Face structure:
% 
% T.image.filename
% T.image.width
% T.image.height
% T.eye(i).open     1 for opened eye, 0 for closed
% T.eye(i).center   eye center coord [row,col]

%
% Author: Kuba Zemlicka, jakub.zemlicka@eyedea.cz
% 2008-02-05 v0.1 new
% 11-07-11 Michal Uricar, corners dataset modification


T = [];

if( strcmp(class(data),'org.apache.xerces.dom.DeferredElementImpl') )
    docNode = data;
else
%% Convert XML-file to the current format if necessary
    % load the file lines
    fid = fopen(data,'r');
    nl = 0; line = [];
    while ~feof(fid)
        nl = nl + 1;
        line{nl} = fgetl(fid);
    end;
    fclose(fid);
    % check the file format 
    % remove the second line "<!DOCTYPE ..." if necessary
    if (length(line)>1) && strcmp(line{2}(1:min(9,length(line{2}))),'<!DOCTYPE')
        fid = fopen(data,'w');
        for jj = 1:nl
            if strcmp(line{jj}(1:min(9,length(line{jj}))),'<!DOCTYPE')
                continue;
            end;
            fprintf(fid,'%s\n',line{jj});
        end;
        fclose(fid);
        fprintf('File %s converted to the current xml-format.\n',data);
    end;
    %% load the file in the current format
    docNode = xmlread( data );
end

%% Fill face data

% --- image
imgList = docNode.getElementsByTagName('image');
imgEl = imgList.item(0);
fnameList = imgEl.getElementsByTagName('filename');
fnameEl = fnameList.item(0);
T.image.filename = char( fnameEl.getTextContent );
T.image.width  = str2num(XML_get_attr_value( imgEl,'width' ));
T.image.height = str2num(XML_get_attr_value( imgEl,'height' ));
try T.image.time   = str2num(XML_get_attr_value( imgEl,'time' ));catch T.image.time = 0; end

% --- eye
faceList = docNode.getElementsByTagName('face');

for i=1:faceList.getLength
    faceEl = faceList.item(i-1);
	if (isempty(faceList)), error('Face_XML_read: Empty face section! Corrupted xml file?'); end;
	    
    gen = XML_get_attr_value(faceEl, 'gen');
    age = XML_get_attr_value(faceEl, 'age');
    tag = XML_get_attr_value(faceEl, 'tag');
    id  = XML_get_attr_value(faceEl, 'id');
    
    if (~isempty(gen) && ~isempty(age) && ~isempty(tag))
        T.face(i).gen = char(gen);
    	T.face(i).age = str2num(age);
        T.face(i).tag = str2num(tag);
    end
	
    if ~isempty(id)
        T.face(i).id = str2num(id);
    end;
    
	eyeList = faceEl.getElementsByTagName('eye');
	if (~isempty(eyeList) && eyeList.getLength == 2)
		eyeEl = eyeList.item(0);
		T.face(i).eyer = [str2num(XML_get_attr_value(eyeEl, 'col')); str2num(XML_get_attr_value(eyeEl, 'row'))];
		eyeEl = eyeList.item(1);
		T.face(i).eyel = [str2num(XML_get_attr_value(eyeEl, 'col')); str2num(XML_get_attr_value(eyeEl, 'row'))];
    end
	
    canthus_rrList = faceEl.getElementsByTagName('canthus_rr');
    if (~isempty(canthus_rrList) && canthus_rrList.getLength > 0)
        canthus_rrE1 = canthus_rrList.item(0);
        T.face(i).canthus_rr = [str2num(XML_get_attr_value(canthus_rrE1, 'col')); str2num(XML_get_attr_value(canthus_rrE1, 'row'))];
    end;
    
    canthus_rlList = faceEl.getElementsByTagName('canthus_rl');
    if (~isempty(canthus_rlList) && canthus_rlList.getLength > 0)
        canthus_rlE1 = canthus_rlList.item(0);
        T.face(i).canthus_rl = [str2num(XML_get_attr_value(canthus_rlE1, 'col')); str2num(XML_get_attr_value(canthus_rlE1, 'row'))];
    end;
    
    canthus_lrList = faceEl.getElementsByTagName('canthus_lr');
    if (~isempty(canthus_lrList) && canthus_lrList.getLength > 0)
        canthus_lrE1 = canthus_lrList.item(0);
        T.face(i).canthus_lr = [str2num(XML_get_attr_value(canthus_lrE1, 'col')); str2num(XML_get_attr_value(canthus_lrE1, 'row'))];
    end;
    
    canthus_llList = faceEl.getElementsByTagName('canthus_ll');
    if (~isempty(canthus_llList) && canthus_llList.getLength > 0)
        canthus_llE1 = canthus_llList.item(0);
        T.face(i).canthus_ll = [str2num(XML_get_attr_value(canthus_llE1, 'col')); str2num(XML_get_attr_value(canthus_llE1, 'row'))];
    end;
    
	mouthList = faceEl.getElementsByTagName('mouth');
	if (~isempty(mouthList) && mouthList.getLength > 0)
		mouthEl = mouthList.item(0);
		T.face(i).mouth = [str2num(XML_get_attr_value(mouthEl, 'col')); str2num(XML_get_attr_value(mouthEl, 'row'))];
    end
    
    mouth_corner_rList = faceEl.getElementsByTagName('mouth_corner_r');
    if (~isempty(mouth_corner_rList) && mouth_corner_rList.getLength > 0)
		mouth_corner_r = mouth_corner_rList.item(0);
		T.face(i).mouth_corner_r = [str2num(XML_get_attr_value(mouth_corner_r, 'col')); str2num(XML_get_attr_value(mouth_corner_r, 'row'))];
    end
    
    mouth_corner_lList = faceEl.getElementsByTagName('mouth_corner_l');
    if (~isempty(mouth_corner_lList) && mouth_corner_lList.getLength > 0)
		mouth_corner_l = mouth_corner_lList.item(0);
		T.face(i).mouth_corner_l = [str2num(XML_get_attr_value(mouth_corner_l, 'col')); str2num(XML_get_attr_value(mouth_corner_l, 'row'))];
    end
	
	noseList = faceEl.getElementsByTagName('nose');
	if (~isempty(noseList) && noseList.getLength > 0)
		noseEl = noseList.item(0);
		T.face(i).nose = [str2num(XML_get_attr_value(noseEl, 'col')); str2num(XML_get_attr_value(noseEl, 'row'))];
    end
    
    detList = faceEl.getElementsByTagName('detection');
 	if (~isempty(detList) && detList.getLength > 0)
		detEl = detList.item(0);
		T.face(i).detection.BoundingBox.TopLeftRow = str2num(XML_get_attr_value(detEl,  'TopLeftRow')); 
		T.face(i).detection.BoundingBox.TopLeftCol = str2num(XML_get_attr_value(detEl,  'TopLeftCol')); 
		T.face(i).detection.BoundingBox.TopRightRow = str2num(XML_get_attr_value(detEl, 'TopRightRow')); 
		T.face(i).detection.BoundingBox.TopRightCol = str2num(XML_get_attr_value(detEl, 'TopRightCol')); 
		T.face(i).detection.BoundingBox.BotLeftRow = str2num(XML_get_attr_value(detEl,  'BotLeftRow')); 
		T.face(i).detection.BoundingBox.BotLeftCol = str2num(XML_get_attr_value(detEl,  'BotLeftCol')); 
		T.face(i).detection.BoundingBox.BotRightRow = str2num(XML_get_attr_value(detEl, 'BotRightRow')); 
		T.face(i).detection.BoundingBox.BotRightCol = str2num(XML_get_attr_value(detEl, 'BotRightCol'));         
    end   
end;

return

% end of face_XML_read.m
