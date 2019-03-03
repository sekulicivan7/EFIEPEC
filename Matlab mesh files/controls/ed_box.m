% Create edit box in current window
% 
% ed_box(h_pos,v_pos,type,var,itx,text,control_dep)
%
% h_pos	= Horizontal position of 1st character of variable name (characters)
%	  For type = 0 or 1, 	optionally  h_pos(2) = Size (characters) of text window
%	  			optionally  h_pos(3) = Size (characters) of edit window
%				for sliders h_pos(2) = Width (characters)
%					    h_pos(3) = Height (characters)
% v_pos	= Vertical position of 1st character of text (lines), from figure top
% type	= control type
%		0 = Numeric variable edit window
%		1 = String  variable edit window.
%		2 = Push button
%		3 = Static text	
%		4 = Check box
%		5 = Frame (h_pos,v_pos,type,var,itx,text)->(x0,y0,5,width,height,[color])
%			color = [R_bg G_bg B_bg R_fg G_fg B_fg]
%			defaults: if lenght(color)<6, fg=[0 0 0]
%				  if lenght(color)<3, bg=get(gcf,'Color')
%		6 = Array of radio buttons: (h_pos, v_pos, itx, text) must be arrays
%		7 = Pop-up menu
%		8 = Matrix input
%		9 = Slider
% var	= variable name (string) for edit windows, or callback (string) for push buttons
% itx	= initial contents of the edit field (string or number), or text color [R G B],
%	  or [initial_value, [list of values for each radio button]]
%	  or [initial_value, min, max] for sliders
%	  For type=7, itx = column array of row-strings. 1st row = name of options file
%	  2ond row = default value
% text	= text
% control_dep = Optional argument indicating controls with visibility dependant from this
%		Only for ckeck(4) and radio(6) controls
%		Control_dep = [handle1 flag1; handle2 flag2; handle3 flag3;...]
%
% version 2.10, Juan M. Rius, Sept 1996

function [ed,var]=ed_box(h_pos,v_pos,type,var,itx,text,control_dep)

screen = get(0, 'ScreenSize');
CHH = 20*(screen(3)/1024); CHW=CHH/2;	% Character size = 20 pixels height
SEP = CHH*1.25;				% Interline separation (pixels)
wpos=get(gcf,'Position'); H_SIZE=wpos(3)/CHW; V_SIZE=wpos(4)/SEP;
back = get(gcf,'Color');

v_pos = V_SIZE-v_pos;	% Vertical position from bottom instead of top        
ED_N=10;	% Size (characters) of number edit window (default)
SED_N=16;	% Size (characters) of string edit window (default)
STX_N=20;	% Size (characters) of text window (default)
NTX_N=20;	% Size (characters) of text window (default)

if type==9,
	ed_pos = [h_pos(1)*CHW v_pos*SEP h_pos(2)*CHW h_pos(3)*SEP-CHH];
	smin = num2str(itx(2)); lmin = (length(smin))*CHW;
	smax = num2str(itx(3)); lmax = (length(smax))*CHW;
	tit = num2str(itx(1));  ltex = (length(tit)+2)*CHW;
				lnom = (length(text))*CHW;
	
	umin = uicontrol(gcf,'Style','Text','BackgroundColor',back,...
	'Position',[ed_pos(1)                ed_pos(2) lmin CHH],'String',smin);
	
	umax = uicontrol(gcf,'Style','Text','BackgroundColor',back,...
	'Position',[ed_pos(1)+ed_pos(3)-lmax ed_pos(2) lmax CHH],'String',smax);

	uit = uicontrol(gcf,'Style','Text','BackgroundColor',back,...
	'Position',[ed_pos(1)+lnom+(ed_pos(3)-lnom-ltex)/2 ed_pos(2)+h_pos(3)*SEP-CHH ltex CHH],'String',tit);

	uit2 = uicontrol(gcf,'Style','Text','BackgroundColor',back,...
	'Position',[ed_pos(1)+(ed_pos(3)-lnom-ltex)/2 ed_pos(2)+h_pos(3)*SEP-CHH lnom CHH],'String',text);
	
	ed_pos(1) = ed_pos(1) + lmin; ed_pos(3) = ed_pos(3) - lmin - lmax;
	call_st=['ed = get(gcf,''CurrentObject''); ' var '= get(ed,''Value'');'...
		'set(get(ed,''Userdata''),''String'',num2str(' var '));'];

	usl = uicontrol(gcf,'Style','Slider','Position',ed_pos,...
	'Min',itx(2),'Max',itx(3),'Value',itx(1),'UserData',uit,...
	'CallBack',call_st);
	ed = [usl umin umax uit uit2];
	return;
end

if type==5,
	bgcolor = back; fgcolor = [0 0 0]; % Defaults
	if     length(text)== 6, bgcolor = text(1:3); fgcolor = text(4:6);
	elseif length(text)== 3, bgcolor = text(1:3);
	end

	fr_pos = [h_pos(1)*CHW v_pos*SEP var*CHW itx*SEP];
	fr = uicontrol(gcf,'Style','frame','Position',fr_pos,...
	'BackGroundColor',bgcolor,'ForeGroundColor',fgcolor);
	ed = fr;								% Handle returned
	return;
end

if type==2, 
	pu_pos = [h_pos(1)*CHW v_pos*SEP (2+length(text))*CHW SEP];
	pu = uicontrol(gcf,'Style','push','Position',pu_pos,'String',text,...
	'Callback',var);
	ed = pu;								% Handle returned
	return
end

if type==3
	tit_pos = [h_pos(1)*CHW v_pos*SEP (1+length(text))*CHW CHH];
	tx = uicontrol(gcf,'Style','Text','String',text,'Position',tit_pos,...
	'BackGroundColor',back,'ForeGroundColor',itx,...
	'HorizontalAlignment','Left');
	ed = tx;								% Handle returned
	return
end

if type==4
	cb_pos = [h_pos(1)*CHW v_pos*SEP (length(text)+2)*CHW SEP];
	ed = uicontrol(gcf,'Style','checkbox','Position',cb_pos,'String',text,...
	'BackGroundColor',back,'ForeGroundColor',[0 0 0],...
	'Value',itx);

	if nargin<7, control_dep = [];
	end
	set(ed,'UserData',control_dep);

	n = 1;
	while n < length(control_dep),
		m = n+1; while rem(control_dep(m),1), m=m+1; end 	% m points to the flag
		if ~xor(control_dep(m),itx), set(control_dep(n:m-1),'Visible','On');
		else set(control_dep(n:m-1),'Visible','Off');
		end;
		n = m+1;	% n points to the next handle
	end; 

	call_st = [ var '=' 'get(get(gcf,''CurrentObject''),''Value'')==1; ',...
		    'UserData = get(get(gcf,''CurrentObject''),''UserData''); ',...
		    'ntmp = 1; ',...
		    'while ntmp < length(UserData), ',...
			    'mtmp = ntmp+1; while rem(UserData(mtmp),1), mtmp=mtmp+1; end; ',...
			    'if ~xor(UserData(mtmp),' var '), set(UserData(ntmp:mtmp-1),''Visible'',''On'');  ',...
			    'else set(UserData(ntmp:mtmp-1),''Visible'',''Off'');  ',...
			    'end; ',... % If
			    'ntmp = mtmp+1; ',...
		    'end; ',... % While
		    'clear ntmp mtmp UserData;' ];
	set(ed,'Callback',call_st);
	return
end

if type==6,
	if nargin<7, control_dep = [];
	end
	N = length(h_pos);
	mat = [zeros(1,N+1) control_dep];
	mat(1) = N;	% mat = [N h1 h2 h3 ... hN, control_dep];

	for i=1:N,
		tx = deblank(text(i,:));
		cb_pos = [h_pos(i)*CHW v_pos(i)*SEP (length(tx)+2)*CHW SEP];
		ed = uicontrol(gcf,'Style','radiobutton','Position',cb_pos,...
		'String',tx,'BackGroundColor',back,'ForeGroundColor',[0 0 0],...
		'Max',itx(i+1),'Min',-1);
		mat(i+1) = ed;		% Save all control numbers
	end

	n = N+2;
	while n < length(mat),
		m1 = n+1;  while  rem(mat(m1),1), m1 = m1+1; end	% m1 -> 1st flag
		m2 = m1+1;
		if m2<=length(mat), while ~rem(mat(m2),1), m2 = m2+1; end	% m2 -> 1st handle
		end
		if any(mat(m1:m2-1)==itx(1)), set(mat(n:m1-1),'Visible','on');
		else set(mat(n:m1-1),'Visible','off');
		end
		n = m2;
	end

	for i=1:N
		set(mat(i+1),'UserData',mat);	% Save mat in all controls
		if itx(1) == itx(i+1), set(mat(i+1),'Value',get(mat(i+1),'Max'));
		else set(mat(i+1),'Value',get(mat(i+1),'Min'));
		end

		call_st = [	'UserData = get(get(gcf,''CurrentObject''),''UserData''); '...
				'for ntmp = 2:UserData(1)+1, '...
				'if get(gcf,''CurrentObject'') == UserData(ntmp), '...
				'set(UserData(ntmp),''Value'',get(UserData(ntmp),''Max'')); '...
				'else '...
				'set(UserData(ntmp),''Value'',get(UserData(ntmp),''Min'')); '...
				'end; '...	% if
				'end; '...	% for
				var '= get(get(gcf,''CurrentObject''),''Value''); '...
				'ntmp = UserData(1)+2; ',...
				'while ntmp < length(UserData), ',...
					'mtmp1 = ntmp+1; while rem(UserData(mtmp1),1), mtmp1 = mtmp1+1; end; ',...
					'mtmp2 = mtmp1+1; ',...
					'if mtmp2<=length(UserData), while ~rem(UserData(mtmp2),1), mtmp2 = mtmp2+1; end; ',...
					'end; ',...
					'if any(UserData(mtmp1:mtmp2-1)==' var '), set(UserData(ntmp:mtmp1-1),''Visible'',''on''); ',...
					'else set(UserData(ntmp:mtmp1-1),''Visible'',''off''); ',...
					'end; ',...
					'ntmp = mtmp2; ',...
				'end; ',...	% while
				'clear UserData ntmp mtmp; '...
			   ];
		set(mat(i+1),'CallBack',call_st);
	end
ed = mat;								% Handles returned
return
end

if type==0,
	if length(h_pos) > 1, ntx_n = h_pos(2); else ntx_n = NTX_N; end
	if length(h_pos) > 2, ed_n = h_pos(3); else ed_n = ED_N; end
	ed_pos = [(h_pos(1)+ntx_n+1)*CHW SEP*v_pos ed_n*CHW CHH];
	tx_pos = [h_pos(1)*CHW SEP*v_pos ntx_n*CHW CHH];
end

if type==1 | type==7 | type==8,
	if length(h_pos) > 1, stx_n = h_pos(2); else stx_n = STX_N; end
	if length(h_pos) > 2, sed_n = h_pos(3); else sed_n = SED_N; end
	ed_pos = [(h_pos(1)+stx_n+1)*CHW SEP*v_pos sed_n*CHW CHH];
	tx_pos = [h_pos(1)*CHW SEP*v_pos stx_n*CHW CHH];
end

ed = uicontrol(gcf,'Position',ed_pos,'HorizontalAlignment','Right',...
		   'BackGroundColor',[0 0 1],'ForeGroundColor',[1 1 1]);

if type==0 | type == 1,
	set(ed,'Style','Edit','String',itx);
end

if type == 8,
	if length(itx) == 0, sitx = '[ ]';
	else,
		sitx = '[';
		for i=1:length(itx), sitx = [sitx num2str(itx(i)) ' '];	end
		sitx = sitx(1:length(sitx)-1);
		sitx = [sitx ']'];
	end
	set(ed,'Style','Edit','String',sitx);
end

if type==7,
	set(ed,'Style','popup');
	[m,n]=size(itx);
	ee_tmp = deblank(itx(3,:));
	if strcmp(itx(3,:),itx(2,:)), tmp_val=1; else tmp_val=0; end
	for i=4:m,
		ee_tmp = [ee_tmp '|' deblank(itx(i,:))];
		if strcmp(itx(i,:),itx(2,:)), tmp_val=i-2; end
	end
	set(ed,'String',ee_tmp);
	if tmp_val > 0, set(ed,'Value',tmp_val); end
end

tx = uicontrol(gcf,'Style','Text','String',text,'Position',tx_pos,...
		   'BackGroundColor',back,'HorizontalAlignment','Right');

if type==7,
	call_st = [ var '=' itx(1,:) ';',...
		    var '=deblank(' var '(get(get(gcf,''CurrentObject''),''Value''),:));'];
end

if type==1,
	call_st = [ var '=get(get(gcf,''CurrentObject''),''String'');'];
end

if type==8,
	call_st = [ var '=get(get(gcf,''CurrentObject''),''String''); '...
			var '=eval(' var ');'];
end

if type==0,
	call_st = [ var '=str2num(get(get(gcf,''CurrentObject''),''String''));'];
end

ed = [tx ed];									% Handles returned
set(ed(2),'Callback',call_st);
