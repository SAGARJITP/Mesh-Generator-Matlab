



 function varargout = MeshGenerator_GUI(varargin)
    % MESHGENERATOR_GUI MATLAB code for MeshGenerator_GUI.fig
    %      MESHGENERATOR_GUI, by itself, creates a new MESHGENERATOR_GUI or raises the existing
    %      singleton*.
    %
    %      H = MESHGENERATOR_GUI returns the handle to a new MESHGENERATOR_GUI or the handle to
    %      the existing singleton*.
    %
    %      MESHGENERATOR_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
    %      function named CALLBACK in MESHGENERATOR_GUI.M with the given input arguments.
    %
    %      MESHGENERATOR_GUI('Property','Value',...) creates a new MESHGENERATOR_GUI or raises the
    %      existing singleton*.  Starting from the left, property value pairs are
    %      applied to the GUI before MeshGenerator_GUI_OpeningFcn gets called.  An
    %      unrecognized property name or invalid value makes property application
    %      stop.  All inputs are passed to MeshGenerator_GUI_OpeningFcn via varargin.
    %
    %      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
    %      instance to run (singleton)".
    %
    % See also: GUIDE, GUIDATA, GUIHANDLES

    % Edit the above text to modify the response to help MeshGenerator_GUI

    % Last Modified by GUIDE v2.5 08-Dec-2017 17:06:31

    % Begin initialization code - DO NOT EDIT
    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
                       'gui_Singleton',  gui_Singleton, ...
                       'gui_OpeningFcn', @MeshGenerator_GUI_OpeningFcn, ...
                       'gui_OutputFcn',  @MeshGenerator_GUI_OutputFcn, ...
                       'gui_LayoutFcn',  [] , ...
                       'gui_Callback',   []);
    if nargin && ischar(varargin{1})
        gui_State.gui_Callback = str2func(varargin{1});
    end

    if nargout
        [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
    else
        gui_mainfcn(gui_State, varargin{:});
    end
    % End initialization code - DO NOT EDIT
end

% --- Executes just before MeshGenerator_GUI is made visible.
function MeshGenerator_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % varargin   command line arguments to MeshGenerator_GUI (see VARARGIN)
    % Choose default command line output for MeshGenerator_GUI
   
    
    load patch_control
    handles.output = hObject;

    % Create the data to plot
    handles.original_patches = patch_control;
    handles.patches = patch_control;
    Q=bfunc(handles.patches);

    % Set the current data value.
    handles.original_data= Q;
    handles.current_data=Q;
   
    xc=handles.patches(:,:,1);
    yc=handles.patches(:,:,2);
    zc=handles.patches(:,:,3);
    datacursormode on
    axes(handles.axes3);
    plot3(xc,yc,zc,'.','MarkerSize',25)
    view(3); box;  view(21,19)
    
    % Update handles structure
    guidata(hObject, handles);

    % UIWAIT makes MeshGenerator_GUI wait for user response (see UIRESUME)
    % uiwait(handles.figure1);
end

    % --- Outputs from this function are returned to the command line.
function varargout = MeshGenerator_GUI_OutputFcn(hObject, eventdata, handles) 
    % varargout  cell array for returning output args (see VARARGOUT);
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Get default command line output from handles structure
    varargout{1} = handles.output;
end

%This button creates mesh
function CreateMesh_Callback(hObject, eventdata, handles)
   
    file = get(handles.GeometryFile,'String');
    elem_size = get(handles.MeshElementSize,'String');
    elem_size = str2double(elem_size);
    set(handles.DisplayStatus,'String', 'Optimizing Mesh...');
    
    %Points, faces and boundaries for the geometry
    [g_points, g_faces, b] = GetInfo(file);
    
    bbox = [min(g_points(:,1)), min(g_points(:,2));max(g_points(:,1)), max(g_points(:,2))];
    ttol = 0.7;
    delta = 0.2;
    
    %CreateMesh(ttol,bbox,elem_size,delta,g_points,b,handles.Mesh);
    %........................................................................................
    dptol = 0.01;
    L0 = 1.2*elem_size;
    
    p = CreateInitialNodes(bbox,elem_size,g_points,b); %Generates initial nodes on the surface
    pfix = FixedNodes(g_points,b,elem_size);
    p = MeshPoints(p,pfix,elem_size);
    axes(handles.Mesh);
    scatter(p(:,1),p(:,2),'x'); %Visualize initial nodes
    N = size(p,1); %No of nodes
    p_old = inf;
    Algorithm = handles.Alg;
    
    iter = 1; %Max number of iterations
    while (iter < 1000)
    
        if ( max(sqrt((p-p_old).^2),2) > ttol)
            p_old = p;
            switch Algorithm
                case 1
                   t = delaunay(p(:,1:2)); %Delaunay triangulation MATLAB                  
                case 2
                   t = createdelaunaytriangles(p(:,1:3)); %Delaunay triangulation using our method
            end
            p_mid = ((p(t(:,1),:)) + (p(t(:,2),:)) + (p(t(:,3),:)))/3;
            t = t(Distance(p_mid,g_points,b) < 0,:);
            edges = [t(:,[1,2]); t(:,[1,3]); t(:,[2,3]);]; %get all the edges of triangular mesh
            edges=unique(edges,'rows'); % Remove duplicate edges
            trimesh(t,p(:,1),p(:,2),zeros(N,1))
            %triplot(handles.Mesh,t,p(:,1),p(:,2))
            view(2),axis equal,axis on,drawnow
        end 
        
        edge_vec = (p(edges(:,1),:) - p(edges(:,2),:)); %edge_vec stors the vector for each edge
        L = sqrt(sum(edge_vec.^2,2)); %calculate length of each edges
        F = (L0-L);
        %F = max(L0-L,0); %Calculte the magnitude of force on each bar
        %F(1:size(pfix,1),:)=0;
        F_vec = F.*edge_vec./(L);
        
        
        %update the position of p
        for i = 1:size(edge_vec,1)  
            if (edges(i,1) > size(pfix,1))
                 p(edges(i,1),:) =  p(edges(i,1),:) + delta*F_vec(i,:);  
            end
                  
        end    
        
        if (max(sqrt(sum(delta*F_vec.^2,2)))) <0.035
            break;
        end   
        
        if strcmp(get(handles.DisplayStatus,'String'),'Mesh Optimized')
            break;
        end
        
        iter = iter + 1;
    end
    %.......................................................................................
       
    %can use handles.p and handles.t to access p and t resp after this lines of codes
    handles.p = p;
    handles.t = t;
    handles.file=file;
    handles.g_points=g_points;
    
    guidata(hObject,handles);
    
    %Will print when iterations are over
    set(handles.DisplayStatus,'String', 'Mesh Optimized');
    
end

% This function exports the mesh file in obj format
function ExportMesh_Callback(hObject, eventdata, handles)
    % hObject    handle to ExportMesh (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
     
    
    p = handles.p;
    t = handles.t;
    file = get(handles.GeometryFile,'String');
    g_points = handles.g_points;
    
    P=retriangulate(file,p,g_points);
    
    p=P;
    
    filename = get(handles.ExportFile,'String');
    fileID = fopen(filename,'w+');
    for i = 1:size(p,1)
        fprintf(fileID,'%s %f %f %f\r\n','v',p(i,1),p(i,2),p(i,3));
    end
    for i = 1:size(t,1)
        fprintf(fileID,'%s %d %d %d\r\n','f',t(i,1),t(i,2),t(i,3));
    end
    
    fclose(fileID);
    
end

%text editer for mesh element size
function MeshElementSize_Callback(hObject, eventdata, handles)
    
    % Hints: get(hObject,'String') returns contents of MeshElementSize as text
    %        str2double(get(hObject,'String')) returns contents of MeshElementSize as a double
end

% Text for inputing mesh element size
function MeshElementSize_CreateFcn(hObject, eventdata, handles)
  
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end

%This function imports geometry file 
function ImportGeometry_Callback(hObject, eventdata, handles)
    [filename, filepath] = uigetfile;
    set(handles.GeometryFile,'String',filename);
end

% This function import the geometry file 
function GeometryFile_CreateFcn(hObject, eventdata, handles)
   
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

end


% Display the status of mesh generation
function DisplayStatus_Callback(hObject, eventdata, handles)

end

% Display the status of mesh generation
function DisplayStatus_CreateFcn(hObject, eventdata, handles)
    
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end

% This button will stop optimizing the mesh
function StopOpt_Callback(hObject, eventdata, handles)
    set(handles.DisplayStatus,'String', 'Mesh Optimized');
end

%The text editer for name of the file to be exported for mesh
function ExportFile_Callback(hObject, eventdata, handles)
    
end

% --- Executes during object creation, after setting all properties.
function ExportFile_CreateFcn(hObject, eventdata, handles)
   
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end

% Decides which algorithm to use
function Alg = Algorithm_Callback(hObject, eventdata, handles)
    % hObject    handle to Algorithm (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Determine the selected data set.
     str = get(hObject, 'String');
     val = get(hObject,'Value');
     Alg = 1;
     switch str{val}  
         case 'MATLAB Delaunay Algorithm'
             %handles.Algorithm = 1;
             Alg = 1;
         case 'Our Delaunay Algorithm'
             %handles.Algorithm = 2;
             Alg = 2;
     end
     handles.Alg = Alg;
     guidata(hObject,handles);
    
    % Hints: contents = cellstr(get(hObject,'String')) returns Algorithm contents as cell array
    %        contents{get(hObject,'Value')} returns selected item from Algorithm
end

% Algorithm to use
function Algorithm_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to Algorithm (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: popupmenu controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end


% --- Executes on button press in Remesh.
function Remesh_Callback(hObject, eventdata, handles)
    % hObject    handle to Remesh (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    p = handles.p;
    t = handles.t;
    xc = str2double(get(handles.remeshx,'string'));
    yc = str2double(get(handles.remeshy,'string'));
    radius = str2double(get(handles.remesh_radius,'string'));
    elem_size = str2double(get(handles.MeshElementSize,'String'));
    
    centre = [xc,yc];
    
    idx = 1;
    for i = 1:size(p,1)      
        d = EuclideanDistance(centre,p(i,:));
        if (d <= radius)
            p_inner(idx,:) = p(i,:);
            indices(idx)= i;
            idx = idx + 1;
        end      
    end
    
   p(indices,:) = [];
   %scatter(p_inner(:,1),p_inner(:,2),'.');
   
   
   %t = delaunay(p(:,1:2));
   
   %trimesh(t,p(:,1),p(:,2),zeros(size(p,1),1));
   %view(2),axis equal,axis on,drawnow
    
    bbox = [min(p_inner(:,1)), min(p_inner(:,2));max(p_inner(:,1)), max(p_inner(:,2))];
    b = boundary(p_inner(:,1),p_inner(:,2),0); %boundary vector
    %plot(p_inner(b,1),p_inner(b,2));
    p_small = CreateInitialNodes(bbox,elem_size,p_inner,b);
   
    p = [p;p_small];
    %scatter(p(:,1),p(:,2));
    t = delaunay(p(:,1:2));
    trimesh(t,p(:,1),p(:,2),zeros(size(p,1),1));
    view(2),axis equal,axis on,drawnow;      
    
    handles.p = p;
    handles.t = t;      
    
    guidata(hObject,handles);
end

% point x cordinaty --remesh
function remeshx_Callback(hObject, eventdata, handles)


end

% point x cordinaty --remesh
function remeshx_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end

% point y cordinaty --remesh
function remeshy_Callback(hObject, eventdata, handles)

end

% point y cordinaty --remesh
function remeshy_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end

%Edit text for remesh radius
function remesh_radius_Callback(hObject, eventdata, handles)

end

%Edit text for remesh radius
function remesh_radius_CreateFcn(hObject, eventdata, handles)
    
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%surface plots



% --- Executes on button press in Rotate.
function Rotate_Callback(hObject, eventdata, handles)
% hObject    handle to Rotate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rotate3d on;
end

% --- Executes on button press in Zoom.
function Zoom_Callback(hObject, eventdata, handles)
% hObject    handle to Zoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
zoom on;
end

% --- Executes on button press in Pan.
function Pan_Callback(hObject, eventdata, handles)
% hObject    handle to Pan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

end

function xField_Callback(hObject, eventdata, handles)
% hObject    handle to xField (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xField as text
%        str2double(get(hObject,'String')) returns contents of xField as a double
end

% --- Executes during object creation, after setting all properties.
function xField_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xField (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

end

function yField_Callback(hObject, eventdata, handles)
% hObject    handle to yField (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yField as text
%        str2double(get(hObject,'String')) returns contents of yField as a double
end

% --- Executes during object creation, after setting all properties.
function yField_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to yField (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: edit controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end


function zField_Callback(hObject, eventdata, handles)
% hObject    handle to zField (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zField as text
%        str2double(get(hObject,'String')) returns contents of zField as a double
end

% --- Executes during object creation, after setting all properties.
function zField_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zField (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on button press in setPoints.
function setPoints_Callback(hObject, eventdata, handles)
% hObject    handle to setPoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dcm_obj = datacursormode(handles.figure1);
c_info = getCursorInfo(dcm_obj);
patch=handles.patches;
row=size(patch,1);
col=size(patch,2);
% counter=1;
for i=1:row
    for j=1:col
            check1=patch(i,j,1)==c_info.Position(1);
            check2=patch(i,j,2)==c_info.Position(2);
            check3=patch(i,j,3)==c_info.Position(3);
            total = check1+check2+check3;
        if (total==3)
            ind=[i,j];
            break
        end
%         counter = counter +1;
    end
    if (total==3)
        break
    end
end
% char(get(handles.xField,'String'))
handles.patches(ind(1),ind(2), 1) = str2double(get(handles.xField,'String')); %handles.newpoint(1,1);
handles.patches(ind(1),ind(2), 2) = str2double(get(handles.yField,'String'));%handles.newpoint(1,2);
handles.patches(ind(1),ind(2), 3) = str2double(get(handles.zField,'String'));%handles.newpoint(1,3);
handles.current_data = bfunc(handles.patches);
guidata(hObject,handles)
end


% --- Executes on button press in RePlot.
function RePlot_Callback(hObject, eventdata, handles)
% hObject    handle to RePlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% hObject    handle to RePlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla(handles.axes3);
% x=handles.current_data(:,:,1);
% y=handles.current_data(:,:,2);
% z=handles.current_data(:,:,3);
% surface(x,y,z);

xc=handles.patches(:,:,1);
yc=handles.patches(:,:,2);
zc=handles.patches(:,:,3);
axes(handles.axes3);
plot3(xc,yc,zc,'.','MarkerSize',25)

end

% --- Executes on button press in Reset.
function Reset_Callback(hObject, eventdata, handles)
    
    
% hObject    handle to Reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.patches = handles.original_patches;
handles.current_data = handles.original_data;
cla

xc=handles.patches(:,:,1);
yc=handles.patches(:,:,2);
zc=handles.patches(:,:,3);
axes(handles.axes3);
plot3(handles.axes3,xc,yc,zc,'.','MarkerSize',25)
view(3); box;  view(21,19)
guidata(hObject,handles)

end

% --- Executes on button press in Toggle.
function Toggle_Callback(hObject, eventdata, handles)
% hObject    handle to Toggle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Toggle
button_state = get(hObject,'Value');
if button_state == get(hObject,'Max')
    datacursormode off
    cla(handles.axes3);
    
    
    x=handles.current_data(:,:,1);
    y=handles.current_data(:,:,2);
    z=handles.current_data(:,:,3);
    axes(handles.axes3);
    surface(x,y,z);
    shading interp
    hold on
    
    xc=handles.patches(:,:,1);
    yc=handles.patches(:,:,2);
    zc=handles.patches(:,:,3);
    axes(handles.axes3);
    plot3(xc,yc,zc,'.','MarkerSize',25)
elseif button_state == get(hObject,'Min')  
    datacursormode on
    cla(handles.axes3);
    xc=handles.patches(:,:,1);
    yc=handles.patches(:,:,2);
    zc=handles.patches(:,:,3);
    axes(handles.axes3);
    plot3(xc,yc,zc,'.','MarkerSize',25)
	
end
end


% --- Executes on button press in Export.
function Export_Callback(hObject, eventdata, handles)
% hObject    handle to Export (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
x=handles.current_data(:,:,1);
y=handles.current_data(:,:,2);
z=handles.current_data(:,:,3);
stlwrite('test.stl',x,y,z,'mode','ascii');

%set end condition to to pull stop reading from stl file
fid=fopen('test.stl','a');
fprintf(fid,'Theend\n');
fclose(fid);

%prepare for obj
fid=fopen('test.stl');
line = fgetl(fid);
i=1;
while 0==strcmp(line,'Theend')

    line = fgetl(fid);
    words = textscan(line,'%s');
    
    if strcmp(words{1,1}{1,1}, 'vertex')
      
      vx=str2double(words{1,1}{2,1});
      vy=str2double(words{1,1}{3,1});
      vz=str2double(words{1,1}{4,1});
      vertexes(i,:) = [vx, vy, vz];
      i=i+1;
    end
end
fclose(fid);

[vertex_Set,ia,ic] = unique(vertexes(:,1:3),'rows');

j=1;
for i=1:length(vertexes)
    if rem(i,3)==0
        [~,indx1]=ismember(vertexes(i-2,:),vertex_Set,'rows');
        [~,indx2]=ismember(vertexes(i-1,:),vertex_Set,'rows');
        [~,indx3]=ismember(vertexes(i,:),vertex_Set,'rows');
        faces(j,:) = [indx1, indx2, indx3];
        j=j+1;
    end
end

[rFaces,~]=size(faces);
[rVertexes,~]=size(vertex_Set);

fid2=fopen('test.obj','w');
k=1;
for i=1:(rFaces+rVertexes)
    if i <= rVertexes
        fprintf(fid2,'v ');
        fprintf(fid2,num2str(vertex_Set(i,:)));
        fprintf(fid2, '\n');
    else
        fprintf(fid2,'f ');
        fprintf(fid2,num2str(faces(k,:)));
        fprintf(fid2, '\n');
        k=k+1;
    end      
end
fclose(fid2);
end
