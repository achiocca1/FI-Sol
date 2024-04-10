%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FileName = 'Results.csv'; % File name of Ansys results
FileNameCoord = 'COORD.csv'; % File name of nodal coordinates
LoadSteps = [2,1]; % Bending Blocks
R = -1; % Load ratio
kFin = 0.668; % Findley material parameter
OmegaSimm = 2*pi; % Symmetry on omega
%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% START DATA READING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RESULTS = importdata(FileName);
Files = find(RESULTS(:,2) == 1);
count = diff(RESULTS(:,2)) == 0;
blocklength = length(find(count == 0)) + 1;
% Variable PREALLOCATION FI
RESFI = zeros(length(Files), 12); % Matrix containing results
% Read file of coordinates to plot the critical plane factors %
Coord = importdata(FileNameCoord);
% Import data from files
index = 1; % Necessary for good preallocation of matrix seizes
% Stress and Strain matrix creation trough cell array
% Preallocate cell memory
E = cell(1, blocklength);
S = cell(1, blocklength);
% Filter for simulations with just on load step
if blocklength == 1
    E{1, LoadSteps(2)} = R*E{1, LoadSteps(1)};
    S{1, LoadSteps(2)} = R*S{1, LoadSteps(1)};
else
end

f = waitbar(0,'Simulation in progress...'); % PROGRESS BAR Option

for q = 1 : length(Files)

    file = RESULTS(Files(q):length(Files):(blocklength - 1)*length(Files) + Files(q), 3:20);
    nodenumber = RESULTS(Files(q), 1); % get the node number
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%% FINISH DATA READING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for j = 1 : blocklength
        E0 = [file(j,1) + file(j,7), (file(j,4) + file(j,10))/2, (file(j,6) + file(j,12))/2;
            (file(j,4) + file(j,10))/2, file(j,2) + file(j,8), (file(j,5) + file(j,11))/2;
            (file(j,6) + file(j,12))/2, (file(j,5) + file(j,11))/2, file(j,3) + file(j,9)];
        E{j} = E0; % Strain Tensor
    end

    for j = 1 : blocklength
        S0 = [file(j,13), file(j,16), file(j,18);
            file(j,16), file(j,14), file(j,17);
            file(j,18), file(j,17), file(j,15)];
        S{j} = S0; % Strain Tensor
    end


    % Small values to zero
    S{LoadSteps(1)}(abs(S{LoadSteps(1)}) < 0.005*max(max(abs(S{LoadSteps(1)})))) = 0;
    S{LoadSteps(2)}(abs(S{LoadSteps(2)}) < 0.005*max(max(abs(S{LoadSteps(2)})))) = 0;

    DeltaS = S{LoadSteps(1)} - S{LoadSteps(2)};


    % Calculate and sort principal strains
    [V, D] = eig(DeltaS);
    D00 = diag(D);
    [D0, ind] = sort(D00, 'descend');
    V0 = V(:,ind);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!!!!!!!!!!!!
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!!!!!!!!!!!!
    %%%%%%%%%%%%%%%%%%%%%%%  ATTENTION  %%%%%%%%%%%%%%%%%%%%%%%%%%%%!!!!!!!!!!!!
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!!!!!!!!!!!!
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!!!!!!!!!!!!
    if abs(abs(D0(1)) - abs(D0(3)))/(abs(max(D0))) > 1e-4 && det(V0) < 0
    elseif abs(abs(D0(1)) - abs(D0(3)))/(abs(max(D0))) < 1e-4 && det(V0) > 0 % PURE TORSION
    else
        DeltaS = S{LoadSteps(2)} - S{LoadSteps(1)};
        [V0, D] = eig(DeltaS);
        D0 = diag(D);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!!!!!!!!!!!!
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!!!!!!!!!!!!
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!!!!!!!!!!!!
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!!!!!!!!!!!!
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!!!!!!!!!!!!

    % Calculate and sort principal stresses S1
    [W1, G1] = eig(S{LoadSteps(1)});
    G01 = diag(G1);
    % Calculate and sort principal stresses S2
    [W2, G2] = eig(S{LoadSteps(2)});
    G02 = diag(G2);

    EigenvalS = sort(D0, 'descend');
    DS1 = max(D0);
    DS3 = min(D0);
    S11 = max(G01);
    S13 = min(G01);
    S21 = max(G02);
    S23 = min(G02);

    tolerance = 0.01*max(abs(D0));
    %% Warning and error messages %%
    if abs(EigenvalS(1) - EigenvalS(2)) < tolerance && abs(EigenvalS(2) - EigenvalS(3)) < tolerance && abs(EigenvalS(1) - EigenvalS(3)) < tolerance
        error('In this case FS = 0! No critical plane exists')
    elseif abs(EigenvalS(1) - EigenvalS(2)) < tolerance || abs(EigenvalS(2) - EigenvalS(3)) < tolerance || abs(EigenvalS(1) - EigenvalS(3)) < tolerance
        warning('The value of FS found is correct! However, infinite critical plane orientations may exist in addition to those found. Use the standard scanning plane method to look for all the possible plane orientations.')
    end

    a = (DS1 - DS3)/4;

    %% S1
    b = 0.5*(S11 + S13);
    c = 0.5*(S11 - S13);
    k = kFin;

    FI_MAX_1 = b*k + sqrt(a^2 + (c^2*k^2));

    x = (c*k)/(sqrt(a^2 + (c^2.*k^2)));
    y = a/(sqrt(a^2 + (c^2.*k^2)));
    OMEGA_MAX = 0.5*atan(y/x) + OmegaSimm;

    % + OMEGA max angle
    % Calculate the angles base on the rotation R = Rotz(Psi)*Roty(Theta)
    RotY = [cos(OMEGA_MAX)   0   sin(OMEGA_MAX);
        0        1        0;
        -sin(OMEGA_MAX)   0   cos(OMEGA_MAX)];


    Matr = V0*RotY;
    %Matr(:,3) = -Matr(:,3);

    Psi = [atan(Matr(2,3)/Matr(1,3)), atan(Matr(2,3)/Matr(1,3)) + pi, atan(Matr(2,3)/Matr(1,3)) + 2*pi];
    idx = ( Psi>= 0) & ( Psi<= 2*pi);
    Psi = Psi(idx);

    for i = 1:length(Psi)
        SINTHETA = Matr(2,3)/sin(Psi(i));
        Theta = atan2(SINTHETA, Matr(3,3));
        if abs(cos(Psi(i))*SINTHETA - Matr(1,3)) < 1e-10
            Theta_1 = Theta;
            Psi_1 = Psi(i);
            break
        else
        end
        Theta_1 = atan2(sqrt(Matr(1,3)^2+Matr(2,3)^2), Matr(3,3));
        Psi_1 = atan2(Matr(2,3), Matr(1,3));
    end


    % - OMEGA max angle
    OMEGA_MAX = - OMEGA_MAX;
    % Calculate the angles base on the rotation R = Rotz(Psi)*Roty(Theta)
    RotY = [cos(OMEGA_MAX)   0   sin(OMEGA_MAX);
        0        1        0;
        -sin(OMEGA_MAX)   0   cos(OMEGA_MAX)];

    Matr = V0*RotY;

    Psi = [atan(Matr(2,3)/Matr(1,3)), atan(Matr(2,3)/Matr(1,3)) + pi, atan(Matr(2,3)/Matr(1,3)) + 2*pi];
    idx = ( Psi>= 0) & ( Psi<= 2*pi);
    Psi = Psi(idx);

    for i = 1:length(Psi)
        SINTHETA = Matr(2,3)/sin(Psi(i));
        Theta = atan2(SINTHETA, Matr(3,3));
        if abs(cos(Psi(i))*SINTHETA - Matr(1,3)) < 1e-10
            Theta_2 = Theta;
            Psi_2 = Psi(i);
            break
        else
        end
        Theta_2 = atan2(sqrt(Matr(1,3)^2+Matr(2,3)^2), Matr(3,3));
        Psi_2 = atan2(Matr(2,3), Matr(1,3));
    end


    %% S2
    b = 0.5*(S21 + S23);
    c = 0.5*(S21 - S23);

    FI_MAX_2 = b*k + sqrt(a^2 + (c^2*k^2));
    x = c*k/(sqrt(a^2 + (c^2.*k^2)));
    y = a/(sqrt(a^2 + (c^2.*k^2)));
    OMEGA_MAX = 0.5*atan(y/x) + OmegaSimm;

    % + OMEGA max angle
    % Calculate the angles base on the rotation R = Rotz(Psi)*Roty(Theta)
    RotY = [cos(OMEGA_MAX)   0   sin(OMEGA_MAX);
        0        1        0;
        -sin(OMEGA_MAX)   0   cos(OMEGA_MAX)];

    Matr = V0*RotY;

    Psi = [atan(Matr(2,3)/Matr(1,3)), atan(Matr(2,3)/Matr(1,3)) + pi, atan(Matr(2,3)/Matr(1,3)) + 2*pi];
    idx = ( Psi>= 0) & ( Psi<= 2*pi);
    Psi = Psi(idx);

    for i = 1:length(Psi)
        SINTHETA = Matr(2,3)/sin(Psi(i));
        Theta = atan2(SINTHETA, Matr(3,3));
        if abs(cos(Psi(i))*SINTHETA - Matr(1,3)) < 1e-10
            Theta_3 = Theta;
            Psi_3 = Psi(i);
            break
        else
        end
        Theta_3 = atan2(sqrt(Matr(1,3)^2+Matr(2,3)^2), Matr(3,3));
        Psi_3 = atan2(Matr(2,3), Matr(1,3));
    end

    % - OMEGA max angle
    OMEGA_MAX = - OMEGA_MAX;
    % Calculate the angles base on the rotation R = Rotz(Psi)*Roty(Theta)
    RotY = [cos(OMEGA_MAX)   0   sin(OMEGA_MAX);
        0        1        0;
        -sin(OMEGA_MAX)   0   cos(OMEGA_MAX)];

    Matr = V0*RotY;

    Psi = [atan(Matr(2,3)/Matr(1,3)), atan(Matr(2,3)/Matr(1,3)) + pi, atan(Matr(2,3)/Matr(1,3)) + 2*pi];
    idx = ( Psi>= 0) & ( Psi<= 2*pi);
    Psi = Psi(idx);

    for i = 1:length(Psi)
        SINTHETA = Matr(2,3)/sin(Psi(i));
        Theta = atan2(SINTHETA, Matr(3,3));
        if abs(cos(Psi(i))*SINTHETA - Matr(1,3)) < 1e-10
            Theta_4 = Theta;
            Psi_4 = Psi(i);
            break
        else
        end
        Theta_4 = atan2(sqrt(Matr(1,3)^2+Matr(2,3)^2), Matr(3,3));
        Psi_4 = atan2(Matr(2,3), Matr(1,3));
    end


    %% RESULTS FI %%
    RESFI(q, :) = [nodenumber max(FI_MAX_1,FI_MAX_2) FI_MAX_1 FI_MAX_2 Theta_1 Theta_2 Theta_3 Theta_4 Psi_1 Psi_2 Psi_3 Psi_4];
    index = index + 1;

    %% PROGRESS BAR
    waitbar(q/length(Files))
    
end

[FI, IndexFI] = max(RESFI(:,2));

X = ['The critical plane factor for the node ',  num2str(RESFI(IndexFI,1)), ' is FI = ', num2str(RESFI(IndexFI,2)), ' MPa '];
disp(X)
if RESFI(IndexFI,3) > RESFI(IndexFI,4)
    X = ['The first critical plane orientation is given by the angles θ = ', num2str(RESFI(IndexFI,5)), ' and Ψ = ', num2str(RESFI(IndexFI,9))];
    Y = ['The second critical plane orientation is given by the angles θ = ', num2str(RESFI(IndexFI,6)), ' and Ψ = ', num2str(RESFI(IndexFI,10))];
else
    X = ['The first critical plane orientation is given by the angles θ = ', num2str(RESFI(IndexFI,7)), ' and Ψ = ', num2str(RESFI(IndexFI,11))];
    Y = ['The second critical plane orientation is given by the angles θ = ', num2str(RESFI(IndexFI,8)), ' and Ψ = ', num2str(RESFI(IndexFI,12))];
end
disp(X)
disp(Y)
delete(f) % wait bar deleted