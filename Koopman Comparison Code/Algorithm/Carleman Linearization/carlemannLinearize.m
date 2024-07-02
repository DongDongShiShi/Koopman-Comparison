function [Acoeffs, truncatedVars, truncatedNumVars] = carlemannLinearize(order,f,omega)
%Carlemann Linearization

%Choose the order of truncation
%order = 2;

%Initialize the starting variables to match those of LB
oldVars = f;

%Initialize the current order of the variables
orderVars = max(polynomialDegree(f));

%Initialize the operator to the collision term of LB
C = omega;

%Initialize the number of variables Q
Q = length(omega);
numVars = Q;

%Retrieve number of variables up to desired order
nidx = length(find(polynomialDegree(oldVars,f) <= order));

%Compute the number of desired variables
ndes = factorial(order+Q+1-1)/(factorial(order)*factorial(Q+1-1));

%Derive the Carlemann linearization matrix
while orderVars <= order
    %Find the higher order variables
    allVars = [];
   
    for i = 1:1:numVars
        [~,vars] = coeffs(C(i),f,'All');
        vars = reshape(vars,[numel(vars) 1]);
        allVars = [allVars; vars];
    end
    
    vars = unique(allVars);
    
    clear allVars;
    
    %Determine the new variables
    newVars = setdiff(vars,oldVars);
    
    clear vars;
    
    %Stop if there are no new variables
    if isempty(newVars)
        break;
    end
    
    %Compute the number of new variables
    numNewVars = length(newVars);
    
    %Expand the current matrix
    %The ODE describing the new variables is obtained from:
    %partial(var)/partial(t) = Sum_i partial(var)/partial(f_i) df_i/dt
    for i = 1:1:numNewVars
        C = [C; sum(gradient(newVars(i),f).*omega)];
    end
    
    clear numNewVars;
    
    C = expand(C);
    
    %Update the existing variables
    oldVars = [oldVars; newVars];
    
    %Update the number of existing variables
    numVars =  length(oldVars);
    
    %Update the order of existing variables
    orderVars = max(polynomialDegree(oldVars,f));
    
    %Retrieve number of variables up to desired order
    nidx = length(find(polynomialDegree(oldVars,f) <= order));
end

clear nidx ndes orderVars;

%Determine the variables of the system
vars = oldVars;

clear oldVars newVars;

%Find the index of "truncated" variables up to the desired order
idx = find(polynomialDegree(vars,f) <= order);

%Determine the "truncated" variables
truncatedVars = vars(idx);

%and their corresponding differential equations
truncatedC =  C(idx);

clear idx C numVars;

%Determine the number of truncated variables
truncatedNumVars = length(truncatedVars);

%Drop the higher order terms from the expressions of the chosen variables
for i = 1:1:truncatedNumVars
    [ms,vs] = coeffs(truncatedC(i),f,'All');
    vs = reshape(vs,[numel(vs),1]);
    ms = reshape(ms,[numel(ms),1]);
    poly = polynomialDegree(vs,f);
    idx = find(poly > order);
    ms(idx) = 0.*ms(idx);
    truncatedC(i) = sum(expand(ms.*vs));
end

clear i ms vs poly idx;

%Rearrange the system of equations in order of variables
tDegree = polynomialDegree(truncatedVars,f);
[tDegree,tSort] = sort(tDegree,'descend');

truncatedVars = truncatedVars(tSort);
truncatedC = truncatedC(tSort);

clear tDegree tSort;

%Replace the variables in terms of the discrete density
%with a new set of variables to obtain the coefficient matrix
%of the Carlemann linearization

%The last variable is not replaced since it corresponds to unity
%needed for constant terms of the corresponding polynomials
%Replace the variables with new variable set

%Declare new set of variables phi
phi = sym('phi',[truncatedNumVars-1,1]);

%Perform the symbolic substitution
for i = 1:truncatedNumVars-1
    
    truncatedC = subs(truncatedC,truncatedVars(i),phi(i));
    
end

%Initialize the coefficients matrix
Acoeffs = sym(zeros(truncatedNumVars,truncatedNumVars));

for i = 1:1:truncatedNumVars-1

    [Acoeff,Avar] = coeffs(truncatedC(i),phi);
    for j = 1:1:length(Avar)
        idx = find([phi;1]==Avar(j));
        Acoeffs(i,idx) = Acoeff(j);
    end
    
end

clear i truncatedC Acoeff Avar j idx phi;
end