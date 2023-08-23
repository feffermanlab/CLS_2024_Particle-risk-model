%Meant to be used with Community.m, so the the buliding and community 
%models are connected

%Community.m gives the porportion of S, I, R people in the buliding

%Sc, Ic, Rc are the susceptible, infected, and recovered people in the
%community-used to set up IC in the buliding
    % use S, I, R values captured from running Community.m for Sc, Ic, Rc

%day tells us which time step from Community.m to base the IC off of 
    %i.e. track particles in the buliding on the 15th day of the outbreak

%To run this code: first run Community.m, [t,S,I,R] = Community, then run
% Multi_Room_Model(S,I,R, 15) (days= 15 is just an example-selected desired
%time step)

%Graphical outputs: 4 plots: one for susceptible, infected, recovered, and
%one for particles. The lines on each plot represent the number of people
%or particles in each of the rooms (number of lines on each plot = number
%of rooms)

%author: Anna Sisk
%email: annahsisk@gmail.com
%Last update: 08/29/2022

function [t] = Multi_Room_Model(Sc, Ic, Rc, day)
close all
N=3; %number of room

m=4;   %number of equations for each room


MaxTime = 8;  %length of time we are tracking particles (hours of a work day?)
delta = 0.15; %proportionality constant-currently arbitrary

%distribution of S, I, R people in the buliding based on the distribution
%in the community on 'day'
Sb= delta*Sc(day);
Ib= delta*Ic(day);
Rb= delta*Rc(day);

%distributing the Sb, Ib, Rb into the rooms randomly
Sr= rand(N);
Ir = rand(N);
Rr = rand(N);
 
Sr= (Sr/sum(Sr))*Sb;
Ir = (Ir/sum(Ir))*Ib;
Rr = (Rr/sum(Rr))*Rb;

%initial conditions in each room
X0 = [];
for i = 1:N
       % [S0,      I0,     R0,   P0]
    ini =[Sr(i);  Ir(i);  Rr(i);  0];
    X0  =[X0; ini];
end 


%create random graph of rooms
netw = graph(createRandRegGraph(N, 2)); %random graph with N patches and fixed connectivity
plot(netw)
%calculate distance for each combination of rooms
dis = distances(netw);    %minus identity matrix (eye(N)), to set diagonal to 0

mov = rand(N,N); %transition rates for people-currently random
                   % mov(i,j) is the fraction of people that moves from i to j
mov(1:N+1:end) = 0;   % to set the diagonal to 0

p = 0.15;   % fraction of people that moves to other patches; this 
            % parameter determines the disease spreading from one location 
            % to others-currently arbitrary
for j=1:N
    mov(j,:)= mov(j,:)/sum(mov(j,:));  %normalized by rows. (people leaving can't be more than 100%)
end 

mov = p* mov; % p is the fraction of people that leave the room


for i=1:N  %disallows transitions between rooms that are not connected
    for j=1:N
        if dis(i,j)>1
            mov(i,j)=0;
        end
    end
end

tvec = linspace(0, MaxTime, MaxTime+1);

%solving the equations
options = odeset('RelTol',1e-8,'AbsTol',1e-10);
[t, pop]=ode45(@Room_eq, tvec, X0, options, mov, N, dis, m)

% separate the results in different variables
S=pop(:,1:m:m*N);
I=pop(:,2:m:m*N);
R=pop(:,3:m:m*N);
P=pop(:,4:m:m*N);

%plot everything
figure;
subplot(2,2,1)
plot(t,S(:,1:N));
xlabel 'Time';
ylabel 'Susceptible';

subplot(2,2,2) 
plot(t,I(:,1:N));
xlabel 'Time';
ylabel 'Infected';

subplot(2,2,3)
plot(t,R(:,1:N));
xlabel 'Time';
ylabel 'Recovered';

subplot(2,2,4) 
plot(t,P(:,1:N));
xlabel 'Time';
ylabel 'Particles';

if nargout<1
   clear t 
end

end

function dPop=Room_eq(t, pop, mov, N, dis, m)
%parameters-all currently arbitrary
%see bottom of page 3 in the write for more information on these parameters
s = 10000; %shedding rate
a = .002*s; %absorption rate
theta = .01*s; %transition rates for particles
d = .001*s; %natural death/decay rates for particles   
X=pop(1:m*N);
dPop=zeros(m*N,1);

% References:
% X(1+m*j) = susceptible
% X(2+m*j) = infected 
% X(3+m*j) = recovered 
% X(4+m*j) = particles


for j = 0:N-1
    %equations for S, I, R people in each room-equation (3) page 3 of the
    %write up 
    
    flux = flux_of_people(X(1:m:end), j+1, mov, N);
    dPop(1+j*m) = flux;
        
    flux = flux_of_people(X(2:m:end), j+1, mov, N);
    dPop(2+j*m) = flux;
        
    flux = flux_of_people(X(3:m:end), j+1, mov, N);
    dPop(3+j*m) = flux;
    
    %equation for particles in each room-equation (4) page 3 of the write up    
    dPop(4+j*m) = s*X(2+m*j) - a*(X(1+m*j) + X(2+m*j) + X(3+m*j))+(theta - d)*X(4+m*j);
end
end

function flux = flux_of_people(X, current_patch, mov, N)
% calculate the term for the in and out flux of people in a room

    j = current_patch;

    in_flux = transpose(X)*mov(:,j);
    
    out_flux = sum(mov(j,:)) * X(j) ;
    
    flux = in_flux - out_flux ;
end

%Two MATLAB functions below-more updated version of MATLAB have these bulit
%in
function A = createRandRegGraph(vertNum, deg)
% createRegularGraph - creates a simple d-regular undirected graph
% simple = without loops or double edges
% d-reglar = each vertex is adjecent to d edges
%
% input arguments :
%   vertNum - number of vertices
%   deg - the degree of each vertex
%
% output arguments :
%   A - A sparse matrix representation of the graph
%
% algorithm :
% "The pairing model" : create n*d 'half edges'.
% repeat as long as possible: pick a pair of half edges 
%   and if it's legal (doesn't creat a loop nor a double edge)
%   add it to the graph
% reference: http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.67.7957&rep=rep1&type=pdf

n = vertNum;
d = deg;
matIter = 10;

%check parameters
if mod(n*d,2)==1   
    disp('createRandRegGraph input err: n*d must be even!');
    A=[];
    return;
end

%a list of open half-edges
U = repmat(1:n,1,d);

%the graphs adajency matrix
A=sparse(n,n);

edgesTested=0; 
repetition=1;

%continue until a proper graph is formed
while ~isempty(U) && repetition < matIter
    
    edgesTested = edgesTested + 1;

    %print progress
    if mod(edgesTested, 5000)==0 
        fprintf('createRandRegGraph() progress: edges=%d/%d\n', edgesTested, n*d);    
    end

    %chose at random 2 half edges
    i1 = ceil(rand*length(U));
    i2 = ceil(rand*length(U));
    v1 = U(i1);
    v2 = U(i2);

    %check that there are no loops nor parallel edges
    if (v1 == v2) || (A(v1,v2) == 1)
        
        %restart process if needed
        if (edgesTested == n*d)           
            repetition=repetition+1;            
            edgesTested = 0;
            U = repmat(1:n,1,d);
            A = sparse(n,n);
        end
    else
        %add edge to graph
        A(v1, v2)=1;
        A(v2, v1)=1;
        
        %remove used half-edges
        v = sort([i1,i2]);
        U = [U(1:v(1)-1), U(v(1)+1:v(2)-1), U(v(2)+1:end)];
    end
end

%check that A is indeed simple regular graph
msg=isRegularGraph(A);
if ~isempty(msg)    
    disp(msg);
end
end%-------------------------------------------------

function msg=isRegularGraph(G)
%is G a simple d-regular graph the function returns []
%otherwise it returns a message describing the problem in G

msg=[];

%check symmetry
if (norm(G-G','fro')>0)
    msg=[msg,' is not symmetric, '];
end

%check parallel edged
if (max(G(:))>1)
    msg=[msg,sprintf(' has %d parallel edges, ',length(find(G(:)>1)) )];
end

%check that d is d-regular
d_vec=sum(G);
if min(d_vec)<d_vec(1) || max(d_vec)>d_vec(1)
    msg=[msg,' not d-regular, '];
end

%check that g doesn't contain any loops
if (norm(diag(G))>0)
    msg=[msg,sprintf(' has %d self loops, ',length(find(diag(G)>0)) )];
end

end   
