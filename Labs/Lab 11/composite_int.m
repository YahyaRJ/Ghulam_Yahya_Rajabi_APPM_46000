% testing composite trap and simpson's rule


a = 0;
b = pi;


f = @(x) sin(x);
Actual = 2;

ntest = 2:2:20;

errors= zeros(3,length(ntest));

for j = 1:length(ntest)
    n = ntest(j);
    h = (b-a)/n;
    errors(1,j) = h;
    qnode = a+[0:n]*h;
    
    
    % trapezoidal rule
    
    Itrap = h/2*(f(a)+f(b)+2*sum(f(qnode(2:end-1))));
    
    errors(2,j) = abs(Itrap-Actual);
    
    % Simpson
    
    Simp = h/3*(f(a)+f(b) + 2*sum(f(qnode(3:2:end-1)))+4*sum(f(qnode(2:2:end-1))));

    errors(3,j) = abs(Simp-Actual);
end

loglog(errors(1,:),errors(2,:),'bo-',errors(1,:),errors(3,:),'rx-')

x = log(errors(1,:));
y = log(errors(2,:));
m = (y(2:end)-y(1:end-1))./(x(2:end)-x(1:end-1));
m

x = log(errors(1,:));
y = log(errors(3,:));
m = (y(2:end)-y(1:end-1))./(x(2:end)-x(1:end-1));
m

keyboard