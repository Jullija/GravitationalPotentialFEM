function gravitationalPotential(n)
a = 0;
b = 3;

h = (b-a)/n;
funP1 =  1./(h.*h);
funP2 = -1./(h.*h);
P1 = funP1 * 2*h;
P2 = funP2 * h;


%główna macierz
M=sparse(n-1,n-1);

%macierz L(ej)
RHS = zeros(n-1, 1);


%uzupełnianie głównej macierzy
M(1, 1) = P1;
M(1, 2) = P2;

for k = 2 : n-2
   M( k , k-1) = P2;
   M ( k , k ) = P1 ;
   M ( k , k+1) = P2;
end

M (n-1,n-2) = P2;
M (n-1,n-1) = P1 ;


%macierz L(ej)
for k = 1 : n-1 
    RHS ( k , 1 ) = f(a+(k-1)*h)*h/6 + f(a+k*h)*2*h/3 + f(a+(k+1)*h)*h/6;
end

result = M\RHS;
points=[a : h : b];
values = zeros ( 1 , n+1 ) ;

for k = 2 : n
   values(k)=result(k-1);
end

for k = 1: n+1
    values(k) = values(k) + shift(a + (k-1)*h);
end

x=points ;
y=values ;

plot(x, y);

end

function y=f(x)
    G = 6.673e-11;
    const = 4 * pi * G;
   if (x > 1 && x <=2)
       y = const * 1;
   else
       y = 0;
   end
end

function y=shift(x)
    y=5 - x/3;
end
