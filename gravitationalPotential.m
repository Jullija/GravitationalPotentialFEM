function gravitationalPotential(n)
a = 0;
b = 3;

h = (b-a)/n;
funP1 = @(h) 1./(h.*h);
funP2 = @(h) -1./(h.*h);
P1 = integral( funP1, 0, 2*h);
G = 6.673e-11;
const = 4 * pi * G;

%główna macierz
M=sparse(n-1,n-1);

%macierz L(ej)
RHS = zeros(n-1, 1);

%uzupełnianie macierzy
for r = 1 : n-1 
    for c = 1 : n-1  
        
        if r == c
            M(r, c) = P1;
        end

        if r == c + 1
            M(r, c) = integral(funP2, c*h, r*h);
        end

        if r + 1 == c 
            M(r, c) = integral(funP2, r*h, c*h);
        end
    end
end

%macierz L(ej)
for i = 2 : n 
    sum = 0;

    %część z całką pierwszą
    sum = sum + 1/3 * integral(@(h) 1./h, h*(i-1), h*i);
    sum = sum + 1/3 * integral(@(h) 1./h, h*i, h*(i+1));

    %część z całką drugą 
    %(tutaj całkujemy na przedziale 1-2, stąd min oraz max)
    %korzystamy ze wzoru na ei
    temp = i;
    if max(h*(i - 1), 1) < min(h*i, 2)
        sum = sum + const * integral(@(h) 1./h.*h - temp + 1, max(h*(i - 1), 1), min(h*i, 2));
    end

     if max(h*i, 1) < min(h*(i + 1), 2) 
        sum = sum + const * integral(@(h) temp - 1./h.*h  + 1, max(h*(i - 1), 1), min(h*i, 2));
    end

    RHS(i-1 , 1) = sum;
end

result = M\RHS;
points=[a : h : b];
values = zeros ( 1 , n+1 ) ;

for k = 2 : n
   values(k)=result(k-1);
end

x=points ;
y=values ;

plot(x, y);



end