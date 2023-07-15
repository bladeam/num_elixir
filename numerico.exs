defmodule Analisis do
# método de newton para ecuaciones
def f(x) do
    x*x-4
end
def df(x) do
    2*x
end
def newton(x,n) do
    if n == 0 do
      x
    else
      newton(x-f(x)/df(x),n-1)
    end
end
# método de Newton para ecuaciones
def newton2(x,0), do: x
def newton2(x,n), do: [x] ++ newton2(x-f(x)/df(x),n-1)
# Método del punto fijo: este metodo depende de la funcion
def g(x), do: x+x*x
def punto_fijo(x,0), do: x
def punto_fijo(x,n), do: [x] ++ punto_fijo(g(x),n-1)
# método de Jacobi
def x(y,z), do: (2*y-3*z+5)/8
def y(x,z), do: (3*x+5*z-2)/9
def z(x,y), do: (2*x-3*y-3)/12
def jacobi(a,b,c,0), do: {a,b,c}
def jacobi(a,b,c,n), do: [{a,b,c}] ++ jacobi(x(b,c),y(a,c),z(a,b),n-1)


#----------------------------------------------------------------------
# Minimos cuadrados
#----------------------------------------------------------------------
# Ejecutar con : Analisis.min_cua([1,2,3],[4,5,6])
def min_cua(l1,l2) do
  sx = Enum.sum(l1)
  sy = Enum.sum(l2)
  n = length(l1) #Enum.map(l1,fn _ -> 1 end)|>Enum.sum
  xy = Enum.sum(for {x,y}<-Enum.zip(l1,l2), do: x*y)
  #x2 = Enum.map(l1,fn x -> x*x end)|>Enum.sum
  x2 = (for x <- l1, do: x*x) |> Enum.sum
  m = (sx*sy-n*xy)/(sx*sx-n*x2)
  b = (sy-m*sx)/n
  [m,b]
end
#----------------------------------------------------------------------
# Runge-Kutta para ED de primer orden
#----------------------------------------------------------------------
# Resolver: y'=e^t    y(0)=1    Hallar y(1)
# Ejecutar con:  Analisis.runge_ed1(5,0,1,1)
def g1(t, x), do: (:math.exp(t))


def runge_ed1(n,t,x,a), do: runge_ed1_aux(n,t,x,(a-t)/n)
def runge_ed1_aux(0,_,x,_), do: x
def runge_ed1_aux(n,t,x,h) do
   k11 = h*g1(t, x)


   k21 = h*g1((t+h/2), (x+k11/2))

   k31 = h*g1((t+h/2), (x+k21/2))

   k41 = h*g1((t+h), (x+k31))

   runge_ed1_aux(n-1,t+h, (x+(1/6)*(k11+2*k21+2*k31+k41)),h)
end
#----------------------------------------------------------------------
# Runge-Kutta para ED de segundo orden
#----------------------------------------------------------------------
# Resolver: y''+2y'+y=e^t    y(0)=1   y'(0)=0   Hallar y(1)
# Ejecutar con: Analisis.runge_ed2(5,0,0,1,1)
def h1(t, x, y), do: (:math.exp(t))-y-2*x
def h2(t, x, y), do: x

def runge_ed2(n,t,x,y,a), do: runge_ed2_aux(n,t,x,y,(a-t)/n)
def runge_ed2_aux(0,_,x,y,_), do: [{x,y}]
def runge_ed2_aux(n,t,x,y,h) do
   k11 = h*h1(t, x, y)
   k12 = h*h2(t, x, y)

   k21 = h*h1((t+h/2), (x+k11/2), (y+k12/2))
   k22 = h*h2((t+h/2), (x+k11/2), (y+k12/2))

   k31 = h*h1((t+h/2), (x+k21/2), (y+k22/2))
   k32 = h*h2((t+h/2), (x+k21/2), (y+k22/2))

   k41 = h*h1((t+h), (x+k31), (y+k32))
   k42 = h*h2((t+h), (x+k31), (y+k32))

   [{x,y}] ++ runge_ed2_aux(n-1,t+h, (x+(1/6)*(k11+2*k21+2*k31+k41)), (y+(1/6)*(k12+2*k22+2*k32+k42)), h)
end
#----------------------------------------------------------------------
# Runge-Kutta para ED de tercer orden
#----------------------------------------------------------------------
# Resolver: y'''+ 4y''+ y'-6y=e^t    y(0)=1   y'(0)=0   y''(0)=2    Hallar y(1)
# ejecutar con: Analisis.runge_ed3(5,0,0,1,2,1)
def p1(t, x, y, z), do: z
def p2(t, x, y, z), do: x
def p3(t, x, y, z), do: (:math.exp(t))+6*y-x-4*z
def runge_ed3(n,t,x,y,z,a), do: runge_ed3_aux(n,t,x,y,z,(a-t)/n)
def runge_ed3_aux(0,_,x,y,z,_), do: [{x,y,z}]
def runge_ed3_aux(n,t,x,y,z,h) do
   k11 = h*p1(t, x, y, z)
   k12 = h*p2(t, x, y, z)
   k13 = h*p3(t, x, y, z)
   k21 = h*p1((t+h/2), (x+k11/2), (y+k12/2), (z+k13/2))
   k22 = h*p2((t+h/2), (x+k11/2), (y+k12/2), (z+k13/2))
   k23 = h*p3((t+h/2), (x+k11/2), (y+k12/2), (z+k13/2))
   k31 = h*p1((t+h/2), (x+k21/2), (y+k22/2), (z+k23/2))
   k32 = h*p2((t+h/2), (x+k21/2), (y+k22/2), (z+k23/2))
   k33 = h*p3((t+h/2), (x+k21/2), (y+k22/2), (z+k23/2))
   k41 = h*p1((t+h), (x+k31), (y+k32), (z+k33))
   k42 = h*p2((t+h), (x+k31), (y+k32), (z+k33))
   k43 = h*p3((t+h), (x+k31), (y+k32), (z+k33))
   [{x,y,z}] ++ runge_ed3_aux(n-1,t+h, (x+(1/6)*(k11+2*k21+2*k31+k41)), (y+(1/6)*(k12+2*k22+2*k32+k42)), (z+(1/6)*(k13+2*k23+2*k33+k43)), h)
end
#----------------------------------------------------------------------
# Runge-Kutta para sistemas de ED 2x2
#----------------------------------------------------------------------
# Resolver: x'= 3*x+2*y+2*t
#           y'= x+4*y-7
# Con el valor inicial: x(0)=0    y(0)=1

# ejecutar con: Analisis.kutta_2x2(5, 0, 0, 1, 1)

def f1(t, x, y), do: 3*x+2*y+2*t
def f2(t, x, y), do: x+4*y-7

def kutta_2x2(n,t,x,y,a), do: kutta_2x2_a(n,t,x,y,(a-t)/n)
def kutta_2x2_a(0,_,x,y,_), do: [{x,y}]
def kutta_2x2_a(n,t,x,y,h) do
   k11 = h*f1(t, x, y)
   k12 = h*f2(t, x, y)

   k21 = h*f1((t+h/2), (x+k11/2), (y+k12/2))
   k22 = h*f2((t+h/2), (x+k11/2), (y+k12/2))

   k31 = h*f1((t+h/2), (x+k21/2), (y+k22/2))
   k32 = h*f2((t+h/2), (x+k21/2), (y+k22/2))

   k41 = h*f1((t+h), (x+k31), (y+k32))
   k42 = h*f2((t+h), (x+k31), (y+k32))

   [{x,y}] ++ kutta_2x2_a(n-1,t+h, (x+(1/6)*(k11+2*k21+2*k31+k41)), (y+(1/6)*(k12+2*k22+2*k32+k42)), h)
end
#----------------------------------------------------------------------
# Runge-Kutta para sistemas de ED 3x3
#----------------------------------------------------------------------
# Ejecutar con : Analisis.kutta_3x3(10,0,1,0,0,1)
def q1(t, x, y, z), do: z+y-x
def q2(t, x, y, z), do: z+x-y
def q3(t, x, y, z), do: x+y+z
def kutta_3x3(n,t,x,y,z,a), do: kutta_3x3_a(n,t,x,y,z,(a-t)/n)
def kutta_3x3_a(0,_,x,y,z,_), do: [{x,y,z}]
def kutta_3x3_a(n,t,x,y,z,h) do
   k11 = h*q1(t, x, y, z)
   k12 = h*q2(t, x, y, z)
   k13 = h*q3(t, x, y, z)
   k21 = h*q1((t+h/2), (x+k11/2), (y+k12/2), (z+k13/2))
   k22 = h*q2((t+h/2), (x+k11/2), (y+k12/2), (z+k13/2))
   k23 = h*q3((t+h/2), (x+k11/2), (y+k12/2), (z+k13/2))
   k31 = h*q1((t+h/2), (x+k21/2), (y+k22/2), (z+k23/2))
   k32 = h*q2((t+h/2), (x+k21/2), (y+k22/2), (z+k23/2))
   k33 = h*q3((t+h/2), (x+k21/2), (y+k22/2), (z+k23/2))
   k41 = h*q1((t+h), (x+k31), (y+k32), (z+k33))
   k42 = h*q2((t+h), (x+k31), (y+k32), (z+k33))
   k43 = h*q3((t+h), (x+k31), (y+k32), (z+k33))
   [{x,y,z}] ++ kutta_3x3_a(n-1,t+h, (x+(1/6)*(k11+2*k21+2*k31+k41)), (y+(1/6)*(k12+2*k22+2*k32+k42)), (z+(1/6)*(k13+2*k23+2*k33+k43)), h)
end
# fin del módulo
end
