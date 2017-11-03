
__precompile__()

module herramientas

export derivacion_simetrica, metodo_newton, integral_rectangulo, integral_trapecio, integral_simpson, interpolador_lagrange, metodo_euler, metodo_euler_implicito, runge_kutta_o4

"""Función para aproximar la derivada numéricamente por método simétrico 
``f`` -> la función a usar, 
 ``x0`` -> condicion inicial,
 ``h`` -> incremento para la aproximación
"""
function derivacion_simetrica(f,x0,h)
    df=(f(x0+h)-f(x0-h))/(2*h)
    return df #Regresa el valor de la derivada simétrica
end

function metodo_newton_ei(g,x0,h)
    dg(z)=derivacion_simetrica(g,z,h) #Definimos a dg(z) como la derivación simétrica, no usamos la variable x porque x será un valor bien definido
    x=x0 # Condición inicial
#Dado que la derivación simétrica tiene una excelente aproximación al valor real de la derivada en el punto deseado y también el método de Newton resulta ser eficiente en la aproximación se predetermina que la iteración se realice 20 veces para fines prácticos.
    for i in 1:20 #Ciclo for para 20 iteraciones
        x=x-g(x)/dg(x) #Método de Newton con una mejor aproximación al valor de la raíz
    end #Fin de la iteración
    return x #Regresa una buena aproximación al valor de la raíz de g
end

"""Método para encontrar las raíces de una función, las variables a declarar son:
 ``f`` -> la función a usar, 
 ``x0`` -> condicion inicial (la cual entre más cercana sea a una raíz de ``f`` garantiza una convergencia más rápida), 
 ``it`` la cantidad de veces a iterar el proceso (número natural), 
 Utiliza derivación simétrica para aproximar el valor de la derivada de f con lo que nos evitamos usar paqueterias ajenas.
 Redacción: `metodo_newton(f,x0,it)`"""
function metodo_newton(f,x0,it)
    df(z)=derivacion_simetrica(f,z,0.001) #Como converge rápido la derivación para evitar pedir otro dato se hace h=0.001
    x=x0
    for i in 1:it
        x=x-f(x)/df(x) #Método de Newton con una mejor aproximación al valor de la raíz
    end #Fin de la iteración
    return x #Regresa una buena aproximación al valor de la raíz de g
end

""" Aproximación a la integral de f(x) por Método de Riemann (Rectángulos) usando particiones homogéneas. Se consideran 4
 parámetros para dar una mayor libertad al usuario al momento de generar el valor de la integral:
 ``f`` -> Se da la función a integrar, 
 ``a`` -> Punto inicial del intervalo de integración, 
 ``b`` -> Punto final del intervalo de integración, por convención a < b, 
 ``n`` -> Cantidad de subintervalos inducidos deseados del intervalo [a,b] n natural, 
 Redacción: `integral_rectangulo(f,a,b,n)`"""
function integral_rectangulo(f,a,b,n)
    longitud=abs((b-a)/n) #La longitud de cada subintervalo dada la partición
    g(x)=f(x) #La función f definida con la variable x
    st=0 #Acumulador del valor de la suma del área de los rectángulos en cada subintervalo
    for i in 1:n #Ciclo for, se detiene al completar los n subintervalos inducidos
        punto=((2*a)+longitud)/2 #Definimos el punto a evaluar en la función, este valor (para ejemplificar tomamos el punto medio de cada subintervalo)
        si=g(punto)*longitud #Hacemos el área del rectangulo inducido en el i-ésimo subintervalo
        st=st+si #Acumulamos el valor en st
        a=a+longitud #Para irnos al siguiente intervalo redefinimos el valor de a
    end
    return st #Muestra el valor final de la suma de Riemann finita
end

""" Aproximación a la integral de f(x) por Método de Trapecios usando particiones homogéneas. Se consideran 4
 parámetros para dar una mayor libertad al usuario al momento de generar el valor de la integral:
 ``f`` -> Se da la función a integrar, 
 ``a`` -> Punto inicial del intervalo de integración, 
 ``b`` -> Punto final del intervalo de integración, por convención a < b, 
 ``n`` -> Cantidad de subintervalos inducidos deseados del intervalo [a,b] n natural, 
 Redacción: `integral_trapecio(f,a,b,n)`"""
function integral_trapecio(f,a,b,n)
    longitud=abs((b-a)/n) #La longitud de cada subintervalo dada la partición
    g(x)=f(x) #La función f definida con la variable x
    factor=(f(a)+f(b))/2 #Factor que se suma a la suma parcial
    sp=0 #Acumulador representando la suma parcial
    for i in 1:(n-1) #Ciclo for para realizar la suma parcial que tiene como límite n-1
        argumento=a+(i*longitud) #Argumento a evaluar en f(x)
        sp=sp+g(argumento) #Suma parcial hasta la i-ésima iteración
    end
    st=longitud*(factor+sp) #Suma total con todos los valores definidos
    return st #Muestra el valor final de la iteración para trapecios
end

""" Aproximación a la integral de f(x) por Método de Simpson usando particiones homogéneas. Se consideran 4
 parámetros para dar una mayor libertad al usuario al momento de generar el valor de la integral:
 ``f`` -> Se da la función a integrar, 
 ``a`` -> Punto inicial del intervalo de integración, 
 ``b`` -> Punto final del intervalo de integración, por convención a < b, 
 ``n`` -> Cantidad de subintervalos inducidos deseados del intervalo [a,b] n natural, 
 Redacción: `integral_simpson(f,a,b,n)`"""
function integral_simpson(f,a,b,n)
    longitud=abs((b-a)/n) #La longitud de cada subintervalo dada la partición
    g(x)=f(x) #La función f definida con la variable x
    st=0 #Acumulador de los valores de las integrales en cada subintervalo
    c=a+longitud #Variable que nos define el supremo (punto final) de cada subintevralo
    for i in 1:n #Ciclo for que corre de 1 a n para evaluar los n intervalos
        aint=((c-a)/6)*(f(a)+4f((a+c)/2)+f(c)) #Aplicando Método de Simpson en cada subintervalo
        st=st+aint #Representación de la suma parcial a la iésima iteración
        a=c #Redefinimos a para analizar el proximo intervalo
        c=c+longitud #Redefinimos c para analizar el proximo intervalo
    end #Fin de for
    return st #Arrojar el resultado de st
end #Fin

"""Dados dos vectores ``x``,``y`` en ``\\mathbb{R}^{n}`` la función `interpolador_lagrange` da el valor del polinomio
 de Lagrange ``L`` que pasa por ``(x_{k},y_{k})\\:\\forall k=1,\\dots,n`` evaluado en ``x_{0}``, i.e. nos da ``L(x_{0})``
 x, y son listas que deben tener la misma cantidad de elementos, y es la lista de x evaluada en f.
 Redacción: `interpolador_lagrange(x,y,x0)`"""
function interpolador_lagrange(x,y,x0)
    lx=length(x)
    l=[]
    for j in 1:lx
        z=1
        for m in 1:lx
            if m != j
                z*=(x0-x[m])/(x[j]-x[m])
            end
        end
        push!(l,z)
    end
    return vecdot(y,l)
end

"""Función que implementa el método de Euler, recibe:
 ``f`` -> función de la ecuación, 
 ``x0`` -> Condición inicial, 
 -OBS: Al ser independiente de la dimensión se pueden evaluar distintas condiciones iniciales a la vez introduciendo en x0
 [condición 1, condición 2, ..., condición n]
 ``listat`` -> Lista la cual aglomera la condición de tiempo inicial (primer valor) y el paso h que pide el método
(incremento la cantidad de veces necesarias), 
 Redacción: metodo_euler(f,listat,x0)"""
function metodo_euler(f,listat,x0)
    longitud=length(listat)
    x=x0
    h=listat[2]-listat[1]
    listax=[] #Solo declaro a listax como un arreglo vacío
    push!(listax,x) #Aplico push! para agregar un valor a la listax el valor x0 que es la condición inicial
    for i in 2:longitud #El Ciclo for ahora corre de i=2 hasta que haya usado todos los valores de la listat
        t=listat[i] #Para usar el elemento t correspondiente a la iteración
        x=x+f(x,t)*h #Aplicando el método de Euler
        push!(listax,x) #Agrego el valor encontrado a la listax el valor de x que encontré
    end
    return listax #Regresó el arreglo
end

"""Implementación del método de Euler implícito, recibe:
 ``f`` -> función de la relación de recurrencia, 
 ``x0`` -> Condición inicial, 
 -OBS: Al ser independiente de la dimensión se pueden evaluar distintas condiciones iniciales a la vez introduciendo en x0
 [condición 1, condición 2, ..., condición n], 
 ``listat`` -> Lista la cual aglomera la condición de tiempo inicial (primer valor) y el paso h que pide el método
(incremento la cantidad de veces necesarias), 
 Redacción: metodo_euler_implicito(f,listat,x0)"""
function metodo_euler_implicito(f,listat,x0)
    longitud=length(listat) #Variable para obtener el número de elementos de la lista, sirve para controlar la iteración
    listax=[] #Declaramos la lista 'listax'
    push!(listax,x0)  #En listax el primer valor es la condición inicial x0
    h=abs(listat[2]-listat[1]) #Para obtener el incremento h en nuestra lista h sacamos la diferencia entre los primeros dos elementos (se crea un linspace o análogos con particiones homogéneas, por eso es útil esto)
    xk=x0
    for i in 1:longitud-1 #Ciclo for, se hace tantas veces como la longitud de la lista menos un valor
        xk=listax[i] #Definimos a xk como el k-ésimo elemento de listax
        t=listat[i+1] #De la misma manera, proponemos a t como el (k+1)-ésimo elemento de listat
        g(z)=z-xk-h*f(z,t) #g es la función con raíz x_(k+1) a la que le aplicaremos el Método de Newton
        push!(listax,metodo_newton_ei(g,xk,h))
        xk=metodo_newton_ei(g,xk,h)
    end #Fin del Ciclo for
    return listax #El resultado es la listax con todos los valores salvados de x_k
end

"""Función que implementa el método de Runge-Kutta de orden 4, recibe:
 ``f`` -> función de la ecuación, 
 ``x0`` -> Condición inicial, 
 -OBS: Al ser independiente de la  dimensión se pueden evaluar distintas condiciones iniciales a la vez introduciendo en x0
 [condición 1, condición 2, ..., condición n], 
 ``listat`` -> Lista la cual aglomera la condición de tiempo inicial (primer valor) y el paso h que pide el método
(incremento la cantidad de veces necesarias), 
 Redacción: Runge_Kutta_O4(f,listat,x0)"""
function runge_kutta_o4(f,listat,x0)
    longitud=length(listat)
    x=x0
    l=(listat[2]-listat[1])/2
    listax=[] #Solo declaro a listax como un arreglo vacío
    push!(listax,x) #Aplico push! para agregar un valor a la listax el valor x0 que es la condición inicial
    for i in 2:length(listat) #El Ciclo for ahora corre de i=2 hasta que haya usado todos los valores de la listat
        t=listat[i] #Para usar el elemento t correspondiente a la iteración
        k1=f(x,t) #Encontramos k1 ahora usando x,t en vez de explicitamente indicar los elementos de una lista
        k2=f(x+l*k1,t+l)
        k3=f(x+l*k2,t+l)
        k4=f(x+2l*k3,t+2l)
        x=x+(l/3)*(k1+2k2+2k3+k4) #Aplicamos el método de Runge-Kutta
        push!(listax,x) #Salvo el valor de x en listax
    end
    return listax #Regreso el arreglo
end

end
