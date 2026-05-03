"""Para ejecutar el código necesitamos instalar las librerías:

        pip install streamlit sympy

Y ejecutar el comando:

        py -m streamlit run app.py
        
o si no funciona escribir python completo (a veces presenta fallas al ejecutar usando uno u otro comando):
        
        python -m streamlit run app.py
        
Podemos seleccionar el localhost o dirección ip que se presenta estando en la misma red. Se abrirá en una ventana de navegador.

"""

import streamlit as st
import sympy as sp
import io
import contextlib

# ==========================================
# TUS MÉTODOS NUMÉRICOS (Intactos)
# ==========================================
def biseccion(f, xl, xu, tol):
    print("\n" + "="*45)
    print(" 1. MÉ MÉTODO DE BISECCIÓN ")
    print("="*45)
    
    if f(xl) * f(xu) >= 0:
        print("Error: f(xl) y f(xu) no tienen signos opuestos.")
        return None, None

    xr_old = 0
    error = float('inf') 
    iteracion = 1

    print(f"{'Iter':<5} | {'xl':<10} | {'xu':<10} | {'xr':<10} | {'Error Abs':<10}")
    print("-" * 55)

    while True:
        xr = (xl + xu) / 2.0
        if iteracion > 1:
            error = abs(xr - xr_old)

        if iteracion == 1:
            print(f"{iteracion:<5} | {xl:<10.5f} | {xu:<10.5f} | {xr:<10.5f} | {'-':<10}")
        else:
            print(f"{iteracion:<5} | {xl:<10.5f} | {xu:<10.5f} | {xr:<10.5f} | {error:<10.5f}")

        if error < tol and iteracion > 1:
            break

        if f(xl) * f(xr) < 0:
            xu = xr
        elif f(xl) * f(xr) > 0:
            xl = xr
        else:
            break 

        xr_old = xr
        iteracion += 1

    print(f"\nRaíz por Bisección: {xr:.7f}")
    return xr, iteracion

def falsa_posicion(f, xl, xu, tol):
    print("\n" + "="*45)
    print(" 2. MÉTODO DE LA FALSA POSICIÓN ")
    print("="*45)
    
    if f(xl) * f(xu) >= 0:
        print("Error: f(xl) y f(xu) no tienen signos opuestos.")
        return None, None

    xr_old = 0
    error = float('inf') 
    iteracion = 1

    print(f"{'Iter':<5} | {'xl':<10} | {'xu':<10} | {'xr':<10} | {'Error Abs':<10}")
    print("-" * 55)

    while True:
        xr = xu - (f(xu) * (xl - xu)) / (f(xl) - f(xu))
        if iteracion > 1:
            error = abs(xr - xr_old)

        if iteracion == 1:
            print(f"{iteracion:<5} | {xl:<10.5f} | {xu:<10.5f} | {xr:<10.5f} | {'-':<10}")
        else:
            print(f"{iteracion:<5} | {xl:<10.5f} | {xu:<10.5f} | {xr:<10.5f} | {error:<10.5f}")

        if error < tol and iteracion > 1:
            break

        if f(xl) * f(xr) < 0:
            xu = xr
        elif f(xl) * f(xr) > 0:
            xl = xr
        else:
            break

        xr_old = xr
        iteracion += 1

    print(f"\nRaíz por Falsa Posición: {xr:.7f}")
    return xr, iteracion

def newton_raphson(f, df, x0, tol):
    print("\n" + "="*45)
    print(" 3. MÉTODO DE NEWTON-RAPHSON ")
    print("="*45)

    iteracion = 1
    error = float('inf')

    print(f"{'Iter':<5} | {'xi':<10} | {'Error Abs':<10}")
    print("-" * 40)

    while True:
        df_val = df(x0)
        if df_val == 0:
            print(f"[!] Derivada cero en x={x0}. Desplazando +0.001...")
            x0 += 0.001  
            df_val = df(x0) 

        x1 = x0 - (f(x0) / df_val)
        error = abs(x1 - x0)

        if iteracion == 1:
            print(f"{iteracion:<5} | {x1:<10.5f} | {'-':<10}")
        else:
            print(f"{iteracion:<5} | {x1:<10.5f} | {error:<10.5f}")

        if error < tol:
            break
            
        if iteracion > 100:
             print("\nEl método no converge después de 100 iteraciones.")
             return None, None

        x0 = x1
        iteracion += 1

    print(f"\nRaíz por Newton-Raphson: {x1:.7f}")
    return x1, iteracion

def secante(f, x_ant, x_act, tol):
    print("\n" + "="*45)
    print(" 4. MÉTODO DE LA SECANTE ")
    print("="*45)

    error = float('inf') 
    iteracion = 1

    print(f"{'Iter':<5} | {'x_i+1':<15} | {'Error Abs':<15}")
    print("-" * 40)

    while error > tol:
        denominador = f(x_ant) - f(x_act)
        if denominador == 0:
            print("División por cero. El método falla.")
            return None, None
            
        numerador = f(x_act) * (x_ant - x_act)
        x_sig = x_act - (numerador / denominador)
        error = abs(x_sig - x_act)
        
        print(f"{iteracion:<5} | {x_sig:<15.5f} | {error:<15.5f}")
        
        x_ant = x_act
        x_act = x_sig
        iteracion += 1
        
        if iteracion > 100:
             print("\nEl método no converge después de 100 iteraciones.")
             return None, None

    print(f"\nRaíz por la Secante: {x_act:.7f}")
    return x_act, iteracion

def punto_fijo(g, xi, tol):
    print("\n" + "="*45)
    print(" 5. MÉTODO DEL PUNTO FIJO ")
    print("="*45)

    error = float('inf')
    iteracion = 1

    print(f"{'Iter':<5} | {'x_i+1':<15} | {'Error Abs':<15}")
    print("-" * 40)

    while error > tol:
        try:
            xi_next = g(xi)
        except Exception as e:
            print(f"\nError al evaluar: {e}")
            return None, None
            
        error = abs(xi_next - xi)
        print(f"{iteracion:<5} | {xi_next:<15.5f} | {error:<15.5f}")
        
        xi = xi_next
        iteracion += 1
        
        if iteracion > 100:
            print("\nEl método diverge con este despeje o punto inicial.")
            return None, None

    print(f"\nRaíz por Punto Fijo: {xi:.7f}")
    return xi, iteracion

# ==========================================
# INTERFAZ WEB CON STREAMLIT
# ==========================================

# Configuración de la página
st.set_page_config(page_title="Calculadora de Raíces", layout="wide")

st.title("Calculadora de Raíces - Métodos Numéricos")
st.markdown("Ingresa tus funciones y parámetros en el menú lateral para ver los resultados.")

# Barra lateral para las entradas (Sidebar)
with st.sidebar:
    st.header("Configuración de Datos")
    
    f_str = st.text_input("1. Función f(x):", value="x**3 - x - 1")
    g_str = st.text_input("2. Función g(x) [Punto Fijo]:", value="(x + 1)**(1/3)")
    
    st.markdown("---")
    st.subheader("Parámetros")
    xl = st.number_input("Límite inferior (xl) / x_ant:", value=1.0, format="%.5f")
    xu = st.number_input("Límite superior (xu) / x_act:", value=2.0, format="%.5f")
    x0 = st.number_input("Punto inicial (x0):", value=1.0, format="%.5f")
    tol = st.number_input("Tolerancia (Ej. 0.001):", value=0.001, format="%.5f")
    
    ejecutar = st.button("Calcular Raíces", use_container_width=True)

# Lógica de ejecución
if ejecutar:
    x = sp.Symbol('x')
    
    # Redirigir la salida de los print() a una variable de texto
    captura_consola = io.StringIO()
    
    with contextlib.redirect_stdout(captura_consola):
        try:
            # Parsear funciones
            f_expr = sp.sympify(f_str)
            f = sp.lambdify(x, f_expr, 'math')
            
            g_expr = sp.sympify(g_str)
            g = sp.lambdify(x, g_expr, 'math')
            
            df_expr = sp.diff(f_expr, x)
            df = sp.lambdify(x, df_expr, 'math')

            print("="*55)
            print(" INICIANDO CÁLCULOS ")
            print("="*55)
            print(f"[+] Función f(x): {f_expr}")
            print(f"[+] Función g(x): {g_expr}")
            print(f"[+] Derivada automática: {df_expr}")
            print(f"[+] Tolerancia: {tol}\n")

            # Ejecutar métodos
            res_bis = biseccion(f, xl, xu, tol)
            res_fp = falsa_posicion(f, xl, xu, tol)
            res_nr = newton_raphson(f, df, x0, tol)
            res_sec = secante(f, xl, xu, tol) 
            res_pf = punto_fijo(g, x0, tol)

            # Resumen final
            print("\n\n" + "="*55)
            print(f"{'RESUMEN FINAL DE RESULTADOS':^55}")
            print("="*55)
            print(f"{'Método':<20} | {'Raíz Aproximada':<15} | {'Iteraciones'}")
            print("-" * 55)
            
            metodos = [
                ("Bisección", res_bis),
                ("Falsa Posición", res_fp),
                ("Newton-Raphson", res_nr),
                ("Secante", res_sec),
                ("Punto Fijo", res_pf)
            ]
            
            for nombre, resultado in metodos:
                raiz, iters = resultado
                if raiz is not None:
                    print(f"{nombre:<20} | {raiz:<15.7f} | {iters}")
                else:
                    print(f"{nombre:<20} | {'Falló/Diverge':<15} | {'-'}")
                    
            print("="*55)
            
        except Exception as e:
            print(f"\n[!] ERROR: Verifica la sintaxis de tus funciones o los valores ingresados.")
            print(f"Detalle técnico: {e}")

    # Mostrar la salida capturada en la página web usando un bloque de código
    st.subheader("Resultados de Consola")
    st.code(captura_consola.getvalue(), language="text")