import math

# Parámetros físicos
epsilon = 0.0002  # rugosidad absoluta (m)
D = 0.05          # diámetro de la tubería (m)
Re = 1e5          # número de Reynolds

# Función de la ecuación de Colebrook
def colebrook(f, Re, epsilon, D):
    return 1 / math.sqrt(f) + 2 * math.log10((epsilon / (3.7 * D)) + (2.51 / (Re * math.sqrt(f))))

# Derivada de la ecuación de Colebrook
def d_colebrook(f, Re, epsilon, D):
    term = (epsilon / (3.7 * D)) + (2.51 / (Re * math.sqrt(f)))
    return (-0.5) * f ** (-1.5) - (2.51 / (math.log(10) * Re * term)) * (-0.5) * f ** (-1.5)

# Estimación inicial usando Haaland
def haaland_approx(Re, epsilon, D):
    return 1 / (-1.8 * math.log10((6.9 / Re) + ((epsilon / (3.7 * D)) ** 1.11)))**2

# Método de Newton-Raphson
def newton_raphson(f0, Re, epsilon, D, tol=1e-6, max_iter=100):
    f = f0
    for i in range(max_iter):
        f_new = f - colebrook(f, Re, epsilon, D) / d_colebrook(f, Re, epsilon, D)
        if abs(f - f_new) < tol:
            return f_new, i+1
        f = f_new
    raise RuntimeError("No converge")

# Aplicación del método
f0 = haaland_approx(Re, epsilon, D)
f_resultado, iteraciones = newton_raphson(f0, Re, epsilon, D)

print(f"Factor de fricción f: {f_resultado:.6f} (convergió en {iteraciones} iteraciones)")

