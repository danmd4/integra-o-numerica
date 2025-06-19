# metodo-de-integra-o-numerico

import numpy as np
import matplotlib.pyplot as plt

k = 8.9875517923e9  # N·m²/C²

def densidade_linear_uniforme(x):
    return 1e-6  # C/m (uniforme)

def densidade_linear_gaussiana(x, mu=0.0, sigma=0.2):
    return 1e-6 * np.exp(-((x - mu)*2) / (2 * sigma*2))

def calcular_potencial_linha(x0, x1, n, r_obs, densidade_func):
    """
    Calcula o potencial elétrico no ponto r_obs devido a uma distribuição linear de carga
    de x0 a x1 com n divisões, usando integração numérica (regra dos trapézios).
    """
    x = np.linspace(x0, x1, n)
    dx = (x1 - x0) / (n - 1)

    V_total = 0
    for i in range(n - 1):
        xi = x[i]
        xi1 = x[i + 1]
        
        # Ponto médio entre xi e xi+1
        xm = 0.5 * (xi + xi1)
        
        # Distância do elemento ao ponto de observação
        R = np.linalg.norm([xm - r_obs[0], -r_obs[1]])  # linha ao longo do eixo x
        if R == 0:
            continue  # evitar singularidade
        
        # Carga aproximada usando densidade média
        lambda_i = densidade_func(xi)
        lambda_i1 = densidade_func(xi1)
        dq = 0.5 * (lambda_i + lambda_i1) * dx

        V_total += k * dq / R

    return V_total

def exemplo_multiplas_distribuicoes():
    # Configuração de distribuições de carga
    distribuicoes = [
        {"x0": -1, "x1": 1, "func": densidade_linear_uniforme},
        {"x0": 2, "x1": 4, "func": lambda x: densidade_linear_gaussiana(x, mu=3, sigma=0.3)}
    ]

    r_obs = [0, 1]  # Ponto de observação (x, y)
    n = 1000  # número de divisões

    V_total = 0
    for dist in distribuicoes:
        V_total += calcular_potencial_linha(dist["x0"], dist["x1"], n, r_obs, dist["func"])

    print(f"Potencial elétrico no ponto {r_obs}: {V_total:.3e} V")

if __name__ == "__main__":
    exemplo_multiplas_distribuicoes()
