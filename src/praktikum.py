import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


# функция решения ДУ
def solve_DU(deltaM1z2Q_, deltaM3z2Q_, deltaf1z1Q_, deltaf3z1Q_):
    """
    :param deltaM1z2Q_:
    :param deltaM3z2Q_:
    :param deltaf1z1Q_:
    :param deltaf3z1Q_:
    :return: two solutions
    """
    print("\nsolving DU for data:", deltaM1z2Q_, deltaM3z2Q_, deltaf1z1Q_, deltaf3z1Q_)
    # константы
    m = 225
    I22 = 12.5
    T = 1
    L = 1
    f1z1 = 20
    f2z1 = 0
    f3z1 = -20
    fz3max = 4
    # обезразмеривание
    u1z1 = f1z1 * T * T / (m * L)
    u2z1 = f2z1 * T * T / (m * L)
    u3z1 = f3z1 * T * T / (m * L)
    v1z1 = deltaf1z1Q_ * T * T / (m * L)
    v1z2 = deltaM1z2Q_ * T * T / I22
    v3z1 = deltaf3z1Q_ * T * T / (m * L)
    v3z2 = deltaM3z2Q_ * T * T / I22
    uz3max = fz3max * T * T / (m * L)

    # управляемая подсистема
    # 1) разгон (н.у = 0)
    t = np.linspace(0, 20, 100)
    y0 = [0, 0, 0, 0, 0, 0]

    def f1(y, t):
        q1, q2, q3, q4, q5, q6 = y
        return [q2, 0, q4, -u1z1 * q5 - uz3max, q6, 0]

    solution1q = odeint(f1, y0, t)
    print("shape of solution:", len(solution1q), len(solution1q[0]))
    # 2) дрейф (н.у = конечная точка для разгона)
    t = np.linspace(20, 40, 100)
    y0 = solution1q[-1]

    def f1(y, t):
        q1, q2, q3, q4, q5, q6 = y
        return [q2, 0, q4, -u2z1 * q5 - uz3max, q6, 0]

    solution2q = odeint(f1, y0, t)
    # 3) торможение (н.у = конечная точка для дрейфа)
    t = np.linspace(40, 60, 100)
    y0 = solution2q[-1]

    def f1(y, t):
        q1, q2, q3, q4, q5, q6 = y
        return [q2, 0, q4, -u3z1 * q5 - uz3max, q6, 0]

    solution3q = odeint(f1, y0, t)
    # проверка того, что интегрирование уравнения по частям и полностью совпадают
    solution4q = odeint(f1, [0, 0, 0, 0, 0, 0], np.linspace(0, 60, 300))
    print("solution1q[20] =", solution1q[-1])
    print("solution2q[40] =", solution2q[-1])
    print("solution3q[60] =", solution3q[-1])
    print("solution4q[60] =", solution4q[-1])

    # возмущаемая подсистема
    # 1) разгон (н.у = 15 по х3)
    t = np.linspace(0, 20, 100)
    y0 = [0, 0, 15, 0, 0, 0]

    def f2(y, t):
        s1, s2, s3, s4, s5, s6 = y
        return [s2, v1z1, s4, -u1z1 * s5, s6, v1z2]

    solution1s = odeint(f2, y0, t)
    # 2) дрейф (н.у = конечная точка для разгона)
    t = np.linspace(20, 40, 100)
    y0 = solution1s[-1]

    def f2(y, t):
        s1, s2, s3, s4, s5, s6 = y
        return [s2, 0, s4, -u2z1 * s5, s6, 0]

    solution2s = odeint(f2, y0, t)
    # 3) торможение (н.у = конечная точка для дрейфа)
    t = np.linspace(40, 60, 100)
    y0 = solution2s[-1]

    def f2(y, t):
        s1, s2, s3, s4, s5, s6 = y
        return [s2, v3z1, s4, -u3z1 * s5, s6, v3z2]

    solution3s = odeint(f2, y0, t)
    print("solution1s[1] =", solution1s[-1])
    print("solution2s[2] =", solution2s[-1])
    print("solution3s[3] =", solution3s[-1])

    # возвращаем s1, s3; q1, q3
    res1 = [solution3s[-1][0], solution3s[-1][2]]
    res2 = [solution3q[-1][0], solution3q[-1][2]]

    return [res1, res2]


# функция поиска самой дальней точки от отрезка
def find_far_dot(x1, y1, x2, y2, d1, d2, d3, d4):
    dists = []
    mx = 0
    for dot in [d1, d2, d3, d4]:
        if dot[1] > x2:
            dst = ((dot[1] - x2) ** 2 + (dot[0] - y2) ** 2) ** 0.5
        elif dot[1] < x1:
            dst = ((dot[1] - x1) ** 2 + (dot[0] - y1) ** 2) ** 0.5
        else:
            dst = abs(dot[0] - y1)
        dists.append(dst)
        if dst > mx:
            mxdot = dot
            mx = dst
    print("max dist = ", max(dists), "=", mx)
    return mxdot


# поиск параметров круга
def find_circle_params(x1=0, y1=0, x2=0, y2=0, x3=0, y3=0):
    if x3 > x2:
        dot = [0, x2]
        dist = ((x3 - x2) ** 2 + (y3 - y2) ** 2) ** 0.5
        center = [x2, y2]
    elif x3 < x1:
        dot = [0, x1]
        dist = ((x3 - x1) ** 2 + (y3 - y1) ** 2) ** 0.5
        center = [x1, y1]
    else:
        dist = abs(y3 - y1)
        center = (x3, y1)
    print("dist =", dist, "center =", center)
    return dist, center


if __name__ == '__main__':
    # набор дискретных возмущений, вариант 1
    P1 = solve_DU(0.005, 0, 0, -0.3)
    P2 = solve_DU(-0.005, 0, 0, -0.3)
    P3 = solve_DU(0.003, -0.04, -0.01, 0)
    P4 = solve_DU(0.003, 0, 0, -0.2)
    print("\nP1 =", P1, "\nP2 =", P2, "\nP3 =", P3, "\nP4 =", P4)
    # x_1 = s_1 - q_1
    x1 = P4[0][0] - P4[1][0]
    print("x1 =", x1)
    # x_3 = s_3 - q_3
    x3 = P4[0][1] - P4[1][1]
    print("x3 =", x3)
    # функционал J
    J = np.linalg.norm([x1, x3])
    print("J =", J)

    # первая картинка
    # рисование точек
    dotx = [P[0][1] for P in [P1, P2, P3, P4]]
    doty = [P[0][0] for P in [P1, P2, P3, P4]]
    plt.plot(dotx, doty, 'ro')
    i = 1
    for P in [P1, P2, P3, P4]:
        plt.annotate(str(i), (P[0][1] + 0.05, P[0][0]))
        i += 1
    # рисование отрезка
    right_dot = max(abs(P1[1][1]), abs(P2[1][1]), abs(P3[1][1]), abs(P4[1][1]))
    left_dot = - right_dot
    plt.plot([left_dot, right_dot], [0, 0], color='green')
    ax = plt.gca()
    ax.grid()
    plt.title('Дискретные возмущения')
    plt.savefig('pic1.png')
    plt.cla()

    # вторая картинка
    # рисование точек
    plt.plot(dotx, doty, 'ro')
    i = 1
    for P in [P1, P2, P3, P4]:
        plt.annotate(str(i), (P[0][1] + 0.05, P[0][0]))
        i += 1
    # рисование отрезка
    right_dot = max(abs(P1[1][1]), abs(P2[1][1]), abs(P3[1][1]), abs(P4[1][1]))
    left_dot = - right_dot
    plt.plot([left_dot, right_dot], [0, 0], color='green')
    # поиск параметров круга
    maxdot = find_far_dot(x1=left_dot, y1=0, x2=right_dot, y2=0, d1=P1[0], d2=P2[0], d3=P3[0], d4=P4[0])
    R, center = find_circle_params(x1=-right_dot, y1=0, x2=right_dot, y2=0, x3=maxdot[1], y3=maxdot[0])
    # рисование круга
    plt.annotate('.O', center)
    circle = plt.Circle(xy=center, radius=R, color='b', fill=False)
    ax = plt.gca()
    ax.grid()
    ax.add_patch(circle)
    plt.title('Поиск седловой точки')
    plt.savefig('pic2.png')
