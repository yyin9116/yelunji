import math
from sko.GA import GA
from sko.DE import DE
import pandas as pd
import matplotlib.pyplot as plt


class Solver:
    def __init__(
        self, 
        T01: int = 1200,
        P01: int = 5,
        expansion_ratio: float = 2.5,
        mass_flow: float = 13,
        rpm: int = 2e4,
        tip_diameter: float =  0.45,
        hub_ratio: float = 0
    ) -> None:
        # 气体常数
        self.R = 287.0  # J/(kg·K)
        self.gamma = 1.4  # 比热比
        self.Cp = self.gamma * self.R / (self.gamma - 1)  # 定压比热
        self.T01 = T01
        self.expansion_ratio = expansion_ratio
        self.T02 = T01 * (1 / self.expansion_ratio)**((self.gamma-1)/self.gamma)
        self.P01 = P01 * 1e5
        self.tip_diameter = tip_diameter
        self.L_u = self.Cp * (self.T01 - self.T02)
        self.u_limit = self.tip_diameter * rpm

        
    

    def solve_inverse_problem(
        self,
        psi: float,
        varphi: float,
        omega: float,
        K: float,
        D_ratio: float,
    ) -> str:
        
        """求解反问题
        
        已知无量纲参数反求速度三角形

        Return: 
        
        完整速度参数
        str: str

        -------
        参数：

        psi: 载荷系数
        varphi: 流量系数
        omega: 反力度
        K: 轴向速比
        D_ratio: 进出口中径比

        """

        u = math.sqrt(self.L_u / psi)
        self.u = u

        v_2_a = varphi * u
        v_1_a = K * v_2_a

        v_1_u = u * (psi / 2 + (1 - omega))
        v_2_u = - u * (psi / 2 - (1 - omega))

        w_1_u = v_1_u - u
        w_2_u = v_2_u - u
        w_1_a = v_1_a
        w_2_a = v_2_a

        v_1 = math.sqrt(v_1_u**2 + v_1_a**2)
        v_2 = math.sqrt(v_2_u**2 + v_2_a**2)
        
        w_1 = math.sqrt(w_1_u**2 + w_1_a**2)
        w_2 = math.sqrt(w_2_u**2 + w_2_a**2)

        alpha_1 = math.degrees(math.atan2(v_1_u, v_1_a))
        alpha_2 = math.degrees(math.atan2(v_2_u, v_2_a))
        beta_1 = math.degrees(math.atan2(w_1_u, w_1_a))
        beta_2 = math.degrees(math.atan2(w_2_u, w_2_a))
        
        return {
            'v_1_u': v_1_u, 'v_1_a': v_1_a, 'v_1': v_1, 'alpha_1': alpha_1,
            'v_2_u': v_2_u, 'v_2_a': v_2_a, 'v_2': v_2, 'alpha_2': alpha_2,
            'w_1_u': w_1_u, 'w_1_a': w_1_a, 'w_1': w_1, 'beta_1': beta_1,
            'w_2_u': w_2_u, 'w_2_a': w_2_a, 'w_2': w_2, 'beta_2': beta_2,
            'u': u
        }

    def solve_direct_problem(
        self,
        v_1_a: float,
        v_1_u: float,
        v_2_a: float,
        v_2_u: float,
        u: float
    ) -> str:
        """求解正问题

        已知速度参数求无量纲参数

        Return: 

        str: str

        -------
        参数:

        v, u
        
        """
        psi = (v_1_u - v_2_u) / u
        varphi = v_2_a / u
        K = v_1_a / v_2_a
        omega = 1 - (v_1_u + v_2_u) / (2 * u)

        v_1 = math.sqrt(v_1_a**2 + v_1_u**2)
        v_2 = math.sqrt(v_2_a**2 + v_2_u**2)
        w_1_u = v_1_u - u
        w_1_a = v_1_a
        w_1 = math.sqrt(w_1_u**2 + w_1_a**2)
        w_2_u = v_2_u -u
        w_2_a = v_2_a
        w_2 = math.sqrt(w_2_u**2 + w_2_a**2)

        # 计算角度
        alpha_1 = math.degrees(math.atan2(v_1_u, v_1_a))
        alpha_2 = math.degrees(math.atan2(v_2_u, v_2_a))
        beta_1 = math.degrees(math.atan2(w_1_u, w_1_a))
        beta_2 = math.degrees(math.atan2(w_2_u, w_2_a))

        return {
            'psi': psi, 'varphi': varphi, 'K': K, 'omega': omega,
            'v_1': v_1, 'v_2': v_2, 'w_1': w_1, 'w_2': w_2,
            'alpha_1': alpha_1, 'alpha_2': alpha_2,
            'beta_1': beta_1, 'beta_2': beta_2
        }

    def calc_turbine_performance(
            self,
            T01: int = 1200,
            P01: int = 5,
            expansion_ratio: float = 2.5,
            mass_flow: float = 13,
            rpm: int = 2e4,
            tip_diameter: float =  0.45,
            hub_ratio: float = 0
    ) -> str:
        """计算涡轮性能参数

        Return: None

        -------
        参数：

        T01: 进口总温 (K)
        P01: 进口总压 (bar)
        expansion_ratio: 膨胀比
        mass_flow: 流量 (kg/s)
        rpm: 转速 (rpm)
        tip_diameter: 叶尖直径 (m)
        hub_ratio: 轮毂比
        """

        # bar -> Pa
        P01 *= 1e5

        u_limit = tip_diameter * rpm
        T02 = T01 * (1/expansion_ratio)**((self.gamma-1)/self.gamma)
        L_u = self.Cp * (self.T01 - T02)

        # 计算中径
        hub_diameter = tip_diameter * hub_ratio
        mean_diameter = (tip_diameter + hub_diameter) / 2

        omega = 2 * math.pi * rpm / 60  # rad/s
        u = omega * mean_diameter / 2

        P02 = P01 / expansion_ratio

        # 计算理想膨胀功
        cp = self.gamma * self.R / (self.gamma - 1)
        T02s = T01 * (P02 / P01) ** ((self.gamma - 1) / self.gamma)
        ideal_work = cp * (T01 - T02s)

        return {
            'mean_diameter': mean_diameter,
            'u': u,
            'ideal_work': ideal_work / 1000  # kJ/kg
        }
    
    def calc_efficiency(
        self, 
        triangle: str, 
        aspect_ratio: float,
        parmas: str,
        is_search: bool = False,
    ) -> float:
        """效率计算
        
        通过半经验损失模型计算效率

        \begin{align*}
        \eta &= \frac{2\psi\eta_{\delta}}{\varphi^2\left[\left(\frac{1}{\zeta_r^2} - 1\right) + \left(\frac{V_{1a}}{V_{2a}}\right)^2\left(\frac{1}{\zeta_s^2} - 1\right)\right] + \frac{\psi^2}{4}\left(\frac{1}{\zeta_r^2} + \frac{1}{\zeta_s^2} - 2\right) + \psi\left\{(1 - \Omega)\left(\frac{1}{\zeta_s^2} - 1\right) + \left[\bar{D}_{2m} - (1 - \Omega)\right]\left(\frac{1}{\zeta_r^2} - 1\right) + 2\right\}} \\
        &\quad + \left[\bar{D}_{2m} - (1 - \Omega)\right]^2 + \left(\frac{1}{\zeta_r^2} - 1\right) + (1 - \Omega)^2\left(\frac{1}{\zeta_s^2} - 1\right)
        \end{align*}

        Return:

        \yta: 效率

        ------
        参数：

        triangle: 反计算出的速度三角形参数

        """
        def cot(value: float) -> float:
            if value != 0:
                return 1 / math.tan(value)
            else:
                raise ZeroDivisionError

        beta_1, beta_2, v_1_a, v_2_a= triangle['beta_1'], abs(triangle['beta_2']), triangle['v_1_a'], triangle['v_2_a']
        # print('beta_1: ', beta_1, '\nbeta_2: ', beta_2)
        # print('r_beta_1: ', math.radians(beta_1), '\nr_beta_2: ', math.radians(beta_2))
        varphi, psi, omega, d_ratio, K = parmas['varphi'], parmas['psi'], parmas['omega'], parmas['d_ratio'], parmas['k']
        # print(f'正在尝试的参数组合: \nvarphi: {varphi:.2f}\npsi: {psi:.2f}\nomega: {omega:.2f}\nK: {K:.2f}')

        # 计算 Kp1 与 Kp2
        if beta_1 + beta_2 < 90:
            K_p_1 = math.sin(math.radians(beta_1)) * math.sin(math.radians(beta_1 + beta_2))
            K_p_2 = math.sin(math.radians(beta_1))
        elif beta_1 + beta_2 >= 90 and beta_1 <= 90:
            K_p_1 = math.sin(math.radians(beta_1)) / math.sin(math.radians(beta_1 + beta_2))
            K_p_2 = math.sin(math.radians(beta_1))
        elif beta_1 > 90:
            K_p_1 = 1 / (math.sin(math.radians(beta_1)) * math.sin(math.radians(beta_1 + beta_2)))
            K_p_2 = 1 / math.sin(math.radians(beta_1))
        else:
            raise ValueError

        # print('K_p_1: ', K_p_1, '\nK_p_2: ', K_p_2)
        # 计算叶型能量损失系数 \xi_p

        xi_denominator = (0.09 * K_p_1 / math.sin(math.radians(beta_2)) + 0.46) * (K_p_2 - math.sin(math.radians(beta_2))) + 0.085
        
        if xi_denominator != 0:
            xi_p = 0.003 / xi_denominator + 0.017 
        else:
            raise ZeroDivisionError
        # print('xi_denominator: ', xi_denominator)
        # 判断展弦比修正系数 K_s, 展弦比(Aspect ratio) $(\frac{h}{b})_2$
        if aspect_ratio >= 2:
            K_s = 1
        else:
            K_s = 1 - 0.25 * math.sqrt(2 - aspect_ratio)

        # 计算二次流能量损失系数 \xi_s
        xi_s = 0.0474 * (math.sin(math.radians(beta_2)) / math.sin(math.radians(beta_1))) * \
            ((cot(math.radians(beta_1)) + cot(math.radians(beta_2))) / aspect_ratio) * K_s + 0.0118
        
        # print('xi_p: ', xi_p, '\n', 'xi_s: ', xi_s)
        # 计算速度损失系数 \zeta_r, \zeta_s

        total_loss = xi_p + xi_s
        total_loss = max(0.05, min(0.95, total_loss))
        zeta_r = math.sqrt(1 - (total_loss))
        zeta_s = zeta_r
        
        # print('zeta_r: ', zeta_r)
        # 计算效率（半经验损失模型）\eta

        term_1 = varphi**2 * ((1 / zeta_r**2 - 1) + ((v_1_a/v_2_a)**2) * (1 / zeta_s**2 - 1))
        term_2 = (psi**2 / 4) * (1 / zeta_r**2 + 1 / zeta_s**2 - 2)
        term_3 = psi * ((1 - omega) * (1 / zeta_s**2 - 1) + (d_ratio - (1 - omega)) * (1 / zeta_r**2 - 1) + 2)
        term_4 = (d_ratio - (1 - omega))**2 + (1 / zeta_r**2 -1) + (1 - omega)**2 * (1 / zeta_r**2 - 1)
        eta_denominator = term_1 + term_2 + term_3 + term_4
        # print(term_1, term_2, term_3, term_4)
        
        # 静叶速度系数 \eta_\delta 理想无损失时为 1
        eta_delta = 1

        if eta_denominator == 0:
            return 0.0
        efficiency = 2 * psi * eta_delta / eta_denominator
        # efficiency = 2 * psi * eta_delta / max(eta_denominator, 1e-6)

        # # 确保效率在[0,1]范围内
        # efficiency = max(0.0, min(1.0, efficiency))
        if is_search:
            return -efficiency
        else:
            return efficiency
        
    
    def problem4search(self, parameter):
        psi, varphi, omega, K = parameter
        parmas = {
            'psi': psi,
            'varphi': varphi,
            'omega': omega,
            'k': K,
            'd_ratio': 1
        }
        
        triangle = self.solve_inverse_problem(psi, varphi, omega, K, 1)
        print('current efficiency: ', self.calc_efficiency(triangle=triangle, aspect_ratio=1, parmas=parmas, is_search=True))

        return self.calc_efficiency(triangle=triangle, aspect_ratio=2, parmas=parmas, is_search=True)
    
    def problem4search4torch(self, parameter):
        psi, varphi, omega, K = parameter
        parmas = {
            'psi': psi,
            'varphi': varphi,
            'omega': omega,
            'k': K,
            'd_ratio': 1
        }
        
        triangle = self.solve_inverse_problem(psi, varphi, omega, K, 1)
        # print('current efficiency: ', self.calc_efficiency(triangle=triangle, aspect_ratio=1, parmas=parmas, is_search=True))

        return self.calc_efficiency(triangle=triangle, aspect_ratio=2, parmas=parmas, is_search=True)

    def search_params_GA(self) -> list:
        """使用遗传算法搜索最优参数

        遗传进化算法

        GA 函数参数：
            func:       问题函数
            n_dim:      搜索参数的数量
            size_pop:   种群数量
            max_iter:   最大迭代次数
            prob_mut:   变异率
            lb:         搜索参数下界
            ub:         搜索参数上界
            precision:  精确度

        Return:

            best_x: list -> 最优参数组合 List
            best_y: list -> 使用效率计算函数计算的出的 efficiency List

        """
        ga = GA(
            func=self.problem4search,
            n_dim=4,
            size_pop=50,
            max_iter=800,
            prob_mut=0.001,
            # psi, varphi, omega, K
            lb=[1.5, 0.4, 0.3, 0.9],
            ub=[2, 0.8, 0.5, 1.1],
            constraint_ueq=[lambda _:self.u - self.u_limit],
            precision=1e-7
        )

        best_x, best_y = ga.run()
        Y_history = pd.DataFrame(ga.all_history_Y)
        fig, ax = plt.subplots(2, 1)
        ax[0].plot(Y_history.index, Y_history.values, '.', color='red')
        Y_history.min(axis=1).cummin().plot(kind='line')
        plt.show()

        print(f'最优参数组合: {best_x}\n最优效率: {-best_y[0]*100:.2f}%')
        return best_x, best_y
    
    def search_params_DE(self) -> list:

        de = DE(
            func=self.problem4search,
            n_dim=4,
            size_pop=50,
            max_iter=800,
            # psi, varphi, omega, K
            lb=[1.5, 0.4, 0.3, 0.9],
            ub=[2, 0.8, 0.5, 1.1],
        )

        best_x, best_y = de.run()

        print(f'最优参数组合: {best_x}\n最优效率: {-best_y[0]*100:.2f}%')
        return best_x, best_y
    
    def search_params_PSO(self) -> list:
        from sko.PSO import PSO
        pso = PSO(
            func=self.problem4search,
            n_dim=4,
            pop=40,
            max_iter=150,
            # psi, varphi, omega, K
            lb=[1.5, 0.4, 0.3, 0.9],
            ub=[2, 0.8, 0.5, 1.1],
        )

        best_x, best_y = pso.run()
        plt.plot(pso.gbest_y_hist)
        plt.show()
        # print(f'最优参数组合: {best_x}\n最优效率: {-best_y[0]*100:.2f}%')
        return best_x, best_y

    def search_params_torch(self):
        import torch
        from torch import nn
        from torch.autograd import Variable

        ranges = torch.tensor([
            [1.5, 2],      # psi 范围
            [0.4, 0.8],    # varphi 范围
            [0.3, 0.5],    # omega 范围
            [0.9, 1.1]     # K 范围
        ])
        params = nn.Parameter(torch.randn(4))
        optimizer = torch.optim.Adam([params], lr=0.1)

        n_iter = 800
        for i in range(n_iter):
            optimizer.zero_grad()
    
            # 将参数映射到目标范围
            sigmoid_params = torch.sigmoid(params)
            scaled_params = ranges[:, 0] + sigmoid_params * (ranges[:, 1] - ranges[:, 0])
            
            # 计算目标函数值
            obj_value = self.problem4search4torch(scaled_params)
            loss = obj_value

            loss.backward()
            optimizer.step()

            # if i % 100 == 0:
                # 打印实际的目标函数值
            print(f"Iteration {i}: objective = {-obj_value.item():.4f}")

    def search_params(self, func):
        best_x, best_y = func()
        best_psi, best_varphi, best_omega, best_k = best_x
        print(best_psi, best_varphi, best_omega, best_k)

   