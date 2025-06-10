import flet as ft
import math
import flet.canvas as cv
import numpy as np
from flet_math import Math
from scipy.optimize import fsolve
from sko.GA import GA


class Solver:
    def __init__(self):
        # 气体常数
        self.R = 287.0  # J/(kg·K)
        self.gamma = 1.4  # 比热比

    def solve_inverse_problem(
        self,
        psi: float,
        varphi: float,
        omega: float,
        K: float,
        D_ratio: float,
        u: float,
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
        u: 轮缘速度 (m/s)

        """

        v_2_a = varphi * u
        v_1_a = K * v_2_a

        v_1_u = u * (varphi / 2 + (1 - omega))
        v_2_u = - u * (varphi / 2 - (1 - omega))

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
                return

        beta_1, beta_2, v_1_a, v_2_a= triangle['beta_1'], triangle['beta_2'], triangle['v_1_a'], triangle['v_2_a']
        print(beta_1, beta_2)
        varphi, psi, omega, d_ratio = parmas['varphi'], parmas['psi'], parmas['omega'], parmas['d_ratio']

        # 计算 Kp1 与 Kp2
        if beta_1 + beta_2 < 90:
            K_p_1 = math.sin(math.radians(beta_1)) * math.sin(math.radians(beta_1 + beta_2))
            K_p_2 = math.sin(math.radians(beta_1))
        elif beta_1 <= 90:
            K_p_1 = math.sin(math.radians(beta_1)) / math.sin(math.radians(beta_1 + beta_2))
            K_p_2 = math.sin(math.radians(beta_1))
        else:
            K_p_1 = 1 / (math.sin(math.radians(beta_1)) * math.sin(math.radians(beta_1 + beta_2)))
            K_p_2 = 1 / math.sin(math.radians(beta_1))

        # 计算叶型能量损失系数 \xi_p

        xi_denominator = (0.09 * K_p_1 / math.sin(math.radians(beta_2)) + 0.46) * (K_p_2 - math.sin(math.radians(beta_2))) + 0.085
        xi_p = 0.003 / xi_denominator + 0.017 if xi_denominator != 0 else 0.0

        # 判断展弦比修正系数 K_s, 展弦比(Aspect ratio) $(\frac{h}{b})_2$
        if aspect_ratio >= 2:
            K_s = 1
        else:
            K_s = 1 - 0.25 * math.sqrt(2 - aspect_ratio)

        # 计算二次流能量损失系数 \xi_s
        xi_s = 0.0474 * (math.sin(math.radians(beta_2)) / math.sin(math.radians(beta_1))) * \
            ((cot(math.radians(beta_1)) + cot(math.radians(beta_2))) / aspect_ratio) * K_s + 0.0118
        
        # 计算速度损失系数 \zeta_r, \zeta_s
        zeta_r = math.sqrt(1 - (xi_p + xi_s))
        zeta_s = zeta_r

        # 计算效率（半经验损失模型）\eta

        term_1 = varphi**2 * ((1 / zeta_r**2 - 1) + ((v_1_a/v_2_a)**2) * (1 / zeta_s**2 - 1))
        term_2 = (psi**2 / 4) * (1 / zeta_r**2 + 1 / zeta_s**2 - 2)
        term_3 = psi * ((1 - omega) * (1 / zeta_s**2 - 1) + (d_ratio - (1 - omega)) * (1 / zeta_r**2 - 1) + 2)
        term_4 = (d_ratio - (1 - omega))**2 + (1 / zeta_r**2 -1) + (1 - omega)**2 * (1 / zeta_r**2 - 1)
        eta_denominator = term_1 + term_2 + term_3 + term_4
        print(term_1, term_2, term_3, term_4)
        # 静叶速度系数 \eta_\delta 理想无损失时为 1
        eta_delta = 1

        if eta_denominator == 0:
            return 0.0
        efficiency = (2 * psi) * eta_delta / eta_denominator

        return efficiency
    
    def problem4search(self, parameter):
        psi, varphi, omega, K = parameter
        parmas = {
            'psi': psi,
            'varphi': varphi,
            'omega': omega,
            'K': K,
            'd_ratio': 1
        }
        u = 200
        triangle = self.solve_inverse_problem(psi, varphi, omega, K, 1, u)
        return -self.calc_efficiency(triangle=triangle, aspect_ratio=1, parmas=parmas)

    def search_params(self) -> str:
        
        ga = GA(
            func=self.problem4search,
            n_dim=4,
            size_pop=50,
            max_iter=800,
            prob_mut=0.001,
            # psi, varphi, omega, K, D_ratio
            lb=[1.5, 0.4, 0.3, 0.9],
            ub=[2, 0.8, 0.5, 1.1],
            precision=1e-7
        )

        best_x, best_y = ga.run()
        # best_psi, best_varphi, best_omega, best_K, best_D_ratio = ga.run()
        # best_params = {                
        #     'psi': best_psi,
        #     'varphi': best_varphi,
        #     'omega': best_omega,
        #     'd_ratio': best_D_ratio,
        # }
        # triangle = self.solve_inverse_problem(best_psi, best_varphi, best_omega, best_K, best_D_ratio)
        # effcience = self.calc_efficiency(triangle=triangle, aspect_ratio=1, parmas=best_params)

        # print(best_params, effcience)

        print(best_x, '\n', best_y)


    
        

def main(page: ft.Page):
    page.title = "单级涡轮一维优化"
    page.window.full_screen = True
    page.padding = 20
    page.theme_mode = ft.ThemeMode.LIGHT
    page.scroll = ft.ScrollMode.AUTO
    
    # 初始参数值 初始化滑块
    u = 350.0
    v_1_a = 250.0
    v_1_u = 150.0
    v_2_a = 250.0
    v_2_u = -50.0
    
    # 创建涡轮设计对象
    solver = Solver()
    
    # 默认参数值
    default_params = {
        'psi': 1.8, 'phi': 0.9, 'omega': 0.4, 'K': 1.0, 'D_ratio': 1.0, 'u': 350,
        'T01': 1200, 'P01': 5, 'expansion_ratio': 2.5, 'rpm': 20000, 'mass_flow': 13, 'tip_diameter': 0.45
    }
    
    # 创建输入控件
    psi_input = ft.TextField(label="载荷系数 (ψ)", value=str(default_params['psi']), width=180)
    phi_input = ft.TextField(label="流量系数 (φ)", value=str(default_params['phi']), width=180)
    omega_input = ft.TextField(label="反力度 (ω)", value=str(default_params['omega']), width=180)
    K_input = ft.TextField(label="轴向速比 (K)", value=str(default_params['K']), width=180)
    D_ratio_input = ft.TextField(label="中径比 (D₂/D₁)", value=str(default_params['D_ratio']), width=180)
    u_input = ft.TextField(label="轮缘速度 (U, m/s)", value=str(default_params['u']), width=180)
    
    # 性能参数输入
    T01_input = ft.TextField(label="进口总温 (K)", value=str(default_params['T01']), width=180)
    P01_input = ft.TextField(label="进口总压 (bar)", value=str(default_params['P01']), width=180)
    expansion_input = ft.TextField(label="膨胀比", value=str(default_params['expansion_ratio']), width=180)
    rpm_input = ft.TextField(label="转速 (rpm)", value=str(default_params['rpm']), width=180)
    mass_flow_input = ft.TextField(label="流量 (kg/s)", value=str(default_params['mass_flow']), width=180)
    tip_diameter_input = ft.TextField(label="叶尖直径 (m)", value=str(default_params['tip_diameter']), width=180)
        
    # 创建绘图控件
    canvas = ft.Container(
        width=1000,
        height=800,
        bgcolor=ft.Colors.GREY_100,
        border_radius=10,
        padding=20,
        content=cv.Canvas(),
    )

    draw = canvas.content
    
    
    # 参数滑块
    u_slider = ft.Slider(
        min=100, 
        max=600, 
        divisions=50, 
        label="轮缘速度 (U): {value} m/s",
        value=u,
        expand=True
    )
    
    v_1_a_slider = ft.Slider(
        min=100, 
        max=500, 
        divisions=40, 
        label="入口轴向速度 (v_1_a): {value} m/s",
        value=v_1_a,
        expand=True
    )
    
    v_1_u_slider = ft.Slider(
        min=-200, 
        max=400, 
        divisions=60, 
        label="入口切向速度 (v_1_u): {value} m/s",
        value=v_1_u,
        expand=True
    )
    
    v_2_a_slider = ft.Slider(
        min=100, 
        max=500, 
        divisions=40, 
        label="出口轴向速度 (v_2_a): {value} m/s",
        value=v_2_a,
        expand=True
    )
    
    v_2_u_slider = ft.Slider(
        min=-400, 
        max=200, 
        divisions=60, 
        label="出口切向速度 (v_2_u): {value} m/s",
        value=v_2_u,
        expand=True
    )
    
    # 计算结果文本
    results = ft.Column(expand=True, spacing=10, width=1000, height=800)
    performance_results = ft.Column(expand=True, spacing=8)
    
    def draw_triangle(
        draw: cv.Canvas,
        triangle, 
        x, 
        y, 
        scale, 
        label, 
        is_inlet=True
    ) -> None:
        """在画布上绘制速度三角形"""
        # 调整缩放比例，使三角形大小合适
        visual_scale = scale * 0.5
        
        draw.shapes.extend(
            [
                # 绘制轮缘速度 (U)
                cv.Line(
                    x, y, 
                    x + triangle['u'] * visual_scale, y, 
                    ft.Paint(color=ft.Colors.RED, stroke_width=3)
                ),
            
                
                # 绘制绝对速度 (C)
                cv.Line(
                    x, y, 
                    x + triangle['v_1_u' if is_inlet else 'v_2_u'] * visual_scale, 
                    y - triangle['v_1_a' if is_inlet else 'v_2_a'] * visual_scale, 
                    ft.Paint(color=ft.Colors.BLUE, stroke_width=3)
                ),
                
                # 绘制相对速度 (W)
                cv.Line(
                    x + triangle['u'] * visual_scale, y, 
                    x + triangle['v_1_u' if is_inlet else 'v_2_u'] * visual_scale, 
                    y - triangle['v_1_a' if is_inlet else 'v_2_a'] * visual_scale, 
                    ft.Paint(color=ft.Colors.GREEN, stroke_width=3)
                ),
                
                # 添加标签
                cv.Text(
                    x + triangle['u'] * visual_scale + 5, y + 10, 
                    f"U = {triangle['u']:.1f} m/s", 
                    ft.TextStyle(size=14, color=ft.Colors.RED)
                ),
                
                cv.Text(
                    x + triangle['v_1_u' if is_inlet else 'v_2_u'] * visual_scale + 5, 
                    y - triangle['v_1_a' if is_inlet else 'v_2_a'] * visual_scale - 20, 
                    f"V = {triangle['v_1' if is_inlet else 'v_2']:.1f} m/s\nα = {triangle['alpha_1' if is_inlet else 'alpha_2']:.1f}°", 
                    ft.TextStyle(size=14, color=ft.Colors.BLUE)
                ),
                
                cv.Text(
                    (x + triangle['u'] * visual_scale + x + triangle['v_1_u' if is_inlet else 'v_2_u'] * visual_scale) / 2,
                    (y + y - triangle['v_1_a' if is_inlet else 'v_2_a'] * visual_scale) / 2 - 10,
                    f"W = {triangle['w_1' if is_inlet else 'w_2']:.1f} m/s\nβ = {triangle['beta_1' if is_inlet else 'beta_2']:.1f}°", 
                    ft.TextStyle(size=14, color=ft.Colors.GREEN)
                ),
                
                # 添加位置标签
                cv.Text(
                    x, y + 30, 
                    label, 
                    ft.TextStyle(size=16, weight=ft.FontWeight.BOLD)
                ),
                
# """
#                 # 添加旋转方向指示器
#                 cv.Circle(x + 700, 50 if is_inlet else 250, 20, 
#                                         ft.Paint(color=ft.Colors.AMBER, style=ft.PaintingStyle.STROKE, stroke_width=2)),
#                 cv.Line(
#                     x + 700, 50 if is_inlet else 250, 
#                     x + 720, 50 if is_inlet else 250, 
#                     ft.Paint(color=ft.Colors.AMBER, stroke_width=2)
#                 ),
#                 cv.Line(
#                     x + 700, 50 if is_inlet else 250, 
#                     x + 695, 40 if is_inlet else 240, 
#                     ft.Paint(color=ft.Colors.AMBER, stroke_width=2)
#                 ),
#                 cv.Line(
#                     x + 700, 50 if is_inlet else 250, 
#                     x + 695, 60 if is_inlet else 260, 
#                     ft.Paint(color=ft.Colors.AMBER, stroke_width=2)
#                 )
# """
            ]
        )
    
    def draw_common_elements(draw: cv.Canvas, u: float):
        """绘制公共元素"""
        # 添加标题
        draw.shapes.extend(
            [
                cv.Text(
                    300, 20, 
                    "涡轮速度三角形可视化", 
                    ft.TextStyle(size=24, weight=ft.FontWeight.BOLD, color=ft.Colors.BLUE_800)
                ),
                
                # 添加坐标系说明
                cv.Line(50, 50, 50, 400, ft.Paint(color=ft.Colors.BLACK, stroke_width=1)),
                cv.Line(50, 400, 700, 400, ft.Paint(color=ft.Colors.BLACK, stroke_width=1)),
                cv.Text(30, 220, "轴向", ft.TextStyle(size=12, color=ft.Colors.BLACK)),
                cv.Text(370, 420, "切向", ft.TextStyle(size=12, color=ft.Colors.BLACK)),
                
                # 添加图例
                # cv.Rect(500, 30, 200, 100, ft.Paint(color=ft.Colors.WHITE)),
                cv.Text(520, 40, "图例", ft.TextStyle(size=16, weight=ft.FontWeight.BOLD)),
                cv.Line(520, 65, 550, 65, ft.Paint(color=ft.Colors.RED, stroke_width=3)),
                cv.Text(560, 60, "轮缘速度 (U)", ft.TextStyle(size=12)),
                cv.Line(520, 85, 550, 85, ft.Paint(color=ft.Colors.BLUE, stroke_width=3)),
                cv.Text(560, 80, "绝对速度 (C)", ft.TextStyle(size=12)),
                cv.Line(520, 105, 550, 105, ft.Paint(color=ft.Colors.GREEN, stroke_width=3)),
                cv.Text(560, 100, "相对速度 (W)", ft.TextStyle(size=12)),
                
                # 添加涡轮叶片示意图
                cv.Line(100, 150, 150, 130, ft.Paint(color=ft.Colors.BROWN, stroke_width=4)),
                cv.Line(150, 130, 180, 150, ft.Paint(color=ft.Colors.BROWN, stroke_width=4)),
                cv.Line(100, 250, 150, 270, ft.Paint(color=ft.Colors.BROWN, stroke_width=4)),
                cv.Line(150, 270, 180, 250, ft.Paint(color=ft.Colors.BROWN, stroke_width=4)),
                cv.Text(140, 290, "涡轮叶片", ft.TextStyle(size=12, color=ft.Colors.BROWN)),
                
                # 添加流动方向
                cv.Line(200, 200, 230, 200, ft.Paint(color=ft.Colors.DEEP_PURPLE, stroke_width=2)),
                # canvas.content.draw_triangle(230, 200, 225, 197, 225, 203, ft.Paint(color=ft.Colors.DEEP_PURPLE)),
                cv.Text(240, 195, "气流方向", ft.TextStyle(size=12, color=ft.Colors.DEEP_PURPLE)),

                # 添加旋转方向
                cv.Arc(650, 200, 40, 40, 0, 2 * math.pi, 
                                    ft.Paint(color=ft.Colors.AMBER, style=ft.PaintingStyle.STROKE, stroke_width=2)),
                cv.Line(690, 200, 685, 195, ft.Paint(color=ft.Colors.AMBER, stroke_width=2)),
                cv.Line(690, 200, 685, 205, ft.Paint(color=ft.Colors.AMBER, stroke_width=2)),
                cv.Text(645, 250, "旋转方向", ft.TextStyle(size=12, color=ft.Colors.AMBER))
            ]
        )
        
    
    def update_canvas(e):
        """更新画布和计算结果"""
        # 更新入口速度三角形

        try:
            # 获取输入参数
            psi = float(psi_input.value)
            varphi = float(phi_input.value)
            omega = float(omega_input.value)
            K = float(K_input.value)
            D_ratio = float(D_ratio_input.value)
            u = float(u_input.value)
            
            # 反问题求解
            triangle = solver.solve_inverse_problem(psi, varphi, omega, K, D_ratio, u)
            
            triangle['u'] = u_slider.value
            triangle['v_1_a'] = v_1_a_slider.value
            triangle['v_2_a'] = v_2_a_slider.value
            triangle['v_1_u'] = v_1_u_slider.value
            triangle['v_2_u'] = v_2_u_slider.value

            params = {
                'psi': psi,
                'varphi': varphi,
                'omega': omega,
                'd_ratio': D_ratio,
            }

            efficiency = solver.calc_efficiency(triangle=triangle, aspect_ratio=1, parmas=params)



            # 清除画布
            draw.shapes.clear()

            draw_common_elements(draw, u)
            
            # 绘制入口速度三角形
            draw_triangle(draw, triangle, 300, 200, 0.5, "入口速度三角形", True)
            
            # 绘制出口速度三角形
            draw_triangle(draw, triangle, 300, 350, 0.5, "出口速度三角形", False)

            search_params_btn =ft.ElevatedButton("开始搜索参数", on_click=lambda _: solver.search_params())
    
            # 更新计算结果
            results.controls.clear()
            results.controls.append(ft.Text("计算结果:", size=18, weight=ft.FontWeight.BOLD, color=ft.Colors.BLUE_700))
            
            results.controls.append(ft.Row([
                ft.Text("入口绝对速度 (V1):", width=200, weight=ft.FontWeight.BOLD),
                ft.Text(f"{triangle['v_1']:.1f} m/s")
            ]))
            
            results.controls.append(ft.Row([
                ft.Text("入口绝对气流角 (α1):", width=200, weight=ft.FontWeight.BOLD),
                ft.Text(f"{triangle['alpha_1']:.1f}°")
            ]))
            
            results.controls.append(ft.Row([
                ft.Text("入口相对速度 (W1):", width=200, weight=ft.FontWeight.BOLD),
                ft.Text(f"{triangle['w_1']:.1f} m/s")
            ]))
            
            results.controls.append(ft.Row([
                ft.Text("入口相对气流角 (β1):", width=200, weight=ft.FontWeight.BOLD),
                ft.Text(f"{triangle['beta_1']:.1f}°")
            ]))
            
            results.controls.append(ft.Divider(height=10))
            
            results.controls.append(ft.Row([
                ft.Text("出口绝对速度 (V2):", width=200, weight=ft.FontWeight.BOLD),
                ft.Text(f"{triangle['v_2']:.1f} m/s")
            ]))
            
            results.controls.append(ft.Row([
                ft.Text("出口绝对气流角 (α2):", width=200, weight=ft.FontWeight.BOLD),
                ft.Text(f"{triangle['alpha_2']:.1f}°")
            ]))
        
            results.controls.append(ft.Row([
                ft.Text("出口相对速度 (W2):", width=200, weight=ft.FontWeight.BOLD),
                ft.Text(f"{triangle['w_2']:.1f} m/s")
            ]))
            
            results.controls.append(ft.Row([
                ft.Text("出口相对气流角 (β2):", width=200, weight=ft.FontWeight.BOLD),
                ft.Text(f"{triangle['beta_2']:.1f}°")
            ]))

            results.controls.append(ft.Row([
                ft.Text("效率( \eta):", width=200, weight=ft.FontWeight.BOLD),
                ft.Text(f"{efficiency:.1f}%")
            ]))

            results.controls.append(search_params_btn)
            
            results.controls.append(ft.Divider(height=10))

            page.update()
        
            

        except Exception as ex:
            print(f"Error: {ex}")

    def calculate_performance(e):
        """计算涡轮性能"""
        try:
            # 获取输入参数
            T01 = float(T01_input.value)
            P01 = float(P01_input.value)
            expansion_ratio = float(expansion_input.value)
            rpm = float(rpm_input.value)
            mass_flow = float(mass_flow_input.value)
            tip_diameter = float(tip_diameter_input.value)
            
            # 计算性能
            perf = solver.calc_turbine_performance(
                T01, P01, expansion_ratio, mass_flow, rpm, tip_diameter
            )
            
            # 更新性能结果
            performance_results.controls.clear()
            performance_results.controls.append(ft.Text("涡轮性能参数:", size=18, weight=ft.FontWeight.BOLD, color=ft.Colors.BLUE_700))
            
            performance_results.controls.append(ft.Row([
                ft.Text("中径:", width=200, weight=ft.FontWeight.BOLD),
                ft.Text(f"{perf['mean_diameter']*1000:.1f} mm")
            ]))
            
            performance_results.controls.append(ft.Row([
                ft.Text("轮缘速度:", width=200, weight=ft.FontWeight.BOLD),
                ft.Text(f"{perf['u']:.1f} m/s")
            ]))
            
            performance_results.controls.append(ft.Row([
                ft.Text("理想膨胀功:", width=200, weight=ft.FontWeight.BOLD),
                ft.Text(f"{perf['ideal_work']:.1f} kJ/kg")
            ]))
            
            # 建议无量纲参数
            performance_results.controls.append(ft.Divider(height=20))
            performance_results.controls.append(ft.Text("建议无量纲参数范围:", size=16, weight=ft.FontWeight.BOLD))
            
            performance_results.controls.append(ft.Row([
                ft.Text("载荷系数 (ψ):", width=200),
                ft.Text("1.5 - 2.0")
            ]))
            
            performance_results.controls.append(ft.Row([
                ft.Text("流量系数 (φ):", width=200),
                ft.Text("0.8 - 1.2")
            ]))
            
            performance_results.controls.append(ft.Row([
                ft.Text("反力度 (ω):", width=200),
                ft.Text("0.3 - 0.5")
            ]))
            
            performance_results.controls.append(ft.Row([
                ft.Text("轴向速比 (K):", width=200),
                ft.Text("0.9 - 1.1")
            ]))
            
            page.update()
        except Exception as ex:
            print(f"Error: {ex}")
    
    # 创建按钮
    calculate_btn = ft.ElevatedButton("计算速度三角形", on_click=update_canvas)
    performance_btn = ft.ElevatedButton("计算性能参数", on_click=calculate_performance)
    
    # 添加滑块事件处理
    u_slider.on_change = update_canvas
    
    v_1_a_slider.on_change = update_canvas
    v_1_u_slider.on_change = update_canvas
    v_2_a_slider.on_change = update_canvas
    v_2_u_slider.on_change = update_canvas
    
    # 程序说明 - 使用Markdown格式
    explanation = ft.Markdown("""
## 航空发动机涡轮一维设计

本程序实现涡轮机械的速度三角形正问题和反问题求解：

### 正问题
给定速度分量计算无量纲参数：
- $ ψ = (V_{1u} - V_{2u}) / U $ (载荷系数)
- $φ = V_{2a} / U$ (流量系数)
- $K = V_{1a} / V_{2a}$ (轴向速比)
- $ω = 1 - (V_{1u} + V_{2u}) / (2U)$ (反力度)
- $D_{2m}/D_{1m}$ (进出口中径比)

### 反问题
给定无量纲参数计算速度分量：
- $V_{2a} = φ \times U$
- $V_{1a} = K × V_{2a}$
- $V_{1u} = U × [2(1-ω) + ψ] / 2$
- $V_{2u} = U × [2(1-ω) - ψ] / 2$

### 设计案例
输入参数：
- 进口总温: 1200 K
- 进口总压: 5 bar
- 膨胀比: 2.5
- 转速: 20000 rpm
- 流量: 13 kg/s
- 叶尖直径: ≤ 0.45 m
- 等中径设计

设计目标：
1. 基于一维设计获得单级涡轮速度三角形
2. 基于半经验损失模型获得级效率
3. 优化获得最优效率设计方案
""", 
    extension_set="gitHubWeb", selectable=True
    )
    
    
    # 初始绘制
    update_canvas(None)
    calculate_performance(None)
    
    # 创建布局
    page.add(
        ft.Row(
            controls=[
                ft.Column(
                    width=400,
                    controls=[
                        explanation,
                        ft.Divider(height=20),
                        ft.Text("无量纲参数输入:", size=18, weight=ft.FontWeight.BOLD),
                        ft.Row([psi_input, phi_input]),
                        ft.Row([omega_input, K_input]),
                        ft.Row([D_ratio_input, u_input]),
                        calculate_btn,
                        
                        ft.Divider(height=20),
                        ft.Text("性能参数输入:", size=18, weight=ft.FontWeight.BOLD),
                        ft.Row([T01_input, P01_input]),
                        ft.Row([expansion_input, rpm_input]),
                        ft.Row([mass_flow_input, tip_diameter_input]),
                        performance_btn,
                        performance_results
                    ],
                    scroll=ft.ScrollMode.AUTO
                ),
                ft.Column(
                    width=400,
                    controls=[
                        
                        ft.Container(
                            content=results,
                            padding=15,
                            bgcolor=ft.Colors.GREY_100,
                            border_radius=10,
                            height=400,
                        ),
                        ft.Divider(),
                        ft.Column(
                            expand=True,
                            controls=[
                                u_slider,
                                v_1_a_slider,
                                v_1_u_slider,
                                v_2_a_slider,
                                v_2_u_slider
                            ]
                        )
                    ],
                    scroll=ft.ScrollMode.AUTO
                ),
                ft.Column(
                    
                    controls=[
                        draw,
                    ]
                )
            
            ],
            spacing=20,
            
        )
    )

ft.app(target=main)