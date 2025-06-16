import flet as ft
import flet.canvas as cv
import math
from solver import Solver    


class TurbineDesignApp:
    def __init__(self, page: ft.Page):
        page.title = "单级涡轮一维优化"
        page.window.full_screen = False
        page.padding = 20

        page.theme_mode = ft.ThemeMode.LIGHT
        page.scroll = ft.ScrollMode.AUTO
        
        self.page = page
        self.page.window_width = 1200
        self.page.window_height = 800

        u = 392.74
        v_1_a = 353.46
        v_1_u = 589.11
        v_2_a = 353.46
        v_2_u = -117.82
        # 创建绘图控件
        self.canvas = ft.Container(
            width=1000,
            height=800,
            bgcolor=ft.Colors.GREY_100,
            border_radius=10,
            padding=20,
            content=cv.Canvas(),
        )

        self.draw = self.canvas.content
        
        self.v_1_a_slider = ft.Slider(
            min=v_1_a-500, 
            max=v_1_a+500, 
            divisions=40, 
            label="入口轴向速度 (v_1_a): {value} m/s",
            value=v_1_a,
            expand=True
        )
        
        self.v_1_u_slider = ft.Slider(
            min=v_1_u-500, 
            max=v_1_u+500, 
            divisions=60, 
            label="入口切向速度 (v_1_u): {value} m/s",
            value=v_1_u,
            expand=True
        )
        
        self.v_2_a_slider = ft.Slider(
            min=v_2_a-500, 
            max=v_2_a+500, 
            divisions=40, 
            label="出口轴向速度 (v_2_a): {value} m/s",
            value=v_2_a,
            expand=True
        )
        
        self.v_2_u_slider = ft.Slider(
            min=v_2_u-500, 
            max=v_2_u+500, 
            divisions=60, 
            label="出口切向速度 (v_2_u): {value} m/s",
            value=v_2_u,
            expand=True
        )
        
        page.add(self.build())
        page.update()

    def build(self):
        # 创建涡轮设计对象
        self.solver = Solver()
        
        # 默认参数值
        default_params = {
            'psi': 1.8, 'phi': 0.9, 'omega': 0.4, 'K': 1.0, 'D_ratio': 1.0,
            'T01': 1200, 'P01': 5, 'expansion_ratio': 2.5, 'rpm': 20000, 'mass_flow': 13, 'tip_diameter': 0.45
        }
        
        # 创建输入控件
        psi_input = ft.TextField(label="载荷系数 (ψ)", value=str(default_params['psi']), width=120)
        phi_input = ft.TextField(label="流量系数 (φ)", value=str(default_params['phi']), width=120)
        omega_input = ft.TextField(label="反力度 (ω)", value=str(default_params['omega']), width=120)
        K_input = ft.TextField(label="轴向速比 (K)", value=str(default_params['K']), width=120)
        D_ratio_input = ft.TextField(label="中径比 (D₂/D₁)", value=str(default_params['D_ratio']), width=120)
        # u_input = ft.TextField(label="轮缘速度 (U, m/s)", value=str(default_params['u']), width=180)

        # 性能参数输入
        T01_input = ft.TextField(label="进口总温 (K)", value=str(default_params['T01']), width=120)
        P01_input = ft.TextField(label="进口总压 (bar)", value=str(default_params['P01']), width=120)
        expansion_input = ft.TextField(label="膨胀比", value=str(default_params['expansion_ratio']), width=120)
        rpm_input = ft.TextField(label="转速 (rpm)", value=str(default_params['rpm']), width=120)
        mass_flow_input = ft.TextField(label="流量 (kg/s)", value=str(default_params['mass_flow']), width=120)
        tip_diameter_input = ft.TextField(label="叶尖直径 (m)", value=str(default_params['tip_diameter']), width=120)
            
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
                
                params = {
                    'psi': psi,
                    'varphi': varphi,
                    'omega': omega,
                    'd_ratio': D_ratio,
                    'k': K
                }

                # 反问题求解
                triangle = self.solver.solve_inverse_problem(params=params)
                
                # 参数滑块
                triangle['v_1_a'] = self.v_1_a_slider.value
                triangle['v_2_a'] = self.v_2_a_slider.value
                triangle['v_1_u'] = self.v_1_u_slider.value
                triangle['v_2_u'] = self.v_2_u_slider.value

                # 在滑块更新上面的4个速度分量后，triangle的其他参数会自动更新
                triangle = self.update_triangle(triangle, params)

                efficiency = self.solver.calc_efficiency(triangle=triangle, aspect_ratio=1/psi, parmas=params)

                # 清除画布
                self.draw.shapes.clear()
               
                self.draw_common_elements(self.draw)
                
                # 绘制入口速度三角形
                self.draw_triangle(self.draw, triangle, 300, 100, 0.5, "入口速度三角形", True)
                
                # 绘制出口速度三角形
                self.draw_triangle(self.draw, triangle, 300, 100, 0.5, "出口速度三角形", False)
                self.get_calc_result(triangle, efficiency)
                search_params_btn =ft.ElevatedButton("开始搜索参数", on_click=lambda _: self.get_search_result())
                self.results.controls.append(search_params_btn)

                self.page.update()
            
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
                perf = self.solver.calc_turbine_performance(
                    T01, P01, expansion_ratio, mass_flow, rpm, tip_diameter
                )
                
                # 更新性能结果
                self.performance_results.controls.clear()
                self.performance_results.controls.extend([
                    ft.Text("涡轮性能参数:", size=18, weight=ft.FontWeight.BOLD, color=ft.Colors.BLUE_700),
                    ft.Row([
                        ft.Text("中径:", width=200, weight=ft.FontWeight.BOLD),
                        ft.Text(f"{perf['mean_diameter']*1000:.1f} mm")
                    ]),
                    ft.Row([
                        ft.Text("轮缘速度:", width=200, weight=ft.FontWeight.BOLD),
                        ft.Text(f"{perf['u']:.1f} m/s")
                    ]),
                    ft.Row([
                        ft.Text("理想膨胀功:", width=200, weight=ft.FontWeight.BOLD),
                        ft.Text(f"{perf['ideal_work']:.1f} kJ/kg")
                    ]),
                    ft.Divider(height=20),
                    ft.Text("建议无量纲参数范围:", size=16, weight=ft.FontWeight.BOLD),
                    ft.Row([
                        ft.Text("载荷系数 (ψ):", width=200),
                        ft.Text("1.5 - 2.0")
                    ]),
                    ft.Row([
                        ft.Text("流量系数 (φ):", width=200),
                        ft.Text("0.8 - 1.2")
                    ]),
                    ft.Row([
                        ft.Text("反力度 (ω):", width=200),
                        ft.Text("0.3 - 0.5")
                    ]),
                    ft.Row([
                        ft.Text("轴向速比 (K):", width=200),
                        ft.Text("0.9 - 1.1")
                    ])
                ])
                
                self.page.update()
            except Exception as ex:
                print(f"Error: {ex}")

        # 计算结果文本
        self.sliders = ft.Column(expand=True, spacing=10, controls=[
            self.v_1_a_slider,
            self.v_1_u_slider,
            self.v_2_a_slider,
            self.v_2_u_slider
        ])
        self.results = ft.Column(expand=True, spacing=10, height=800)
        self.search_result = ft.Column(expand=True, spacing=10, height=500)
        self.performance_results = ft.Column(expand=True, spacing=8)

        # 创建按钮
        calculate_btn = ft.ElevatedButton("计算速度三角形", on_click=update_canvas)
        performance_btn = ft.ElevatedButton("计算性能参数", on_click=calculate_performance)
        
        # 添加滑块事件处理

        self.v_1_a_slider.on_change = update_canvas
        self.v_1_u_slider.on_change = update_canvas
        self.v_2_a_slider.on_change = update_canvas
        self.v_2_u_slider.on_change = update_canvas
        
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

        self.page.update()

        widget = ft.Container(
            width=self.page.window_width,
            height=self.page.window_height,
            padding=10,
            border_radius=20,
            bgcolor=ft.Colors.WHITE,
            shadow= ft.BoxShadow(
                color=ft.Colors.BLACK,
                offset=ft.Offset(0, 3)
            ),
            content=ft.Row(
                controls=[
                    # 左侧参数输入区
                    ft.Container(
                        width=300,
                        border_radius=10,
                        bgcolor=ft.Colors.GREY_100,
                        padding=ft.padding.only(top=10),
                        content=ft.Column(
                            
                            expand=True,
                            # scroll=ft.ScrollMode.AUTO,
                            controls=[
                                ft.Text("无量纲参数", size=16, weight=ft.FontWeight.BOLD),
                                # ft.Divider(height=10, thickness=1),
                                ft.Row([psi_input, ft.VerticalDivider(width=10), phi_input]),
                                ft.Row([omega_input, ft.VerticalDivider(width=10), K_input]),
                                ft.Row([D_ratio_input, calculate_btn], alignment=ft.MainAxisAlignment.SPACE_BETWEEN),
                                
                                ft.Divider(height=15),
                                
                                ft.Text("性能参数", size=16, weight=ft.FontWeight.BOLD),
                                # ft.Divider(height=10, thickness=1),
                                ft.Row([T01_input, ft.VerticalDivider(width=10), P01_input]),
                                ft.Row([expansion_input, ft.VerticalDivider(width=10), rpm_input]),
                                ft.Row([mass_flow_input, ft.VerticalDivider(width=10), tip_diameter_input]),
                                performance_btn,
                                # self.performance_results，
                                ft.Divider(height=15),
                                ft.Text("计算结果", size=16, weight=ft.FontWeight.BOLD),
                                ft.Divider(height=10, thickness=1),
                                ft.Container(
                                    content=self.results,
                                    padding=10,
                                    border_radius=8,
                                    bgcolor=ft.Colors.WHITE,
                                    expand=True
                                ),
                            ]
                        )
                    ),
                    
                    # 垂直分隔线
                    ft.VerticalDivider(width=10, color=ft.Colors.TRANSPARENT),
                    
                    # 中间结果展示区
                    ft.Container(
                        width=300,
                        border_radius=10,
                        bgcolor=ft.Colors.GREY_100,
                        padding=15,
                        content=ft.Row(
                            expand=True,
                            scroll=ft.ScrollMode.AUTO,
                            controls=[
                                ft.Container(
                                    content=ft.Column(
                                        controls=[
                                            ft.Container(
                                                width=400,
                                                border_radius=10,
                                                bgcolor=ft.Colors.GREY_100,
                                                padding=15,
                                                content=self.draw
                                            ),
                                        ]
                                    )
                                ),
                                ft.VerticalDivider(width=10, color=ft.Colors.TRANSPARENT),
                                ft.Container(
                                    content=ft.Column(
                                        controls=[
                                            ft.Text("参数搜索", size=16, weight=ft.FontWeight.BOLD),
                                            ft.Divider(height=10, thickness=1),
                                            ft.Container(
                                                content=self.search_result,
                                                padding=10,
                                                border_radius=8,
                                                bgcolor=ft.Colors.WHITE,
                                                height=200
                                            ),
                                            
                                            ft.Divider(height=15),
                                            
                                            ft.Text("参数调节", size=16, weight=ft.FontWeight.BOLD),
                                            ft.Divider(height=10, thickness=1),
                                            ft.Container(
                                                content=self.sliders,
                                                padding=10,
                                                border_radius=8,
                                                bgcolor=ft.Colors.WHITE
                                            )
                                        ]
                                    )
                                )
                                
                            ]
                        )
                    ),
                
                ],
                expand=True,
                spacing=0
            )
        )
        
        return widget
        
    
    def draw_triangle(
            self,
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
                        x - triangle['w_1_u' if is_inlet else 'w_2_u'] * visual_scale, y + triangle['v_1_a' if is_inlet else 'v_2_a'] * visual_scale, 
                        x - triangle['v_1_u' if is_inlet else 'v_2_u'] * visual_scale, y + triangle['v_1_a' if is_inlet else 'v_2_a'] * visual_scale, 
                        ft.Paint(color=ft.Colors.RED, stroke_width=3)
                    ),      
                    
                    # 绘制绝对速度 (C)
                    cv.Line(
                        x, y, 
                        x - triangle['v_1_u' if is_inlet else 'v_2_u'] * visual_scale, 
                        y + triangle['v_1_a' if is_inlet else 'v_2_a'] * visual_scale, 
                        ft.Paint(color=ft.Colors.BLUE, stroke_width=3)
                    ),
                    
                    # 绘制相对速度 (W)
                    cv.Line(
                        x, y, 
                        x - triangle['w_1_u' if is_inlet else 'w_2_u'] * visual_scale, 
                        y + triangle['v_1_a' if is_inlet else 'v_2_a'] * visual_scale, 
                        ft.Paint(color=ft.Colors.GREEN, stroke_width=3)
                    ),
                    
                    cv.Text(
                        x - (triangle['v_1_u' if is_inlet else 'v_2_u'] + triangle['w_1_u' if is_inlet else 'w_2_u']) / 2 * visual_scale, 
                        y + triangle['v_1_a' if is_inlet else 'v_2_a'] * visual_scale + 20,
                        f"U = {triangle['u']:.1f} m/s", 
                        ft.TextStyle(size=14, color=ft.Colors.RED)
                    ),

                    # 绝对速度 
                    cv.Text(
                        x - triangle['v_1_u' if is_inlet else 'v_2_u'] * visual_scale, 
                        y + triangle['v_1_a' if is_inlet else 'v_2_a'] * visual_scale - 25,
                        f"V = {triangle['v_1' if is_inlet else 'v_2']:.1f} m/s\nα = {triangle['alpha_1' if is_inlet else 'alpha_2']:.1f}°", 
                        ft.TextStyle(size=14, color=ft.Colors.BLUE)
                    ),

                    # 相对速度 
                    cv.Text(
                        x - (triangle['w_1_u' if is_inlet else 'w_2_u'] * visual_scale) / 2 - 10, 
                        y + (triangle['v_1_a' if is_inlet else 'v_2_a'] * visual_scale) / 2 - 20,
                        f"W = {triangle['w_1' if is_inlet else 'w_2']:.1f} m/s\nβ = {triangle['beta_1' if is_inlet else 'beta_2']:.1f}°", 
                        ft.TextStyle(size=14, color=ft.Colors.GREEN)
                    ),

                    # 添加位置标签
                    cv.Text(
                        x, y + 30, 
                        label, 
                        ft.TextStyle(size=16, weight=ft.FontWeight.BOLD)
                    ),
                ]
            )
        
    def draw_common_elements(self, draw: cv.Canvas):
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
                cv.Text(140, 2120, "涡轮叶片", ft.TextStyle(size=12, color=ft.Colors.BROWN)),
                
                # 添加流动方向
                cv.Line(200, 200, 230, 200, ft.Paint(color=ft.Colors.DEEP_PURPLE, stroke_width=2)),
                # canvas.content.draw_triangle(230, 200, 225, 197, 225, 203, ft.Paint(color=ft.Colors.DEEP_PURPLE)),
                cv.Text(240, 195, "气流方向", ft.TextStyle(size=12, color=ft.Colors.DEEP_PURPLE)),

                # 添加旋转方向
                cv.Arc(650, 200, 40, 40, 0, 2 * math.pi, 
                                    ft.Paint(color=ft.Colors.AMBER, style=ft.PaintingStyle.STROKE, stroke_width=2)),
                cv.Line(6120, 200, 685, 195, ft.Paint(color=ft.Colors.AMBER, stroke_width=2)),
                cv.Line(6120, 200, 685, 205, ft.Paint(color=ft.Colors.AMBER, stroke_width=2)),
                cv.Text(645, 250, "旋转方向", ft.TextStyle(size=12, color=ft.Colors.AMBER))
            ]
        )
        
    def get_search_result(self):
        best_x, best_y = self.solver.search_params(self.solver.search_params_DE)
        method = "差分进化" if self.solver.search_params_DE else "粒子群优化"
        best_psi, best_varphi, best_omega, best_k = best_x

        self.search_result.controls.clear()
        self.search_result.controls.extend([
            ft.Text("搜索结果:", size=18, weight=ft.FontWeight.BOLD, color=ft.Colors.BLUE_700),
            ft.Row([
                ft.Text("搜索算法:", width=200, weight=ft.FontWeight.BOLD),
                ft.Text(f"{method}")]),
            ft.Row([
                ft.Text("最优载荷系数(psi):", width=200, weight=ft.FontWeight.BOLD),
                ft.Text(f"{best_psi:.3f}")]),
            ft.Row([
                ft.Text("最优流量系数 (varphi):", width=200, weight=ft.FontWeight.BOLD),
                ft.Text(f"{best_varphi:.3f}")]),
            ft.Row([
                ft.Text("最优反力度(omega):", width=200, weight=ft.FontWeight.BOLD),
                ft.Text(f"{best_omega:.3f}")]),
            ft.Row([
                ft.Text("最优轴向速比(K):", width=200, weight=ft.FontWeight.BOLD),
                ft.Text(f"{best_k:.3f}")]),
            ft.Row([
                ft.Text("最优效率:", width=200, weight=ft.FontWeight.BOLD),
                ft.Text(f"{-best_y[0]*100:.2f}%")]),  
            ft.Row([
                ft.TextButton("一键填入当前值", on_click=lambda _: self.get_search_result(), width=200)])
        ])

        # print(f"Best params: psi={best_psi}, varphi={best_varphi}, omega={best_omega}, k={best_k}, efficiency={-best_y[0]*100:.2f}%")
        self.page.update()

    def get_calc_result(self, triangle, efficiency):
        # 更新计算结果
        self.results.controls.clear()
        self.results.controls.extend([
            ft.Text("计算结果:", size=18, weight=ft.FontWeight.BOLD, color=ft.Colors.BLUE_700),
            ft.Row([
                ft.Text("入口绝对速度 (V1):", width=200, weight=ft.FontWeight.BOLD),
                ft.Text(f"{triangle['v_1']:.1f} m/s")
            ]),
            ft.Row([
                ft.Text("入口绝对气流角 (α1):", width=200, weight=ft.FontWeight.BOLD),
                ft.Text(f"{triangle['alpha_1']:.1f}°")
            ]),
            ft.Row([
                ft.Text("入口相对速度 (W1):", width=200, weight=ft.FontWeight.BOLD),
                ft.Text(f"{triangle['w_1']:.1f} m/s")
            ]),
            ft.Row([
                ft.Text("入口相对气流角 (β1):", width=200, weight=ft.FontWeight.BOLD),
                ft.Text(f"{triangle['beta_1']:.1f}°")
            ]),
            ft.Divider(height=10),
            ft.Row([
                ft.Text("出口绝对速度 (V2):", width=200, weight=ft.FontWeight.BOLD),
                ft.Text(f"{triangle['v_2']:.1f} m/s")
            ]),
            ft.Row([
                ft.Text("出口绝对气流角 (α2):", width=200, weight=ft.FontWeight.BOLD),
                ft.Text(f"{triangle['alpha_2']:.1f}°")
            ]),
            ft.Row([
                ft.Text("出口相对速度 (W2):", width=200, weight=ft.FontWeight.BOLD),
                ft.Text(f"{triangle['w_2']:.1f} m/s")
            ]),
            ft.Row([
                ft.Text("出口相对气流角 (β2):", width=200, weight=ft.FontWeight.BOLD),
                ft.Text(f"{triangle['beta_2']:.1f}°")
            ]),
            ft.Row([
                ft.Text("轮缘速度 (U):", width=200, weight=ft.FontWeight.BOLD),
                ft.Text(f"{triangle['u']:.1f} m/s")
            ]),
            ft.Row([
                ft.Text("效率(eta):", width=200, weight=ft.FontWeight.BOLD),
                ft.Text(f"{efficiency*100:.2f}%")
            ]),
        ])

        # self.page.update()

    def update_triangle(self, triangle, params):
        """更新速度三角形参数"""
        # psi = params['psi']
        # varphi = params['varphi']
        # omega = params['omega']
        # K = params['k']
        # d_ratio = params['d_ratio']
        # 计算其他参数
        triangle['v_1'] = math.sqrt(triangle['v_1_a']**2 + triangle['v_1_u']**2)
        triangle['alpha_1'] = math.degrees(math.atan2(triangle['v_1_a'], triangle['v_1_u']))
        
        triangle['w_1_u'] = triangle['v_1_u'] - triangle['u']
        triangle['w_1_a'] = triangle['v_1_a'] 
        triangle['w_1'] = math.sqrt(triangle['w_1_u']**2 + triangle['v_1_a']**2)
        triangle['beta_1'] = math.degrees(math.atan2(triangle['w_1_a'], triangle['w_1_u']))
        
        # 更新出口速度三角形参数
        triangle['v_2'] = math.sqrt(triangle['v_2_a']**2 + triangle['v_2_u']**2)
        triangle['alpha_2'] = 180 - math.degrees(math.atan2(triangle['v_2_a'], triangle['v_2_u']))
        
        triangle['w_2_a'] = triangle['v_2_a']
        triangle['w_2_u'] = triangle['v_2_u'] - triangle['u']
        triangle['w_2'] = math.sqrt(triangle['w_2_u']**2 + triangle['w_2_a']**2)
        triangle['beta_2'] = 180 - math.degrees(math.atan2(triangle['w_2_a'], triangle['w_2_u']))

        # triangle['alpha_1'] = math.degrees(math.atan2(varphi*K, (psi/2) + (1-omega)))
        # triangle['alpha_2'] = math.degrees(math.atan2(varphi, (psi/2) - (1-omega)))
        # triangle['beta_1'] = math.degrees(math.atan2(varphi*K, (psi/2) - omega))
        # triangle['alpha_2'] = math.degrees(math.atan2(varphi, ((psi/2) - (1-omega) + d_ratio)))  
        return triangle


    
        