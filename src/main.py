import flet as ft
import math
import flet.canvas as cv


class VelocityTriangle:
    def __init__(self, u, v_a, v_t):
        """
        初始化速度三角形
        u: 轮缘速度 (m/s)
        v_a: 绝对速度的轴向分量 (m/s)
        v_t: 绝对速度的切向分量 (m/s)
        """
        self.u = u
        self.v_a = v_a
        self.v_t = v_t
        
        # 计算绝对速度
        self.c = math.sqrt(v_a**2 + v_t**2)
        
        # 计算相对速度分量
        self.w_t = v_t - u
        self.w = math.sqrt(v_a**2 + self.w_t**2)
        
        # 计算角度
        self.alpha = math.degrees(math.atan2(v_t, v_a))  # 绝对气流角
        self.beta = math.degrees(math.atan2(self.w_t, v_a))  # 相对气流角


def main(page: ft.Page):
    page.title = "航空发动机涡轮速度三角形可视化"
    page.window_width = 1000
    page.window_height = 800
    page.padding = 20
    page.theme_mode = ft.ThemeMode.LIGHT
    page.scroll = ft.ScrollMode.AUTO
    
    # 初始参数值
    u = 300.0
    v_a_in = 250.0
    v_t_in = 150.0
    v_a_out = 250.0
    v_t_out = -50.0
    
    # 创建入口速度三角形
    inlet_triangle = VelocityTriangle(u, v_a_in, v_t_in)
    outlet_triangle = VelocityTriangle(u, v_a_out, v_t_out)
    
    # 创建绘图控件
    canvas = ft.Container(
        width=800,
        height=500,
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
    
    v_a_in_slider = ft.Slider(
        min=100, 
        max=500, 
        divisions=40, 
        label="入口轴向速度 (v_a1): {value} m/s",
        value=v_a_in,
        expand=True
    )
    
    v_t_in_slider = ft.Slider(
        min=-200, 
        max=400, 
        divisions=60, 
        label="入口切向速度 (v_t1): {value} m/s",
        value=v_t_in,
        expand=True
    )
    
    v_a_out_slider = ft.Slider(
        min=100, 
        max=500, 
        divisions=40, 
        label="出口轴向速度 (v_a2): {value} m/s",
        value=v_a_out,
        expand=True
    )
    
    v_t_out_slider = ft.Slider(
        min=-400, 
        max=200, 
        divisions=60, 
        label="出口切向速度 (v_t2): {value} m/s",
        value=v_t_out,
        expand=True
    )
    
    # 计算结果文本
    results = ft.Column(expand=True, spacing=10)
    
    def draw_triangle(draw, triangle, x, y, scale, label):
        """在画布上绘制速度三角形"""
        # 调整缩放比例，使三角形大小合适
        visual_scale = scale * 0.5
        
        draw.shapes.extend([
        # 绘制轮缘速度 (U)
        cv.Line(
            x, y, 
            x + triangle.u * visual_scale, y, 
            ft.Paint(color=ft.Colors.RED, stroke_width=3)
        ),
        
        # 绘制绝对速度 (C)
        cv.Line(
            x, y, 
            x + triangle.v_t * visual_scale, 
            y - triangle.v_a * visual_scale, 
            ft.Paint(color=ft.Colors.BLUE, stroke_width=3)
        ),
        
        # 绘制相对速度 (W)
        cv.Line(
            x + triangle.u * visual_scale, y, 
            x + triangle.v_t * visual_scale, 
            y - triangle.v_a * visual_scale, 
            ft.Paint(color=ft.Colors.GREEN, stroke_width=3)
        ),
        
        # 添加标签
        cv.Text(
            x + triangle.u * visual_scale + 5, y + 10, 
            f"U = {triangle.u:.1f} m/s", 
            ft.TextStyle(size=14, color=ft.Colors.RED)
        ),
        
        cv.Text(
            x + triangle.v_t * visual_scale + 5, 
            y - triangle.v_a * visual_scale - 20, 
            f"C = {triangle.c:.1f} m/s\nα = {triangle.alpha:.1f}°", 
            ft.TextStyle(size=14, color=ft.Colors.BLUE)
        ),
        
        cv.Text(
            (x + triangle.u * visual_scale + x + triangle.v_t * visual_scale) / 2,
            (y + y - triangle.v_a * visual_scale) / 2 - 10,
            f"W = {triangle.w:.1f} m/s\nβ = {triangle.beta:.1f}°", 
            ft.TextStyle(size=14, color=ft.Colors.GREEN)
        ),
        
        # 添加位置标签
        cv.Text(
            x, y + 20, 
            label, 
            ft.TextStyle(size=16, weight=ft.FontWeight.BOLD)
        )])
    
    def update_canvas(e):
        """更新画布和计算结果"""
        # 更新入口速度三角形
        inlet_triangle.u = u_slider.value
        inlet_triangle.v_a = v_a_in_slider.value
        inlet_triangle.v_t = v_t_in_slider.value
        inlet_triangle.c = math.sqrt(inlet_triangle.v_a**2 + inlet_triangle.v_t**2)
        inlet_triangle.w_t = inlet_triangle.v_t - inlet_triangle.u
        inlet_triangle.w = math.sqrt(inlet_triangle.v_a**2 + inlet_triangle.w_t**2)
        inlet_triangle.alpha = math.degrees(math.atan2(inlet_triangle.v_t, inlet_triangle.v_a))
        inlet_triangle.beta = math.degrees(math.atan2(inlet_triangle.w_t, inlet_triangle.v_a))
        
        # 更新出口速度三角形
        outlet_triangle.u = u_slider.value
        outlet_triangle.v_a = v_a_out_slider.value
        outlet_triangle.v_t = v_t_out_slider.value
        outlet_triangle.c = math.sqrt(outlet_triangle.v_a**2 + outlet_triangle.v_t**2)
        outlet_triangle.w_t = outlet_triangle.v_t - outlet_triangle.u
        outlet_triangle.w = math.sqrt(outlet_triangle.v_a**2 + outlet_triangle.w_t**2)
        outlet_triangle.alpha = math.degrees(math.atan2(outlet_triangle.v_t, outlet_triangle.v_a))
        outlet_triangle.beta = math.degrees(math.atan2(outlet_triangle.w_t, outlet_triangle.v_a))
        
        # 清除画布
        draw.shapes.clear()
        
        # 绘制入口速度三角形
        draw_triangle(draw, inlet_triangle, 150, 300, 0.5, "入口速度三角形")
        
        # 绘制出口速度三角形
        draw_triangle(draw, outlet_triangle, 550, 300, 0.5, "出口速度三角形")
        
        # 添加图例
        cv.Text(
            650, 50, 
            "图例:\n红色: 轮缘速度 (U)\n蓝色: 绝对速度 (C)\n绿色: 相对速度 (W)",
            ft.TextStyle(size=14)
        )
        
        # 更新计算结果
        results.controls.clear()
        results.controls.append(ft.Text("计算结果:", size=18, weight=ft.FontWeight.BOLD))
        results.controls.append(ft.Text(f"入口绝对速度 (C1): {inlet_triangle.c:.1f} m/s"))
        results.controls.append(ft.Text(f"入口绝对气流角 (α1): {inlet_triangle.alpha:.1f}°"))
        results.controls.append(ft.Text(f"入口相对速度 (W1): {inlet_triangle.w:.1f} m/s"))
        results.controls.append(ft.Text(f"入口相对气流角 (β1): {inlet_triangle.beta:.1f}°"))
        results.controls.append(ft.Text(f"出口绝对速度 (C2): {outlet_triangle.c:.1f} m/s"))
        results.controls.append(ft.Text(f"出口绝对气流角 (α2): {outlet_triangle.alpha:.1f}°"))
        results.controls.append(ft.Text(f"出口相对速度 (W2): {outlet_triangle.w:.1f} m/s"))
        results.controls.append(ft.Text(f"出口相对气流角 (β2): {outlet_triangle.beta:.1f}°"))
        
        # 计算功和效率
        work = inlet_triangle.u * (inlet_triangle.v_t - outlet_triangle.v_t) / 1000  # kJ/kg
        results.controls.append(ft.Text(f"涡轮功: {work:.1f} kJ/kg", weight=ft.FontWeight.BOLD))
        
        page.update()
    
    # 添加滑块事件处理
    u_slider.on_change = update_canvas
    v_a_in_slider.on_change = update_canvas
    v_t_in_slider.on_change = update_canvas
    v_a_out_slider.on_change = update_canvas
    v_t_out_slider.on_change = update_canvas
    
    # 程序说明
    explanation = ft.Column(
        expand=True,
        controls=[
            ft.Text("航空发动机涡轮速度三角形解析", size=24, weight=ft.FontWeight.BOLD),
            ft.Text("速度三角形是涡轮机械设计和分析中的基本工具，展示了气流在转子入口和出口的速度矢量关系。"),
            ft.Text("三个关键速度矢量：", weight=ft.FontWeight.BOLD),
            ft.Text("1. 绝对速度 (C) - 相对于静止坐标系的气流速度"),
            ft.Text("2. 轮缘速度 (U) - 转子在半径处的切向速度"),
            ft.Text("3. 相对速度 (W) - 相对于旋转坐标系的气流速度"),
            ft.Text("重要角度：", weight=ft.FontWeight.BOLD),
            ft.Text("• α: 绝对气流角 (绝对速度与轴向的夹角)"),
            ft.Text("• β: 相对气流角 (相对速度与轴向的夹角)"),
            ft.Text("计算公式：", weight=ft.FontWeight.BOLD),
            ft.Text("• 绝对速度: C = √(Cₐ² + Cₜ²)"),
            ft.Text("• 相对速度: W = √(Cₐ² + (Cₜ - U)²)"),
            ft.Text("• 绝对气流角: α = arctan(Cₜ/Cₐ)"),
            ft.Text("• 相对气流角: β = arctan((Cₜ - U)/Cₐ)"),
            ft.Text("• 涡轮功: W = U × (Cₜ₁ - Cₜ₂)"),
        ]
    )
    
    # 初始绘制
    update_canvas(None)
    
    # 创建布局
    page.add(
        explanation,
        ft.Divider(height=20),
        canvas,
        ft.Divider(height=20),
        ft.Row(
            controls=[
                ft.Column(
                    expand=True,
                    controls=[
                        u_slider,
                        v_a_in_slider,
                        v_t_in_slider,
                    ]
                ),
                ft.Column(
                    expand=True,
                    controls=[
                        ft.Container(width=20),
                        results,
                    ]
                )
            ]
        ),
        ft.Row(
            controls=[
                v_a_out_slider,
                v_t_out_slider,
            ]
        )
    )

ft.app(target=main)

