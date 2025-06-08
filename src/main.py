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
        self.v = math.sqrt(v_a**2 + v_t**2)
        
        # 计算相对速度分量
        self.w_t = v_t - u
        self.w = math.sqrt(v_a**2 + self.w_t**2)
        
        # 计算角度
        self.alpha = math.degrees(math.atan2(v_t, v_a))  # 绝对气流角
        self.beta = math.degrees(math.atan2(self.w_t, v_a))  # 相对气流角


def main(page: ft.Page):
    page.title = "航空发动机涡轮速度三角形可视化"
    page.window_width = 1200
    page.window_height = 850
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
        width=750,
        height=450,
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
    
    def draw_triangle(draw, triangle, x, y, scale, label, is_inlet=True):
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
            f"C = {triangle.v:.1f} m/s\nα = {triangle.alpha:.1f}°", 
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
            x, y + 30, 
            label, 
            ft.TextStyle(size=16, weight=ft.FontWeight.BOLD)
        ),

        # 添加旋转方向指示器
        cv.Circle(x + 700, 50 if is_inlet else 250, 20, 
                                  ft.Paint(color=ft.Colors.AMBER, style=ft.PaintingStyle.STROKE, stroke_width=2)),
        cv.Line(
            x + 700, 50 if is_inlet else 250, 
            x + 720, 50 if is_inlet else 250, 
            ft.Paint(color=ft.Colors.AMBER, stroke_width=2)
        ),
        cv.Line(
            x + 700, 50 if is_inlet else 250, 
            x + 695, 40 if is_inlet else 240, 
            ft.Paint(color=ft.Colors.AMBER, stroke_width=2)
        ),
        cv.Line(
            x + 700, 50 if is_inlet else 250, 
            x + 695, 60 if is_inlet else 260, 
            ft.Paint(color=ft.Colors.AMBER, stroke_width=2)
        )])
    
    def draw_common_elements(canvas):
        """绘制公共元素"""
        # 添加标题
        cv.Text(
            300, 20, 
            "涡轮速度三角形可视化", 
            ft.TextStyle(size=24, weight=ft.FontWeight.BOLD, color=ft.Colors.BLUE_800)
        )
        
        # 添加坐标系说明
        cv.Line(50, 50, 50, 400, ft.Paint(color=ft.Colors.BLACK, stroke_width=1))
        cv.Line(50, 400, 700, 400, ft.Paint(color=ft.Colors.BLACK, stroke_width=1))
        cv.Text(30, 220, "轴向", ft.TextStyle(size=12, color=ft.Colors.BLACK))
        cv.Text(370, 420, "切向", ft.TextStyle(size=12, color=ft.Colors.BLACK))
        
        # 添加图例
        cv.Rect(500, 30, 200, 100, ft.Paint(color=ft.Colors.WHITE))
        cv.Text(520, 40, "图例", ft.TextStyle(size=16, weight=ft.FontWeight.BOLD))
        cv.Line(520, 65, 550, 65, ft.Paint(color=ft.Colors.RED, stroke_width=3))
        cv.Text(560, 60, "轮缘速度 (U)", ft.TextStyle(size=12))
        cv.Line(520, 85, 550, 85, ft.Paint(color=ft.Colors.BLUE, stroke_width=3))
        cv.Text(560, 80, "绝对速度 (C)", ft.TextStyle(size=12))
        cv.Line(520, 105, 550, 105, ft.Paint(color=ft.Colors.GREEN, stroke_width=3))
        cv.Text(560, 100, "相对速度 (W)", ft.TextStyle(size=12))
        
        # 添加涡轮叶片示意图
        cv.Line(100, 150, 150, 130, ft.Paint(color=ft.Colors.BROWN, stroke_width=4))
        cv.Line(150, 130, 180, 150, ft.Paint(color=ft.Colors.BROWN, stroke_width=4))
        cv.Line(100, 250, 150, 270, ft.Paint(color=ft.Colors.BROWN, stroke_width=4))
        cv.Line(150, 270, 180, 250, ft.Paint(color=ft.Colors.BROWN, stroke_width=4))
        cv.Text(140, 290, "涡轮叶片", ft.TextStyle(size=12, color=ft.Colors.BROWN))
        
        # 添加流动方向
        cv.Line(200, 200, 230, 200, ft.Paint(color=ft.Colors.DEEP_PURPLE, stroke_width=2))
        # canvas.content.draw_triangle(230, 200, 225, 197, 225, 203, ft.Paint(color=ft.Colors.DEEP_PURPLE))
        cv.Text(240, 195, "气流方向", ft.TextStyle(size=12, color=ft.Colors.DEEP_PURPLE))
    
    
    def update_canvas(e):
        """更新画布和计算结果"""
        # 更新入口速度三角形
        inlet_triangle.u = u_slider.value
        inlet_triangle.v_a = v_a_in_slider.value
        inlet_triangle.v_t = v_t_in_slider.value
        inlet_triangle.v = math.sqrt(inlet_triangle.v_a**2 + inlet_triangle.v_t**2)
        inlet_triangle.w_t = inlet_triangle.v_t - inlet_triangle.u
        inlet_triangle.w = math.sqrt(inlet_triangle.v_a**2 + inlet_triangle.w_t**2)
        inlet_triangle.alpha = math.degrees(math.atan2(inlet_triangle.v_t, inlet_triangle.v_a))
        inlet_triangle.beta = math.degrees(math.atan2(inlet_triangle.w_t, inlet_triangle.v_a))
        
        # 更新出口速度三角形
        outlet_triangle.u = u_slider.value
        outlet_triangle.v_a = v_a_out_slider.value
        outlet_triangle.v_t = v_t_out_slider.value
        outlet_triangle.v = math.sqrt(outlet_triangle.v_a**2 + outlet_triangle.v_t**2)
        outlet_triangle.w_t = outlet_triangle.v_t - outlet_triangle.u
        outlet_triangle.w = math.sqrt(outlet_triangle.v_a**2 + outlet_triangle.w_t**2)
        outlet_triangle.alpha = math.degrees(math.atan2(outlet_triangle.v_t, outlet_triangle.v_a))
        outlet_triangle.beta = math.degrees(math.atan2(outlet_triangle.w_t, outlet_triangle.v_a))
        
        # 清除画布
        draw.shapes.clear()

        draw_common_elements(draw)
        
        # 绘制入口速度三角形
        draw_triangle(draw, inlet_triangle, 300, 200, 0.5, "入口速度三角形", True)
        
        # 绘制出口速度三角形
        draw_triangle(draw, outlet_triangle, 300, 350, 0.5, "出口速度三角形", False)
        
        # 绘制连接线
        cv.Line(
            300 + outlet_triangle.v_t * 0.25, 
            350 - outlet_triangle.v_a * 0.25, 
            300 + inlet_triangle.v_t * 0.25, 
            200 - inlet_triangle.v_a * 0.25, 
            ft.Paint(color=ft.Colors.PURPLE, stroke_width=1, stroke_dash_pattern=[5, 5])
        )
        # 更新计算结果
        results.controls.clear()
        results.controls.append(ft.Text("计算结果:", size=18, weight=ft.FontWeight.BOLD, color=ft.Colors.BLUE_700))
        
        results.controls.append(ft.Row([
            ft.Text("入口绝对速度 (C1):", width=200, weight=ft.FontWeight.BOLD),
            ft.Text(f"{inlet_triangle.v:.1f} m/s")
        ]))
        
        results.controls.append(ft.Row([
            ft.Text("入口绝对气流角 (α1):", width=200, weight=ft.FontWeight.BOLD),
            ft.Text(f"{inlet_triangle.alpha:.1f}°")
        ]))
        
        results.controls.append(ft.Row([
            ft.Text("入口相对速度 (W1):", width=200, weight=ft.FontWeight.BOLD),
            ft.Text(f"{inlet_triangle.w:.1f} m/s")
        ]))
        
        results.controls.append(ft.Row([
            ft.Text("入口相对气流角 (β1):", width=200, weight=ft.FontWeight.BOLD),
            ft.Text(f"{inlet_triangle.beta:.1f}°")
        ]))
        
        results.controls.append(ft.Divider(height=10))
        
        results.controls.append(ft.Row([
            ft.Text("出口绝对速度 (C2):", width=200, weight=ft.FontWeight.BOLD),
            ft.Text(f"{outlet_triangle.v:.1f} m/s")
        ]))
        
        results.controls.append(ft.Row([
            ft.Text("出口绝对气流角 (α2):", width=200, weight=ft.FontWeight.BOLD),
            ft.Text(f"{outlet_triangle.alpha:.1f}°")
        ]))
        
        results.controls.append(ft.Row([
            ft.Text("出口相对速度 (W2):", width=200, weight=ft.FontWeight.BOLD),
            ft.Text(f"{outlet_triangle.w:.1f} m/s")
        ]))
        
        results.controls.append(ft.Row([
            ft.Text("出口相对气流角 (β2):", width=200, weight=ft.FontWeight.BOLD),
            ft.Text(f"{outlet_triangle.beta:.1f}°")
        ]))
        
        results.controls.append(ft.Divider(height=10))
        
        # 计算功和效率
        work = inlet_triangle.u * (inlet_triangle.v_t - outlet_triangle.v_t) / 1000  # kJ/kg
        results.controls.append(ft.Row([
            ft.Text("涡轮功:", width=200, weight=ft.FontWeight.BOLD, color=ft.Colors.RED),
            ft.Text(f"{work:.1f} kJ/kg", weight=ft.FontWeight.BOLD, color=ft.Colors.RED)
        ]))
        
        efficiency = 100 * (work / (inlet_triangle.v**2 / 2000))  # 简化效率计算
        results.controls.append(ft.Row([
            ft.Text("涡轮效率:", width=200, weight=ft.FontWeight.BOLD, color=ft.Colors.GREEN),
            ft.Text(f"{efficiency:.1f}%", weight=ft.FontWeight.BOLD, color=ft.Colors.GREEN)
        ]))

        page.update()
    
    # 添加滑块事件处理
    u_slider.on_change = update_canvas
    v_a_in_slider.on_change = update_canvas
    v_t_in_slider.on_change = update_canvas
    v_a_out_slider.on_change = update_canvas
    v_t_out_slider.on_change = update_canvas
    
    # 程序说明 - 使用Markdown格式
    explanation = ft.Markdown("""
## 航空发动机涡轮速度三角形解析

速度三角形是涡轮机械设计和分析中的基本工具，展示了气流在转子入口和出口的速度矢量关系。

### 关键速度矢量

1. **绝对速度 (C)** - 相对于静止坐标系的气流速度
2. **轮缘速度 (U)** - 转子在半径处的切向速度
3. **相对速度 (W)** - 相对于旋转坐标系的气流速度

### 重要角度

- **α (绝对气流角)**: 绝对速度与轴向的夹角
- **β (相对气流角)**: 相对速度与轴向的夹角

### 计算公式

1. **绝对速度**:
   ```
   C = √(Cₐ² + Cₜ²)
   ```
   
2. **相对速度**:
   ```
   W = √(Cₐ² + (Cₜ - U)²)
   ```
   
3. **绝对气流角**:
   ```
   α = arctan(Cₜ / Cₐ)
   ```
   
4. **相对气流角**:
   ```
   β = arctan((Cₜ - U) / Cₐ)
   ```
   
5. **涡轮功**:
   ```
   Work = U × (Cₜ₁ - Cₜ₂)
   ```

### 设计原则

- 入口角度设计影响冲击损失
- 出口角度设计影响余速损失
- 速度三角形应避免过大的扩散度
- 相对马赫数应控制在合理范围

### 应用领域

- 航空发动机涡轮设计
- 燃气轮机优化
- 蒸汽轮机分析
- 涡轮机械教学演示
""", extension_set="gitHubWeb", selectable=True)
    
    
    # 初始绘制
    update_canvas(None)
    
    # 创建布局
    page.add(
        ft.Row(
            controls=[
                ft.Container(
                    content=explanation,
                    width=400,
                    height=750,
                    padding=15,
                    bgcolor=ft.Colors.GREY_100,
                    border_radius=10,
                    border=ft.border.all(1, ft.Colors.BLUE_100),
                ),
                ft.Column(
                    expand=True,
                    controls=[
                        canvas,
                        ft.Container(
                            content=results,
                            padding=15,
                            bgcolor=ft.Colors.GREY_100,
                            border_radius=10,
                            height=200
                        ),
                        ft.Card(
                            content=ft.Column(
                                controls=[
                                    ft.Text("入口参数", weight=ft.FontWeight.BOLD),
                                    u_slider,
                                    v_a_in_slider,
                                    v_t_in_slider,
                                ]
                            ),
                            elevation=5
                        ),
                        ft.Card(
                            content=ft.Column(
                                controls=[
                                    ft.Text("出口参数", weight=ft.FontWeight.BOLD),
                                    v_a_out_slider,
                                    v_t_out_slider,
                                ]
                            ),
                            elevation=5
                        )
                    ],
                    spacing=15
                )
            ],
            spacing=20
        )
    )

ft.app(target=main)