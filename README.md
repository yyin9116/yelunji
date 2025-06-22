# 基于试验设计/遗传算法的单级涡轮一维优化

## 问题描述：

### 输入参数
进口总温1200K、进口总压5bar、膨胀比2.5、转速20000rpm、流量13kg/s、等中径设计、
叶尖直径不大于0.45m，不考虑间隙。

载荷系数参考值：1.5-2.0，反力度参考值：0.3-0.5。


## 项目概述：

实现了一个基于Flet框架的交互式航空发动机涡轮一维优化设计系统。系统提供涡轮机械的速度三角形正问题和反问题求解功能，支持多种优化算法进行参数搜索，并通过可视化界面展示涡轮设计结果。

### 主要功能
- 速度三角形正问题求解（给定速度分量计算无量纲参数）

- 速度三角形反问题求解（给定无量纲参数计算速度分量）

- 涡轮性能参数计算

- 基于半经验损失模型的效率计算

- 多种优化算法进行参数搜索

- 交互式可视化界面

### 系统架构

#### 主要文件结构
- main.py  **程序入口**： 初始化Flet应用程序，创建主窗口并设置参数，启动涡轮设计应用程序
- solver.py **涡轮设计核心算法**：Solver类实现涡轮设计相关计算，使用遗传算法(GA)，粒子群优化(PSO)，差分进化(DE)等
- widget.py  **用户界面组件**：实现应用界面，速度三角形可视化，参数输入控件，计算结果展示以及优化算法选择与执行

## 安装与运行：

### uv管理器

运行桌面应用程序:

```
uv run flet run
```

运行在浏览器中:

```
uv run flet run --web
```

### 运行程序
```
python main.py
```

> 运行部分优化算法需要安装环境

## 功能说明：

### 效率计算

基于半经验损失模型计算涡轮效率：

$$
\begin{align*}
\eta &= \frac{2\psi\eta_{\delta}}{\varphi^2\left[\left(\frac{1}{\zeta_r^2} - 1\right) + \left(\frac{V_{1a}}{V_{2a}}\right)^2\left(\frac{1}{\zeta_s^2} - 1\right)\right] + \frac{\psi^2}{4}\left(\frac{1}{\zeta_r^2} + \frac{1}{\zeta_s^2} - 2\right) + \psi\left\{(1 - \Omega)\left(\frac{1}{\zeta_s^2} - 1\right) + \left[\bar{D}_{2m} - (1 - \Omega)\right]\left(\frac{1}{\zeta_r^2} - 1\right) + 2\right\}} \\
&\quad + \left[\bar{D}_{2m} - (1 - \Omega)\right]^2 + \left(\frac{1}{\zeta_r^2} - 1\right) + (1 - \Omega)^2\left(\frac{1}{\zeta_s^2} - 1\right)
\end{align*}
$$


### 参数优化

#### 支持四种优化算法搜索最优参数：
- 遗传算法(GA)

- 粒子群优化(PSO)

- 差分进化(DE)

- Pytorch优化器(Adam)

### 界面说明

<img src='/src/assets/readme_1.png' width=1000>


### 操作步骤
1. 在**左侧输入区**输入无量纲参数和性能参数

2. 点击"**计算速度三角形**"按钮生成结果

3. 可使用**滑块**微调速度分量

4. 在右侧选择**优化算法**并点击"**开始搜索参数**"

5. 查看搜索结果并点击"一键填入当前值"应用最优参数

6. 此时点击**计算速度三角形**可更新速度三角形图像
   
> 需要注意更新后左小角显示的效率会根据滑块数据计算，以搜索得到的最优效率为准。

## 示例

### 初始方案
程序预设了初始方案：
| 预设参数 | $\varphi$ | $\psi$ | $\omega$ | $K$ | $D_{ratio}$ |
|------|-----------|--------|----------|------|------------|
| 初始值  | 0.9       | 1.8    | 0.4      | 1.0  | 1.0        |

初始效率为：$\eta=78.02\% $
初始速度三角形：
<img src='/src/assets/readme_3.png' width=400>
### 优化算法搜索


使用遗传算法搜索，得到最优参数组合以及计算的最优效率

<img src='/src/assets/readme_2.png' width=400>

### 最终方案
一键填入后，重绘速度三角形得到：

<img src='/src/assets/readme_4.png' width=400>

## 待改进：

- 几种算法搜索得到的最优参数组合结果接近，考虑陷入局部最优或搜索范围太窄，需要调整
- 优化后的效率较低
- 添加其他参数的相关计算

## *Flet 构建应用程序

### Android

```
flet build apk -v
```

For more details on building and signing `.apk` or `.aab`, refer to the [Android Packaging Guide](https://flet.dev/docs/publish/android/).

### iOS

```
flet build ipa -v
```

For more details on building and signing `.ipa`, refer to the [iOS Packaging Guide](https://flet.dev/docs/publish/ios/).

### macOS

```
flet build macos -v
```

For more details on building macOS package, refer to the [macOS Packaging Guide](https://flet.dev/docs/publish/macos/).

### Linux

```
flet build linux -v
```

For more details on building Linux package, refer to the [Linux Packaging Guide](https://flet.dev/docs/publish/linux/).

### Windows

```
flet build windows -v
```

For more details on building Windows package, refer to the [Windows Packaging Guide](https://flet.dev/docs/publish/windows/).