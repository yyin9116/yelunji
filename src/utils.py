import flet as ft
import time
from functools import wraps


no_hover_style = ft.ButtonStyle(
    overlay_color=ft.Colors.TRANSPARENT,  # 使悬停效果透明
)


def measure_search_time(func):
    """装饰器：测量搜索算法的执行时间并在结果中显示"""
    @wraps(func)
    def wrapper(self, *args, **kwargs):
        start_time = time.time()
        result = func(self, *args, **kwargs)
        end_time = time.time()
        
        # 计算耗时（秒）
        duration = end_time - start_time
        
        # 将耗时信息添加到结果中
        if isinstance(result, tuple) and len(result) == 2:
            # 如果返回的是 (best_x, best_y) 格式
            best_x, best_y = result
            return best_x, best_y, duration
        else:
            # 其他返回格式
            return result, duration
    
    return wrapper