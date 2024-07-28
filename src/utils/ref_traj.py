"""
用于实现参考轨迹的生成 

参考轨迹变量包含
x,y,psi,v,a,delta （任一种轨迹都会提供如上变量)

traj_type 
    0: sin曲线,控制周期 ,const_arg控制周期， 2pi / const_arg
    1: 直线， const_arg控制速度
    2: 常控制量曲线
    -1: 序列控制量轨迹 a,delta 
"""

import numpy as np
import math
import random

class RefTraj:
    def __init__(self,traj_type = 0,const_arg=1,control_type="",a = -1,delta = -1,v=-1,total_time = 10,dt =0.1) -> None:
        self.x = [] 
        self.y = [] 
        self.psi = [] 
        self.v = [] 
        self.a = []
        self.delta = [] 
        self.total_time = total_time
        self.dt = dt 
        self.num = int(self.total_time / self.dt)
        self.L = 2.3 
        # base tra 
        if traj_type == 0:
            self.generate_sin_curve(const_arg)
        elif traj_type == 1:
            self.generate_straight_line(const_arg)
        elif traj_type == 2:
            self.generate_control_curve(control_type,a,delta,v)
        else:
            raise ValueError("Please input correct traj_type or add type ")
    def generate_sin_curve(self,k):
        self.clear()
        for i in range(1000):
            self.x.append(0.1 * i)
            self.y.append(2 * math.sin(self.x[-1] / k))

        for i in range(len(self.x)):
            if i == 0:
                dx = self.x[i+1] - self.x[i]
                dy = self.y[i+1] - self.y[i]
                ddx = self.x[i+2] - 2 * self.x[i+1] + self.x[i]
                ddy = self.y[i+2] - 2 * self.y[i+1] + self.y[i]

            elif i == len(self.x) - 1:
                dx = self.x[i] - self.x[i-1]
                dy = self.y[i] - self.y[i-1]
                ddx = self.x[i] - 2 * self.x[i-1] + self.x[i-2]
                ddy = self.y[i] - 2 * self.y[i-1] + self.y[i-2]

            else:
                dx = self.x[i+1] - self.x[i-1]
                dy = self.y[i+1] - self.y[i-1]
                ddx = self.x[i+1] - 2 * self.x[i] + self.x[i-1]
                ddy = self.y[i+1] - 2 * self.y[i] + self.y[i-1]
        
            self.psi.append(math.atan(dy/dx)) # 车身横摆角
            # 计算曲率
            kappa = (ddy * dx - ddx * dy) / (dx ** 2 + dy ** 2) ** 1.5
            self.delta.append(math.atan(kappa * self.L))
            # 速度信息
            self.v.append((dy*dy + dx*dx)**0.5 / self.dt)  
            self.a.append(0)
        
    def generate_straight_line(self,k):
        self.clear()
        t = 0
        while t < self.total_time:
            self.x.append(t)
            self.y.append(k*t)
            self.psi.append(k)
            self.v.append(1)
            self.a.append(0)
            self.delta.append(0)
            t += self.dt
    
    def generate_control_curve(self,control_type,a,delta,v):
        # 判断控制变量类型
        con_type = control_type.split('-')
        for i in range(len(con_type)):
            if con_type[i] == 'a':
                if isinstance(a,list):
                    self.a = a 
                else:
                    self.a = [a for i in range(self.num)]
            elif con_type[i] == 'delta':
                if isinstance(delta,list):
                    self.delta = delta
                else:
                    self.delta = [delta for i in range(self.num)]
            elif con_type[i] == 'v':
                if isinstance(v,list):
                    self.v = v
                else:
                    self.v = [v for i in range(self.num)]
            else:
                raise ValueError("Please input correct control type")
            
        # 因为是基于控制去生成轨迹，取前 num-1个控制量 
        if len(self.a)>0 and len(self.delta)>0:
            self._gen_by_delta_a(list(self.delta),list(self.a))
        elif len(self.v) and len(self.delta):
            self._gen_by_delta_v(list(self.delta),list(self.v))

    def _gen_by_delta_a(self,delta,a):
        self.clear()
        self.x.append(0)
        self.y.append(0)
        self.psi.append(0)
        self.v.append(0)
        print(delta,a)
        for ind in range(len(a)):
            self.x.append(self.v[-1] * math.cos(self.psi[-1]) * self.dt + self.x[-1] )
            self.y.append(self.v[-1] * math.sin(self.psi[-1]) * self.dt + self.y[-1] )
            self.psi.append(self.v[-1] / self.L * math.tan(delta[ind]) * self.dt + self.psi[-1] )
            self.v.append(self.v[-1] + a[ind] * self.dt )
            self.a.append(a[ind])
            self.delta.append(delta[ind])
        print(self.x,self.y)

    
    def _gen_by_delta_v(self,delta,v):
        self.clear()
        self.x.append(0)
        self.y.append(0)
        self.psi.append(0)
        for ind in range(len(v)):
            self.x.append(v[ind] * math.cos(self.psi[-1]) * self.dt + self.x[-1] )
            self.y.append(v[ind] * math.sin(self.psi[-1]) * self.dt + self.y[-1] )
            self.psi.append(v[ind] / self.L * math.tan(delta[ind]) * self.dt + self.psi[-1])
            self.v.append(v[ind])
            self.a.append(0)
            self.delta.append(delta[ind])


    # valid pass
    def generate_traj_base_seq_control(self,a_seq,delta_seq,L=2.3):
        self.clear()
        self.x.append(0)
        self.y.append(0)
        self.v.append(0)
        self.psi.append(np.pi / 2 )
        for ind in range(len(a_seq)):
            self.x.append(self.v[-1] * math.cos(self.psi[-1]) * self.dt + self.x[-1] )
            self.y.append(self.v[-1] * math.sin(self.psi[-1]) * self.dt + self.y[-1] )
            self.psi.append(self.v[-1] / L * math.tan(delta_seq[ind]) * self.dt + self.psi[-1] )
            self.v.append(self.v[-1] + a_seq[ind] * self.dt )
            self.a.append(a_seq[ind])
            self.delta.append(delta_seq[ind])

    def clear(self):
        self.x.clear()
        self.y.clear()
        self.psi.clear()
        self.v.clear()
        self.a.clear()
        


if __name__ == "__main__":
    import draw
    import pprint
    draw_util = draw.DrawUtil()
    # a-delta 
    # ref_tra = RefTraj(traj_type=2,control_type="a-delta",a=[0.1 for i in range(100)],delta=[0.1 for i in range(100)])
    # X = [[ref_tra.x[ind],ref_tra.y[ind],ref_tra.psi[ind],ref_tra.v[ind]] for ind in range(len(ref_tra.x))]
    # pprint.pprint(X)
    # draw_util.draw_state(X,"x-y-psi-v","ref_traj")
    # draw_util.draw_contour([[ref_tra.a[ind],ref_tra.delta[ind]] for ind in range(len(ref_tra.a))],"a-delta","ref_traj_U")
    # draw_util.draw()

    # v-delta 
    ref_tra = RefTraj(traj_type=2,control_type="v-delta",v=[0.1 for i in range(100)],delta=[0.1 for i in range(100)])
    X = [[ref_tra.x[ind],ref_tra.y[ind],ref_tra.psi[ind]] for ind in range(len(ref_tra.x))]
    pprint.pprint(X)
    draw_util.draw_state(X,"x-y-psi","ref_traj")
    draw_util.draw_contour([[ref_tra.v[ind],ref_tra.delta[ind]] for ind in range(len(ref_tra.v))],"v-delta","ref_traj_U")
    draw_util.draw()

