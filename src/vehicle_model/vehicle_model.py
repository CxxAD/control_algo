"""
车辆常用模型实现文件 
"""
import math 
import matplotlib.pyplot as plt 
import numpy as np 

sin = math.sin
cos = math.cos
tan = math.tan
# 基类
class SingleStackModelBase():
    
    def update(self,state,control):
        # 给定状态量和控制量，更新模型 
        pass 

    def state_space(self,x_r,u_r):
        pass 


"""
X = x,y,\psi
U = \delta 
"""
class SingleVehicleLon(SingleStackModelBase):
    def __init__(self, L,dt=0.1) -> None:
        self.L = L 
        self.dt = 0.1

    """
    state: x,y,v,psi 
    control:delta
    """
    def update(self,state,control):
        x,y,v,psi  = state
        delta = control 
        x_new = x + v * cos(psi) * self.dt
        y_new = y + v * sin(psi) * self.dt 
        psi_new = psi + v / self.L * tan(delta) * self.dt
        return x_new,y_new,v,psi_new

    def state_space(self, x_r, u_r):
        x,y,v,psi = x_r
        delta = u_r 
        A = np.array([
            [1,0,-v*sin(psi)*self.dt],
            [0,1, v*cos(psi) *self.dt],
            [0,0,1]
        ])
        B = np.array([
            [0],
            [0],
            [v * self.dt / self.L / cos(delta) /cos(delta)]
        ])
        return A,B

"""
X:x,y,psi
U:delta,v
"""
class SingleVehicleControlV(SingleStackModelBase):
    def __init__(self,x,y,psi,v, L,dt=0.1) -> None:
        self.L = L 
        self.dt = 0.1
        self.x = x 
        self.y = y 
        self.psi = psi 
        self.v = v 

    def update(self,X,U):
        x,y,psi = X 
        delta,v = U 
        self.x = self.x + v * cos(self.psi) * self.dt 
        self.y = self.y + v * sin(self.psi) * self.dt
        self.psi = self.psi + v / self.L * delta * self.dt
        # self.v = self.v + a * self.dt 
        self.v = v 

        return [self.x,self.y,self.psi]
    def get_state(self):
        return [self.x,self.y,self.psi,self.v]

    def state_space(self, x_r, u_r):
        x,y,psi = x_r
        delta,v = u_r 
        A = np.array([
            [1,0,-v*sin(psi) * self.dt],
            [0,1,v*cos(psi)* self.dt],
            [0,0,1]
        ])

        B = np.array([
            [0,self.dt * cos(psi)],
            [0,self.dt * sin(psi)],
            [v * self.dt / self.L / cos(delta) / cos(delta),self.dt * tan(delta) / self.L]
        ])
        return A,B


    

"""
X:x,y,v,psi
U:delta,a
"""
class SingleVehicleControlA(SingleStackModelBase):
    def __init__(self, L,dt=0.1) -> None:
        self.L = L 
        self.dt = 0.1
    
    def update(self, state, control):
        x,y,v,psi = state
        delta,a = control

        x_new = x + v * cos(psi) * self.dt 
        y_new = y + v * sin(psi) * self.dt 
        psi_new = psi + v / self.L * tan(delta) * self.dt 
        v_new = v + a * self.dt 
        return x_new,y_new,v_new,psi_new
    # 
    def state_space(self, x_r, u_r):
        x,y,v,psi = x_r
        delta,a = u_r

        A = np.array([
            [1,0,-v*sin(psi) * self.dt , cos(psi) * self.dt],
            [0,1,v * cos(psi) * self.dt, sin(psi) * self.dt],
            [0,0,1,tan(delta) / self.L * self.dt],
            [0,0,0,1]
        ])

        B = np.array([
            [0,0],
            [0,0],
            [v * self.dt / self.L / cos(delta) / cos(delta),0],
            [0,self.dt]
        ])

        return A,B 


def draw_model1(state_seq1,state_seq2):
    x = [state[0] for state in state_seq1]
    y = [state[1] for state in state_seq1]
    plt.plot(x,y,label='update_ans')
    x = [state[0] for state in state_seq2]
    y = [state[1] for state in state_seq2]
    plt.plot(x,y,label='state_space_ans')
    plt.legend()
    plt.show()

    





        

        

