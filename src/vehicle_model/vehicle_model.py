"""
车辆常用模型实现文件 
"""
import os 
import sys
from unittest.mock import Base

# 将当前目录的上一级目录加入到系统路径中
current_file_path = os.path.abspath(__file__)
current_dir = os.path.dirname(current_file_path)
parent_dir = os.path.abspath(os.path.join(current_dir, os.pardir))
if parent_dir not in sys.path:
    sys.path.append(parent_dir)

from utils.math_util import normalize_angle

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
        print("x,y,psi",x,y,psi)
        print("self:x,y,psi",self.x,self.y,self.psi)
        delta,v = U 
        # delta = normalize_angle(delta)
        self.x +=  v * cos(self.psi) * self.dt 
        self.y +=  v * sin(self.psi) * self.dt
        self.psi += v / self.L * tan(delta) * self.dt
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
    
    def reset(self,x,y,psi,v):
        self.x = x 
        self.y = y 
        self.psi = psi 
        self.v = v 


    

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

class CILQRModel(SingleStackModelBase):
    """
    A vehicle model with 4 dof. 
    State - [x, y, vel, theta]
    Control - [acc, yaw_rate]
    """
    def __init__(self, L,dt=0.1):
        self.wheelbase = L
        # self.steer_min = args.steer_angle_limits[0]
        # self.steer_max = args.steer_angle_limits[1]
        # self.accel_min = args.acc_limits[0]
        # self.accel_max = args.acc_limits[1]
        # self.max_speed = args.max_speed
        self.Ts = dt
        
    def update(self, state, control):
        """
        Find the next state of the vehicle given the current state and control input
        """
        # Clips the controller values between min and max accel and steer values
        # control[0] = np.clip(control[0], self.accel_min, self.accel_max)
        # control[1] = np.clip(control[1], state[2]*tan(self.steer_min)/self.wheelbase, state[2]*tan(self.steer_max)/self.wheelbase)
        
        next_state = np.array([state[0] + cos(state[3])*(state[2]*self.Ts + (control[0]*self.Ts**2)/2),
                               state[1] + sin(state[3])*(state[2]*self.Ts + (control[0]*self.Ts**2)/2),
                               state[2] + control[0]*self.Ts,
                               state[3] + control[1]*self.Ts])  # wrap angles between 0 and 2*pi - Gave me error
        return next_state

    def state_space(self,x_r,u_r):
        x,y,v,theta = x_r 
        v_dot,delta = u_r 
        A = np.array([[1, 0, cos(theta)*self.Ts, -(v*self.Ts + (v_dot*self.Ts**2)/2)*sin(theta)],
                      [0, 1, sin(theta)*self.Ts,  (v*self.Ts + (v_dot*self.Ts**2)/2)*cos(theta)],
                      [0, 0,                  1,                                        0],
                      [0, 0,                  1/self.wheelbase * tan(delta),                    1]])
        B = np.array([[self.Ts**2*cos(theta)/2,        0],
                      [self.Ts**2*sin(theta)/2,        0],
                      [         self.Ts,        0],
                      [                0, self.Ts * v / self.wheelbase / (cos(delta)*cos(delta))]])
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

    





        

        

