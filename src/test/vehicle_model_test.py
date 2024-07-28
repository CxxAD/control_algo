import os 
import sys 
import numpy as np 
import matplotlib.pyplot as plt 
# print(sys.path)
sys.path.append("D:\per\开发项目\control_algo\src")
# print(sys.path)
from vehicle_model.vehicle_model import SingleVehicleControlV,SingleVehicleControlA,SingleVehicleLon
from utils.ref_traj import RefTraj
from control_algorithm.lqr import LQR
import pprint
import random
import math
import utils.draw as draw

draw_util = draw.DrawUtil()



def get_ref_traj(traj_type = -1,const_arg=1,a = -1,delta = -1):
    return RefTraj(traj_type=traj_type,const_arg=const_arg,a = a,delta =delta)
     

"""
X:x,y,v,psi
U:delta,a
"""
def test_vehicle_model_controlA():
    ref_traj = RefTraj(-1)
    ref_traj.generate_traj_base_seq_control([1 + random.uniform(-0.2,0.2) for i in range(100)],[0.1  for i in range(100)])
    ref = [[ref_traj.x[ind],ref_traj.y[ind],ref_traj.v[ind],ref_traj.psi[ind]]\
                for ind in range(len(ref_traj.x))] # n
    ref_u = np.array([[ref_traj.delta[ind],ref_traj.a[ind]] for ind in range(len(ref_traj.x)-1)]) # n-1

    # base control gener ref 
    pprint.pprint(ref)
    pprint.pprint(ref_u)
    draw_util.draw_state(ref,"x-y-v-psi","ref_traj")
    draw_util.draw_contour(ref_u,"delta-a","ref_traj_U")
    # return 
    # model 
    vehicle_model = SingleVehicleControlA(2.3)

    # algorithm
    lqr = LQR()
    # opt ans 
    opt_traj = [ref[0]] 
    opt_u = []
    # 求最优控制序列 
    step = len(ref_traj.x) 
    # 状态权重
    Q = np.eye(4)
    # 控制权重
    R = np.eye(2)
    for i in range(step - 1):
        # 状态空间
        # ref_ind = -1 
        # old_dis = 1e9
        # for ind, ref_p in enumerate(ref):
        #     if ind == len(ref) - 1:
        #         break
        #     dis = (ref_p[0] - opt_traj[i][0])**2 + (ref_p[1] - opt_traj[i][1])**2
        #     if dis < old_dis:
        #         ref_ind = ind 
        #         old_dis = dis
        ref_ind = i
        A,B = vehicle_model.state_space(ref[ref_ind],ref_u[ref_ind])
        # 求解最优控制序列
        K,P = lqr.solve(A,B,Q,R)
        # 计算控制量
        u = -K @ (np.array(opt_traj[i]) - np.array(ref[ref_ind])) + np.array(ref_u[ref_ind])
        opt_u.append(u)
        # u = [0.1,1]
        # 更新状态量
        opt_traj.append(vehicle_model.update(opt_traj[i],u))
    # pprint.pprint(opt_u)
    draw_util.draw_state(opt_traj,'x-y-v-psi',"opt_traj")
    draw_util.draw_contour(opt_u,'delta-a',"opt_traj_U")

def test_vehicle_model_controlA_sin_curva():
    ref_traj = get_ref_traj(0,3)
    ref_X = [[ref_traj.x[ind],ref_traj.y[ind],ref_traj.psi[ind]]\
              for ind in range(len(ref_traj.x))] # n
    # ref U 
    ref_U = np.array([[ref_traj.delta[ind],ref_traj.v[ind]] for ind in range(len(ref_traj.x)-1)]) # n-1  
    # 参考轨迹
    # ref_traj_base_U = RefTraj() 
    # ref_traj_base_U.generate_traj_base_seq_control(ref_U[:,1],ref_U[:,0])
    # # ref X 
    # ref_X_base_U = [[ref_traj_base_U.x[ind],ref_traj_base_U.y[ind],ref_traj_base_U.v[ind],ref_traj_base_U.psi[ind]]\
    #           for ind in range(len(ref_traj_base_U.x))] # n
    # pprint.pprint(ref_X_base_U)
    # pprint.pprint(ref_U)
    draw_util.draw_state(ref_X,"x-y-psi","ref_traj_base_U") # 
    draw_util.draw_contour(ref_U,"delta-v","ref_traj_base_U")

    # return 

    vehicle_model = SingleVehicleControlV(2.3)

    # algorithm
    lqr = LQR()

    # opt ans 
    opt_traj = [ref_X[0]] 
    opt_u = []
    # 求最优控制序列 
    step = len(ref_traj.x) 
    # 状态权重
    Q = np.array([
        [3,0,0],
        [0,3,0],
        [0,0,1]
    ])
    # 控制权重
    R = np.array([
        [1,0],
        [0,1]
    ])

    # LQR:动规迭代
    Dx = np.array([[1],[1],[1]])
    Du = np.array([[1],[1]])
    A = np.zeros((len(ref_X)-1,len(ref_X[0]),len(ref_X[0])))
    B = np.zeros((len(ref_X)-1,len(ref_X[0]),len(ref_U[0])))
    for ind in range(len(ref_X)-1):
        _A,_B = vehicle_model.state_space(ref_X[ind],ref_U[ind])
        A[ind] = _A
        B[ind] = _B

    opt_traj,opt_u = lqr.solve_iter(np.array(ref_X),np.array(ref_U),A,B,Q,R,Dx,Du)

    # LQR:解黎卡提方程
    # for i in range(step - 1):
    #     # 状态空间
    #     # ref_ind = -1 
    #     # old_dis = 1e9
    #     # for ind, ref_p in enumerate(ref_X):
    #     #     if ind == len(ref_X) - 1:
    #     #         break
    #     #     dis = (ref_p[0] - opt_traj[i][0])**2 + (ref_p[1] - opt_traj[i][1])**2
    #     #     if dis < old_dis:
    #     #         ref_ind = ind 
    #     #         old_dis = dis
    #     ref_ind = i 
    #     A,B = vehicle_model.state_space(ref_X[ref_ind],ref_U[ref_ind])
    #     # 求解最优控制序列
    #     K,P = lqr.solve(A,B,Q,R)
    #     # 计算控制量
    #     u = -K @ (np.array(opt_traj[i]) - np.array(ref_X[ref_ind])) + np.array(ref_U[ref_ind])
    #     opt_u.append(u)
    #     # u = [0.1,1]
    #     # 更新状态量
    #     opt_traj.append(vehicle_model.update(opt_traj[i],u))
    # pprint.pprint(opt_u)
    draw_util.draw_state(opt_traj,"x-y-psi","opt_traj")
    draw_util.draw_contour(opt_u,"delta-v","opt_traj_U")


def test_base_seq_gengerate():
    ref_traj = RefTraj(-1)
    ref_traj.generate_traj_base_seq_control([1 + random.uniform(-0.2,0.2) for i in range(100)],[0.1  for i in range(100)])
    ref_X = [[ref_traj.x[ind],ref_traj.y[ind],ref_traj.v[ind],ref_traj.psi[ind]]\
                for ind in range(len(ref_traj.x))] # n
    ref_U = np.array([[ref_traj.delta[ind],ref_traj.a[ind]] for ind in range(len(ref_traj.x)-1)]) # n-1
    draw_util.draw_state(ref_X,"x-y-v-psi","ref_traj")
    draw_util.draw_contour(ref_U,"delta-a","ref_traj")
    draw()


# test_vehicle_model_controlA()

test_vehicle_model_controlA_sin_curva()
draw_util.draw()
    
# test_base_seq_gengerate()