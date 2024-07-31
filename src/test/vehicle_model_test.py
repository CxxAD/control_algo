import os 
import sys 
import numpy as np 
import matplotlib.pyplot as plt 
# 将当前目录的上一级目录加入到系统路径中
current_file_path = os.path.abspath(__file__)
current_dir = os.path.dirname(current_file_path)
parent_dir = os.path.abspath(os.path.join(current_dir, os.pardir))
if parent_dir not in sys.path:
    sys.path.append(parent_dir)

from vehicle_model.vehicle_model import SingleVehicleControlV,SingleVehicleControlA,SingleVehicleLon,CILQRModel
from utils.ref_traj import RefTraj
from control_algorithm.lqr import LQR
from control_algorithm.iLQR import iLQR
from control_algorithm.PID  import PID
import pprint
import random
import math
import utils.draw as draw

draw_util = draw.DrawUtil()


def get_ref_traj(traj_type = -1,const_arg=1,a = -1,delta = -1,v=-1,control_type=""):
    return RefTraj(traj_type=traj_type,const_arg=const_arg,a = a,delta =delta,v=v,control_type=control_type)
     

def test_vehicle_model_controlA_sin_curva():
    ref_traj = get_ref_traj(traj_type=0,const_arg=3)
    
    # ref_X = [[ref_traj.x[ind],ref_traj.y[ind],ref_traj.psi[ind]]\
    #           for ind in range(len(ref_traj.x))] # n
    # # ref U 
    # ref_U = np.array([[ref_traj.delta[ind],ref_traj.v[ind]] for ind in range(len(ref_traj.x)-1)]) # n-1  
    # 参考轨迹
    # ref_traj_base_U = RefTraj() 
    # ref_traj_base_U.generate_traj_base_seq_control(ref_U[:,1],ref_U[:,0])
    # # ref X 
    # ref_X_base_U = [[ref_traj_base_U.x[ind],ref_traj_base_U.y[ind],ref_traj_base_U.v[ind],ref_traj_base_U.psi[ind]]\
    #           for ind in range(len(ref_traj_base_U.x))] # n
    # pprint.pprint(ref_X_base_U)
    # pprint.pprint(ref_U)
    # ref_traj = get_ref_traj(2,control_type='v-delta',v=[1 + random.uniform(-0.2,0.2) for i in range(100)],delta=[0.1  for i in range(100)])
    # ref_X = [[ref_traj.x[ind],ref_traj.y[ind],ref_traj.psi[ind]]\
    #             for ind in range(len(ref_traj.x))] # n
    # for i in range(len(ref_X)):
        # if i > 40:
            # ref_X[i][0] += 0
            # ref_X[i][1] += 10
    # ref_U = np.array([[ref_traj.delta[ind],ref_traj.v[ind]] for ind in range(len(ref_traj.x)-1)]) # n-1
    # vehicle_model = SingleVehicleControlV(2.3)

    def LQR_solver():
    # algorithm
        vehicle_model = SingleVehicleControlV(*ref_X[0],ref_U[0][1],2.3)
        lqr = LQR()
        # 求最优控制序列 
        x_num = len(ref_traj.x) 
        u_num = len(ref_traj.x) - 1
        xn = len(ref_X[0])
        un = len(ref_U[0])
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
            # 如何选择展开点？以自身轨迹作为展开点？
            _A,_B = vehicle_model.state_space(ref_X[ind+1],ref_U[ind]) # TODO:改为最近点
            A[ind] = _A
            B[ind] = _B
        Lx = np.zeros((u_num,xn,1))
        Lu = np.zeros((u_num,un,1))
        Lxx = np.zeros((u_num,xn,xn))
        Luu = np.zeros((u_num,un,un))
        Lux = np.zeros((u_num,un,xn))
        Lxu = np.zeros((u_num,xn,un))
        for i in range(u_num):
            Lx[i] = Dx
            Lu[i] = Du
            Lxx[i] = Q 
            Lxu[i] = np.zeros((xn,un))
            Lux[i] = np.zeros((un,xn))
            Luu[i] = R 

        
    # self,x_ref,u_ref,A,B,lxx,luu,l_x,l_u,lxu,lux,xn,un,num
        return lqr.solve_iter(ref_X[0],np.array(ref_X),np.array(ref_U),A,B,Lxx,Luu,Lx,Lu,Lxu,Lux,xn,un,u_num,vehicle_model.update)

    # LQR-racati求解，调试完成，关注重点：如何保证参考点正确更新。-->额外控制v，确保自车速度贴近参考速度。
    def LQR_racati_solver():
        vehicle_model = SingleVehicleControlV(*ref_X[0],ref_U[0][1],2.3)
        lqr = LQR()
        pid = PID(1.1,0.5,0.2,2,4,-4)
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
        v_init = ref_U[0][1]
        print("vvvv:",v_init)
        for i in range(700):
            # 状态空
            # ref_ind = -1 
            # old_dis = 1e9
            x = opt_traj[i][0]
            y = opt_traj[i][1]
            # for ind, ref_p in enumerate(ref_X):
            #     if ind == len(ref_X) - 1:
            #         break
            #     dis = (ref_p[0] - opt_traj[i][0])**2 + (ref_p[1] - opt_traj[i][1])**2
            #     if dis < old_dis:
            #         ref_ind = ind 
            #         old_dis = dis
            # ref_ind = i

            # 
            d_x = []
            d_y = []
            d = []
            for k in range(len(ref_X)):
                if k == len(ref_X) - 1:
                    break
                d_x.append(ref_X[k][0]- x)
                d_y.append(ref_X[k][1] -y)
                d.append((d_x[k] ** 2 + d_y[k] ** 2)**0.5)
            # 
            min_d = min(d)
            ref_ind = d.index(min_d)
            # print("ref_ind:",ref_ind,min_index)   
            # ref_ind = i
            A,B = vehicle_model.state_space(ref_X[ref_ind],ref_U[ref_ind])
            # 求解最优控制序列
            K,P = lqr.solve(A,B,Q,R)
            # 计算控制量
            # print("K",K)
            print("ref_ind",ref_ind)
            print("opt_traj - ref_X",np.array(opt_traj[i]) - np.array(ref_X[ref_ind]))
            u = -K @ (np.array(opt_traj[i]) - np.array(ref_X[ref_ind])) + np.array(ref_U[ref_ind])
            # u = [0.1,1]
            pid.set_target(ref_U[ref_ind][1])
            v_u = pid.cal_output(v_init)
            v_init = 0.1 * v_u + v_init
            u = [u[0],v_init]
            opt_u.append(u)
            # 更新状态量
            opt_traj.append(vehicle_model.update(opt_traj[i],u))
        return opt_traj,opt_u

    def iLQR_solver():
        ref_X = [[ref_traj.x[ind],ref_traj.y[ind],ref_traj.v[ind],ref_traj.psi[ind]]\
                for ind in range(len(ref_traj.x))] # n
        # ref U 
        ref_U = np.array([[ref_traj.a[ind],ref_traj.raw_rate[ind]] for ind in range(len(ref_traj.x)-1)]) # n-1 
        draw_util.draw_state(ref_X,"x-y-v-psi","ref_traj_base_U") # 
        draw_util.draw_contour(ref_U,"a-raw_rate","ref_traj_base_U")
        # 第一个点用ref的起点，是ok的
        vehicle_model = CILQRModel(2.3)
        xn = len(ref_X[0])
        un = len(ref_U[0])
        ilqr = iLQR(xn,un)
        return ilqr.opt(vehicle_model.state_space,np.array(ref_X),np.array(ref_U),10,vehicle_model.update,vehicle_model)
    
    opt_traj,opt_u = iLQR_solver()
    # opt_traj,opt_u = LQR_solver()
    # opt_traj,opt_u = LQR_racati_solver()
    draw_util.draw_state(opt_traj,"x-y-psi","opt_traj")
    draw_util.draw_contour(opt_u,"delta-v","opt_traj_U")


def test_base_seq_gengerate():
    # model-v
    ref_traj = get_ref_traj(2,control_type='v-delta',v=[1 + random.uniform(-0.2,0.2) for i in range(100)],delta=[0.1  for i in range(100)])
    ref_X = [[ref_traj.x[ind],ref_traj.y[ind],ref_traj.psi[ind]]\
                for ind in range(len(ref_traj.x))] # n
    ref_U = np.array([[ref_traj.delta[ind],ref_traj.v[ind]] for ind in range(len(ref_traj.x)-1)]) # n-1
    draw_util.draw_state(ref_X,"x-y-psi","ref_traj")
    draw_util.draw_contour(ref_U,"delta-v","ref_traj")
    draw_util.draw()


# test_vehicle_model_controlA()

test_vehicle_model_controlA_sin_curva()
draw_util.draw()
    
# test_base_seq_gengerate()