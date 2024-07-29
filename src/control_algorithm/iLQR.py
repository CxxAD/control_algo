"""
iLQR求解器:假设模型参考量个数 3,控制量个数2
Q = diag([1,1,1])
min c = \sum_0^{N-1} 1/2(x^T @ Q @ x + u^T @ R @ u) + D_x^T @ x + D_u^T @ u
min delta_c = sum(0,N-1) 1/2(dot{x}^T @ \Delta/dx @ dot{x} + dot{u}^T @ \Delta/du @ dot{u}) + \nolba/dx @ dot{x} + \nolba /du @ dot{u}

\Delta/dx = Q
\Delta/du = R
\nolba/dx = Q @ x + Dx
\nolba/du = R @ u + Du

lx = Q @ x + Dx
lu = R @ u + Du
lxx = Q
luu = R
lux = 0
"""
import os 
import sys 
current_file_path = os.path.abspath(__file__)
current_dir = os.path.dirname(current_file_path)
parent_dir = os.path.abspath(os.path.join(current_dir, os.pardir))
if current_dir not in sys.path:
    sys.path.append(current_dir)
    sys.path.append(parent_dir)

import numpy as np
import lqr

from utils.draw import DrawUtil


drwa_util = DrawUtil()
class iLQR(lqr.LQR):
    def __init__(self,xn,un) -> None:
        # 初始化原始代价函数的各参数
        # self.Q = np.zeros((xn,xn))
        # self.R = np.zeros((un,un))
        # self.Dx = np.zeros((xn,1))
        # self.Du = np.zeros((un,1))

        # for i in range(xn-1):
        #     self.Q[i,i] = 3
        #     self.Dx[i] = 1 
        # self.Q[xn-1,xn-1] = 1

        # for i in range(un):
        #     self.R[i,i] = 1
        #     self.Du[i] = 1

        self.Q = np.array([
            [1,0,0],
            [0,1,0],
            [0,0,1]
        ])
        self.R = np.array([
            [1,0],
            [0,1]
        ])
        self.Dx = np.array([
            [1],
            [1],
            [1]
        ])
        self.Du = np.array([
            [1],
            [1]
        ])
        

    def gener_drivatives(self,x_new,x_old,u_new,u_old,xn,un,u_num):
        Lx = np.zeros((u_num,xn,1))
        Lu = np.zeros((u_num,un,1))
        Lxx = np.zeros((u_num,xn,xn))
        Luu = np.zeros((u_num,un,un))
        Lux = np.zeros((u_num,un,xn))
        Lxu = np.zeros((u_num,xn,un))
        for i in range(u_num):
            # Lu[i] = self.R @ (u_new[i] - u_old[i]).reshape((un,1)) + self.Du
            Lu[i] = self.Du
            Luu[i] = self.R
            Lxu[i] = np.zeros((xn,un))
            Lux[i] = np.zeros((un,xn))
            # Lx[i] =  self.Q @ (x_new[i] - x_old[i]).reshape((xn,1)) + self.Dx
            Lx[i]  = self.Dx
            Lxx[i] = self.Q
        return Lx,Lu,Lxx,Luu,Lux,Lxu

    def opt(self,f_AB,ref_X,ref_U,iter_num):
        x_init = ref_X[0]
        x_num = len(ref_X)
        u_num = len(ref_U) 
        xn = len(ref_X[0])
        un = len(ref_U[0])
        A = np.zeros((x_num-1,xn,xn))
        B = np.zeros((x_num-1,xn,un))
        x_old = ref_X
        u_old = ref_U
        x_new = ref_X
        u_new = ref_U
        # 构建误差模型二次型的Q,R矩阵
        for i in range(iter_num):
            print("iter num:",i)
            for ind in range(len(ref_X)-1):
                _A,_B = f_AB(x_new[ind],u_new[ind]) # 参考点怎么选？
                A[ind] = _A
                B[ind] = _B
            Lx,Lu,Lxx,Luu,Lux,Lxu = self.gener_drivatives(x_new,x_old,u_new,u_old,xn,un,u_num) 
            # 输入都是最原始的形式，模型的误差形式在算法内部构建，输出也是标准形式，而非误差量。
            print("ILQR 输出：")
            print("refX",ref_X[:4])
            print("refU",ref_U[:4])
            print("A",A[:4])
            print("B",B[:4])
            print("Lxx",Lxx[:4])
            print("Luu",Luu[:4])
            print("Lx",Lx[:4])
            print("Lu",Lu[:4])
            print("Lxu",Lxu[:4])
            print("Lux",Lux[:4])
            x_tmp,u_tmp = self.solve_iter(x_new[0],x_new[1:],u_new,A,B,Lxx,Luu,Lx,Lu,Lxu,Lux,xn,un,u_num)
            drwa_util.draw_state(x_new,"x-y-psi","opt_traj")
            drwa_util.draw_contour(u_new,"delta-v","opt_traj_U")
            drwa_util.draw()
            # 退出条件: 误差小于阈值
            if np.linalg.norm(x_tmp - x_new) < 1e-5:
                print("Converged！")
                break
            x_old = x_new
            u_old = u_new
            x_new = x_tmp
            u_new = u_tmp 
        return x_tmp,u_tmp
