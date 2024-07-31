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

"""
iLQR依托误差来实现优化

"""

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
            [1,0,0,0],
            [0,1,0,0],
            [0,0,1,0],
            [0,0,0,0]
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

    def _find_cloest_point(self,x,ref_X):
        """
        找到最近点
        """
        d_x = []
        d_y = []
        d = []
        for k in range(len(ref_X)):
            if k == len(ref_X) - 1:
                break
            d_x.append(ref_X[k][0]- x[0])
            d_y.append(ref_X[k][1] -x[1])
            d.append((d_x[k] ** 2 + d_y[k] ** 2)**0.5)
        # 
        min_d = min(d)
        ref_ind = d.index(min_d)
        return ref_ind

    def gener_drivatives(self,X,U,x_ref,u_ref,xn,un,u_num):
        Lx = np.zeros((u_num,xn,1))
        Lu = np.zeros((u_num,un,1))
        Lxx = np.zeros((u_num,xn,xn))
        Luu = np.zeros((u_num,un,un))
        Lux = np.zeros((u_num,un,xn))
        Lxu = np.zeros((u_num,xn,un))
        # 只要求跟好状态量，因此，控制量期望可以为0.
        for i in range(u_num):
            # 控制量期望为0，所以 U[i] - [0,0..]
            _u = np.array([U[i][0] - 0,U[i][1] - 2])
            Lu[i] = 2*self.R @ _u.reshape((un,1))
            Luu[i] = 2*self.R
            Lxu[i] = np.zeros((xn,un))
            Lux[i] = np.zeros((un,xn))
            c_ind = self._find_cloest_point(X[i],x_ref)
            # 以参考为期望，合理的。-->以最近点为期望也行，但是就变成迭代的过程了。
            _x = np.array([X[i][0] - x_ref[c_ind][0],X[i][1] - x_ref[c_ind][1],X[i][2] - 2,0])
            Lx[i] =  2*self.Q @ _x.reshape((xn,1)) 
            Lxx[i] = 2*self.Q
        return Lx,Lu,Lxx,Luu,Lux,Lxu

    def gener_init_traj(self,X_0,U,func):
        X = np.zeros((len(U)+1,len(X_0)))
        X[0] = X_0
        for i in range(len(U)):
            # print(func(X[i],U[i]))
            X[i+1] = func(X[i],U[i])
        return X

    def opt(self,f_AB,ref_X,ref_U,iter_num,func,model):
        x_num = len(ref_X)
        u_num = len(ref_U) 
        xn = len(ref_X[0])
        un = len(ref_U[0])
        A = np.zeros((x_num-1,xn,xn))
        B = np.zeros((x_num-1,xn,un))
        
        X_0 = ref_X[0]
        U = np.zeros((u_num,un))
        # 生成初始解
        U[:,0] = 0.5
        X = self.gener_init_traj(X_0,U,func)  # 

        # 构建误差模型二次型的Q,R矩阵
        for i in range(iter_num):
            print("iter num:",i)
            # 重置运动学模型：
            # model.reset(*ref_X[0],ref_U[0][1])
            drwa_util = DrawUtil()
            for ind in range(len(ref_X)-1):
                _A,_B = f_AB(X[ind+1],U[ind]) # 优化轨迹自身作为展开点,通过不断迭代来逼近参考轨迹。
                A[ind] = _A
                B[ind] = _B
            Lx,Lu,Lxx,Luu,Lux,Lxu = self.gener_drivatives(X[1:],U,ref_X[1:],ref_U,xn,un,u_num) 
            # 输入都是最原始的形式，模型的误差形式在算法内部构建，输出也是标准形式，而非误差量。
            print("ILQR 输出：")
            # 返回优化轨迹
            drwa_util.draw_state(X,"x-y-psi","ref")
            drwa_util.draw_contour(U,"delta-v","opt_traj_U")
            x_tmp,u_tmp = self.solve_iter(X[0],X[1:],U,A,B,Lxx,Luu,Lx,Lu,Lxu,Lux,xn,un,u_num,func)
            drwa_util.draw_state(x_tmp,"x-y-psi","opt_traj")
            drwa_util.draw_contour(u_tmp,"delta-v","opt_traj_U")
            drwa_util.draw()
            # 退出条件: 误差小于阈值
            if np.linalg.norm(x_tmp - X) < 1e-5:
                print("Converged！")
                break
            X = x_tmp
            U = u_tmp
        return x_tmp,u_tmp
