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
import numpy as np

class iLQR():
    def __init__(self,A,B,Q,R,Dx,Du,xn,un) -> None:
        self.A = A  # 参考点香菇n
        self.B = B 
        self.Q = Q
        self.R = R 
        self.Dx = Dx
        self.Du = Du
        self.xn = xn 
        self.un = un 
    # 使用对应索引的点作为参考点？或者使用参考轨迹上最近的点作为参考点？
    def backword(self,ref_X,ref_U):
        # 生成 AB 
        u_num = len(ref_U)
        x_num = len(ref_X)
        if x_num -1 != u_num:
            raise ValueError("状态和控制数目不匹配")
        xn = len(ref_X[0])
        un = len(ref_U[0])
        lx = np.zeros(xn)
        lu = np.zeros(un)
        np.fill_diagonal(lx,1)
        np.fill_diagonal(lu,1)
        lxx = np.eye(xn)
        luu = np.eye(un)
        lxu = np.zeros((xn,un))
        lux = np.zeros((un,xn))

        v_x = lx
        v_xx = lxx 

        k_seq = np.zeros((u_num,un))
        K_seq = np.zeros((u_num,un,xn))
        for ind in range(len(ref_U)-1,-1,-1):
            # 反向传播
            A = self.A(ref_X[ind],ref_U[ind])
            B = self.B(ref_X[ind],ref_U[ind])
            # 计算最优控制率 

            lx = lx + A.T @ v_x
            lu = lu + B.T @ v_x
            lxx = lxx + A.T @ v_xx @ A
            luu = luu + B.T @ v_xx @ B
            lux = lux + B.T @ v_xx @ A
            lxu = lxu + A.T @ v_xx @ B

            k = -np.linalg.inv(luu) @ lu
            K = -np.linalg.inv(luu) @ lux
            v_x = lx + K.T @ luu @ k + K.T @ lu + lxu @ k
            v_xx = lxx + K.T @ luu @ K + K.T @ lux + lxu @ K

            k_seq[ind] = k
            K_seq[ind] = K


    def forward(self,k,K,ref_X,ref_U):
        X = np.zeros((len(ref_X),len(ref_X[0])))
        X[0] = ref_X[0]
        init_X = np.array(ref_X[0])
        for ind in range(len(ref_U)):
            opt_u = ref_U[ind] + k + K @ (X[ind] - ref_X[ind])
            X[ind+1] = self.A(ref_X[ind],opt_u) @ ( X[ind] - ref_X[ind]) + self.B(ref_X[ind],opt_u) @ opt_u + X[ind]
    
    def gener_drivatives(self,ref_X,ref_U,u_num):
        Lx = np.zeros((u_num,self.xn,1))
        Lu = np.zeros((u_num,self.un,1))
        Lxx = np.zeros((u_num,self.xn,self.xn))
        Luu = np.zeros((u_num,self.un,self.un))
        Lux = np.zeros((u_num,self.un,self.xn))
        Lxu = np.zeros((u_num,self.xn,self.un))
        for i in range(u_num):
            Lx[i] =self.Q @ ref_X[i] + self.Dx
            Lu[i] = self.R @ ref_U[i] + self.Du
            Lxx[i] = self.Q
            Luu[i] = self.R
            Lxu[i] = np.zeros((self.xn,self.un))
            Lux[i] = np.zeros((self.un,self.xn))
        return Lx,Lu,Lxx,Luu,Lux,Lxu

    def opt(self,ref_X,ref_U):
        x_num = len(ref_X)
        u_num = len(ref_U) 
        # 构建误差模型二次型的Q,R矩阵
        Lx,Lu,Lxx,Luu,Lux,Lxu = self.gener_drivatives(ref_X,ref_U,u_num)
        # LQR求解器。