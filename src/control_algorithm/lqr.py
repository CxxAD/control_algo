"""
LQR 求解器，黎卡提递归.
"""
import numpy as np 

class LQR:

    def __init__(self) -> None:
        pass

    def solve(self, A, B, Q, R):
        """
        求解离散时间LQR问题
        """
        P = Q
        max_iter = 150
        for i in range(max_iter):
            P_last = P
            P = Q + A.T @ P @ A - A.T @ P @ B @ np.linalg.inv(R + B.T @ P @ B) @ B.T @ P @ A
            if np.linalg.norm(P - P_last) < 1e-6:
                break
        K = np.linalg.inv(R + B.T @ P @ B) @ B.T @ P @ A
        return K, P
    
    def LM(self,C_uu):
        """
        LM算法
        """
        lam = 1e-3
        while True:
            try:
                np.linalg.inv(C_uu + lam * np.eye(C_uu.shape[0]))
                break
            except:
                lam *= 10
        return C_uu + lam * np.eye(C_uu.shape[0])
    
    def _backward_pass(self,A,B,Q_xx,Q_uu,D_x,D_u,Q_xu,Q_ux,xn,un,num):
        """
        动规求解
        """ 
        C_xx = Q_xx[-1]
        C_xu = Q_xu[-1]
        C_ux = Q_ux[-1]
        C_uu = Q_uu[-1] 

        C_x = D_x[-1]
        C_u = D_u[-1]

        V_xx = C_xx 
        V_x = C_x
        k_seq = np.zeros((num,un,1)) # 2 * 1 
        K_seq = np.zeros((num,un,xn)) # 2 * 3 
        for ind in range(num-1,-1,-1):
            C_uu = self.LM(C_uu)
            K = -np.linalg.inv(C_uu) @ C_ux 
            k = -np.linalg.inv(C_uu) @ C_u

            k_seq[ind] = k
            K_seq[ind] = K
            # 更新V 3*3 + 3*3   + 
            V_xx = C_xx + K.T @ C_ux + C_xu @ K + K.T @ C_uu @ K
            V_x = C_xu @ k + K.T @ C_u + C_x + K.T @ C_uu @ k
            # const = const + k.T @ C_u + 0.5 * k.T @ C_uu @ k

            C_xx = Q_xx[ind] + A[ind].T @ V_xx @ A[ind]
            C_xu = Q_xu[ind] + A[ind].T @ V_xx @ B[ind]
            C_ux = Q_ux[ind] + B[ind].T @ V_xx @ A[ind]
            C_uu = Q_uu[ind] + B[ind].T @ V_xx @ B[ind]
            C_x = D_x[ind] + A[ind].T @ V_x
            C_u = D_u[ind] + B[ind].T @ V_x
        return k_seq,K_seq

            
    def _forward_pass(self,x_init,k_seq,K_seq,A,B,x_ref,u_ref,num):
        x_new_seq = np.zeros((num,len(x_init)))
        u_new_seq = np.zeros((num-1,len(k_seq[0])))
        x_new_seq[0] = x_init
        for ind in range(num-1):
            # 拿最合适的做展开点or参考点 相对于展开点的偏差，
            u_new = (k_seq[ind].flatten() + (K_seq[ind] @ (x_new_seq[ind] - x_ref[ind]).T) + u_ref[ind]).T
            x_new_seq[ind+1] = A[ind] @ (x_new_seq[ind] - x_ref[ind]).T + B[ind] @ u_new + x_ref[ind].T
            u_new_seq[ind] = u_new
        return x_new_seq,u_new_seq

    def solve_iter(self,x_init,x_ref,u_ref,A,B,lxx,luu,l_x,l_u,lxu,lux,xn,un,num):
        """
        迭代求解
        """
        k_seq,K_seq = self._backward_pass(A,B,lxx,luu,l_x,l_u,lxu,lux,xn,un,num)
        x_new_seq,u_new_seq = self._forward_pass(x_init,k_seq,K_seq,A,B,x_ref,u_ref,num+1)
        return x_new_seq,u_new_seq

