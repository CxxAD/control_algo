"""
LQR 求解器，黎卡提递归.

理解：
1. 代价导数求的是给定状态下的最优控制量，即，当前状态下，使用控制量U*,可以使得代价函数最小。
于是，当前状态下的期望控制量就得到了，此时，如果假定代价函数为x-x_r 的SME，那么，就可以求的当前状态下的最优控制量，实地的最近接参考。
因此，此时参考点的作用是作为期望。
那么，该怎么选择参考点呢？
方案一：选择最近的点作为参考点，这样，可以保证当前状态下的最优控制量是最接近的。
方案二：选择对应索引的点作为参考点，这样，可以保证和整个轨迹在逐步接近。

Q:
1. 轨迹点间gap过大的话，可能导致连续很多次求解都是同一个参考点，从而规划出一个逐渐减速停止的轨迹，始终无法到达目标点，如何解决？
速度不采用LQR控制，可以保证速度量不会因为参考点的选择而接近0，而是接近参考点的速度，从而可以保证车辆能够突破gap过大导致的停止问题。
2. 动规形式LQR，如果动力学模型的线性化采用泰勒展开，则A,B需要参考点来做展开，导致只能使用参考轨迹ind做参考点，从而误差较大；且如果轨迹有突变，前后两个点无法光滑过渡，则很可能无解。
A：LQR的局限性，能使用的信息有限，只能依靠起点、终点两个点来进行优化，一旦模型本身需要展开点信息，则要求参考轨迹相对合理，能光滑过渡才行。理论上，迭代解法应该可以解决该问题，<待尝试>。
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
            C_uu_inv = np.linalg.inv(self.LM(C_uu))
            K = -C_uu_inv @ C_ux 
            k = -C_uu_inv @ C_u

            k_seq[ind] = k
            K_seq[ind] = K
            # 更新V 3*3 + 3*3   + 
            V_xx = C_xx + K.T @ C_ux + C_xu @ K + K.T @ C_uu @ K
            V_x = C_x + C_xu @ k + K.T @ C_u  + K.T @ C_uu @ k
            # const = const + k.T @ C_u + 0.5 * k.T @ C_uu @ k
            # 这里这个ind，怎么选择？
            C_xx = Q_xx[ind] + A[ind].T @ V_xx @ A[ind]
            C_xu = Q_xu[ind] + A[ind].T @ V_xx @ B[ind]
            C_ux = Q_ux[ind] + B[ind].T @ V_xx @ A[ind]
            C_uu = Q_uu[ind] + B[ind].T @ V_xx @ B[ind]
            C_x = D_x[ind] + A[ind].T @ V_x
            C_u = D_u[ind] + B[ind].T @ V_x
        return k_seq,K_seq

            
    def _forward_pass(self,x_init,k_seq,K_seq,A,B,X,U,num,func):
        x_new_seq = np.zeros((num,len(x_init)))
        u_new_seq = np.zeros((num-1,len(k_seq[0])))
        x_new_seq[0] = x_init
        for ind in range(num-1):
            u_new = (k_seq[ind].flatten() + (K_seq[ind] @ (x_new_seq[ind] - X[ind]).T) + U[ind]).T
            print(x_new_seq[ind],u_new)
            x_new_seq[ind+1] = func(x_new_seq[ind],u_new) # 用原非线性函数更新状态量，减小误差。
            u_new_seq[ind] = u_new
        return x_new_seq,u_new_seq

    def solve_iter(self,x_init,x,u,A,B,lxx,luu,l_x,l_u,lxu,lux,xn,un,num,func):
        """
        迭代求解
        """
        k_seq,K_seq = self._backward_pass(A,B,lxx,luu,l_x,l_u,lxu,lux,xn,un,num)
        # 因为是用的x,u做的展开点，所以这里生成的时候，还是用的x,u 
        x_new_seq,u_new_seq = self._forward_pass(x_init,k_seq,K_seq,A,B,x,u,num+1,func)
        return x_new_seq,u_new_seq

