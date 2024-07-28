
## 车辆模型设计对应说明文档

### 1.1 单车模型 
单车运动学模型描述(非线性)：
$$
\begin{aligned}
\dot{x} &= v \cdot \cos(\theta) \\
\dot{y} &= v \cdot \sin(\theta) \\
\dot{\psi} &= \frac{v}{L} \cdot \tan(\delta) \\
\dot{v} &= a
\end{aligned}
$$

变量：$x,y,v,\psi,\delta,a$

*不同状态量、控制量的选择可以构建不同精度的模型*

### 1.2 非线性模型线性化过程分析

假设状态量$X$,控制量$U$，状态方程为：
非线性模型表示为:
$$
f(X,U)
$$
为了实现线性化，对模型在*参考点*做一阶泰勒展开，模型具有一阶精度：
$$
\begin{aligned}
f(X,U) &\approx f(X_r,U_r) + \frac{\partial f}{\partial X_r} \cdot (X-X_r) + \frac{\partial f}{\partial U_r} \cdot (U-U_r) + O(X^2,U^2)\\
f(X,U) - f(X_r,U_r) &= \frac{\partial f}{\partial X_r} \cdot (X-X_r) + \frac{\partial f}{\partial U_r} \cdot (U-U_r)\\
\end{aligned}
$$
在展开后的近似下，如果给定参考点$X_0,U_0$，则在参考点附近，可以用线性模型近似原模型。又因为导数的加减等于加减的导数，于是，设置新的状态量:$\hat{X} = X - X_r$，控制量U的增量$\hat{U} = U - U_r$，则有：
$$
\begin{aligned}
\dot{\hat{X}} &= \frac{\partial f}{\partial X_r} \cdot \hat{X} + \frac{\partial f}{\partial U_r} \cdot \hat{U} \\
\frac{\hat{X_{k+1}} - \hat{X_{k}}}{T} &= A \cdot \hat{X_k} + B \cdot \hat{U_k} \\
X_{k+1} &= (AT+I) \cdot \hat{X_k} + BT \cdot \hat{U_k}
\end{aligned}
$$


#### 1.1.1 模型1 
状态量：$x,y,\psi$  车身横摆角
控制量: $\delta$ 前轮转角
该模型设计没有纳入速度、加速度的概率，很明显只能用于横向的控制跟踪。

非线性模型：
$$
\begin{aligned}
\dot{x} &= v \cdot \cos(\psi) \\
\dot{y} &= v \cdot \sin(\psi) \\
\dot{\psi} &= \frac{v}{L} \cdot \tan(\delta) \\
\end{aligned}
$$

线性化模型：
$$
\begin{aligned}
AT+I = \begin{bmatrix}
1 & 0 & -v \cdot \sin(\psi)T \\
0 & 1 & v \cdot \cos(\psi)T \\
0 & 0 & 1
\end{bmatrix}  \\
BT = \begin{bmatrix}
0 \\
0 \\
\frac{v}{L \cos(\delta)^2} \cdot T
\end{bmatrix}  
\end{aligned}
$$

#### 1.1.2 模型2
状态量：$x,y,\psi$  车身横摆角
控制量: $\delta,v$ 前轮转角、速度

非线性模型：同上

线性化模型：
$$
\begin{aligned}
AT+I = \begin{bmatrix}
1 & 0 & -v \cdot \sin(\psi)T \\
0 & 1 & v \cdot \cos(\psi)T \\
0 & 0 & 1
\end{bmatrix}  \\
BT = \begin{bmatrix}
0 && T\cos(\psi)\\
0 && T\sin(\psi)\\
\frac{v}{L \cos(\delta)^2} \cdot T && T\frac{\tan(\delta)}{L}
\end{bmatrix}  
\end{aligned}
$$

#### 1.1.3 模型3

状态量：$x,y,\psi,v$  车身横摆角、车速
控制量: $\delta,v$ 前轮转角、速度

非线性模型：
$$
\begin{aligned}
\dot{x} &= v \cdot \cos(\psi) \\
\dot{y} &= v \cdot \sin(\psi) \\
\dot{\psi} &= \frac{v}{L} \cdot \tan(\delta) \\
\dot{v} &= a
\end{aligned}
$$

线性化模型：
$$
\begin{aligned}
AT+I = \begin{bmatrix}
1 & 0 & -v \cdot \sin(\psi)T & \cos(\psi)T \\
0 & 1 & v \cdot \cos(\psi)T & \sin(\psi)T \\
0 & 0 & 1 & \frac{\tan(\delta)}{L}T \\
0 & 0 & 0 & 1
\end{bmatrix}  \\
BT = \begin{bmatrix}
0 && 0\\
0 && 0\\
\frac{v}{L \cos(\delta)^2} \cdot T && 0 \\
0 && T
\end{bmatrix}
\end{aligned}
$$

虽然转成了 $A,B$ 矩阵，但是这个矩阵中存在非线性项三角函数，
模型需要做进一步的简化，即，当参考点附近的$\delta$很小，可以近似$\tan(\delta) \approx \delta$，这样就可以得到一个线性模型。