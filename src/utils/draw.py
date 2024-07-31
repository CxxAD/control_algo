"""
提供需要的绘图工具 
"""

# [x,y,v,psi]
import matplotlib.pyplot as plt
from numpy import size 


class DrawUtil():
    def __init__(self) -> None:
        self.ax_xy = plt.subplot(3,2,1)
        self.ax_xy.set_xlabel("x-y",fontsize=12)
        self.ax_v = plt.subplot(3,2,3)
        self.ax_v.set_xlabel("v",fontsize=12)
        self.ax_psi = plt.subplot(3,2,5)
        self.ax_psi.set_xlabel("psi",fontsize=12)
        self.ax_a = plt.subplot(3,2,2)
        self.ax_a.set_xlabel("a",fontsize=12)
        self.ax_delta = plt.subplot(3,2,4)
        self.ax_delta.set_xlabel("delta",fontsize=12)
        self.ax_raw_rate = plt.subplot(3,2,6)
        self.ax_raw_rate.set_xlabel("raw_rate",fontsize=12)
        plt.subplots_adjust(wspace = 0.5,hspace = 0.5)

    def draw_contour(self,traj,lab,name):
        state_eles = lab.split('-')
        if len(state_eles) == 2:
            if state_eles[1] == 'a':
                a = [traj[i][1] for i in range(len(traj))]
                self.ax_a.plot([i*0.1 for i in range(len(a))],a,label=name)
            elif state_eles[1] == 'delta':
                delta = [traj[i][1] for i in range(len(traj))]
                self.ax_delta.plot([i*0.1 for i in range(len(delta))],delta,label=name)
            elif state_eles[1] == 'v':
                v = [traj[i][1] for i in range(len(traj))]
                self.ax_v.plot([i*0.1 for i in range(len(v))],v,label=name)
            elif state_eles[1] == 'raw_rate':
                raw_rate = [traj[i][1] for i in range(len(traj))]
                self.ax_raw_rate.plot([i*0.1 for i in range(len(raw_rate))],raw_rate,label=name)
            else:
                raise ValueError("Please input correct state")
            if state_eles[0] == 'a':
                a = [traj[i][0] for i in range(len(traj))]
                self.ax_a.plot([i*0.1 for i in range(len(a))],a,label=name)
            elif state_eles[0] == 'delta':
                delta = [traj[i][0] for i in range(len(traj))]
                self.ax_delta.plot([i*0.1 for i in range(len(delta))],delta,label=name)
            elif state_eles[0] == 'v':
                v = [traj[i][0] for i in range(len(traj))]
                self.ax_v.plot([i*0.1 for i in range(len(v))],v,label=name)
            elif state_eles[1] == 'raw_rate':
                raw_rate = [traj[i][1] for i in range(len(traj))]
                self.ax_raw_rate.plot([i*0.1 for i in range(len(raw_rate))],raw_rate,label=name)
            else:
                raise ValueError("Please input correct state")
        
    def draw_state(self,traj,lab,name):
        state_eles = lab.split('-')
        # 必须x-y起头
        x = [traj[i][0] for i in range(len(traj))]
        y = [traj[i][1] for i in range(len(traj))]
        self.ax_xy.scatter(x,y,s=0.1,label=name)
        if state_eles[2] == 'v':
            v = [traj[i][2] for i in range(len(traj))]
            self.ax_v.plot([i*0.1 for i in range(len(v))],v,label=name,marker='o')
        elif state_eles[2] == 'psi':
            psi = [traj[i][2] for i in range(len(traj))]
            self.ax_psi.plot([i*0.1 for i in range(len(psi))],psi,label=name,marker='.')
        else:
            raise ValueError("Please input correct state")
        if len(state_eles) == 4:
            if state_eles[3] == 'v':
                v = [traj[i][3] for i in range(len(traj))]
                self.ax_v.plot([i*0.1 for i in range(len(v))],v,label=name,marker='o')
            elif state_eles[3] == 'psi':
                psi = [traj[i][3] for i in range(len(traj))]
                self.ax_psi.plot([i*0.1 for i in range(len(psi))],psi,label=name,marker='.')
            else:
                raise ValueError("Please input correct state")
            
    def draw(self):
        self.ax_psi.legend()
        self.ax_v.legend()
        self.ax_xy.legend()
        self.ax_a.legend()
        self.ax_delta.legend()
        plt.show()

