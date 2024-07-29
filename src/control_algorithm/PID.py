# This file is used to define the PID controller class
import matplotlib.pyplot as plt 
class PID():
    def __init__(self,Kp,Ki,Kd,target,upper,lower) -> None:
        self.Kp=Kp
        self.Ki = Ki
        self.Kd = Kd
        self.target = target
        self.upper = upper
        self.lower = lower

        self.error = 0.0 
        self.pre_error = 0.0 
        self.sum_error = 0.0 

    def set_target(self,target):
        self.target = target

    def set_k(self,Kp,Ki,Kd):
        self.Kp = Kp
        self.Ki = Ki
        self.Kd = Kd

    def set_bound(self,upper,lower):
        self.upper = upper
        self.lower = lower

    def cal_output(self,state):
        self.error = self.target - state 
        u = self.error * self.Kp + self.sum_error * self.Ki + (self.error - self.pre_error) * self.Kd
        if u > self.upper:
            u = self.upper
        elif u < self.lower:
            u = self.lower
        self.pre_error = self.error
        self.sum_error += self.error
        # print("error >>> ",self.pre_error,self.sum_error)
        return u 

    def reset(self):
        self.error = 0.0
        self.pre_error = 0.0
        self.sum_error = 0.0


    def set_sum_error(self,sum_error):
        self.sum_error = sum_error

