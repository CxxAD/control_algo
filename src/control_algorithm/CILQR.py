from iLQR import iLQR

class CILQR(iLQR):
    def __init__(self) -> None:
        pass 

    """
    依据需求构建对应的导数:
    [
        [lxx, lxu],
        [lux,,luu]
    ]

    [lx, lu].T
    """
    def gener_drivatives(self):
        pass 


    def fix(self,f_AB,ref_X,ref_U,iter_num):
        return self.opt(f_AB,ref_X,ref_U,iter_num)