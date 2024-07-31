"""

"""

PI = 3.1415926

def normalize_angle(angle):
    while angle > 2 * PI:
        angle -= 2 * PI
    while angle < -2*PI:
        angle += 2 * PI
    return angle