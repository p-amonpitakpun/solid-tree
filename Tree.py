from solid import translate, rotate, union
from solid import cylinder, sphere

import numpy as np


def decay(y_0, y_N, k):
    def y(x):
        return y_N - (y_N - y_0) * np.exp(- k * x)
    return y


def percent(y_0, y_N, i, N):
    def y(x):
        return y_0 + (y_N - y_0) * x / N
    return y

class Branch:
    def __init__(self, btype, origin, orientation, i, N, r_0=2):
        self.btype = btype
        self.origin = origin
        self.orientation = orientation
        self.i = i
        
        self.cylinders = []
        self.N = N

        L_0 = 50
        L_N = 0
        self.L = percent(L_0, L_N, i, N)

        self.r_0 = r_0
        r_N = 0
        self.r = percent(r_0, r_N, i, N)

        theta_0 = 0
        theta_N = 45
        self.theta = percent(theta_0, theta_N, i, N)

        self.r1 = r_0 if i == 0 else self.r(i)
        self.L_list = []

        self.Theta = [orientation[1]]
        self.Phi = [orientation[2]]
        self.X = [origin[0]]
        self.Y = [origin[1]]
        self.Z = [origin[2]]
        
        self.config = []
        
    def grow(self, i):
        l = self.L(i)
        self.r2 = self.r(i + 1)
        t = 15 - 30 * np.random.rand()
        phi = 180 - 360 * np.random.rand()
        z = l * np.cos(np.deg2rad(np.sum(self.Theta) + t))
        w = l * np.sin(np.deg2rad(np.sum(self.Theta) + t))
        x = w * np.cos(np.deg2rad(np.sum(self.Phi) + phi))
        y = w * np.sin(np.deg2rad(np.sum(self.Phi) + phi))

#         if i <= self.i + 1:
#             print(f'btype: {self.btype} \t L: {l:.03f} \t r: {self.r1:.03f} , {self.r2:.03f} \t t: {t:.03f}')

        translation = [sum(self.X), sum(self.Y), sum(self.Z)]
        rotation = [0, sum(self.Theta) + t, sum(self.Phi) + phi]
        
        c = translate(translation)(
            rotate(rotation)(
                union()(
                    cylinder(h=l, r1=self.r1, r2=self.r2),
                    sphere(self.r1)
                )
            )
        )
    
#         self.config.append({
#             'translation': translation,
#             'rotation': rotation,
#             'cylinder': { 'h': l, 'r1': self.r1, 'r2': self.r2 },
#             'sphere': { 'radius': self.r1 },
#         })
        self.config.append([
            translation,
            rotation,
            [l, self.r1, self.r2],
            [self.r1],
        ])


        self.cylinders.append(c.copy())
        self.L_list.append(l)
        self.r1 = self.r2

        self.Theta.append(t)
        self.Phi.append(phi)
        self.X.append(x)
        self.Y.append(y)
        self.Z.append(z)
    
    def branching(self, i, N):
            origin = [sum(self.X[: -1]), sum(self.Y[: -1]), sum(self.Z[: -1])]
            orientation = [0, sum(self.Theta[: -1]) + 90 - 45 * np.random.rand(), sum(self.Phi[: -1]) + 90 * np.random.rand()]
#             print(f'new branch {i} {N} from {self.btype}')
            return Branch(self.btype + 1, origin, orientation, i, N, r_0=self.r_0 * 0.75)
    
    
class Tree:
    def __init__(self, origin, orientation):
        self.origin = origin
        self.orientation = orientation
        
        self.branches = []
        
    def render(self, N):
        config = {
            'r': 2
        }
        trunk = Branch(0, self.origin, self.orientation, 0, N)
        self.branches.append(trunk)

        for i in range(N):
            new_branches = []
            for branch in self.branches:
#                 print(i, 'grow existed', branch)
                branch.grow(i)

                if i > 0 and self.is_branching:
                    new_branch = branch.branching(i, N)
#                     print(i, 'grow new', new_branch, 'from', branch)
                    if new_branch is not None:
                        new_branches.append(new_branch)
                        
            for branch in new_branches:
                branch.grow(i)
                self.branches.append(branch)
                
    def getParts(self):
#         print(self.branches)
        parts = []
        for branch in self.branches[: ]:
            for cylinder in branch.cylinders:
                parts.append(cylinder)
        return parts
    
    def getConfig(self):
        all_config = []
        for branch in self.branches:
            for config in branch.config:
                all_config.append(config)
        return all_config

    def is_branching(self, i, N):
#         k = - np.log(0.05) / N
#         chance = 1 - np.exp(- k * i)
        chance = percent(0.125, 0.50, i, N)()
        return np.random.choice([False, True], p=[1 - chance, chance])
    